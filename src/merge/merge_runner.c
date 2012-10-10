#include "merge_runner.h"

#define TREE_LIMIT  10000

KHASH_MAP_INIT_INT64(pos, array_list_t*);

int run_merge(shared_options_data_t *shared_options_data, merge_options_data_t *options_data) {
    if (options_data->num_files == 1) {
        LOG_INFO("Just one VCF file specified, no need to merge");
        return 0;
    }
    
    list_t *read_list[options_data->num_files];
    memset(read_list, 0, options_data->num_files * sizeof(list_t*));
    list_t *output_list = (list_t*) malloc (sizeof(list_t));
    list_init("output", shared_options_data->num_threads, shared_options_data->max_batches * shared_options_data->batch_lines, output_list);
    
    int ret_code = 0;
    double start, stop, total;
    vcf_file_t *files[options_data->num_files];
    memset(files, 0, options_data->num_files * sizeof(vcf_file_t*));
    
    // Initialize variables related to the different files
    for (int i = 0; i < options_data->num_files; i++) {
        files[i] = vcf_open(options_data->input_files[i], shared_options_data->max_batches);
        if (!files[i]) {
            LOG_FATAL_F("VCF file %s does not exist!\n", options_data->input_files[i]);
        }
        
        read_list[i] = (list_t*) malloc(sizeof(list_t));
        list_init("text", 1, shared_options_data->max_batches, read_list[i]);
    }
    
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
#pragma omp parallel sections private(start, stop, total)
    {
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
            // Reading
            start = omp_get_wtime();

            ret_code = vcf_multiread_batches(read_list, shared_options_data->batch_lines, files, options_data->num_files);

            stop = omp_get_wtime();
            total = stop - start;

            if (ret_code) {
                LOG_ERROR_F("Error %d while reading VCF files\n", ret_code);
            }

            LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
        
#pragma omp section
        {
            // Enable nested parallelism and set the number of threads the user has chosen
            omp_set_nested(1);
            omp_set_num_threads(shared_options_data->num_threads);
            
            LOG_DEBUG_F("Thread %d processes data\n", omp_get_thread_num());
            
            int num_eof_found = 0;
            int eof_found[options_data->num_files];
            memset(eof_found, 0, options_data->num_files * sizeof(int));
            
            list_item_t *items[options_data->num_files];
            memset(items, 0, options_data->num_files * sizeof(list_item_t*));
            char *texts[options_data->num_files];
            memset(texts, 0, options_data->num_files * sizeof(char*));
            
            khash_t(pos) *positions_read = kh_init(pos);
            
            long max_position_merged = LONG_MAX;
            char *max_chromosome_merged = NULL;
            
            int num_chromosomes;
            char **chromosome_order = get_chromosome_order(shared_options_data->host_url, shared_options_data->species, 
                                                           shared_options_data->version, &num_chromosomes);
            
            start = omp_get_wtime();

            while (num_eof_found < options_data->num_files) {
                /* Process:
                 * - N threads getting batches of VCF records and inserting them in a data structure. The common minimum 
                 * position of each group of batches will also be stored.
                 * - If the data structure reaches certain size or the end of a chromosome, merge positions prior to the 
                 * last minimum registered.
                 */
                
// #pragma omp parallel num_threads(shared_options_data->num_threads) firstprivate(items, texts)
//                 {
// #pragma omp critical
//                 {
                // Getting text elements in a critical region guarantees that each thread gets variants in positions in the same range
                for (int i = 0; i < options_data->num_files; i++) {
                    if (eof_found[i]) {
                        continue;
                    }
                    
                    items[i] = list_remove_item(read_list[i]);
                    if (items[i] == NULL || !strcmp(items[i]->data_p, "")) {
                        LOG_INFO_F("[%d] EOF found in file %s\n", omp_get_thread_num(), options_data->input_files[i]);
                        eof_found[i] = 1;
                        num_eof_found++;
                        
                        if(items[i] != NULL && !strcmp(items[i]->data_p, "")) {
                            free(items[i]->data_p);
                            list_item_free(items[i]);
                            LOG_DEBUG_F("[%d] Text batch freed\n", omp_get_thread_num());
                        } else {
                            LOG_DEBUG_F("[%d] No need to free text batch\n", omp_get_thread_num());
                        }
                        
                        continue;
                    }
                    
                    assert(items[i]->data_p != NULL);
                    texts[i] = items[i]->data_p;
                    
//                     printf("[%d] text batch from file %d\n", omp_get_thread_num(), i);
//                     printf("[%d] contents = '%s'\n", texts[i]);
                }
//                 }
                
                for (int i = 0; i < options_data->num_files; i++) {
                    if (eof_found[i]) {
                        continue;
                    }
                    
                    char *text_begin = texts[i];
                    char *text_end = text_begin + strlen(text_begin);
                    assert(text_end != NULL);
                    
                    // Get VCF batches from text batches
                    vcf_reader_status *status = vcf_reader_status_new(shared_options_data->batch_lines, 1);
                    ret_code = run_vcf_parser(text_begin, text_end, shared_options_data->batch_lines, files[i], status);
                    
                    vcf_batch_t *batch = fetch_vcf_batch(files[i]);
                    
//                     printf("[%d] vcf batch from file %d\n", omp_get_thread_num(), i);
                    
                    // Insert records into hashtable
                    khiter_t iter;
                    int ret;
                    for (int j = 0; j < batch->records->size; j++) {
                        vcf_record_t *record = vcf_record_copy(array_list_get(j, batch->records));
                        vcf_record_file_link *link = vcf_record_file_link_new(record, files[i]);
                        array_list_t *records_in_position;
//                         merged_node *node;
                        iter = kh_get(pos, positions_read, record->position);
                        if (iter != kh_end(positions_read)) {
                            records_in_position = kh_val(positions_read, iter);
                            ret = array_list_insert(link, records_in_position);
                            assert(ret);
//                             printf("record inserted %s:%ld\n", record->chromosome, record->position);
                        } else {
                            records_in_position = array_list_new(8, 1.5, COLLECTION_MODE_SYNCHRONIZED);
                            ret = array_list_insert(link, records_in_position);
                            assert(ret);
                            
                            iter = kh_put(pos, positions_read, record->position, &ret);
                            assert(ret);
                            if (ret) {
                                kh_value(positions_read, iter) = records_in_position;
//                                 printf("New value = %ld\n", record->position);
                            }
                            
                        }
                    }
                    
                    // Update minimum position being a maximum of these batches
                    vcf_record_t *current_record = (vcf_record_t*) array_list_get(batch->records->size - 1, batch->records);
                    
// #pragma omp critical
//                     {
                    if (max_chromosome_merged == NULL) {
                        // Max merged chrom:position not set, assign without any other consideration
                        max_chromosome_merged = strndup(current_record->chromosome, current_record->chromosome_len);
                        max_position_merged = current_record->position;
                    } else {
                        char *current_chromosome = strndup(current_record->chromosome, current_record->chromosome_len);
                        long unsigned int current_position = current_record->position;
                        
//                         printf("current = %s:%ld\tmax = %s:%ld\n", current_chromosome, current_position, max_chromosome_merged, max_position_merged);
                        int chrom_comparison = compare_chromosomes(current_chromosome, max_chromosome_merged, chromosome_order, num_chromosomes);
                        int position_comparison = compare_positions(&current_position, &max_position_merged);
                        
                        // Max merged chrom:position is posterior to the last one in this batch
                        if (chrom_comparison < 0 || (chrom_comparison == 0 && position_comparison < 0)) {
                            max_chromosome_merged = current_chromosome;
                            max_position_merged = current_position;
                        } else {
                            free(current_chromosome);
                        }
                    }
                    
                    // Free batch and its contents
                    vcf_reader_status_free(status);
                    vcf_batch_free(batch);
                    list_item_free(items[i]);
//                     }
                }
                
                if (num_eof_found == options_data->num_files) {
                    max_chromosome_merged = chromosome_order[num_chromosomes-1];
                    max_position_merged = LONG_MAX;
                }
                
                // If the data structure reaches certain size or the end of a chromosome, 
                // merge positions prior to the last minimum registered TODO
                if (num_eof_found < options_data->num_files && kh_size(positions_read) > TREE_LIMIT) {
                    LOG_INFO_F("Merging until position = %s:%ld\n", max_chromosome_merged, max_position_merged);
                    
                    // Free records in the table
                    khiter_t k;
                    for (k = kh_begin(positions_read); k != kh_end(positions_read); ++k) {
                        if (kh_exist(positions_read, k)) {
                            array_list_t *records_in_position = kh_val(positions_read, k);
                            assert(records_in_position);
//                                     printf("records in position %ld = %zu (before)\n", keys[k], records_in_position->size);
                            vcf_record_file_link *link = NULL;
                            vcf_record_t *record = NULL;
                            
                            int num_links = 0;
                            vcf_record_file_link **links = NULL;
//                             vcf_record_file_link **links = malloc (records_in_position->size * sizeof(vcf_record_file_link*));
                            // Free records in a node, and whose positions are prior to the last chromosome:position to merge
                            for (int i = records_in_position->size - 1; i >= 0; i--) {
                                assert(records_in_position);
//                                 assert(records_in_position->size);
                                link = array_list_get(i, records_in_position);
                                record = link->record;
//                                     printf("record ? %d\n", record != NULL);
                                assert(record);
//                                 printf("in position %s:%ld\n", record->chromosome, record->position);
                                /* Each node is identified by its position, but several chromosomes can coexist at the same time.
                                    * Merged records must be positioned:
                                    * - In a previous chromosome
                                    * - In the same chromosome, in a previous position
                                    */
                                // TODO if there are multiple chromosomes previous to the limit, this will merge all of them together!!!
                                if (compare_chromosomes(record->chromosome, max_chromosome_merged, 
                                                        chromosome_order, num_chromosomes) < 0 ||
                                    (compare_chromosomes(record->chromosome, max_chromosome_merged, 
                                                        chromosome_order, num_chromosomes) == 0 &&
                                    compare_positions(&(record->position), &max_position_merged) <= 0)) {
                                    // TODO added to check if leaks are reduced
                                    if (!links) {
//                                         assert(records_in_position->size);
//                                         assert(sizeof(vcf_record_file_link*));
                                        assert(records_in_position->size * sizeof(vcf_record_file_link*));
                                        links = malloc (records_in_position->size * sizeof(vcf_record_file_link*));
                                        // TODO maybe? 
//                                         links = malloc ((i+1) * sizeof(vcf_record_file_link*));
                                        // That would reduce the number of positions allocated to the minimum neccessary
                                    }
//                                     printf("comparing %s:%ld against %s:%ld\n", record->chromosome, record->position, max_chromosome_merged, max_position_merged);
                                    assert(array_list_remove_at(i, records_in_position));
//                                     vcf_record_free_deep(record);
                                    links[num_links] = link;
                                    num_links++;
                                }
                            }
                            
                            // TODO launch merge
                            if (num_links > 0) {
//                             if (links) {
                                int err_code = 0;
//                                 printf("links[0] = %s:%ld in file %s\n", links[0]->record->chromosome, links[0]->record->position, links[0]->file->filename);
                                vcf_record_t *merged = merge_position(links, num_links, files, options_data->num_files, options_data, &err_code);
                                
                                if (!err_code) {
                                    list_item_t *item = list_item_new(k, MERGED_RECORD, merged);
                                    list_insert_item(item, output_list);
                                }
                                
//                                 printf("finished %s:%ld\n", record->chromosome, record->position);
                                
                                for (int i = 0; i < num_links; i++) {
                                    vcf_record_file_link_free(links[i]);
                                }
                                free(links);
                                
                            }
                            
                            // Free empty nodes (lists of records in the same position)
                            if (array_list_size(records_in_position) == 0) {
//                                     printf("position %ld to be removed!\n", keys[k]);
                                array_list_free(records_in_position, NULL);
                                kh_del(pos, positions_read, k);
                            }
                        } // End kh_exist
                    }
                }
                
                // When reaching EOF for all files, merge the remaining entries
                if (num_eof_found == options_data->num_files && kh_size(positions_read) > 0) {
                    LOG_INFO_F("Cleaning until last position = %s:%ld\n", chromosome_order[num_chromosomes-1], LONG_MAX);
                    khiter_t k;
                    for (k = kh_begin(positions_read); k != kh_end(positions_read); ++k) {
                        if (kh_exist(positions_read, k)) {
                            array_list_t *records_in_position = kh_val(positions_read, k);
                            assert(records_in_position != NULL);
                            // TODO if there are multiple chromosomes previous to the limit, this will merge all of them together!!!
                            int err_code = 0;
                            vcf_record_t *merged = merge_position(records_in_position->items, records_in_position->size, 
                                                                  files, options_data->num_files, options_data, &err_code);
                            if (!err_code) {
                                list_item_t *item = list_item_new(k, MERGED_RECORD, merged);
                                list_insert_item(item, output_list);
                            }
                            
                            array_list_free(records_in_position, vcf_record_file_link_free);
                            kh_del(pos, positions_read, k);
                        }
                    }
                }
                
                // Set variables ready for next iteration of the algorithm
                if (max_chromosome_merged) {
                    free(max_chromosome_merged);
                }
                max_chromosome_merged = NULL;
                max_position_merged = LONG_MAX;
//                 }
            }
            
            kh_destroy(pos, positions_read);
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%d] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%d] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);

            // Decrease list writers count
            for (int i = 0; i < shared_options_data->num_threads; i++) {
                list_decr_writers(output_list);
            }
        }
        
#pragma omp section
        {
            LOG_DEBUG_F("Thread %d writes the output\n", omp_get_thread_num());
    
            start = omp_get_wtime();

            // Create file streams for results
            int dirname_len = strlen(shared_options_data->output_directory);
            
            list_item_t* item = NULL;
//             merge_result_t *merge;
//             FILE *merge_fd = NULL;
//             char merge_filename[1024];
//             char input_filename[256];
//             get_filename_from_path(shared_options_data->vcf_filename, input_filename);
            
            vcf_record_t *record;
            while ((item = list_remove_item(output_list)) != NULL) {
                record = item->data_p;
//                 printf("record to free %s:%ld\n", record->chromosome, record->position);
//                 merge = item->data_p;
//                 
//                 memset(merge_filename, 0, 1024 * sizeof(char));
//                 sprintf(merge_filename, "%s/%s_%s", shared_options_data->output_directory, merge->merge_name, input_filename);
// //                 sprintf(merge_filename, "%s/%s.vcf", shared_options_data->output_directory, merge->merge_name, shared_options_data->vcf_filename);
//                 
// //                 printf("Merge filename = '%s'\n", merge_filename);
//                 
//                 merge_fd = cp_hashtable_get(output_files, merge->merge_name);
//                 if (!merge_fd) {
//                     // TODO If its the first line to write into the file, create file and include the header
//                     merge_fd = fopen(merge_filename, "w");
//                     cp_hashtable_put(output_files, merge->merge_name, merge_fd);
//                     
//                     vcf_write_to_file(file, merge_fd);
//                 }
//                 
//                 // TODO write line into the file
//                 write_record(merge->record, merge_fd);
                vcf_record_free_deep(record);
                list_item_free(item);
            }
            
            stop = omp_get_wtime();

            total = stop - start;

            LOG_INFO_F("[%dW] Time elapsed = %f s\n", omp_get_thread_num(), total);
            LOG_INFO_F("[%dW] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
        }
    }

    // TODO Free variables related to the different files
    for (int i = 0; i < options_data->num_files; i++) {
        if(files[i]) { vcf_close(files[i]); }
        if(read_list[i]) { free(read_list[i]); }
    }
    free(output_list);
    
    return ret_code;
}
