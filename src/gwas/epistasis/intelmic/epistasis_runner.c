/*
 * Copyright (c) 2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2013 Ignacio Medina (ICM-CIPF)
 *
 * This file is part of hpg-variant.
 *
 * hpg-variant is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * hpg-variant is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with hpg-variant. If not, see <http://www.gnu.org/licenses/>.
 */

#include "epistasis_runner.h"


int run_epistasis(shared_options_data_t* shared_options_data, epistasis_options_data_t* options_data) {
    int ret_code = 0;
    
    // Load binary input dataset
    int num_affected, num_unaffected;
    size_t num_variants, file_len, genotypes_offset;
    
    uint8_t *input_file = epistasis_dataset_load(&num_affected, &num_unaffected, &num_variants, &file_len, &genotypes_offset, options_data->dataset_filename);
    if (!input_file) {
        LOG_FATAL_F("File %s does not exist!\n", options_data->dataset_filename);
    }
    
    uint8_t *genotypes = input_file + genotypes_offset;
    
    // Try to create the directory where the output files will be stored
    ret_code = create_directory(shared_options_data->output_directory);
    if (ret_code != 0 && errno != EEXIST) {
        LOG_FATAL_F("Can't create output directory: %s\n", shared_options_data->output_directory);
    }
    
    /*************** Precalculate the rest of variables the algorithm needs  ***************/
    
    int order = options_data->order;
    int stride = options_data->stride;
    int num_folds = options_data->num_folds;
    int num_samples = num_affected + num_unaffected;
    int num_blocks_per_dim = ceil((double) num_variants / stride);
    size_t max_num_block_coords = (size_t) pow(num_blocks_per_dim, order);
    
    LOG_INFO_F("Combinations of order %d, %d variants per block\n", order, stride);
    LOG_INFO_F("%d variants, %d blocks per dimension\n", num_variants, num_blocks_per_dim);
    
    // Precalculate which genotype combinations can be tested for a given order (order 2 -> {(0,0), (0,1), ... , (2,1), (2,2)})
    int num_genotype_permutations;
    uint8_t **genotype_permutations = get_genotype_combinations(order, &num_genotype_permutations);
    
    // Ranking of best models in each repetition
    struct heap *best_models[options_data->num_cv_repetitions];
    
    compare_risky_heap_func heap_max_func = NULL;
    compare_risky_heap_func heap_min_func = NULL;
    if (options_data->eval_mode == CV_A) {
        LOG_INFO("Using CV-a as ranking criteria");
        heap_max_func = compare_risky_heap_accuracy_max;
        heap_min_func = compare_risky_heap_accuracy_min;
    } else if (options_data->eval_mode == CV_C) {
        LOG_INFO("Using CV-c as ranking criteria");
        heap_max_func = compare_risky_heap_count_max;
        heap_min_func = compare_risky_heap_count_min;
    } else {
        LOG_FATAL("Rank criteria not specified! Must be 'count' or 'accu'");
    }
   
    // Variables to offload to the MIC
    int nThreads = shared_options_data->num_threads;
    enum evaluation_mode evalMode = options_data->eval_mode;
 
    /**************************** End of variables precalculus  ****************************/
    
    
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        LOG_INFO_F("Running cross-validation #%d...\n", r+1);
#pragma offload target (mic) in(genotypes : length(file_len - genotypes_offset)) \
                             in(genotype_permutations : length(num_genotype_permutations * order)) \
                             in(nThreads) in(evalMode)
    {
//        omp_set_num_threads(240); 

        // Initialize folds, first block coordinates, genotype combinations and rankings for each repetition
        unsigned int *testing_sizes, *training_sizes;
        int **folds = get_k_folds(num_affected, num_unaffected, num_folds, &testing_sizes);
        uint8_t *fold_masks = get_k_folds_masks(num_affected, num_unaffected, num_folds, folds, testing_sizes);
        
/*
        printf("fold_masks = {\n");
        for (int i = 0; i < num_folds; i++) {
            for (int j = 0; j < info.num_samples_with_padding; j++) {
                printf("%u ", fold_masks[i * info.num_samples_with_padding + j]);
            }
            printf("\n");
        }
        printf("}\n");
*/
        
        // Calculate size of training datasets
        training_sizes = calloc(3 * num_folds, sizeof(unsigned int));
        for (int f = 0; f < num_folds; f++) {
            training_sizes[3 * f] = num_samples - testing_sizes[3 * f];
            training_sizes[3 * f + 1] = num_affected - testing_sizes[3 * f + 1];
            training_sizes[3 * f + 2] = num_unaffected - testing_sizes[3 * f + 2];
        }
        
        // Initialize rankings for each repetition
        struct heap **ranking_risky = malloc(num_folds * sizeof(struct heap*));
        for (int i = 0; i < num_folds; i++) {
            ranking_risky[i] = malloc(sizeof(struct heap));
            heap_init(ranking_risky[i]);
        }
        
        int *block_coords = calloc(max_num_block_coords * order, sizeof(int));
        
        // Calculate all blocks coordinates
        int curr_idx = 0, next_idx = 0;
        size_t num_block_coords = 0;
        do {
            curr_idx = num_block_coords * order;
            next_idx = curr_idx + order;
            memcpy(block_coords + next_idx, block_coords + curr_idx, order * sizeof(int));
            curr_idx = next_idx;
            num_block_coords++;
        } while (get_next_block(num_blocks_per_dim, order, block_coords + curr_idx));
        
        #pragma omp parallel for num_threads(nThreads)
        for (int i = 0; i < num_block_coords; i++) {
            int my_block_coords[order];
            
            // Initialize rankings for each repetition
            struct heap **ranking_risky_local = malloc(num_folds * sizeof(struct heap*));
            for (int i = 0; i < num_folds; i++) {
                ranking_risky_local[i] = malloc(sizeof(struct heap));
                heap_init(ranking_risky_local[i]);
            }
        
            memcpy(my_block_coords, block_coords + i * order, order * sizeof(int));
            // printf("%d) cv %d, block %d %d\n", omp_get_thread_num(), r, my_block_coords[0], my_block_coords[1]);

            // ***************** Variables private to each task (block) *****************

            // Masks information (number (un)affected with padding, buffers, and so on)
            masks_info info; masks_info_init(order, COMBINATIONS_ROW_SSE, num_affected, num_unaffected, &info);

            // Scratchpad for block genotypes
            uint8_t *scratchpad[order];
            for (int s = 0; s < order; s++) {
                scratchpad[s] = _mm_malloc(stride * info.num_samples_with_padding * sizeof(uint8_t), 16);
            }
            // Genotypes for the current block
            uint8_t *block_genotypes[order];

            // Counts per genotype combination
            // Grouped by fold, then combination, then permutation, so there is spatial locality when getting confusion matrix
            int max_num_counts = 16 * (int) ceil(((double) info.num_cell_counts_per_combination * info.num_combinations_in_a_row * num_folds) / 16);
            int *counts_aff = _mm_malloc(max_num_counts * sizeof(int), 16);
            int *counts_unaff = _mm_malloc(max_num_counts * sizeof(int), 16);

            // Confusion matrix
            unsigned int conf_matrix[4];

            // Buffer for genotypes masks
            uint8_t *masks = _mm_malloc(info.num_combinations_in_a_row * info.num_masks * sizeof(uint8_t), 16);

            // Functions for ranking combinations
            compare_risky_heap_func heap_max_func_local = NULL;
            compare_risky_heap_func heap_min_func_local = NULL;
            if (evalMode == CV_A) {
                heap_max_func_local = compare_risky_heap_accuracy_max;
                heap_min_func_local = compare_risky_heap_accuracy_min;
            } else {
                heap_max_func_local = compare_risky_heap_count_max;
                heap_min_func_local = compare_risky_heap_count_min;
            }
    
            // *************** Variables private to each task (block) (end) ***************



            // -------------------- Get genotypes of block --------------------

            uint8_t *block_starts[order];
            for (int s = 0; s < order; s++) {
                block_starts[s] = genotypes + my_block_coords[s] * stride * num_samples;
            }

/*
            // Initialize first coordinate (only if it's different from the previous)
            block_genotypes[0] = get_genotypes_of_block_coord(num_variants, num_samples, info, stride,
                                                              my_block_coords[0], block_starts[0], scratchpad[0]);

            // Initialize the rest of coordinates
            for (int m = 1; m < order; m++) {
                bool already_present = false;
                // If any coordinate is the same as a previous one, don't copy, but reference directly
                for (int n = 0; n < m; n++) {
                    if (my_block_coords[m] == my_block_coords[n]) {
//                             printf("taking %d -> %d\n", n, m);
                        block_genotypes[m] = block_genotypes[n];
                        already_present = true;
                        break;
                    }
                }

                if (!already_present) {
                    // If not equals to a previous one, retrieve data
                    block_genotypes[m] = get_genotypes_of_block_coord(num_variants, num_samples, info, stride,
                                                                      my_block_coords[m], block_starts[m], scratchpad[m]);
                }
            }
*/


//            printf("padded block (%d*%d) = {\n", stride, info.num_samples_with_padding);
//            for (int m = 0; m < MIN(stride, num_variants); m++) {
//                for (int n = 0; n < info.num_samples_with_padding; n++) {
//                    printf("%d ", block_genotypes[0][m * info.num_samples_with_padding + n]);
//                }
//                printf("\n");
//            }
//            printf("}\n");
//
//            printf("padded block (%d*%d) = {\n", stride, info.num_samples_with_padding);
//            for (int m = 0; m < MIN(stride, num_variants); m++) {
//                for (int n = 0; n < info.num_samples_with_padding; n++) {
//                    printf("%d ", block_genotypes[1][m * info.num_samples_with_padding + n]);
//                }
//                printf("\n");
//            }
//            printf("}\n-------------------------\n");

            // -------------------- Get genotypes of block (end) --------------------

            // Combination of variants being tested
            int comb[order];
            // Array of combinations to process in a row
            int combs[info.num_combinations_in_a_row * order];
            int cur_comb_idx = 0;

            // Test first combination in the block
            get_first_combination_in_block(order, comb, my_block_coords, stride);

/*
            do {
                memcpy(combs + cur_comb_idx * order, comb, order * sizeof(int));
                cur_comb_idx++;

                if (cur_comb_idx < info.num_combinations_in_a_row) {
                    continue; // Nothing to do until we have an amount (COMBINATIONS_ROW_SSE) of combinations ready
                }

                process_set_of_combinations(info.num_combinations_in_a_row, combs, order, stride, num_folds, fold_masks,
                                            training_sizes, testing_sizes,  block_genotypes, genotype_permutations,
                                            masks, options_data->eval_subset, info, 
                                            heap_min_func_local, counts_aff, counts_unaff, conf_matrix, 
                                            options_data->max_ranking_size, ranking_risky_local);
                
                cur_comb_idx = 0;
            } while (get_next_combination_in_block(order, comb, my_block_coords, stride, num_variants));

            
            // Process combinations out of a full set
            process_set_of_combinations(cur_comb_idx, combs, order, stride, num_folds, fold_masks,
                                        training_sizes, testing_sizes, block_genotypes, genotype_permutations,
                                        masks, options_data->eval_subset, info, 
                                        heap_min_func_local, counts_aff, counts_unaff, conf_matrix, 
                                        options_data->max_ranking_size, ranking_risky_local);

            // Insert best model of this block in the global ranking (critical section)
            for (int i = 0; i < num_folds; i++) {
                while (!heap_empty(ranking_risky_local[i])) {
                    struct heap_node *hn = heap_take(heap_min_func_local, ranking_risky_local[i]);
                    risky_combination *risky_comb = (risky_combination*) hn->value;
                    
                    int position = -1;
                    
#pragma omp critical
                    {
                        position = add_to_model_ranking(risky_comb, options_data->max_ranking_size, ranking_risky[i], heap_min_func_local);
                    }
                    if (position < 0) {
                        risky_combination_free(risky_comb);
                    }
                    free(hn);
                }
                free(ranking_risky_local[i]);
            }
            free(ranking_risky_local);
            
            _mm_free(masks);
            for (int s = 0; s < order; s++) {
                _mm_free(scratchpad[s]);
            }
            _mm_free(counts_aff);
            _mm_free(counts_unaff);
*/

            // Notify a block has been processed
            char end_block_msg[256]; memset(end_block_msg, 0, 256 * sizeof(char));
            sprintf(end_block_msg, "[%d] Block finished: (", omp_get_thread_num());
            for (int s = 0; s < order; s++) {
                sprintf(end_block_msg + strlen(end_block_msg), "%d,", my_block_coords[s] + 1);
            }
            size_t end_block_msg_len = strlen(end_block_msg);
            end_block_msg[end_block_msg_len - 1] = ')';
            // TODO Add when merged with 'next' branch because new logging system does not add newline at the end
            end_block_msg[end_block_msg_len] = '\n';
            
            printf(end_block_msg);
        }
        
/*
        for (int f = 0; f < num_folds; f++) {
            printf("Ranking fold %d = {\n", f);
            risky_combination *element = NULL;
            linked_list_iterator_t* iter = linked_list_iterator_new(ranking_risky[f]);
            while(element = linked_list_iterator_curr(iter)) {
                printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
                linked_list_iterator_next(iter);
            }
            linked_list_iterator_free(iter);
            printf("\n\n");
        }
*/
/*        
        // Merge all rankings in one
        best_models[r] = merge_rankings(num_folds, ranking_risky, heap_min_func, heap_max_func);
        
        // Show the best model of this repetition
        char *path, default_path[32];
        sprintf(default_path, "hpg-variant.cv%d.epi", r+1);
        FILE *fd = get_output_file(shared_options_data, default_path, &path);
        epistasis_report(order, r, options_data->eval_mode, options_data->eval_subset, best_models[r], options_data->max_ranking_size, heap_max_func, fd);
        fclose(fd);
*/        
        // Free data por this repetition
        for (int i = 0; i < num_folds; i++) {
            free(folds[i]);
            free(ranking_risky[i]);
        }
        free(ranking_risky);

        free(folds);
        free(testing_sizes);
        free(training_sizes);
        _mm_free(fold_masks);
    }
    }
    
    // Free data for the whole epistasis check
    for (int i = 0; i < num_genotype_permutations; i++) {
        free(genotype_permutations[i]);
    }
    free(genotype_permutations);
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        struct heap_node *hn;
        risky_combination *element = NULL;

        while (!heap_empty(best_models[r])) {
            hn = heap_take(heap_max_func, best_models[r]);
            risky_combination_free((risky_combination*) hn->value);
            free(hn);
        }
    }
    epistasis_dataset_close(input_file, file_len);
    
    return ret_code;
}
