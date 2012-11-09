/*
 ============================================================================
 Name        : haploviewc.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include "main_haplo.h"

/**
 * argv[1] = input file full path
 * argv[2] = max_batches
 * argv[3] = number of threads
 * argv[4] = entries per thread
 * argv[5] = the size of the batch
 * argv[6] = how to determine the main allele; 1 - use Haploview approach and set the main as the one
 * 			having the most appearances or 2 - use the reference and alternates from the VCF file
 */
int main( int argc, char **argv) {
    
	return haplo_main(argc, argv, "haplo.conf");
}

int haplo_main(int argc, char *argv[], const char *configuration_file)
{
    array_list_t *all_markers, *result;
	// This is the return list
	unsigned int num_samples = 0;
        
        
    /* ******************************
     *       Modifiable options     *
     * ******************************/

    shared_options_t *shared_opts = new_shared_cli_options();
    haplo_options_t *haplo_opts = new_haplo_cli_options();
    
    // Step 1: read options from configuration file
    int config_errors = read_global_configuration(configuration_file, shared_opts);
    config_errors &= read_haplo_configuration(configuration_file, haplo_opts, shared_opts);
    LOG_INFO_F("Config read with errors = %d\n", config_errors);
    
    if (config_errors) {
        return CANT_READ_CONFIG_FILE;
    }
    
     // Step 2: parse command-line options
    void **argtable = parse_haplo_options(argc, argv, haplo_opts, shared_opts);
    
    
    // Step 4: Create XXX_options_data_t structures from valid XXX_options_t
    shared_options_data_t *shared_opts_data = new_shared_options_data(shared_opts);
    haplo_options_data_t *opts_data = new_haplo_options_data(haplo_opts);
 
//    conf_params cparams = {.file_path = argv[1],
//    		.num_threads = (int) strtol(argv[3], (char **)NULL, 10),
//    		.entries_per_thread = (int) strtol(argv[4], (char **)NULL, 10),
//    		.batch_size = 20,//(int) strtol(argv[5], (char **)NULL, 10),
//    		.max_simultaneous_batches = INT_MAX};//(int) strtol(argv[2], (char **)NULL, 10)};
	// Alloc memo for all markers; this way each thread which process a block will be able to write its
	// share on unique positions
	all_markers = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

	get_markers_array(all_markers, shared_opts_data, opts_data, &num_samples);

	//clean markers
	size_t idx =0, len = all_markers->size;
	while (idx<len)
		if (((marker *) array_list_get(idx, all_markers))->rating <= 0)
		{
			array_list_remove_at(idx, all_markers);--len;
		} else
			idx++;


	result = exec_gabriel(all_markers, num_samples);
        
        printf("\nPrint blocks:\n");
	for (unsigned int idx = 0; idx< result->size; idx++) {
		array_list_t *a = (array_list_t *)array_list_get(idx, result);
		for (unsigned int idx2 = 0; idx2< a->size; idx2++) {
			int *e = (int *)array_list_get(idx2, a);
			printf("%d ", *e);
		}
		printf("\n");
	}

	array_list_free(all_markers, NULL);
	array_list_free(result, NULL);
	//free_marker_array(all_markers, samples_names->size);
        free(haplo_opts);
        //free_shared_options_data(shared_opts);
        arg_freetable(argtable, haplo_opts->num_opts + shared_opts->num_options);
        return EXIT_SUCCESS;
}

int read_haplo_configuration(const char *filename, haplo_options_t *haplo_options, shared_options_t *shared_options) {
    if (filename == NULL || haplo_options == NULL || shared_options == NULL) {
        return -1;
    }

    config_t *config = (config_t*) calloc (1, sizeof(config_t));
    int ret_code = config_read_file(config, filename);
    if (ret_code == CONFIG_FALSE) {
        LOG_ERROR_F("Configuration file error: %s\n", config_error_text(config));
        return CANT_READ_CONFIG_FILE;
    }

    const char *tmp_string;
    
    // Read number of threads that will make request to the web service
    ret_code = config_lookup_int(config, "haplo.num-threads", shared_options->num_threads->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Number of threads not found in config file, must be set via command-line");
    } else {
        LOG_DEBUG_F("num-threads = %ld\n", *(shared_options->num_threads->ival));
    }

    // Read maximum number of batches that can be stored at certain moment
    ret_code = config_lookup_int(config, "haplo.max-batches", shared_options->max_batches->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Maximum number of batches not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("max-batches = %ld\n", *(shared_options->max_batches->ival));
    }
    
    // Read size of a batch (in lines or bytes)
    ret_code = config_lookup_int(config, "haplo.batch-lines", shared_options->batch_lines->ival);
    ret_code |= config_lookup_int(config, "haplo.batch-bytes", shared_options->batch_bytes->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Neither batch lines nor bytes found in configuration file, must be set via command-line");
    } 
//     else {
//         LOG_INFO_F("batch-lines = %ld\n", *(shared_options->batch_lines->ival));
//     }
    
    // Read host URL
    ret_code = config_lookup_string(config, "haplo.url", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Web services URL not found in configuration file, must be set via command-line");
    } else {
        *(shared_options->host_url->sval) = strdup(tmp_string);
        LOG_DEBUG_F("web services host URL = %s (%zu chars)\n",
                   *(shared_options->host_url->sval), strlen(*(shared_options->host_url->sval)));
    }
    
    // Read version
    ret_code = config_lookup_string(config, "haplo.version", &tmp_string);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        *(shared_options->version->sval) = strdup(tmp_string);
        LOG_DEBUG_F("version = %s (%zu chars)\n",
                   *(shared_options->version->sval), strlen(*(shared_options->version->sval)));
    }
    
    // Read method for determining the ref and minor alleles 
    ret_code = config_lookup_int(config, "haplo.alleleref", haplo_options->alleleref->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("alleleref = %ld\n", *(haplo_options->alleleref->ival));
    }
    
    // Read HWCUT 
    ret_code = config_lookup_float(config, "haplo.hwcut", haplo_options->hwcut->dval);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("hwcut = %lf\n", *(haplo_options->hwcut->dval));
    }
    
    // Read Failed genotype cut
    ret_code = config_lookup_int(config, "haplo.fgcut", haplo_options->fgcut->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("fgcut = %lf\n", *(haplo_options->fgcut->ival));
    }
    
    // Read Number of mendelian errors for cutting
    ret_code = config_lookup_int(config, "haplo.mendcut", haplo_options->mendcut->ival);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("mendcut = %lf\n", *(haplo_options->mendcut->ival));
    }
    
    // Read Limit for the minor allele frequency when cutting
    ret_code = config_lookup_float(config, "haplo.mafcut", haplo_options->mafcut->dval);
    if (ret_code == CONFIG_FALSE) {
        LOG_WARN("Version not found in configuration file, must be set via command-line");
    } else {
        LOG_DEBUG_F("mafcut = %lf\n", *(haplo_options->mafcut->dval));
    }

    config_destroy(config);
    free(config);

    return 0;
}

void **parse_haplo_options(int argc, char *argv[], haplo_options_t *haplo_options, shared_options_t *shared_options) {
    struct arg_end *end = arg_end(haplo_options->num_opts + shared_options->num_options);
    void **argtable = merge_haplo_options(haplo_options, shared_options, end);
    
    int num_errors = arg_parse(argc, argv, argtable);
    if (num_errors > 0) {
        arg_print_errors(stdout, end, "hpg-variant");
    }
    
    return argtable;
}

void **merge_haplo_options(haplo_options_t *haplo_opts, shared_options_t *shared_options, struct arg_end *arg_end) {
    size_t opts_size = haplo_opts->num_opts + shared_options->num_options + 1;
    void **tool_options = malloc (opts_size * sizeof(void*));
    tool_options[0] = shared_options->vcf_filename;
    tool_options[1] = shared_options->ped_filename;
    tool_options[2] = shared_options->output_filename;
    tool_options[3] = shared_options->output_directory;
    
    tool_options[4] = shared_options->host_url;
    tool_options[5] = shared_options->version;
    tool_options[6] = shared_options->species;
    
    tool_options[7] = shared_options->max_batches;
    tool_options[8] = shared_options->batch_lines;
    tool_options[9] = shared_options->batch_bytes;
    tool_options[10] = shared_options->num_threads;
    tool_options[11] = shared_options->entries_per_thread;
    
    tool_options[12] = shared_options->num_alleles;
    tool_options[13] = shared_options->coverage;
    tool_options[14] = shared_options->quality;
    tool_options[15] = shared_options->region;
    tool_options[16] = shared_options->region_file;
    tool_options[17] = shared_options->snp;
    
    tool_options[18] = shared_options->config_file;
    tool_options[19] = shared_options->mmap_vcf_files;
    
    tool_options[20] = haplo_opts->alleleref;
    tool_options[21] = haplo_opts->hwcut;
    tool_options[22] = haplo_opts->fgcut;
    tool_options[23] = haplo_opts->mendcut;
    tool_options[24] = haplo_opts->mafcut;
               
    tool_options[25] = arg_end;
    
    return tool_options;
}


haplo_options_t *new_haplo_cli_options(void) {
    haplo_options_t *options = (haplo_options_t*) malloc (sizeof(haplo_options_t));
    options->alleleref = arg_int0(NULL, "alleleref", NULL, "Determine which algorithm will be used to count the alleles");
    options->hwcut = arg_dbl0(NULL, "hwcut", NULL, "Hw cut");
    options->fgcut = arg_int0(NULL, "fgcut", NULL, "Failed genotype cut");
    options->mendcut = arg_int0(NULL, "mendcut", NULL, " Number of mendelian errors for cutting");
    options->mafcut = arg_dbl0(NULL, "mafcut", NULL, "Limit for the minor allele frequency when cutting");
    options->num_opts = NUM_HAPLO_OPTS;
    return options;
}

haplo_options_data_t *new_haplo_options_data(haplo_options_t *options) {
    haplo_options_data_t *options_data = (haplo_options_data_t*) calloc (1, sizeof(haplo_options_data_t));
    options_data->alleleref = *(options->alleleref->ival);
    options_data->hwcut = *(options->hwcut->dval);
    options_data->fgcut = *(options->fgcut->ival);
    options_data->mendcut = *(options->mendcut->ival);
    options_data->mafcut = *(options->mafcut->dval);
    return options_data;
}


void free_haplo_options(haplo_options_t *options_data) {
    free(options_data);
}