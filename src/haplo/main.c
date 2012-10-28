/*
 ============================================================================
 Name        : haploviewc.c
 Author      : 
 Version     :
 Copyright   : Your copyright notice
 Description : Hello World in C, Ansi-style
 ============================================================================
 */

#include "globals.h"
#include "ld.h"
#include "file_handling.h"

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
    array_list_t *all_markers, *result;
	// This is the return list
	unsigned int num_samples = 0;
    user_params uparams = {.ref_allele_method = ALLELE_REF_HAPLOVIEW,
    		.hwCut = 0.001,
    			.failedGenoCut = 75,
    			.numMendErrCut = 1,
    			.mafCut = 0.001};
    conf_params cparams = {.file_path = argv[1],
    		.num_threads = (int) strtol(argv[3], (char **)NULL, 10),
    		.entries_per_thread = (int) strtol(argv[4], (char **)NULL, 10),
    		.batch_size = 20,//(int) strtol(argv[5], (char **)NULL, 10),
    		.max_simultaneous_batches = INT_MAX};//(int) strtol(argv[2], (char **)NULL, 10)};
	// Alloc memo for all markers; this way each thread which process a block will be able to write its
	// share on unique positions
	all_markers = array_list_new(10, 1.5, COLLECTION_MODE_SYNCHRONIZED);

	get_markers_array(all_markers, &cparams, &uparams, &num_samples);

	//clean markers
	size_t idx =0, len = all_markers->size;
	while (idx<len)
		if (((marker *) array_list_get(idx, all_markers))->rating <= 0)
		{
			array_list_remove_at(idx, all_markers);--len;
		} else
			idx++;


	result = exec_gabriel(all_markers, num_samples);
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
	return EXIT_SUCCESS;
}
