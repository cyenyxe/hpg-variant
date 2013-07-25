/*
 * ld.h
 *
 *  Created on: Jul 4, 2012
 *      Author: Andrei Alic
 */
#ifndef LD_H_
#define LD_H_


/** SYSTEM HEADERS - SHOULD BE AVAILABLE ON ANY SYSTEM*/
#include <string.h>
#include <float.h>
#include <limits.h>
#include <stdbool.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <unistd.h>
#include <stdint.h>

/** CUSTOM HEADERS */
//#include <ped_file.h>
#include <bioformats/vcf/vcf_file.h>
//#include <family.h>
//include <stats.h>
#include <bioformats/vcf/vcf_stats.h>
//#include <bioformats/vcf/vcf_batch.h>
#include <containers/array_list.h>
#include <containers/list.h>
#include <commons/file_utils.h>

#include <bioformats/vcf/vcf_util.h>
#include <bioformats/vcf/vcf_file_structure.h>

//#include "hpg_vcf_tools_utils.h"
#include "globals.h"

/** Different information about a pair of SNPs, kept in a table having these properties for all
 * considered SNPs */
typedef struct pairwise_linkage {
	/** linkage disequilibrum measure for D' */
	double dprime;
	/** LOD is the log of the likelihood odds ratio, a measure of confidence in the value of D' */
	double lod;
	/** Confidence low */
	double ci_low;
	/** Confidence high*/
	double ci_high;
} pairwise_linkage;

/**********************************
* DECLARE ALL FUNCTION PROTOTYPES *
**********************************/

/**
 * Check if 2 refs array have the same elements; their length must be the same
 * @param **char alt1 The set of alternates for the first ref allele, as an array of strings
 * @param **char alt2 The set of alternates for the second ref allele as an array of strings
 * @param int numAlts1 The num of alternates for the first allele
 * @param int numAlts2 The num of alternates for the second allele
 * @return int The result of the comparison of the 2 char arrays
 */
int is_same_ref(char **alt1, char **alt2, int numAlts1, int numAlts2);
/**
 * Check whether a certain ref SNP appears in the list of alternates
 */
int is_alternate(char *ref, char **alt, int numAlts);
pairwise_linkage *compute_pairwise_linkage(array_list_t *markers_arr, int pos1, int pos2, int num_samples);
double compute_dprime(double *prob_haps, double *num_haps, double *known, double pA2, double pB2);
array_list_t *exec_gabriel(array_list_t *markers_arr, int num_samples);

/**
 * Free the markers array and their elements (deep free)
 */
void free_marker_array(marker *array, int len);

void count_haps(int em_round, double *num_haps, double *known, double *prob_haps,
		double doublehet, double const_prob);

void estimate_p(double *num_haps, double *prob_haps, double const_prob);

pairwise_linkage **generate_pairwise_linkage_tbl(array_list_t *markers_arr, int num_samples);

/**
 * Provides a method to compare 2 elements used in quicksort; It is using the sep field
 * of the marker_info struct; The order is descending
 */
inline static int compare_markers(const void *markeri1, const void *markeri2);

inline static int compare_markers_asc(const void *markeri1, const void *markeri2) ;

/**
 * Check if a marker's rating is above 0, meaning it is valid and it can be used in the processing;
 * @param marker_check The marker which will be checked to see if it is right
 * @return If the rating of the marker is above 0 or not
 */
inline static bool filter_rating(marker *marker_check);

#endif /* LD_H_ */
