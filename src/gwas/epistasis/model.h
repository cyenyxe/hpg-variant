#ifndef EPISTASIS_MODEL
#define EPISTASIS_MODEL

#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <commons/log.h>
#include <containers/array_list.h>
#include <containers/linked_list.h>

#include "cross_validation.h"
#include "dataset.h"
#include "mdr.h"


typedef struct {
    double accuracy;
    int order;
    int *combination;
    uint8_t *genotypes;
    int num_risky_genotypes;
    void *auxiliary_info;
} risky_combination;


enum eval_mode { TESTING, TRAINING };

/**
 * @brief Functions available for evaluating the results of a MDR model
 * @details Functions available for evaluating the results of a MDR model. These are:
 * - CA: accuracy
 * - BA: balanced accuracy, used by the original MDR
 * - wBA: weighted balanced accuracy, Namkung et al. (2008) TODO
 * - Gamma: Goodman and Kruskal (1954)
 * - Tau-b: Kendall rank correlation coefficient
 **/
enum eval_function { CA, BA, wBA, GAMMA, TAU_B };


/* **************************
 *       Main pipeline      *
 * **************************/

risky_combination *get_model_from_combination_in_fold(int order, int comb[order], uint8_t *val, unsigned int num_affected_in_training, unsigned int num_unaffected_in_training,
                                                      int num_genotype_combinations, uint8_t **genotype_combinations);

double test_model(int order, risky_combination *risky_comb, uint8_t *val, 
                  unsigned int num_affected_in_testing, unsigned int num_unaffected_in_testing);

int add_to_model_ranking(risky_combination *risky_comb, int max_ranking_size, linked_list_t *ranking_risky);


/* **************************
 *          Counts          *
 * **************************/

/**
 * @brief Gets the number of ocurrences of each genotype both in affected and unaffected groups.
 * @details Gets the number of ocurrences of each genotype both in affected and unaffected groups. For using 
 * these values in a contingency table, number of not-occurrences can be calculated like the following:
 * not_occur_affected = num_affected - occur_affected
 *
 * @param order Number of SNPs combined
 * @param genotypes Genotypes of all SNPs (grouped by SNP)
 * @param genotype_combinations Possible genotype combinations for the given order
 * @param num_genotype_combinations Number of genotype combinations for the given order
 * @param num_affected Number of affected samples
 * @param num_unaffected Number of unaffected samples
 * @param num_counts Number of counts returned
 * @return List of counts, paired in (affected,unaffected)
 **/
int* get_counts(int order, uint8_t* genotypes, uint8_t **genotype_combinations, int num_genotype_combinations, int num_affected, int num_unaffected, int *num_counts);

uint8_t* get_masks(int order, uint8_t *genotypes, int num_samples, int *num_masks);


/* **************************
 *         High risk        *
 * **************************/

int* get_high_risk_combinations(unsigned int *counts, unsigned int num_counts, unsigned int num_affected, unsigned int num_unaffected, 
                                unsigned int *num_risky, void** aux_ret,
                                bool (*test_func)(unsigned int, unsigned int, unsigned int, unsigned int, void **));


risky_combination *risky_combination_new(int order, int comb[order], uint8_t **possible_genotypes_combinations, int num_risky, int *risky_idx, void *aux_info);

void risky_combination_free(risky_combination *combination);



/* **************************
 *  Evaluation and ranking  *
 * **************************/

unsigned int *get_confusion_matrix(int order, risky_combination *combination, int num_affected_in_fold, int num_unaffected_in_fold, uint8_t *genotypes);

double evaluate_model(unsigned int *confusion_matrix, enum eval_function function);

#endif