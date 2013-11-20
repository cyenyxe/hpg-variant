/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012-2013 Ignacio Medina (ICM-CIPF)
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

#ifndef EPISTASIS_H
#define EPISTASIS_H

/** 
 * @file epistasis.h
 * @brief Structures and functions associated to the options for the epistasis tool
 * 
 * This file defines the structures which store the options for the epistasis tool, and also 
 * functions to read their value from a configuration file or the command line.
 */ 

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef _USE_MPI
#include <mpi.h>
#endif

#include <commons/argtable/argtable2.h>
#include <commons/config/libconfig.h>
#include <commons/log.h>
#include <math_utils.h>

#include "error.h"
#include "shared_options.h"

#include "model.h"

/**
 * Number of options applicable to the epistasis tool.
 */
#define NUM_EPISTASIS_OPTIONS  7

KHASH_MAP_INIT_STR(cvc, int);


typedef struct epistasis_options {
    int num_options;
    
    struct arg_file *dataset_filename;    /**< Binary file used as input. */
    struct arg_int *order;
    struct arg_int *stride;
    struct arg_int *num_folds;
    struct arg_int *num_cv_repetitions;
    struct arg_int *max_ranking_size;
    struct arg_str *evaluation_subset;
    struct arg_str *evaluation_mode;
} epistasis_options_t;

/**
 * @brief Values for the options of the epistasis tool.
 * 
 * This struct contains the options specific to the epistasis tool, such as the 
 * order of the SNP combinations or the maximum number of cross-validation runs.
 */
typedef struct epistasis_options_data {
    char *dataset_filename;    /**< Binary file used as input. */
    int order;
    int stride;
    int num_folds;
    int num_cv_repetitions;
    int max_ranking_size;
    enum evaluation_subset eval_subset;
    enum evaluation_mode eval_mode;
} epistasis_options_data_t;


static epistasis_options_t *new_epistasis_cli_options(void);

/**
 * @brief Initializes an epistasis_options_data_t structure mandatory members.
 * @return A new epistasis_options_data_t structure.
 * 
 * Initializes a new epistasis_options_data_t structure mandatory members, which are the buffers for 
 * the URL parts, as well as its numerical ones.
 */
static epistasis_options_data_t *new_epistasis_options_data(epistasis_options_t *options);

/**
 * @brief Free memory associated to a epistasis_options_data_t structure.
 * @param options_data the structure to be freed
 * 
 * Free memory associated to a epistasis_options_data_t structure, including its text buffers.
 */
static void free_epistasis_options_data(epistasis_options_data_t *options_data);


/* **********************************************
 *                Options parsing               *
 * **********************************************/

/**
 * @brief Reads the configuration parameters of the epistasis tool.
 * @param filename file the options data are read from
 * @param options_data local options values (host URL, species, num-threads...)
 * @return Zero if the configuration has been successfully read, non-zero otherwise
 * 
 * Reads the basic configuration parameters of the tool. If the configuration
 * file can't be read, these parameters should be provided via the command-line
 * interface.
 */
int read_epistasis_configuration(const char *filename, epistasis_options_t *epistasis_options, shared_options_t *shared_options);

/**
 * @brief Parses the tool options from the command-line.
 * @param argc Number of arguments from the command-line
 * @param argv List of arguments from the command line
 * @param[out] options_data Struct where the tool-specific options are stored in
 * @param[out] global_options_data Struct where the application options are stored in
 * 
 * Reads the arguments from the command-line, checking they correspond to an option for the 
 * epistasis tool, and stores them in the local or global structure, depending on their scope.
 */
void **parse_epistasis_options(int argc, char *argv[], epistasis_options_t *epistasis_options, shared_options_t *shared_options);

void **merge_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options, struct arg_end *arg_end);

/**
 * @brief Checks semantic dependencies among the tool options.
 * @param global_options_data Application-wide options to check
 * @param options_data Tool-wide options to check
 * @return Zero (0) if the options are correct, non-zero otherwise
 * 
 * Checks that all dependencies among options are satisfied, i.e.: option A is mandatory, 
 * option B can't be provided at the same time as option C, and so on.
 */
int verify_epistasis_options(epistasis_options_t *epistasis_options, shared_options_t *shared_options);


/* **********************************************
 *               Epistasis execution            *
 * **********************************************/

__attribute__ (( target(mic) ))
void process_set_of_combinations(int num_combinations, int *combs, int order, int stride, 
                                 int num_folds, uint8_t *fold_masks, int *training_sizes, int *testing_sizes,
                                 uint8_t **block_genotypes, uint8_t *genotype_permutations,
                                 uint8_t *masks, enum evaluation_subset subset, masks_info info, 
                                 compare_risky_heap_func cmp_heap_func,
                                 int *counts_aff, int *counts_unaff, unsigned int conf_matrix[4], 
                                 int max_ranking_size, struct heap **ranking_risky_local);

__attribute__ (( target(mic) ))
struct heap* merge_rankings(int num_folds, struct heap **ranking_risky, compare_risky_heap_func heap_min_func, compare_risky_heap_func heap_max_func);

__attribute__ (( target(mic) ))
int compare_risky(const void *risky_1, const void *risky_2);


/* **********************************************
 *                  Final report                *
 * **********************************************/

__attribute__ (( target(mic) ))
void epistasis_report(int order, int cv_repetition, enum evaluation_mode mode, enum evaluation_subset subset, 
                      struct heap *best_models, int max_ranking_size, compare_risky_heap_func cmp_heap_max, FILE *fd);

#endif
