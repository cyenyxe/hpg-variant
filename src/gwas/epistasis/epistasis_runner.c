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


static void show_best_models_per_repetition(int order, int num_cv_repetitions, struct heap *best_models[]);
static khash_t(cvc) *prepare_models_for_cvc(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[]);
static int choose_best_model(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[], khash_t(cvc) *models_for_cvc, char **bestkey);


int run_epistasis(shared_options_data_t* shared_options_data, epistasis_options_data_t* options_data) {
    int ret_code = 0;
    
    // Load binary input dataset
    int num_affected, num_unaffected;
    size_t num_variants, file_len, genotypes_offset;
    
    uint8_t *input_file = epistasis_dataset_load(&num_affected, &num_unaffected, &num_variants, &file_len, &genotypes_offset, options_data->dataset_filename);
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
    
    LOG_INFO_F("Combinations of order %d, %d variants per block\n", order, stride);
    LOG_INFO_F("%d variants, %d blocks per dimension\n", num_variants, num_blocks_per_dim);
    
    // Precalculate which genotype combinations can be tested for a given order (order 2 -> {(0,0), (0,1), ... , (2,1), (2,2)})
    int num_genotype_permutations;
    uint8_t **genotype_permutations = get_genotype_combinations(order, &num_genotype_permutations);
    
    // Ranking of best models in each repetition
//    linked_list_t *best_models[options_data->num_cv_repetitions];
    struct heap *best_models[options_data->num_cv_repetitions];
    
    /**************************** End of variables precalculus  ****************************/
    
    
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        LOG_INFO_F("Running cross-validation #%d...\n", r+1);
        
        // Initialize folds, first block coordinates, genotype combinations and rankings for each repetition
        unsigned int *sizes, *training_sizes;
        int **folds = get_k_folds(num_affected, num_unaffected, num_folds, &sizes);
        uint8_t *fold_masks = get_k_folds_masks(num_affected, num_unaffected, num_folds, folds, sizes);
        
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
        for (int i = 0; i < num_folds; i++) {
            training_sizes[3 * i] = num_samples - sizes[3 * i];
            training_sizes[3 * i + 1] = num_affected - sizes[3 * i + 1];
            training_sizes[3 * i + 2] = num_unaffected - sizes[3 * i + 2];
        }
        
        // Initialize rankings for each repetition
        struct heap **ranking_risky = malloc(num_folds * sizeof(struct heap*));
        for (int i = 0; i < num_folds; i++) {
            ranking_risky[i] = malloc(sizeof(struct heap));
            heap_init(ranking_risky[i]);
        }
        
        // Coordinates of the block being tested
        int block_coords[order]; memset(block_coords, 0, order * sizeof(int));

#pragma omp parallel num_threads(shared_options_data->num_threads) shared(ranking_risky)
{
#pragma omp single
    {
        do {

            // OpenMP parallelization: Each block will be run in a separate thread
#pragma omp task
            {

            int tid = omp_get_thread_num();
            int my_block_coords[order];

#pragma omp critical
            {
                memcpy(my_block_coords, block_coords, order * sizeof(int));
            }
//            printf("%d) cv %d, block %d %d\n", omp_get_thread_num(), r, my_block_coords[0], my_block_coords[1]);

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

            // *************** Variables private to each task (block) (end) ***************



            // -------------------- Get genotypes of block --------------------

            uint8_t *block_starts[order];
//             printf("block starts = { ");
            for (int s = 0; s < order; s++) {
                block_starts[s] = genotypes + my_block_coords[s] * stride * num_samples;
//                 printf("%d ", block_coords[s] * stride);
            }
//             printf("}\n");

            // Initialize first coordinate (only if it's different from the previous)
            block_genotypes[0] = get_genotypes_for_block(num_variants, num_samples, info, stride,
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
                    block_genotypes[m] = get_genotypes_for_block(num_variants, num_samples, info, stride,
                            my_block_coords[m], block_starts[m], scratchpad[m]);
                }
            }


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
            int cur_comb_idx = -1;

            // Test first combination in the block
            get_first_combination_in_block(order, comb, my_block_coords, stride);

            do {

                // print_combination(comb, 0, order);

                cur_comb_idx = (cur_comb_idx < info.num_combinations_in_a_row - 1) ? cur_comb_idx+1 : 0;
                memcpy(combs + cur_comb_idx * order, comb, order * sizeof(int));

                if (cur_comb_idx < info.num_combinations_in_a_row - 1) {
                    continue; // Nothing to do until we have an amount (COMBINATIONS_ROW_SSE) of combinations ready
                }

                // Get genotypes of a row of combinations
                uint8_t *combination_genotypes[info.num_combinations_in_a_row * order];
                for (int c = 0; c < info.num_combinations_in_a_row; c++) {
                    for (int s = 0; s < order; s++) {
                        // Derive combination address from block
                        combination_genotypes[c * order + s] = block_genotypes[s] +
                                                               (combs[c * order + s] % stride) * info.num_samples_with_padding;
                    }
                }

                uint8_t *masks = set_genotypes_masks(order, combination_genotypes, info.num_combinations_in_a_row, info); // Grouped by SNP

                // Get counts for the provided genotypes
                combination_counts_all_folds(order, fold_masks, num_folds, genotype_permutations, info, counts_aff, counts_unaff);

                /* Right now the rest of the pipeline will be executed as it previously was, but for the sake of parallelization
                 * it could be better to make the choose_high_risk_combinations function work over
                 */
                for (int f = 0; f < num_folds; f++) {
                    // Get high risk pairs for those counts
                    void *aux_info;
                    unsigned int num_risky[info.num_combinations_in_a_row];
                    memset(num_risky, 0, info.num_combinations_in_a_row * sizeof(int));

                    int *risky_idx = choose_high_risk_combinations2(counts_aff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                                    counts_unaff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                                    info.num_combinations_in_a_row, info.num_cell_counts_per_combination,
                                                                    info.num_affected, info.num_unaffected,
                                                                    num_risky, &aux_info, mdr_high_risk_combinations2);

                    int risky_begin_idx = 0;
                    for (int rc = 0; rc < info.num_combinations_in_a_row; rc++) {
                        int *comb = combs + rc * order;
                        uint8_t **my_genotypes = combination_genotypes + rc * order;

                        risky_combination *risky_comb = NULL;

                        // Filter non-risky SNP combinations
                        if (num_risky > 0) {
                            // Put together the info about the SNP combination and its genotype combinations
                            risky_comb = risky_combination_new(order, comb, genotype_permutations,
                                                               num_risky[rc], risky_idx + risky_begin_idx,
                                                               aux_info);
                        }

                        risky_begin_idx += num_risky[rc];

                        if (risky_comb) {
                            // Check the model against the testing dataset
                            double accuracy = 0.0f;

    //                        if (options_data->evaluation_mode == TESTING) {
    //                             uint8_t *testing_genotypes = get_genotypes_for_combination_and_fold(order, risky_comb->combination,
    //                                                                                                 num_samples, sizes[3 * i + 1] + sizes[3 * i + 2],
    //                                                                                                 folds[i], stride, block_starts);
    //                             accuracy = test_model(order, risky_comb, testing_genotypes, sizes[3 * i + 1], sizes[3 * i + 2], &confusion_time);
    //                             free(testing_genotypes);
    //                        } else {
                                accuracy = test_model(order, risky_comb, my_genotypes, info, conf_matrix);
    //                        }
    //                         printf("*  Balanced accuracy: %.3f\n", accuracy);

                            int position = -1;

#pragma omp critical
                            {
                            position = add_to_model_ranking(risky_comb, options_data->max_ranking_size, ranking_risky[f], compare_risky_heap_min);
                            }

                            // If not inserted it means it is not among the most risky combinations, so free it
                            if (position < 0) {
                                risky_combination_free(risky_comb);
                            }
                        }

                    }

                    free(risky_idx);


//                     free(reference);
//
//                     for (int c = 0; c < num_samples; c++) {
//                         free(genotypes_for_testing[c]);
//                     }
//                     free(genotypes_for_testing);

                }

            } while (get_next_combination_in_block(order, comb, my_block_coords, stride, num_variants));

            // -------------------- Process combinations out of a full set --------------------

            // Get genotypes of a row of combinations
            uint8_t *combination_genotypes[info.num_combinations_in_a_row * order];
            for (int c = 0; c < cur_comb_idx + 1; c++) {
                for (int s = 0; s < order; s++) {
                    // Derive combination address from block
                    combination_genotypes[c * order + s] = block_genotypes[s] +
                                                            (combs[c * order + s] % stride) * info.num_samples_with_padding;
                }
            }

            uint8_t *masks = set_genotypes_masks(order, combination_genotypes, cur_comb_idx + 1, info); // Grouped by SNP

            // Get counts for the provided genotypes
            combination_counts_all_folds(order, fold_masks, num_folds, genotype_permutations, info, counts_aff, counts_unaff);

            /* Right now the rest of the pipeline will be executed as it previously was, but for the sake of parallelization
             * it could be better to make the choose_high_risk_combinations function work over
             */
            for (int f = 0; f < num_folds; f++) {
                // Get high risk pairs for those counts
                void *aux_info;
                unsigned int num_risky[info.num_combinations_in_a_row];
                memset(num_risky, 0, info.num_combinations_in_a_row * sizeof(int));

                int *risky_idx = choose_high_risk_combinations2(counts_aff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                                counts_unaff + f * info.num_combinations_in_a_row * info.num_cell_counts_per_combination,
                                                                info.num_combinations_in_a_row, info.num_cell_counts_per_combination,
                                                                info.num_affected, info.num_unaffected,
                                                                num_risky, &aux_info, mdr_high_risk_combinations2);

/*
                printf("num risky = { ");
                for (int rc = 0; rc < info.num_combinations_in_a_row; rc++) {
                    printf("%d ", num_risky[rc]);
                }
                printf("}\n");

                printf("risky gts = { ");
                for (int rc = 0; rc < info.num_combinations_in_a_row * info.num_counts_per_combination; rc++) {
                    printf("%d ", risky_idx[rc]);
                }
                printf("}\n");
*/

                int risky_begin_idx = 0;
                for (int rc = 0; rc < cur_comb_idx + 1; rc++) {
                    int *comb = combs + rc * order;
                    uint8_t **my_genotypes = combination_genotypes + rc * order;

                    // ------------------- BEGIN get_model_from_combination_in_fold -----------------------

                    risky_combination *risky_comb = NULL;

                    // Filter non-risky SNP combinations
                    if (num_risky > 0) {
                        // Put together the info about the SNP combination and its genotype combinations
                        risky_comb = risky_combination_new(order, comb, genotype_permutations,
                                                        num_risky[rc], risky_idx + risky_begin_idx,
                                                        aux_info);
                    }

                    risky_begin_idx += num_risky[rc];

                    // ------------------- END get_model_from_combination_in_fold -----------------------

                    if (risky_comb) {
                        // Check the model against the testing dataset
                        double accuracy = 0.0f;

//                        if (options_data->evaluation_mode == TESTING) {
//                             uint8_t *testing_genotypes = get_genotypes_for_combination_and_fold(order, risky_comb->combination,
//                                                                                                 num_samples, sizes[3 * i + 1] + sizes[3 * i + 2],
//                                                                                                 folds[i], stride, block_starts);
//                             accuracy = test_model(order, risky_comb, testing_genotypes, sizes[3 * i + 1], sizes[3 * i + 2], &confusion_time);
//                             free(testing_genotypes);
//                        } else {
                            accuracy = test_model(order, risky_comb, my_genotypes, info, conf_matrix);
//                        }
//                         printf("*  Balanced accuracy: %.3f\n", accuracy);

                        int position = -1;

#pragma omp critical
                        {
                        position = add_to_model_ranking(risky_comb, options_data->max_ranking_size, ranking_risky[f], compare_risky_heap_min);
                        }

                        // If not inserted it means it is not among the most risky combinations, so free it
                        if (position < 0) {
                            risky_combination_free(risky_comb);
                        }
                    }

                }

                free(risky_idx);

//                     free(reference);
//
//                     for (int c = 0; c < num_samples; c++) {
//                         free(genotypes_for_testing[c]);
//                     }
//                     free(genotypes_for_testing);

            }

            _mm_free(info.masks);
            for (int s = 0; s < order; s++) {
                _mm_free(scratchpad[s]);
            }
            _mm_free(counts_aff);
            _mm_free(counts_unaff);

            }
        } while (get_next_block(num_blocks_per_dim, order, block_coords));

    }
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
        
        // Merge all rankings in one
        size_t repetition_ranking_size = 0;
        for (int i = 0; i < num_folds; i++) {
            repetition_ranking_size += ranking_risky[i]->size;
        }
        risky_combination *repetition_ranking[repetition_ranking_size];
        size_t current_index = 0;
        
        for (int i = 0; i < num_folds; i++) {
            struct heap_node *hn;
            risky_combination *element = NULL;

//            printf("Ranking fold %d = {\n", i);
            while (!heap_empty(ranking_risky[i])) {
                hn = heap_take(compare_risky_heap_min, ranking_risky[i]);
                element = (risky_combination*) hn->value;
                repetition_ranking[current_index] = element;
                current_index++;
                free(hn);

//                printf("(%d ", element->combination[0]);
//                for (int s = 1; s < order; s++) {
//                    printf("%d ", element->combination[s]);
//                }
//                printf("- %.3f) ", element->accuracy);
            }
//            printf("}\n\n");
        }

        assert(current_index == repetition_ranking_size);
        
        // qsort by coordinates
        qsort(repetition_ranking, repetition_ranking_size, sizeof(risky_combination*), compare_risky);
        
/*
        printf("Ranking sorted = {\n");
        for (int i = 0; i < repetition_ranking_size; i++) {
            risky_combination *element = repetition_ranking[i];

            printf("(%2d ", element->combination[0]);
            for (int s = 1; s < order; s++) {
                printf("%2d ", element->combination[s]);
            }
            printf("- %.3f) ", element->accuracy);
        }
        printf("\n");
*/
        
        // Sum all values of each position and get the mean of accuracies
        struct heap *sorted_repetition_ranking = malloc(sizeof(struct heap)); heap_init(sorted_repetition_ranking);
        risky_combination *current = repetition_ranking[0];
        
        for (int i = 1; i < repetition_ranking_size; i++) {
            risky_combination *element = repetition_ranking[i];
/*
            printf("current = (%d %d - %.3f)\telement = (%d %d - %.3f)\n", 
                    current->combination[0], current->combination[1], current->accuracy, 
                    element->combination[0], element->combination[1], element->accuracy);
*/
            if (!compare_risky(&current, &element)) {
                assert(current != element);
                current->accuracy += element->accuracy;
                risky_combination_free(element);
            } else {
                current->accuracy /= num_folds;
                
/*
                printf("final acc = (%d %d - %.3f)\n", 
                        current->combination[0], current->combination[1], current->accuracy);
*/
                
                add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, compare_risky_heap_max);
                current = element;
            }
        }
        // Don't leave last element out!
        current->accuracy /= num_folds;
        add_to_model_ranking(current, repetition_ranking_size, sorted_repetition_ranking, compare_risky_heap_max);
        
        // Show full ranking
/*
        linked_list_iterator_t* iter = linked_list_iterator_new(sorted_repetition_ranking);
        risky_combination *element = NULL;
        printf("Final ranking = {\n");
        while(element = linked_list_iterator_curr(iter)) {
            printf("(%d %d - %.3f) ", element->combination[0], element->combination[1], element->accuracy);
            linked_list_iterator_next(iter);
        }
        printf("}\n");
        linked_list_iterator_free(iter);
*/

        // Save the models ranking
        best_models[r] = sorted_repetition_ranking;
        
        // Free data por this repetition
        for (int i = 0; i < num_folds; i++) {
            free(folds[i]);
            free(ranking_risky[i]);
        }
        free(ranking_risky);
        free(folds);
        free(sizes);
        free(training_sizes);
        _mm_free(fold_masks);
    }
    
    
    // Show the best model of each repetition
    show_best_models_per_repetition(order, options_data->num_cv_repetitions, best_models);
    
    // CVC (get the model that appears more times in the first ranking position)
    int max_val_len = log10f(num_variants);
    khash_t(cvc) *models_for_cvc = prepare_models_for_cvc(order, options_data->num_cv_repetitions, max_val_len, best_models);
    
    char *bestkey;
    int bestvalue = choose_best_model(order, options_data->num_cv_repetitions, max_val_len, best_models, models_for_cvc, &bestkey);
    
    assert(bestkey);
    LOG_INFO_F("Best model is %s with a CVC of %d/%d\n", bestkey, bestvalue, options_data->num_cv_repetitions);
    
    kh_destroy(cvc, models_for_cvc);

    // Free data for the whole epistasis check
    for (int i = 0; i < num_genotype_permutations; i++) {
        free(genotype_permutations[i]);
    }
    free(genotype_permutations);
    for (int r = 0; r < options_data->num_cv_repetitions; r++) {
        struct heap_node *hn;
        risky_combination *element = NULL;

        while (!heap_empty(best_models[r])) {
            hn = heap_take(compare_risky_heap_max, best_models[r]);
            risky_combination_free((risky_combination*) hn->value);
            free(hn);
        }
    }
    epistasis_dataset_close(input_file, file_len);
    
    return ret_code;
}


/* ******************************
 *      Auxiliary functions     *
 * ******************************/

int compare_risky(const void *risky_1, const void *risky_2) {
    risky_combination *r1 = *((risky_combination**) risky_1);
    risky_combination *r2 = *((risky_combination**) risky_2);
    
    for (int i = 0; i < r1->order; i++) {
        if (r1->combination[i] < r2->combination[i]) {
            return -1;
        }
        if (r1->combination[i] > r2->combination[i]) {
            return 1;
        }
    }
    
    return 0;
}

static void show_best_models_per_repetition(int order, int num_cv_repetitions, struct heap *best_models[]) {
    for (int r = 0; r < num_cv_repetitions; r++) {
        struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
        risky_combination *element = (risky_combination*) hn->value;
        assert(element);
        assert(element->combination);
        printf("CV %d\t(", r);
        for (int i = 0; i < order; i++) {
            printf(" %d ", element->combination[i]);
        }
        printf(") - %.3f)\n", element->accuracy);
    }
}

static khash_t(cvc) *prepare_models_for_cvc(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[]) {
    khash_t(cvc) *models_for_cvc = kh_init(cvc);
    for (int r = 0; r < num_cv_repetitions; r++) {
        struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
        risky_combination *risky = (risky_combination*) hn->value;
        
        // key = snp1_snp2_..._snpN
        char *key = calloc(order * (max_val_len + 1), sizeof(char));

        for (int i = 0; i < order-1; i++) {
            sprintf(key + strlen(key), "%d_", risky->combination[i]);
        }
        sprintf(key + strlen(key), "%d", risky->combination[order-1]);

        int ret;
        khiter_t iter = kh_get(cvc, models_for_cvc, key);
        if (iter != kh_end(models_for_cvc)) {
            (kh_value(models_for_cvc, iter))++; // Increment number of occurrences
        } else {
            iter = kh_put(cvc, models_for_cvc, key, &ret);
            if (ret) {
                kh_value(models_for_cvc, iter) = 1;
            }
        }
    }
    
    return models_for_cvc;
}

static int choose_best_model(int order, int num_cv_repetitions, int max_val_len, struct heap *best_models[], khash_t(cvc) *models_for_cvc, char **bestkey) {
    int bestvalue = 0;
    for (int k = kh_begin(models_for_cvc); k < kh_end(models_for_cvc); k++) {
        if (kh_exist(models_for_cvc, k)) {
            char *key = kh_key(models_for_cvc, k);
            int value = kh_value(models_for_cvc, k);
//             printf("%s -> %d\n", key, value);
            if (value > bestvalue) {
                *bestkey = key;
                bestvalue = value;
            } else if (value == bestvalue) {
                // If CVC(best) == CVC(candidate) ---> use CV-a
                double acc_best = 0.0f;
                double acc_candidate = 0.0f;
                
                // Sum all accuracies for the best and the candidate
                for (int r = 0; r < num_cv_repetitions; r++) {
                    struct heap_node *hn = heap_peek(compare_risky_heap_max, best_models[r]);
                    risky_combination *element = (risky_combination*) hn->value;
                    
                    // maybe_key = snp1_snp2_..._snpN
                    char *maybe_key = calloc(order * (max_val_len + 1), sizeof(char));
                    
                    for (int i = 0; i < order-1; i++) {
                        sprintf(maybe_key + strlen(maybe_key), "%d_", element->combination[i]);
                    }
                    sprintf(maybe_key + strlen(maybe_key), "%d", element->combination[order-1]);
                    
                    if (!strcmp(maybe_key, key)) {
                        acc_candidate += element->accuracy;
                    } else if (!strcmp(maybe_key, *bestkey)) {
                        acc_best += element->accuracy;
                    }
                }
                
                // Check which one is greater
                if (acc_candidate > acc_best) {
                    *bestkey = key;
                    bestvalue = value;
                }
            }
        }
    }
    
    return bestvalue;
}
