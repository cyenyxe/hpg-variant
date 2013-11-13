#ifndef MIC_ARRAY_UTILS_H
#define MIC_ARRAY_UTILS_H

#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "math_utils.h"


/**
 * @brief Shuffles an array of doubles
 * @details Shuffles an array of doubles using the Fisher–Yates shuffle algorithm, which guarantees that 
 * every permutation is equally likely (unbiased)
 * @see http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * 
 * @param[in,out] values The array to shuffle
 * @param n Number of elements in the array
 **/
__attribute__ (( target(mic) ))
void array_shuffle(double *values, size_t n);

/**
 * @brief Shuffles an array of integers
 * @details Shuffles an array of integers using the Fisher–Yates shuffle algorithm, which guarantees that 
 * every permutation is equally likely (unbiased)
 * @see http://en.wikipedia.org/wiki/Fisher%E2%80%93Yates_shuffle
 * 
 * @param[in,out] values The array to shuffle
 * @param n Number of elements in the array
 **/
__attribute__ (( target(mic) ))
void array_shuffle_int(int *values, size_t n);


#endif
