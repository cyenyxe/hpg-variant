/*
 * Copyright (c) 2012-2013 Cristina Yenyxe Gonzalez Garcia (ICM-CIPF)
 * Copyright (c) 2012 Ignacio Medina (ICM-CIPF)
 * Copyright (c) 2012 Andrei Alic
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

#ifndef FILE_HANDLING_H_
#define FILE_HANDLING_H_

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
#include <inttypes.h>

/** CUSTOM HEADERS */
#include <bioformats/vcf/vcf_file.h>
#include <bioformats/vcf/vcf_stats.h>
#include <bioformats/ped/ped_read.h>
#include <bioformats/ped/ped_file.h>
#include <bioformats/ped/ped_file_structure.h>
#include <containers/array_list.h>
#include <containers/list.h>
#include <commons/file_utils.h>
#include <commons/log.h>

#include <bioformats/vcf/vcf_util.h>
#include <bioformats/vcf/vcf_file_structure.h>

//#include "hpg_vcf_tools_utils.h"
#include "haplo.h"
#include "error.h"

/**********************************
* DECLARE ALL FUNCTION PROTOTYPES *
**********************************/

//bool get_markers_array_new(array_list_t *all_markers, const conf_params *cparams,
//		const user_params *uparams, unsigned int *num_samples);

bool get_markers_array(array_list_t *all_markers, const shared_options_data_t *share_data,
		const haplo_options_data_t *haplo_data, unsigned int  *num_samples);

/** @return char* As many array of chars as the number of elements in a block. The representation is
 * actually a matrix encoded as an array.
 * This is not really a real char array. We need its small dimension (8b) for storing the data associated with a sample
	 * First 4 bits are the first base, next 4 are the second base of a SNP */
marker **get_markers(array_list_t *variants, const unsigned int num_samples,
		const haplo_options_data_t *params);

/**
 * VCF files won't tell you the base for each allele in the SNP, instead they are giving a position
 * which if it is 0 then means it the ref allele otherwise is the position in the alternates
 * @param unsigned char allele_file The value found in the file for an individual
 * @param unsigned char* alternates Array with all the alternates as read from the VCF
 * @param unsigned char reference The reference allele, used in case of 0 for allele_file
 * @return unsigned char The code of the base as encoded in either A_BASE 1, T_BASE 2, G_BASE 3 or C_BASE 4
 */
static unsigned char allele_translation(int allele_file, unsigned char *alternates, unsigned char reference);

/**
 * Converts a char representation of a base (nucleic acid) to it's int correspondent to save
 * space for addressing
 * @param char allele The value as a char of the base
 * @return unsigned char One of the following values representing a base: A_BASE, T_BASE, G_BASE, C_BASE
 */
static unsigned char base_to_int(char allele);

/**
 * Calculate the rating of a marker (chromosome SNP for the input set of samples) using different params
 * and user set limits for each one of them. If the limit is broken then the rating is subtracted with a certain value
 * @return The overall accumulated result of the rating after it's params were compared with the user's ones
 */
inline static int calc_rating(double genopct, double pval, int menderr, double maf, const haplo_options_data_t *params);

/**
 * Check if a chromosome is X; It checks the input against a predefined list of values representing X
 * @param chromosome The char array containing the name of the Chromosome as read from the input file
 * @return True if the read chromosome is equal with some predefined values for X, false otherwise
 */
inline static bool is_x(const char *chromosome);


inline static double get_geno_percent(double called, double missing);

static double get_pvalue(const uint32_t *parent_hom, size_t parent_hom_len, uint32_t parent_het);

static double hw_calculate(uint32_t obsAA, uint32_t obsAB, uint32_t obsBB);

#endif /* FILE_HANDLING_H_ */
