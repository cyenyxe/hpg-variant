/*
 * globals.h
 *
 *  Created on: Oct 11, 2012
 *      Author: andrei
 */

#ifndef GLOBALS_H_
#define GLOBALS_H_

#include <string.h>
#include <stdbool.h>

#define A_BASE (unsigned char)1
#define T_BASE (unsigned char)4
#define G_BASE (unsigned char)3
#define C_BASE (unsigned char)2
#define UNK_BASE (unsigned char)0
/** use this when you need to distinguish between a real base and a special case which needs another starting code*/
#define NOT_A_BASE (unsigned char)0x5
/** Total number of representations available for the bases (known/unk/not_a_base)
 * in an input file; used to alloc memory when getting some statistics about the input, to avoid the need for a branch*/
#define NUM_KINDS_BASES_F 6

#define UNK_ALLELE -1

/** There are 2 options when it comes to determine the ref and alt allele:
 * first is to use the Haploview approach and to consider the ref the one having the largest number of appearances in the samples and the alt the other one
 * and the second is to actually use what the VCF file considers as being ref and alt*/
#define ALLELE_REF_HAPLOVIEW 1
#define ALLELE_REF_VCF 2

/** In VCF, the reference allele is always 0, but for a better understanding it is defined*/
#define REF_ALLELE_IDX 0

/** The idea of the mask is to obtain the desired number by annulling a certain part of the input;
 * Here we encoded 2 number , 4 bits each, in a 8 bit type, namely char, and we need to obtaib the
 * first 4 bits and then the next 4; we will use AND op to let the part we need unchanged and to
 * convert to 0 the other part. For reader's convenience the notations are in hex format*/
#define MASK_VAR1 0xF0 //11110000
#define MASK_VAR2 0x0F //00001111
/** The number of bits to be shifted to obtain the second encoded number */
#define NUM_BITS_SHIFT 4
/** Right now there are 4 bases in the DNA, when the science will discover more we will increase the number*/
#define NUM_BASES 4

/************************************************************************************
 * Global structures
 ************************************************************************************/
/** Configuration related to inner workings of the programs like
 * number of threads used to read the data from the disk */
typedef struct conf_params
{
	/** @var The full path including the name of the file */
	char *file_path;
	size_t max_simultaneous_batches;
	size_t num_threads;
	size_t batch_size;
	size_t entries_per_thread;
} conf_params;

/** The definition of a marker which is actually a line from VCF having
 * multiple samples from different individuals */
typedef struct marker {
	/** main allele*/
	unsigned char reference;
	int num_alternates;
	/** Array of alternates, now only one letter supported */
	unsigned char *alternates;
	/** The actual position in the genome */
	unsigned long int position;
	/** All the samples, representing individuals from the VCF columns
	 * It is not really a char array, but a method to store the first allele
	 * in the 1st 4 bits of the char and the second allele from the 2nd chr
	 * in the next 4 bits of the char */
	unsigned char *samples;
	/** minor allele frequency */
	double maf;
	/** rating for the marker calculated using mend_errors, maf, geno percentage and
	 *
	 */
	int rating;
	/** Special field stating if this marker is from X chromosome */
	bool is_x;
} marker;

/** To keep the program as expandable as possible, keep all the settings that can be changed by the user
 * in one place */
typedef struct user_params {
	/** method to use when determining the reference/alternate allele; Can be ALLELE_REF_HAPLOVIEW or ALLELE_REF_VCF*/
	int ref_allele_method;

	double hwCut;// = 0.001;
	int failedGenoCut;// = 75;
	int numMendErrCut;// = 1;
	double mafCut;// = 0.001;

} user_params;

/** Information about 2 markers used when creating the strong pairs */
typedef struct marker_info {
	/** position of the 1st marker in the strong array*/
	unsigned int marker_p1;
	/** position of the 2nd marker in the strong array*/
	unsigned int marker_p2;
	/** the distance between the 2 as extracted from VCF*/
	unsigned long sep;
} marker_info;

#endif /* GLOBALS_H_ */
