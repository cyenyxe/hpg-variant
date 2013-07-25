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
#include <commons/argtable/argtable2.h>

#include "global_options.h"

#define UNK_BASE (unsigned char)0
#define A_BASE (unsigned char)1
#define C_BASE (unsigned char)2
#define G_BASE (unsigned char)3
#define T_BASE (unsigned char)4
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

/** Multiple types of input options for the program; It is able to read from different
 types of data sources and thus they must be handled accordingly*/
typedef enum 
{
    /**Uses a VCF file to get the data set with the genes and all other options
     and it also loads a PED file to get the sex and the parents which aren't
     defined in the VCF standard*/
   IN_VCF_PED,
   /** A PED file together with a info file*/
   IN_PED_INFO        
} INPUT_TYPE;

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


#define NUM_HAPLO_OPTS 5
/** Options used when parsing arguments from the cmd line */
typedef struct haplo_options {
        int num_opts;
    
	struct arg_int *alleleref;
	struct arg_dbl *hwcut;// = 0.001;
	struct arg_int *fgcut;// = 75;
	struct arg_int *mendcut;// = 1;
	struct arg_dbl *mafcut;// = 0.001;

} haplo_options_t;

/** To keep the program as easy to expand as possible, keep all the settings that can be changed by the user
 * in one place */
typedef struct haplo_options_data {
	/** method to use when determining the reference/alternate allele; Can be ALLELE_REF_HAPLOVIEW or ALLELE_REF_VCF*/
	int alleleref;

	double hwcut;// = 0.001;
        /* Failed genotype cut*/
	int fgcut;// = 75;
        /* Number of mendelian errors for cutting */
	int mendcut;// = 1;
        /* Limit for the minor allele frequency when cutting */
	double mafcut;// = 0.001;

} haplo_options_data_t;

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
