/*
 * file_handling.c
 *
 *  Created on: Oct 11, 2012
 *      Author: andrei
 */

#include "file_handling.h"

bool get_markers_array(array_list_t *all_markers, const conf_params *cparams,
		const user_params *uparams, unsigned int  *num_samples)
{
	int ret_code;

    vcf_file_t* vcf_file_p = vcf_open(cparams->file_path, cparams->max_simultaneous_batches);
    array_list_t *samples_names = vcf_file_p->samples_names;
//#pragma omp parallel sections //private(start, stop, total)
//{
//    #pragma omp section
//    {
        LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
        // Reading
        //start = omp_get_wtime();

        ret_code = vcf_parse_batches(cparams->batch_size, vcf_file_p, 1);
        notify_end_reading(vcf_file_p);

       // stop = omp_get_wtime();
        //total = (stop - start);

        if (ret_code) { LOG_FATAL_F("[%dR] Error code = %d\n", omp_get_thread_num(), ret_code); }
        //LOG_INFO_F("[%dR] Time elapsed = %f s\n", omp_get_thread_num(), total);
        //LOG_INFO_F("[%dR] Time elapsed = %e ms\n", omp_get_thread_num(), total*1000);
//    }

//    #pragma omp section
//    {

        //start = omp_get_wtime();

        int header_written = 0;

        int i = 0;
        vcf_batch_t *batch;
        while ( (batch = fetch_vcf_batch(vcf_file_p)) != NULL ) {
        	marker **segment_markers = get_markers(batch->records, samples_names->size, uparams);
        	// Heck knows why it's impossible to get array's size directly with the indirection operator. level 2
        	array_list_t *temp = batch->records;
        	size_t num_rec = temp->size;

        	array_list_insert_all((void **) segment_markers, num_rec, all_markers);

            vcf_batch_free(batch);

            i++;
        }

//    }
//}
    *num_samples = samples_names->size;

    vcf_close(vcf_file_p);
    return true;
}

static marker **get_markers(array_list_t *variants, const int num_samples,
		const user_params *params) {
	//int num_samples;
    char *copy_buf, *copy_buf2, *token, *sample;
    char *save_strtok, *str_dup_buf;
    // allele 1 & 2 definition to be packed further in a unsigned char
    unsigned char a1, a2, auxuc;

    int num_alternates, gt_pos, cur_pos, numa1, numa2;
    int allele1, allele2, alleles_code;

    // Temporary variables for file stats updating
    int variants_count = 0, samples_count = 0, snps_count = 0, indels_count = 0, pass_count = 0;
    int transitions_count = 0, transversions_count = 0, biallelics_count = 0, multiallelics_count = 0;
    double accum_quality = 0, minor_allele_num, maf, auxf, pvalue;
    // Counter for the bases in a file
    int64_t alleles_count[NUM_KINDS_BASES_F], sum=0, num=0, mincount = -1,
    		total_alleles_count = 0, total_genotypes_count = 0, num_hets, num_alleles, female_count = 0;
    uint32_t founderHetCount = 0, *founderHomCount;
    double missing_data, called;

    // Variant stats management
    vcf_record_t *record;
    size_t num_variants = variants->size;
    //num_samples = variants[0]->samples->size;
    marker **result = malloc(num_variants * sizeof(marker*));
    bool is_chr_x;
    founderHomCount = malloc(sizeof(*founderHomCount) * NUM_KINDS_BASES_F);
    for (size_t i = 0; i < num_variants; i++) {
    	//init the alleles' counter, better than a for because the compiler can decide how to expand it
    	memset(alleles_count, 0, sizeof(*alleles_count)*NUM_KINDS_BASES_F);

        record = (vcf_record_t *) array_list_get(i, variants);
        result[i] = malloc(sizeof(**result));
        // Reset counters
        total_alleles_count = 0;
        total_genotypes_count = 0;
        numa1 = 0;
        numa2 = 0;

        // Create list of alternates
        copy_buf = (char*) calloc (strlen(record->alternate)+1, sizeof(char));
        strcat(copy_buf, record->alternate);
        //for now use this, even if it is possible to have more than 1 base as an alternate
        char **alt = split(copy_buf, ",", &num_alternates);
        (result[i])->alternates =  malloc(sizeof(*((result[i])->alternates))*num_alternates);
        for (uint32_t idxAlt=0;idxAlt<num_alternates; idxAlt++) {
        	(result[i])->alternates[idxAlt] = base_to_int(alt[idxAlt][0]);
        	//printf("alt %u\n", (result[i])->alternates[idxAlt]);
        }
        (result[i])->reference = base_to_int(record->reference[0]);//strtol(record->reference, NULL, 10);

        (result[i])->num_alternates = num_alternates;
//         if (!strncmp(stats->alternates[0], ".", 1)) {
//             stats->num_alleles = 1;
//         } else {
//             stats->num_alleles = num_alternates + 1;
//         }

        // Create lists of allele and genotypes counters and frequencies
        /*stats->alleles_count = (int*) calloc (stats->num_alleles, sizeof(int));
        stats->genotypes_count = (int*) calloc (stats->num_alleles * stats->num_alleles, sizeof(int));
        stats->alleles_freq = (double*) calloc (stats->num_alleles, sizeof(double));
        stats->genotypes_freq = (double*) calloc (stats->num_alleles * stats->num_alleles, sizeof(double));*/

        // Get position where GT is in sample
        str_dup_buf = (char *) strdup(record->format);
        gt_pos = get_field_position_in_format("GT", str_dup_buf);
        LOG_DEBUG_F("Genotype position = %d\n", gt_pos);
        if (gt_pos < 0) { continue; }   // This variant has no GT field
        (result[i])->samples = malloc(num_samples * sizeof(*((result[i])->samples)));
        (result[i])->position = record->position;
        minor_allele_num = 0;
        called = 0.0;
        missing_data = 0.0;
        female_count = 0;
        // check if the current chromosome record is the X one ( we have some special treatment for it)
        is_chr_x = is_x(record->chromosome);
        char id[70];
         strncpy(id, record->id, 70);
        printf("%s\n", id);
        if (strncmp(id, "rs142725888", 11) == 0)
        	printf("stop");
        (result[i])->is_x = is_chr_x;
        founderHetCount = 0;
        memset(founderHomCount, 0, sizeof(*founderHomCount) * NUM_KINDS_BASES_F);

        // Traverse samples and find the existing and missing alleles
        for(uint32_t j = 0; j < num_samples; j++) {
        	sample = (char *) array_list_get(j, record->samples);
            // Get to GT position
            alleles_code = get_alleles(strdup(sample), gt_pos, &allele1, &allele2);

            // We need to determine how many times the reference and the alternates appear in the file
            // We consider for now that ONLY ONE alt EXISTS
            a1 == REF_ALLELE_IDX ? numa1++ : numa2++;
            a2 == REF_ALLELE_IDX ? numa1++ : numa2++;

            // Translate allele to a base to be easier to handle later; we don't need the idx, we need the
            // base representation
            a1 = allele_translation(allele1, (result[i])->alternates,
            		(result[i])->reference);
            a2 = allele_translation(allele2, (result[i])->alternates,
                        		(result[i])->reference);
            // before we actually prepare the data for D' and any further execution, we need to determine
            // some statistics from the original file to be used when determining the MAF for instance
            // DISABLED BECAUSE OF THE HAPLOVIEW CODE
//            alleles_count[a1]++;
//            alleles_count[a2]++;
            if (a1 >= 6 || a2 >= 6)
            	printf("cute");
            //CODE FROM HAPLOVIEW
                    	if(a1 > 0 && a2 >0){

            					//if (a1 != 9){  //value of 9 means an 'h' allele for haps files...
            						alleles_count[a1]++;
            //					}else{
            //						alleles_count[5]++;
            //					}
            					// TODO NORMAL - There is no gender column in the VCF format; please advise
            					if (!is_chr_x) {// || currentInd.getGender() != 1) {
            						if(a1 != a2) {// || a1 == 9 || a2 == 9) {
            							founderHetCount++;
            						}else{
            							founderHomCount[a1]++;
            						}
            						//if(a2 != 9){
            							alleles_count[a2]++;
//            						}else{
//            							alleles_count[5]++;
//            						}
            						// TODO  There is no sex in the VCF
            						female_count++;
            					}


            				called++;
            			}
            			//missing data
            			else missing_data++;

            if (a1 == a2 || a1 == 0 || a2 == 0) {
            				/*a1 = NOT_A_BASE + a1;
            				a2 = NOT_A_BASE + a2;*/
            			} else {
			//if (a1 != a2 && a1 >= 0 && a2 >= 0) {
				/*a1 = NOT_A_BASE + a1;
				a2 = NOT_A_BASE + a2;*/
				a1 = 4 + a1;
				a2 = 4 + a2;
			}

            // create the magic variable with the 2 alleles (first one on the first 4 bits, second one
            // on the next 4 bits)
        	(result[i])->samples[j] = (a1 << NUM_BITS_SHIFT) + a2 ;
        }

        // Check the type of ref allele determination required and switch if necessary
		if (params->ref_allele_method == ALLELE_REF_HAPLOVIEW
				&& numa2 >= numa1) {//we just exchange the first position
			auxuc = (result[i])->reference;
			(result[i])->reference = (result[i])->alternates[0];
			(result[i])->alternates[0] = auxuc;
		}// Else don't change anything

        // CODE FROM HAPLOVIEW
        num_hets = alleles_count[5];
        alleles_count[5] = 0;
        if (num_hets > 0){
			num_alleles = 0;
			for (int i = 1; i < NUM_KINDS_BASES_F-1; i++){
				if (alleles_count[i] > 0){
					num_alleles++;
				}
			}

			if (num_alleles == 0){
				alleles_count[1] += num_hets/2;
				alleles_count[3] += num_hets/2;
			}else if (num_alleles ==  1){
				for (int i = 1; i < NUM_KINDS_BASES_F-1; i++){
					if (alleles_count[i] > 0){
						alleles_count[i] += num_hets/2;
						if (i == 4){
							alleles_count[3] += num_hets/2;
						}else{
							alleles_count[i+1] += num_hets/2;
						}
						break;
					}
				}
			}else if (num_alleles == 2){
				for (int i = 1; i < NUM_KINDS_BASES_F -1; i++){
					if (alleles_count[i] > 0){
						alleles_count[i] += num_hets/2;
					}
				}
			}
		}
        maf = 0;
        //sumsq=0;
        sum=0; num=0; mincount = -1;
		//int num_alleles = 0;
        //zero the missing allele
        alleles_count[0] = 0;
		for(int i=0;i<NUM_KINDS_BASES_F;i++){
			if(alleles_count[i] != 0){
				//numberOfAlleles++;
				num = alleles_count[i];
				//sumsq += num*num;
				sum += num;
				if (mincount < 0 || mincount > num){
					mincount = num;
				}
			}
		}

		/*if (numberOfAlleles > 2){
			throw new PedFileException("More than two alleles!");
		}*/

		if (sum == 0){
			maf = 0.0f;
		}else{
			if ((auxf = mincount/(double)sum) == 1.0f){
				maf = 0.0f;//numa2/(numa1+numa2);
			}else{
				maf = auxf;
			}
		}

		(result[i])->maf = maf;

		// TODO NORMAL - female count can't be done in VCF; Undefined behavior
		//This will cause the values to show up as NA since there aren't enough females to calculate
//		if(female_count < 10 && is_chr_x){
//			pvalue = DBL_MAX;
//		}

		// TODO MINOR - There are no families thus the mendelian errors can be calculated correctly
		//this param is calculated when there are families (in a ped file for instance)
		int mend_error_num = 0;
		double genopct = get_geno_percent(called, missing_data);
		pvalue = get_pvalue(founderHomCount, NUM_KINDS_BASES_F, founderHetCount);
		(result[i])->rating = calc_rating(genopct, pvalue, mend_error_num, maf, params);

        free(copy_buf);
    }
    free(founderHomCount);
    return result;
}/*get_markers*/

inline static bool is_x(const char *chromosome)
{
	if (strcmp(chromosome, "X") == 0 || strcmp(chromosome, "x") == 0)
		return true;
	else
		return false;
}

inline static double get_geno_percent(double called,
		double missing)
{
	if (called == 0){
		return 0;
	} else{
		return 100.0*(called/(called+missing));
	}
} /*get_geno_percent*/

inline static int calc_rating(double genopct, double pval, int menderr, double maf, const user_params *params)
{
	int rating = 0;
	if (genopct < params->failedGenoCut){
		rating -= 2;
	}
	if (pval < params->hwCut){
		rating -= 4;
	}
	if (menderr > params->numMendErrCut){
		rating -= 8;
	}
	if (maf < params->mafCut){
		rating -= 16;
	}
	if (rating == 0){
		rating = 1;
	}
	return rating;
} /*getRating*/

static double get_pvalue(const uint32_t *parent_hom, size_t parent_hom_len, uint32_t parent_het)
{
	//ie: 11 13 31 33 -> homA =1 homB = 1 parentHet=2
	uint32_t homA=0, homB=0;
	double pvalue=0;
	for(size_t i=0;i<parent_hom_len;i++){
		if(parent_hom[i] !=0){
			if(homA>0) homB = parent_hom[i];
			else homA = parent_hom[i];
		}
	}
	//caculate p value from homA, parentHet and homB
	if (homA + parent_het + homB <= 0){
		pvalue=0;
	}else{
		pvalue = hwCalculate(homA, parent_het, homB);
	}
	return pvalue;
} /*get_pvalue*/

 static double hwCalculate(uint32_t obsAA, uint32_t obsAB, uint32_t obsBB)
 {
        //Calculates exact two-sided hardy-weinberg p-value. Parameters
        //are number of genotypes, number of rare alleles observed and
        //number of heterozygotes observed.
        //
        // (c) 2003 Jan Wigginton, Goncalo Abecasis
	 double result = 0.0;
	  uint64_t diplotypes =  obsAA + obsAB + obsBB;
	  uint64_t aux = (obsAA*2) + obsAB;
	 //make sure "rare" allele is really the rare allele
	 int64_t rare =  aux > diplotypes ? 2*diplotypes-aux : aux;
	  int32_t hets = obsAB;


//
//        if (rare > diplotypes){
//            rare = 2*diplotypes-rare;
//        }
        //make sure numbers aren't screwy
        if (hets > rare){
            throw_error("Exception: HW test: %" PRIu32"heterozygotes but only %" PRIu64"rare alleles.", hets, rare);
        }

        size_t tail_probs_len = rare + 1;
        double *tail_probs = malloc(sizeof(*tail_probs) * tail_probs_len);
        for (int z = 0; z < tail_probs_len; z++){
        	tail_probs[z] = 0;
        }

        //start at midpoint
        //all the casting is to make sure we don't overflow ints if there are 10's of 1000's of inds
        int mid = (int)((double)rare * (double)(2 * diplotypes - rare) / (double)(2 * diplotypes));

        //check to ensure that midpoint and rare alleles have same parity
        if (((rare & 1) ^ (mid & 1)) != 0){
            mid++;
        }
        int het = mid;
        int hom_r = (rare - mid) / 2;
        int hom_c = diplotypes - het - hom_r;

        //Calculate probability for each possible observed heterozygote
        //count up to a scaling constant, to avoid underflow and overflow
        tail_probs[mid] = 1.0;
        double sum = tail_probs[mid];
        for (het = mid; het > 1; het -=2){
        	tail_probs[het-2] = (tail_probs[het] * het * (het-1.0))/(4.0*(hom_r + 1.0) * (hom_c + 1.0));
            sum += tail_probs[het-2];
            //2 fewer hets for next iteration -> add one rare and one common homozygote
            hom_r++;
            hom_c++;
        }

        het = mid;
        hom_r = (rare - mid) / 2;
        hom_c = diplotypes - het - hom_r;
        for (het = mid; het <= rare - 2; het += 2){
        	tail_probs[het+2] = (tail_probs[het] * 4.0 * hom_r * hom_c) / ((het+2.0)*(het+1.0));
            sum += tail_probs[het+2];
            //2 more hets for next iteration -> subtract one rare and one common homozygote
            hom_r--;
            hom_c--;
        }

        for (int z = 0; z < tail_probs_len; z++){
        	tail_probs[z] /= sum;
        }

        double top = tail_probs[hets];
        for (int i = hets+1; i <= rare; i++){
            top += tail_probs[i];
        }
        double otherSide = tail_probs[hets];
        for (int i = hets-1; i >= 0; i--){
            otherSide += tail_probs[i];
        }

        if (top > 0.5 && otherSide > 0.5){
        	result = 1.0;
        }else{
            if (top < otherSide){
            	result =  top * 2;
            }else{
            	result =  otherSide * 2;
            }
        }
        free(tail_probs);
        return result;
    }

static inline unsigned char base_to_int(char allele) {
	allele = toupper(allele);
	switch (allele) {
		case 'A':
			return A_BASE;
		case 'C':
			return C_BASE;
		case 'G':
			return G_BASE;
		case 'T':
			return T_BASE;
		case '.':
			return UNK_BASE;
		default://we can't detect the representation, send an error message
			return NOT_A_BASE;
	}
}/* base_to_int*/

static inline unsigned char allele_translation(int allele_file, unsigned char *alternates, unsigned char reference) {
	//printf("%d\n", allele_file);
	if (allele_file == 0) {//if the pos is 0 then it is the ref
		return reference;
	} else if (allele_file == UNK_ALLELE) {//if the base is unknown just set as UNK VALUE
		return UNK_BASE;
	}else {//else search in the alt array and get the code
		//char representation = alternates[allele_file-1][0];
		return //base_to_int(alternates[allele_file-1][0]);
				alternates[allele_file-1];
	}
}/*allele_translation*/

/**************************************************************
 * OLD CODE
 * ***********************************************************/

//static marker **get_markers_old(vcf_record_t **variants, const int num_variants, const int num_samples,
//		const user_params *params) {
//	//int num_samples;
//    char *copy_buf, *copy_buf2, *token, *sample;
//    char *save_strtok, *str_dup_buf;
//    // allele 1 & 2 definition to be packed further in a unsigned char
//    unsigned char a1, a2, auxuc;
//
//    int num_alternates, gt_pos, cur_pos, numa1, numa2;
//    int allele1, allele2, alleles_code;
//
//    // Temporary variables for file stats updating
//    int variants_count = 0, samples_count = 0, snps_count = 0, indels_count = 0, pass_count = 0;
//    int transitions_count = 0, transversions_count = 0, biallelics_count = 0, multiallelics_count = 0;
//    double accum_quality = 0, minor_allele_num, maf, auxf;
//    // Temporary variables for variant stats updating
//    int total_alleles_count = 0, total_genotypes_count = 0, num_hets, num_alleles;
//    // Counter for the bases in a file
//    long alleles_count[NUM_KINDS_BASES_F], sum=0, num=0, mincount = -1;
//
//    // Variant stats management
//    vcf_record_t *record;
//    //num_samples = variants[0]->samples->size;
//    marker **result = malloc(num_variants * sizeof(marker*));
//    for (int i = 0; i < num_variants; i++) {
//    	//init the alleles' counter, better than a for because the compiler can decide how to expand it
//    	memset(alleles_count, 0, sizeof(*alleles_count)*NUM_KINDS_BASES_F);
//
//        record = variants[i];
//        result[i] = malloc(sizeof(**result));
//        // Reset counters
//        total_alleles_count = 0;
//        total_genotypes_count = 0;
//        numa1 = 0;
//        numa2 = 0;
//
//        // Create list of alternates
//        copy_buf = (char*) calloc (strlen(record->alternate)+1, sizeof(char));
//        strcat(copy_buf, record->alternate);
//        //for now use this, even if it is possible to have more than 1 base as an alternate
//        char **alt = split(copy_buf, ",", &num_alternates);
//        (result[i])->alternates =  malloc(sizeof(*((result[i])->alternates))*num_alternates);
//        for (int idxAlt=0;idxAlt<num_alternates; idxAlt++) {
//        	(result[i])->alternates[idxAlt] = base_to_int(alt[idxAlt][0]);
//        	//printf("alt %u\n", (result[i])->alternates[idxAlt]);
//        }
//        (result[i])->reference = base_to_int(record->reference[0]);//strtol(record->reference, NULL, 10);
//
//        (result[i])->num_alternates = num_alternates;
////         if (!strncmp(stats->alternates[0], ".", 1)) {
////             stats->num_alleles = 1;
////         } else {
////             stats->num_alleles = num_alternates + 1;
////         }
//        //LOG_DEBUG_F("num alternates = %d\tnum_alleles = %d\n", num_alternates, stats->num_alleles);
//
//        // Create lists of allele and genotypes counters and frequencies
//        /*stats->alleles_count = (int*) calloc (stats->num_alleles, sizeof(int));
//        stats->genotypes_count = (int*) calloc (stats->num_alleles * stats->num_alleles, sizeof(int));
//        stats->alleles_freq = (double*) calloc (stats->num_alleles, sizeof(double));
//        stats->genotypes_freq = (double*) calloc (stats->num_alleles * stats->num_alleles, sizeof(double));*/
//
//        // Get position where GT is in sample
//        str_dup_buf = strdup(record->format);
//        gt_pos = get_field_position_in_format("GT", str_dup_buf);
//        LOG_DEBUG_F("Genotype position = %d\n", gt_pos);
//        if (gt_pos < 0) { continue; }   // This variant has no GT field
//        (result[i])->samples = malloc(num_samples * sizeof(*((result[i])->samples)));
//        (result[i])->position = record->position;
//        minor_allele_num = 0;
//        // Traverse samples and find the existing and missing alleles
//        for(int j = 0; j < num_samples; j++) {
//        	sample = (char *) array_list_get(j, record->samples);
//            // Get to GT position
//            alleles_code = get_alleles(strdup(sample), gt_pos, &allele1, &allele2);
//            LOG_DEBUG_F("sample = %s, alleles = %d/%d\n", sample, allele1, allele2);
//            // create the magic variable with the 2 alleles (first one on the first 4 bits, second one
//            // on the next 4 bits)
//
//            /*if (allele1 < 0 || allele2 < 0){
//                                    thisMarkerA = 0;
//                                    thisMarkerB = 0;
//                                }else{
//                                    thisMarkerA = currentInd.getAllele(i,0);
//                                    thisMarkerB = currentInd.getAllele(i,1);
//                                }*/
//            // We need to determine how many times the reference and the alternates appear in the file
//            // We consider for now that ONLY ONE alt EXISTS
//            a1 == REF_ALLELE_IDX ? numa1++ : numa2++;
//            a2 == REF_ALLELE_IDX ? numa1++ : numa2++;
//
//            // Translate allele to a base to be easier to handle later; we don't need the idx, we need the
//            // base representation
//            a1 = allele_translation(allele1, (result[i])->alternates,
//            		(result[i])->reference);
//            a2 = allele_translation(allele2, (result[i])->alternates,
//                        		(result[i])->reference);
//            // before we actually prepare the data for D' and any further execution, we need to determine
//            // some statistics from the original file to be used when determining the MAF for instance
//            alleles_count[a1]++;
//            alleles_count[a2]++;
//
//            if (a1 == a2 || a1 == 0 || a2 == 0) {
//            				/*a1 = NOT_A_BASE + a1;
//            				a2 = NOT_A_BASE + a2;*/
//            			} else {
//			//if (a1 != a2 && a1 >= 0 && a2 >= 0) {
//				/*a1 = NOT_A_BASE + a1;
//				a2 = NOT_A_BASE + a2;*/
//				a1 = 4 + a1;
//				a2 = 4 + a2;
//			}
//
//        	(result[i])->samples[j] = (a1 << NUM_BITS_SHIFT) + a2 ;
//
//        }
//
//        // Check the type of ref allele determination required and switch if necessary
//		if (params->ref_allele_method == ALLELE_REF_HAPLOVIEW
//				&& numa2 >= numa1) {//we just exchange the first position
//			auxuc = (result[i])->reference;
//			(result[i])->reference = (result[i])->alternates[0];
//			(result[i])->alternates[0] = auxuc;
//		}// Else don't change anything
//
//        // CODE FROM HAPLOVIEW
//        num_hets = alleles_count[5];
//        alleles_count[5] = 0;
//        if (num_hets > 0){
//			num_alleles = 0;
//			for (int i = 1; i < NUM_KINDS_BASES_F-1; i++){
//				if (alleles_count[i] > 0){
//					num_alleles++;
//				}
//			}
//
//			if (num_alleles == 0){
//				alleles_count[1] += num_hets/2;
//				alleles_count[3] += num_hets/2;
//			}else if (num_alleles ==  1){
//				for (int i = 1; i < NUM_KINDS_BASES_F-1; i++){
//					if (alleles_count[i] > 0){
//						alleles_count[i] += num_hets/2;
//						if (i == 4){
//							alleles_count[3] += num_hets/2;
//						}else{
//							alleles_count[i+1] += num_hets/2;
//						}
//						break;
//					}
//				}
//			}else if (num_alleles == 2){
//				for (int i = 1; i < NUM_KINDS_BASES_F -1; i++){
//					if (alleles_count[i] > 0){
//						alleles_count[i] += num_hets/2;
//					}
//				}
//			}
//		}
//        maf = 0;
//        //sumsq=0;
//        sum=0; num=0; mincount = -1;
//		//int num_alleles = 0;
//        //zero the missing allele
//        alleles_count[0] = 0;
//		for(int i=0;i<NUM_KINDS_BASES_F;i++){
//			if(alleles_count[i] != 0){
//				//numberOfAlleles++;
//				num = alleles_count[i];
//				//sumsq += num*num;
//				sum += num;
//				if (mincount < 0 || mincount > num){
//					mincount = num;
//				}
//			}
//		}
//
//		/*if (numberOfAlleles > 2){
//			throw new PedFileException("More than two alleles!");
//		}*/
//
//		if (sum == 0){
//			maf = 0.0f;
//		}else{
//			if ((auxf = mincount/(double)sum) == 1.0f){
//				maf = 0.0f;//numa2/(numa1+numa2);
//			}else{
//				maf = auxf;
//			}
//		}
//
//		(result[i])->maf = maf;
//        free(copy_buf);
//    }
//    return result;
//}/*get_markers*/


// bool get_markers_array(array_list_t *all_markers, const conf_params *cparams,
// 		const user_params *uparams, unsigned int  *num_samples)
// {
// 	vcf_file_t *vcf_file_p;
// 	int ret_code;
// 	list_t *read_list, *output_list;
// 	list_item_t *item;
//     marker **segment_markers;
//     file_stats_t *file_stats;
//     int *chunk_sizes, *chunk_starts, num_chunks = 0;
//
//     //lists inits
// 	read_list = (list_t *) malloc(sizeof(list_t));
// 	list_init("batches", 1, INT_MAX, read_list);
// 	output_list = (list_t*) malloc (sizeof(list_t));
// 	list_init("output", cparams->num_threads, INT_MAX, output_list);
// 	file_stats = file_stats_new();
//
//
// 	/*for(count = 0; count < argc; count++ )
// 		        printf( "  argv[%d]   %s\n", count, argv[count] );*/
// 		//printf("%s\n", getcwd(directory, sizeof(directory)));
//
// 		vcf_file_p = vcf_open(cparams->file_path, cparams->max_simultaneous_batches);
//
// 		// READ THE FILE AND INIT THE STATS
// 		// taken from bioinfo-c/libs/bioformats/vcf/vcf_stats.c & bioinfo-c/hpg-vcf-tools/stats/stats_runner.c
// 	/*#pragma omp parallel sections
// 	    {
// 	#pragma omp section
// 	        {*/
//
// 		ret_code = vcf_parse_batches(1000000, vcf_file_p, 1);//(read_list, 1000, vcf_file_p, 100);
//
// 		//unsuccessfully opened/parsed the file
// 		if (ret_code)
// 			return false;
// 		list_decr_writers(read_list);
// 	 /*       }//end read batch #pragma omp section*/
//
// 		array_list_t *samples_names = vcf_file_p->samples_names;
// 				//(marker *) malloc(sizeof(*all_markers) * vcf_file_p->num_records);
//
// 	/*#pragma omp section
// 	        {
// 		// Enable nested parallelism and set the number of threads the user has chosen
// 		omp_set_nested(1);
// 		omp_set_num_threads(num_threads);*/
//
// 				size_t i = 0;
// 				//list_item_t* item = NULL;
// 				item = NULL;
// 				while ((item = list_remove_item(read_list)) != NULL) {
// 					vcf_batch_t *batch = (vcf_batch_t*) item->data_p;
// 					array_list_t *input_records = batch;
//
// 					/*if (i % 50 == 0) {
// 						LOG_INFO_F("Batch %d reached by thread %d - %zu/%zu records \n",
// 									i, omp_get_thread_num(),
// 									batch->length, batch->max_length);
// 					}*/
//
// 					// Divide the list of passed records in ranges of size defined in config file
// 					//int max_chunk_size = entries_per_thread;
//
//
// 					/*list_item_t **chunk_starts = create_chunks(input_records, max_chunk_size, &num_chunks);
//
// 					// OpenMP: Launch a thread for each range
// 					//#pragma omp parallel for
// 					for (int j = 0; j < num_chunks; j++) {
// 						//LOG_DEBUG_F("[%d] Stats invocation\n", omp_get_thread_num());
// 						ret_code = get_variants_stats(chunk_starts[j], max_chunk_size, output_list, file_stats);
// 					}*/
//
// 					chunk_starts = create_chunks(input_records->size, cparams->entries_per_thread, &num_chunks, &chunk_sizes);
// 					//                 int max_chunk_size = shared_options_data->entries_per_thread;
// 					//                 int num_chunks;
// 					//					list_item_t **chunk_starts = create_chunks(input_records, max_chunk_size, &num_chunks);
//
// 					                // OpenMP: Launch a thread for each range
// 					//#pragma omp parallel for
// 					for (int j = 0; j < num_chunks; j++) {
// 						//LOG_INFO_F("[%d] Stats invocation\n", omp_get_thread_num());
// 	//                     ret_code = get_variants_stats(chunk_starts[j], max_chunk_size, output_list, file_stats);
// 						/*ret_code = get_variants_stats((vcf_record_t**) (input_records->items + chunk_starts[j]),
// 													  chunk_sizes[j],
// 													  output_list,
// 													  file_stats);*/
// 						segment_markers = get_markers((vcf_record_t**) (input_records->items + chunk_starts[j]),
// 										  chunk_sizes[j], samples_names->size, uparams);
//
// 						/*printf("segment markers[0] = %ld, segment markers[1] = %ld\n",
// 								segment_markers[0].position, segment_markers[1].position);*/
//
// 						/*for (idxsm=0;idxsm<chunk_sizes[j];idxsm++) {
// 							printf("segment markers[%d] before insertion = %ld\n",
// 									idxsm, (&(segment_markers[idxsm]))->position);
// 							array_list_insert((void*) &(segment_markers[idxsm]), all_markers);
// 						}*/
//
// 						array_list_insert_all((void **) segment_markers, chunk_sizes[j], all_markers);
// 						//for (int k = 0; k < all_markers->size; k++) {
// 						//	marker *m = all_markers->items[k];
// 						//	printf("marker[%d] in %ld\n", k, m->position);
// 						//}
// 						//memcpy(all_markers + j*chunk_sizes[j], segment_markers, chunk_sizes[j] * sizeof(*segment_markers));
// 					}
// 					if (i % 25 == 0) { LOG_INFO_F("*** %dth stats invocation finished\n", i); }
//
// 					free(chunk_starts);
// 					vcf_batch_free(item->data_p);
// 					list_item_free(item);
//
// 					i++;
// 				}
// 			    free(chunk_sizes);
// 				// Decrease list writers count
// 				for (i = 0; i < cparams->num_threads; i++) {
// 					list_decr_writers(output_list);
// 				}
// 	/*
//
// 	        }//end get variants #pragma omp section
//
// 	    }//end omp parallel sections
// 	*/
//
// 	vcf_close(vcf_file_p);
//     free(read_list);
//     free(output_list);
//     free(file_stats);
//
//     *num_samples = samples_names->size;
//
// 	return true;
// }
