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

#include "file_handling.h"

bool get_markers_array(array_list_t *all_markers, const shared_options_data_t *share_data,
		const haplo_options_data_t *haplo_data, unsigned int  *num_samples)
{
    int ret_code;
    list_t *read_list = (list_t*) malloc (sizeof(*read_list));
    int header_written = 0;

    int i = 0;
/*
    int add_ret_code = 0;
    list_item_t* item = NULL;
    ped_batch_t *batch = NULL;
    list_item_t *batch_item = NULL;
    
    // STEP 1 - parse PED file and get the individuals with their sex and parents
    ped_file_t* file = ped_open(share_data->ped_filename);
    
    printf("markers 1.1\n");
    
    // ret_code = ped_read(file);
    ret_code = ped_read_batches(read_list, share_data->batch_lines, file);
    
    printf("markers 1.2\n");
    while ( (item = list_remove_item(read_list)) != NULL ) {
        batch = (ped_batch_t*) item->data_p;
        if (i % 200 == 0) 
        {
            int debug = 1;
            LOG_DEBUG_F("Batch %d reached by thread %d - %zu/%zu records \n", i, omp_get_thread_num(), 
                ((ped_batch_t*) item->data_p)->length, ((ped_batch_t*) item->data_p)->max_length);
        }

        while ( (batch_item = list_remove_item_async(batch)) != NULL) {
            add_ret_code = add_ped_record(batch_item->data_p, file);
            if (add_ret_code > 0) {
                LOG_ERROR_F("%s - %s\n", ((ped_record_t*) batch_item->data_p)->family_id, get_ped_semantic_error_msg(add_ret_code));
            }
        }

        ped_batch_free(item->data_p);
        list_item_free(item);
        i++;
    }
        
    printf("markers 2\n");
       
     list_decr_writers(read_list);
*/
    
    //STEP 2 - Resolve families and approve the list of individuals to be used 
    //         in computation, the rest of them being discarded
    
    
    // STEP 3 - Parse the VCF and get the samples (take into consideration the
    //          previous list of allowed samples after creating the families
    
    vcf_file_t* vcf_file_p = vcf_open(share_data->vcf_filename, share_data->max_batches);
    array_list_t *samples_names = vcf_file_p->samples_names;
//#pragma omp parallel sections //private(start, stop, total)
//{
//    #pragma omp section
//    {
        LOG_DEBUG_F("Thread %d reads the VCF file\n", omp_get_thread_num());
        // Reading
        //start = omp_get_wtime();

/*
        ret_code = vcf_parse_batches(share_data->batch_lines, vcf_file_p, 1);
*/
        ret_code = vcf_read(vcf_file_p, 1, 1000, 1);
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

        vcf_batch_t *vcf_batch = NULL;
        vcf_batch = fetch_vcf_batch(vcf_file_p);
        printf("markers 4.1\n");
        marker **segment_markers = get_markers(vcf_batch->records, samples_names->size, haplo_data);
        printf("markers 4.2\n");
        // Heck knows why it's impossible to get array's size directly with the indirection operator. level 2
        array_list_t *temp = vcf_batch->records;
        size_t num_rec = temp->size;

        array_list_insert_all((void **) segment_markers, num_rec, all_markers);
        printf("markers 4.3\n");

        vcf_batch_free(vcf_batch);
        printf("markers 4.4\n");

        i++;

        notify_end_parsing(vcf_file_p);
//    }
//}
    *num_samples = samples_names->size;

    vcf_close(vcf_file_p);
    
    return true;
}

marker **get_markers(array_list_t *variants, const unsigned int num_samples,
		const haplo_options_data_t *params) {
    char *copy_buf, *str_dup_buf;
    // allele 1 & 2 definition to be packed further in a unsigned char
    unsigned char a1 = 0, a2 = 0, auxuc;

    int num_alternates, gt_pos, cur_pos, numa1, numa2;
    int allele1, allele2, alleles_code;

    // Counter for the bases in a file
    int64_t alleles_count[NUM_KINDS_BASES_F], sum=0, num=0, mincount = -1;
    int female_count = 0;
    uint32_t founderHetCount = 0, *founderHomCount;
    double missing_data, called;

    // Variant stats management
    vcf_record_t *record;
    size_t num_variants = variants->size;
    //num_samples = variants[0]->samples->size;
    marker **result = malloc(num_variants * sizeof(marker*));
    list_t *output_list = malloc (sizeof(list_t));
    list_init("stats", 1, 10, output_list);
    
    bool is_chr_x;
    founderHomCount = malloc(sizeof(*founderHomCount) * NUM_KINDS_BASES_F);
    for (size_t i = 0; i < num_variants; i++) {
    	//init the alleles' counter, better than a for because the compiler can decide how to expand it
    	memset(alleles_count, 0, sizeof(int64_t)*NUM_KINDS_BASES_F);

        record = (vcf_record_t *) array_list_get(i, variants);
        result[i] = malloc(sizeof(**result));
        // Reset counters
        numa1 = 0;
        numa2 = 0;

        // Create list of alternates
        copy_buf = strndup(record->alternate, record->alternate_len);
        //for now use this, even if it is possible to have more than 1 base as an alternate
        char **alt = split(copy_buf, ",", &num_alternates);
        result[i]->alternates =  malloc(sizeof(*(result[i]->alternates))*num_alternates);
        for (uint32_t idxAlt=0;idxAlt<num_alternates; idxAlt++) {
        	result[i]->alternates[idxAlt] = base_to_int(alt[idxAlt][0]);
        	//printf("alt %u\n", result[i]->alternates[idxAlt]);
        }
        for (int j = 0; j < num_alternates; j++) {
            free(alt[j]);
        }
        free(alt);
        free(copy_buf);
        
        result[i]->reference = base_to_int(record->reference[0]);//strtol(record->reference, NULL, 10);
        result[i]->num_alternates = num_alternates;

        // Get position where GT is in sample
        copy_buf = strndup(record->format, record->format_len);
        gt_pos = get_field_position_in_format("GT", copy_buf);
        free(copy_buf);
        
        if (gt_pos < 0) { continue; }   // This variant has no GT field
        result[i]->samples = malloc(num_samples * sizeof(*(result[i]->samples)));
        result[i]->position = record->position;
        called = 0.0;
        missing_data = 0.0;
        female_count = 0;
        // check if the current chromosome record is the X one ( we have some special treatment for it)
        is_chr_x = is_x(record->chromosome);
        result[i]->is_x = is_chr_x;
        founderHetCount = 0;
        memset(founderHomCount, 0, sizeof(*founderHomCount) * NUM_KINDS_BASES_F);

        // Traverse samples and find the existing and missing alleles
        for(uint32_t j = 0; j < num_samples; j++) {
            char *sample = strdup((char *) array_list_get(j, record->samples));
            // Get to GT position
            alleles_code = get_alleles(sample, gt_pos, &allele1, &allele2);
            free(sample);

            // We need to determine how many times the reference and the alternates appear in the file
            // We consider for now that ONLY ONE alt EXISTS
            a1 == REF_ALLELE_IDX ? numa1++ : numa2++;
            a2 == REF_ALLELE_IDX ? numa1++ : numa2++;

            // Translate allele to a base to be easier to handle later; we don't need the idx, we need the
            // base representation
            a1 = allele_translation(allele1, result[i]->alternates, result[i]->reference);
            a2 = allele_translation(allele2, result[i]->alternates, result[i]->reference);
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
            result[i]->samples[j] = (a1 << NUM_BITS_SHIFT) + a2 ;
        }
        
        // Check the type of ref allele determination required and switch if necessary (HaploView mode)
        // TODO after testing, remove this Haploview mode! Is it sooo incorrect
        if (alleles_count[1] >= alleles_count[0]) {
                auxuc = result[i]->reference;
                result[i]->reference = result[i]->alternates[0];
                result[i]->alternates[0] = auxuc;
        }// Else don't change anything

        file_stats_t *file_stats = file_stats_new();
        get_variants_stats(&record, 1, NULL, NULL, 0, output_list, file_stats);
        list_item_t* item = list_remove_item(output_list);
        variant_stats_t *variant_stats = item->data_p;
        
        float maf = 1.0f;
        for (int i = 0; i < variant_stats->num_alleles; i++) {
            if (variant_stats->alleles_freq[i] < maf) {
                maf = variant_stats->alleles_freq[i];
            }
        }
        result[i]->maf = maf;
        double missing_genotypes_percent = 1 - (double) variant_stats->missing_genotypes / num_samples ;
        
        double pvalue = get_pvalue(founderHomCount, NUM_KINDS_BASES_F, founderHetCount);
        result[i]->rating = calc_rating(missing_genotypes_percent, pvalue, variant_stats->mendelian_errors, maf, params);
        LOG_DEBUG_F("%.*s) MAF = %.3f\tmendel_err = %d\trating = %d\n", record->id_len, record->id, 
                    result[i]->maf, variant_stats->mendelian_errors, result[i]->rating);

        list_decr_writers(output_list);
        list_item_free(item);
        variant_stats_free(variant_stats);
        file_stats_free(file_stats);
    }
    free(founderHomCount);
    return result;
}


inline static bool is_x(const char *chromosome) {
    return strcmp(chromosome, "X") == 0 || strcmp(chromosome, "x") == 0;
}

inline static int calc_rating(double genopct, double pval, int menderr, double maf, const haplo_options_data_t *options) {
    int rating = 0;
    if (genopct < options->fgcut){
        rating -= 2;
    }
    if (pval < options->hwcut){
        rating -= 4;
    }
    if (menderr > options->mendcut){
        rating -= 8;
    }
    if (maf < options->mafcut){
        rating -= 16;
    }
    if (rating == 0){
        rating = 1;
    }
    return rating;
} /*calc_rating*/

static double get_pvalue(const uint32_t *parent_hom, size_t parent_hom_len, uint32_t parent_het) {
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
		pvalue = hw_calculate(homA, parent_het, homB);
	}
	return pvalue;
} /*get_pvalue*/

 static double hw_calculate(uint32_t obsAA, uint32_t obsAB, uint32_t obsBB) {
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
            LOG_FATAL_F("Exception: HW test: %d heterozygotes but only %d rare alleles.", hets, rare);
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
