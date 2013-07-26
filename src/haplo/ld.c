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

#include "ld.h"

/** Easier to read, define the square for macro expansion */
#define square(x) ((x)*(x))

#define AA  0
#define AB  1
#define BA  2
#define BB  3

/** Pre-compute the val of ln(10) to avoid the calculus repetition (each time we need the log) */
#define LN10 2.302585

/** All the values under this threshold are considered to be 0 */
const double LIMIT_0 = 0.00000001;
#define TOLERANCE 1e-8 //in JAVA is 0.00000001

/** max distance allowed between 2 SNP when computing the table of D'*/
#define MAX_DISTANCE    500000

// Settings for the block finding algorithm
double cutHighCI = 0.98;
double cutLowCI = 0.70;
double cutLowCIVar[] = {0,0,0.80,0.50,0.50};
double maxDist[] = {0,0,20000,30000,1000000};
double recHighCI = 0.90;
double informFrac = 0.95;


inline static bool filter_rating(marker *marker_check) {
    return marker_check->rating > 0;
}

array_list_t *exec_gabriel(array_list_t *markers_arr, int num_samples, haplo_options_data_t *options) { //pairwise_linkage *pairs_table) {
	array_list_t *blocks, *this_block;
	pairwise_linkage **pairs_table, *apair;
	size_t markers_num, *cpyP, numStrong = 0, count = 0, blocks_num = 0, block_idx = 0, numRec = 0, num_in_group = 0, strong_pairs_num = 0;

	//Vector blocks = new Vector();
	//first set up a filter of markers which fail the MAF threshold
	bool *skip_marker, *used_in_block;
	// a block array comprising all the elements between 2 positions from strong_pairs
	//int *this_block;
	// the pairs filtered from the distances table which are strong
	marker_info *strong_pairs;

	//generate table
	pairs_table = generate_pairwise_linkage_tbl(markers_arr, num_samples);
	//stats_mat_result = generate_pairwise_linkage_tbl(stats_num, stats_arr);

	markers_num = markers_arr->size; /*Here is not marker_num - is the filtered markers num with objects having a positive rating */

	strong_pairs = (marker_info *) malloc(sizeof(*strong_pairs) * ((markers_num-1)*markers_num/2));
	skip_marker = (bool *) malloc(sizeof(*skip_marker) * markers_num);
//#pragma omp parallel shared(markers_num, markers_arr, skip_marker)
//{
//#pragma omp for
	size_t x = 0;
	while (x < markers_arr->size){
            marker *temp = (marker *)array_list_get(x, markers_arr);
            if (filter_rating(temp)) {
                //checks minor allele frequency
                skip_marker[x] = temp->maf < options->mafcut;
            }
            
            LOG_DEBUG_F("%d = %.3f -> %d\n", temp->position, temp->maf, skip_marker[x]);
//		else
//			skip_marker[x] = true;
            x++;
//		else
//			array_list_remove_at(x, markers_arr);
	}
//}// End #pragma omp parallel
        
	markers_num = markers_arr->size;/*Here is not marker_num - is the filtered markers num with objects having a positive rating */
uint32_t countn = 0;
        //#pragma omp parallel shared(markers_num, markers_arr, skip_marker)
//{
//#pragma omp for
	//next make a list of marker pairs in "strong LD", sorted by distance apart
//	for (size_t x = 0; x < markers_num-1; x++){
//		for (size_t y = x+1; y < markers_num; y++){
        
	for (int64_t x = markers_num-2; x >=0; x--){
		for (int64_t y =markers_num-1 ; y > x; y--){

/*
        for (int x = 0; x < markers_num-1; x++){
            for (int y = x+1; y < markers_num; y++){
*/
                apair = pairs_table[x*markers_num + y];
                
                if (apair != NULL &&
                    (!skip_marker[x] && !skip_marker[y]) &&
                    apair->lod >= -90 &&
                    (apair->ci_high >= cutHighCI && apair->ci_low >= cutLowCI) //&&//must pass "strong LD" test
//					filter_rating(temp1) &&
//					filter_rating(temp2)
                    )
                {
                    long sep;
                    //compute actual separation
                    sep = labs(((marker *)array_list_get(y, markers_arr))->position -
                                    ((marker *)array_list_get(x, markers_arr))->position);
                    strong_pairs[strong_pairs_num].marker_p1 = x;
                    strong_pairs[strong_pairs_num].marker_p2 = y;
                    strong_pairs[strong_pairs_num++].sep = sep;
                    LOG_DEBUG_F("+ (%d, %d, %ld)\n", x, y, sep);
                }
            }
	}
//}// End #pragma omp parallel
        //printf("\nNUm nulls:%u\n", countn);

	strong_pairs = realloc(strong_pairs, sizeof(*strong_pairs) *strong_pairs_num);

/*
	for (int i=0;i<strong_pairs_num;i++)
            printf("m1=%d m2=%d sep=%ld\n", strong_pairs[i].marker_p1, strong_pairs[i].marker_p2, strong_pairs[i].sep);
*/
        
	// Sort descending the strong pairs
	qsort(strong_pairs, strong_pairs_num, sizeof(*strong_pairs), compare_markers);

/*
	for (int i=0;i<strong_pairs_num;i++)
            printf("m1=%d m2=%d sep=%ld\n", strong_pairs[i].marker_p1, strong_pairs[i].marker_p2, strong_pairs[i].sep);
*/

	//@ Now take this list of pairs with "strong LD" and construct blocks
	used_in_block = (bool *) calloc(markers_num + 1, sizeof(*used_in_block));
//	printf("\nPairs:\n");
//	for(size_t v=0; v<strong_pairs_num;v++)
//	{
//		printf("%d %d, ", strong_pairs[v].marker_p1, strong_pairs[v].marker_p2);
//	}
//	printf("\n");
	//@ Init the final array of blocks
	blocks = array_list_new(10, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
//#pragma omp parallel private(numStrong, numRec, numInGroup, this_block) /
//	shared(strong_pairs)
//	{
//#pragma omp parallel for nowait
	for(size_t v=0; v<strong_pairs_num;v++) {
		numStrong = 0; numRec = 0; num_in_group = 0;
		uint32_t first = strong_pairs[v].marker_p1;
		uint32_t last =  strong_pairs[v].marker_p2;
		uint64_t sep = strong_pairs[v].sep;
		// We know the positions in the array, we can calculate how many elements we have between the 2
		// and to allocate the necessary number of elements
		this_block = array_list_new(last - first + 1, 1.5, COLLECTION_MODE_ASYNCHRONIZED);
				//(int *) malloc(sizeof(*this_block) * (last - first + 1));

		//first see if this block overlaps with another:
		if (used_in_block[first] || used_in_block[last]) continue;

		//next, count the number of markers in the block.
		for (size_t x = first; x <=last ; x++){
			if(!skip_marker[x]) num_in_group++;
		}

		//skip it if it is too long in bases for it's size in markers
		if (num_in_group < 4 && sep > maxDist[num_in_group]) continue;

		//this_block[0] = first;
		cpyP = (size_t *) malloc(sizeof(*cpyP));
		*cpyP = first;
		array_list_insert(cpyP, this_block);
		count = 1;
		//test this block. requires 95% of informative markers to be "strong"
		for (size_t y = first+1; y <= last; y++)
		{
			if (skip_marker[y]) continue;
			cpyP = (size_t *) malloc(sizeof(*cpyP));
			*cpyP = y;
			array_list_insert(cpyP, this_block);
			//count++;
			//this_block[count++] = y;
			//loop over columns in row y
			for (size_t x = first; x < y; x++){
				if (!skip_marker[x]) { //continue;
					apair = pairs_table[x*markers_num + y];
					if (apair != NULL){//	continue;
						//get the right bits
						double lod = apair->lod;
						double lowCI = apair->ci_low;
						double highCI = apair->ci_high;
						if (lod < -90) continue;   //monomorphic marker error
						if (lod == 0 && lowCI == 0 && highCI == 0) continue; //skip bad markers

						//for small blocks use different CI cutoffs
						if (num_in_group < 5){
							if (lowCI > cutLowCIVar[num_in_group] && highCI >= cutHighCI) numStrong++;
						}else{
							if (lowCI > cutLowCI && highCI >= cutHighCI) numStrong++; //strong LD
						}
						if (highCI < recHighCI) numRec++; //recombination
					}
				}
			}
		}
#ifdef DEBUG
		//printf("size %d\n", this_block->size);
#endif
		//change the definition somewhat for small blocks
		if (num_in_group > 3){
			if (numStrong + numRec < 6) continue;
		}else if (num_in_group > 2){
			if (numStrong + numRec < 3) continue;
		}else{
			if (numStrong + numRec < 1) continue;
		}


		//System.out.println(first + " " + last + " " + numStrong + " " + numRec);
		if ((double)numStrong/(double)(numStrong + numRec) > informFrac){ //this qualifies as a block
			//#pragma omp critical
			//{
				//@ List_insert_item(list_item_new(0, 0, this_block), blocks);
				//add to the block list, but in order by first marker number:
				if (blocks->size == 0) {
					array_list_insert(this_block, blocks);
				} else {
					bool placed = false;
					for (block_idx = 0; block_idx < blocks->size; block_idx++){

						if (first < *((size_t *) array_list_get(0, ((array_list_t *)array_list_get(block_idx, blocks)))) ){
							array_list_insert_at(block_idx, this_block, blocks);
							placed = true;
							break;
						}
					}
					//make sure to put in blocks which fall on the tail end
					if (!placed) 
                                            array_list_insert(this_block, blocks);
				}
				//mark the SNPs as used in block to avoid overlapping
				for (size_t used = first; used <= last; used++){
					used_in_block[used] = true;
				}
				blocks_num++;
			//} // end omp critical
		} else {
			//cleanup
			 array_list_free(this_block, free);
		}
	}
	//}// End parallel section

	//list_init("output", num_threads, INT_MAX, result);
	free(pairs_table);
	free(used_in_block);
	free(strong_pairs);
	free(skip_marker);
	return blocks;
}


pairwise_linkage **generate_pairwise_linkage_tbl(array_list_t *markers_arr, int num_samples) {
    pairwise_linkage **stats_mat_result, *aux;
    size_t markers_num = markers_arr->size;
    int idxl_real = 0;
    int idxc_real = 0;
    // Alloc the necessary memo for all linkages
    stats_mat_result = (pairwise_linkage **) malloc(sizeof(pairwise_linkage*) * square(markers_num));
    printf("pairwise linkage table: \n");

    for (size_t idxL = 0; idxL < markers_num-1; idxL++) {
        for (size_t idxC = idxL+1; idxC < markers_num; idxC++) {
            LOG_DEBUG_F("[%d,%d] %ld - %ld = %d (<= %d?)\n", idxC, idxL, 
                        ((marker *)array_list_get(idxC, markers_arr))->position,
                        ((marker *)array_list_get(idxL, markers_arr))->position,
                        ((marker *)array_list_get(idxC, markers_arr))->position - 
                        ((marker *)array_list_get(idxL, markers_arr))->position,
                        MAX_DISTANCE);
            if (MAX_DISTANCE > 0) {
                if (((marker *)array_list_get(idxC, markers_arr))->position -
                        ((marker *)array_list_get(idxL, markers_arr))->position <= MAX_DISTANCE) {
                    aux = compute_pairwise_linkage(markers_arr, idxL, idxC, num_samples);
                    stats_mat_result[idxL * markers_num + idxC] = aux;
                    if (aux) {
                        printf("%04.3f\t%04.2f\t%04.2f\t%04.2f\n", aux->dprime, aux->lod, aux->ci_low, aux->ci_high);
                    }
                } else {
                        stats_mat_result[idxL * markers_num + idxC] = NULL;
                }
            } else {//we need all the possible combinations when the distance between SNPs is not specified
                stats_mat_result[idxL * markers_num + idxC] = compute_pairwise_linkage(markers_arr, idxL, idxC, num_samples);
            }
        }
    }

    // We return the matrix with the linkage between each 2 elements from the stats array
    // It is up to the main mechanism to free the space alloc'ed by this function
    return stats_mat_result;
}


pairwise_linkage *compute_pairwise_linkage( array_list_t *markers_arr, int pos1, int pos2, int num_samples) {
    marker *marker1 = ((marker *)array_list_get(pos1, markers_arr));
    marker *marker2 = ((marker *)array_list_get(pos2, markers_arr));
    //printf("%2.4f/%2.4f ", marker1->maf, marker2->maf);
    
    //check for non-polymorphic markers
    if (marker1->maf < LIMIT_0 || marker2->maf < LIMIT_0){
        return NULL;
    }

    pairwise_linkage *aux_pair;
    // the denominator used to calc D'
    double min, min1, min2, ci_low, ci_high;
    // Get the value of D with D = fAA*fBB - fAB*fBA, where f stands for frequency
    double d;
    // We need 4 elements since there are 4 combinations between 2 SNPs, A1A2 and B1B2
    double *known = (double *) calloc(NUM_BASES, sizeof(*known));
    int *bases_map_two_m1 = (int *) calloc(NUM_KINDS_BASES_F, sizeof(*bases_map_two_m1));
    int *bases_map_two_m2 = (int *) calloc(NUM_KINDS_BASES_F, sizeof(*bases_map_two_m2));
    marker aux;
    unsigned char a1, a2, b1, b2;
    int total_chroms = 0, count = 0;
    double pA1, pA2, pB1, pB2, const_prob, loglike, oldloglike, doublehet = 0.0f,
    		loglike1, loglike0, dpr, total_prob, sum_prob, auxF;
    // Probability for each of the AA, AB, BA and BB combination
    double *prob_haps = (double *) calloc(NUM_BASES, sizeof(*prob_haps));
    // Number of haploids; keep them in double since we need them in divisions and other operations
    //which are not exactly fit for an int
    double *num_haps = (double *) calloc(NUM_BASES, sizeof(*num_haps));

    double tmp[NUM_BASES];
    double lsurface[101];

    aux_pair = malloc(sizeof(*aux_pair));

   //aux = markers_arr[pos1];
    bases_map_two_m1[marker1->reference] = 1;
    bases_map_two_m1[marker1->alternates[0]] = 0;
    //aux = markers_arr[pos2];
    bases_map_two_m2[marker2->reference] = 1;
    bases_map_two_m2[marker2->alternates[0]] = 0;

    for (int idx=0;idx<num_samples;idx++) {
    	//aux = markers_arr[pos1];
    	//get allele 1 from the first 4 bits and shift them to obtain the actual number
    	a1 = (((marker *)array_list_get(pos1, markers_arr))->samples[idx] & MASK_VAR1) >> NUM_BITS_SHIFT;
    	//get allele 2 by masking the first 4 bits and retaining the last 2
    	b1 = (((marker *)array_list_get(pos1, markers_arr))->samples[idx] & MASK_VAR2);
    	//aux = markers_arr[pos2];
    	//do the same with the second marker
    	a2 = (((marker *)array_list_get(pos2, markers_arr))->samples[idx] & MASK_VAR1) >> NUM_BITS_SHIFT;
    	b2 = (((marker *)array_list_get(pos2, markers_arr))->samples[idx] & MASK_VAR2);
    	//printf("%u %u %u %u\n", a1, b1, a2, b2);

    	 if (a1 == 0 || a2 == 0 || b1 == 0 || b2 == 0){
    		 continue;
    	                    //skip missing data
		} else if (((a1 >= NOT_A_BASE || b1 >= NOT_A_BASE) && (a2 >= NOT_A_BASE || b2 >= NOT_A_BASE)) ||
				(a1 >= NOT_A_BASE && !(a2 == b2)) || (a2 >= NOT_A_BASE && !(a1 == b1))){
			doublehet++;
			//find doublehets and resolved haplotypes
		} else if (a1 >= NOT_A_BASE || b1 >= NOT_A_BASE){
			known[bases_map_two_m2[a2]]++;
			known[2 + bases_map_two_m2[a2]]++;
			total_chroms+=2;
		} else if (a2 >= NOT_A_BASE || b2 >= NOT_A_BASE){

			known[2*bases_map_two_m1[a1]]++;
			known[2*bases_map_two_m1[a1] + 1]++;
			total_chroms+=2;
		} else {

			known[bases_map_two_m1[a1]*2 + bases_map_two_m2[a2]]++;
			known[bases_map_two_m1[b1]*2 + bases_map_two_m2[b2]]++;
			total_chroms+=2;
		}
    }

    //another monomorphic marker check
   if ( (known[0] + known[1] == 0 //r1
		   || known[2] + known[3] == 0 //r2
		   || known[0] + known[2] == 0 //c1
		   || known[1] + known[3] == 0) //c2
		   && doublehet == 0){
	   aux_pair->dprime = 1.0;
	   aux_pair->ci_high = 0.0;
	   aux_pair->ci_low = 0.0;
	   aux_pair->lod = 0.0;
	   return aux_pair;
   }

   // Calculate the total number of chromosomes and then the probabilities for each allele
    total_chroms= total_chroms + (2*doublehet);
    pA1 = (known[AA]+known[AB]+doublehet) / (double) total_chroms;
    pB1 = 1.0-pA1;
    pA2 = (known[AA]+known[BA]+doublehet) / (double) total_chroms;
    pB2 = 1.0-pA2;
    const_prob = 0.1f;

    prob_haps[AA] = const_prob;
    prob_haps[AB] = const_prob;
    prob_haps[BA] = const_prob;
    prob_haps[BB] = const_prob;

	/* so that the first count step will produce an
	initial estimate without inferences (this should
	be closer and therefore speedier than assuming
	they are all at equal frequency) */

	count_haps(0, num_haps, known, prob_haps, doublehet, const_prob);
	estimate_p(num_haps, prob_haps, const_prob);

	/* now we have an initial reasonable guess at p we can
	        start the EM - let the fun begin */

	const_prob=0.0f;
	count=1;
	loglike = -999999999.0;

	do {
		oldloglike = loglike;
		count_haps(count, num_haps, known, prob_haps, doublehet, const_prob);
		loglike = (known[AA]*log(prob_haps[AA]) + known[AB]*log(prob_haps[AB]) + known[BA]*log(prob_haps[BA]) +
				known[BB]*log(prob_haps[BB]))/LN10 +
				(doublehet*log(prob_haps[AA]*prob_haps[BB] + prob_haps[AB]*prob_haps[BA]))/LN10;
		if (fabs(loglike-oldloglike) < TOLERANCE) break;
		estimate_p(num_haps, prob_haps, const_prob);
	} while(++count < 1000);
	/* in reality I've never seen it need more than 10 or so iterations
			to converge so this is really here just to keep it from running off into eternity */

	loglike1 = (known[AA]*log(prob_haps[AA]) + known[AB]*log(prob_haps[AB]) +
			known[BA]*log(prob_haps[BA]) + known[BB]*log(prob_haps[BB]) +
			doublehet*log(prob_haps[AA]*prob_haps[BB] + prob_haps[AB]*prob_haps[BA]))/LN10;
	loglike0 = (known[AA]*log(pA1*pA2) + known[AB]*log(pA1*pB2) + known[BA]*log(pB1*pA2) +
			known[BB]*log(pB1*pB2) + doublehet*log(2*pA1*pA2*pB1*pB2))/LN10;
	/*printf("test=%f\n", pA1*pB2);
	printf("test=%f\n", known[BA]*log(pB1*pA2));
	printf("test=%f\n", known[BB]*log(pB1*pB2));
	printf("test=%f\n", doublehet*log(2*pA1*pA2*pB1*pB2));
	printf("loglike1=%f loglike0=%f\n", loglike1, loglike0);*/

	//LOD is the log of the likelihood odds ratio, a measure of confidence in the value of D'
	aux_pair->lod = loglike1 - loglike0;
	/*for (int it=0;it<4;it++) {
		printf("%f ", known[it]);
		//printf("%f ", prob_haps[it]);
	}
	printf("\n%f %f\n", loglike1, loglike0);*/

	d =prob_haps[AA]*prob_haps[BB] - prob_haps[AB]*prob_haps[BA];;
	if (d < 0) {
		/* flip matrix so we get the positive D' */
		/* flip AA with AB and BA with BB */
		auxF=prob_haps[AA]; prob_haps[AA]=prob_haps[AB]; prob_haps[AB]=auxF;
		auxF=prob_haps[BB]; prob_haps[BB]=prob_haps[BA]; prob_haps[BA]=auxF;
		/* flip frequency of second allele */
		pA2 = pA2 + pB2;
		pB2 = pA2 - pB2;
		pA2 = pA2 - pB2;
		//pA2=pB2;pB2=temp;
		/* flip counts in the same fashion as p's */
		auxF=num_haps[AA]; num_haps[AA]=num_haps[AB]; num_haps[AB]=auxF;
		auxF=num_haps[BB]; num_haps[BB]=num_haps[BA]; num_haps[BA]=auxF;
		/* num has now undergone a sign change */
		d = prob_haps[AA]*prob_haps[BB] - prob_haps[AB]*prob_haps[BA];
		/* flip known array for likelihood computation */
		auxF=known[AA]; known[AA]=known[AB]; known[AB]=auxF;
		auxF=known[BB]; known[BB]=known[BA]; known[BA]=auxF;
	}

	min1 = (prob_haps[AA]+prob_haps[BA])*(prob_haps[BA]+prob_haps[BB]);
	min2 = (prob_haps[AA]+prob_haps[AB])*(prob_haps[AB]+prob_haps[BB]);
	if (min1 < min2) { min = min1; }
	else { min = min2; }

	aux_pair->dprime = d / min;//compute_dprime(prob_haps, num_haps, known, pA2, pB2);


	//calculate ci_high and ci_low
	 for (int i=0; i<=100; i++) {
		dpr = (double)i*0.01f;
		tmp[AA] = dpr*min + pA1*pA2;
		tmp[AB] = pA1-tmp[AA];
		tmp[BA] = pA2-tmp[AA];
		tmp[BB] = pB1-tmp[BA];
		if (i==100) {
			/* one value will be 0 */
			if (tmp[AA] < 1e-10) tmp[AA]=1e-10;
			if (tmp[AB] < 1e-10) tmp[AB]=1e-10;
			if (tmp[BA] < 1e-10) tmp[BA]=1e-10;
			if (tmp[BB] < 1e-10) tmp[BB]=1e-10;
		}
		lsurface[i] = (known[AA]*log(tmp[AA]) + known[AB]*log(tmp[AB]) + known[BA]*log(tmp[BA]) +
				known[BB]*log(tmp[BB]) + doublehet*log(tmp[AA]*tmp[BB] + tmp[AB]*tmp[BA]))/LN10;
	}

	/* Confidence bounds #2 - used in Gabriel et al (2002) - translate into posterior dist of D' -
	assumes a flat prior dist. of D' - someday we may be able to make
	this even more clever by adjusting given the distribution of observed
	D' values for any given distance after some large scale studies are complete */

	total_prob=sum_prob=0.0;

	for (int i=0; i<=100; i++) {
		lsurface[i] -= loglike1;
		lsurface[i] = pow(10.0,lsurface[i]);
		total_prob += lsurface[i];
	}

	for (int i=0; i<=100; i++) {
		sum_prob += lsurface[i];
		if (sum_prob > 0.05*total_prob &&
				sum_prob-lsurface[i] < 0.05*total_prob) {
			ci_low = i-1;
			break;
		}
	}

	sum_prob=0.0;
	for (int i=100; i>=0; i--) {
		sum_prob += lsurface[i];
		if (sum_prob > 0.05*total_prob &&
				sum_prob-lsurface[i] < 0.05*total_prob) {
			ci_high = i+1;
			break;
		}
	}
	if (ci_high > 100){ ci_high = 100; }

	aux_pair->ci_high = ci_high /100.0;
	aux_pair->ci_low = ci_low / 100.0;

	//printf("D' = %G\n", aux_pair.dprime);

	//it's time to clear the garbage and move on
	free(known);
	free(prob_haps);
	free(num_haps);
	free(bases_map_two_m1);
	free(bases_map_two_m2);
	return aux_pair;
}/*compute_pairwise_linkage*/


__inline__ int is_same_ref(char **alt1, char **alt2, int numAlts1, int numAlts2) {
	if (numAlts1 != numAlts2)
		return 1;//false if it doesn't exist
	else {
		int found;
		// right now check if the elements are equal or not
		return strcmp(alt1[0], alt2[0]);
	}
}/*is_same_ref*/

__inline__ int is_alternate(char *ref, char **alt, int numAlts) {
	for (int i=0;i<numAlts;i++)
		if (strcmp(ref, alt[i]) == 0)
			return 0;//true if we found the ref in alts
	return 1;//false if it doesn't exist
}/*is_alternate*/

void count_haps(int em_round, double *num_haps, double *known, double *prob_haps,
		double doublehet, double const_prob){
        /* only the double heterozygote [AB][AB] results in
        ambiguous reconstruction, so we'll count the obligates
        then tack on the [AB][AB] for clarity */

	num_haps[AA] = known[AA];
	num_haps[AB] = known[AB];
	num_haps[BA] = known[BA];
	num_haps[BB] = known[BB];
	if (em_round) {
		num_haps[AA] += doublehet* (prob_haps[AA]*prob_haps[BB])/((prob_haps[AA]*prob_haps[BB])+(prob_haps[AB]*prob_haps[BA]));
		num_haps[BB] += doublehet* (prob_haps[AA]*prob_haps[BB])/((prob_haps[AA]*prob_haps[BB])+(prob_haps[AB]*prob_haps[BA]));
		num_haps[AB] += doublehet* (prob_haps[AB]*prob_haps[BA])/((prob_haps[AA]*prob_haps[BB])+(prob_haps[AB]*prob_haps[BA]));
		num_haps[BA] += doublehet* (prob_haps[AB]*prob_haps[BA])/((prob_haps[AA]*prob_haps[BB])+(prob_haps[AB]*prob_haps[BA]));
	}
}/*count_haps*/

void estimate_p(double *num_haps, double *prob_haps, double const_prob) {
	double total= num_haps[AA]+num_haps[AB]+num_haps[BA]+num_haps[BB]+(4.0*const_prob);
	prob_haps[AA]=(num_haps[AA]+const_prob)/total; if (prob_haps[AA] < 1e-10) prob_haps[AA]=1e-10;
	prob_haps[AB]=(num_haps[AB]+const_prob)/total; if (prob_haps[AB] < 1e-10) prob_haps[AB]=1e-10;
	prob_haps[BA]=(num_haps[BA]+const_prob)/total; if (prob_haps[BA] < 1e-10) prob_haps[BA]=1e-10;
	prob_haps[BB]=(num_haps[BB]+const_prob)/total; if (prob_haps[BB] < 1e-10) prob_haps[BB]=1e-10;
}/*estimate_p*/

void free_marker_array(marker *array, int len) {
	//iterate through the items and eliminate them (deep freeing)
	for (int idx=0;idx<len;idx++) {
		if ((array[idx]).samples != NULL)
			free((array[idx]).samples);
		if ((array[idx]).alternates != NULL)
			free((array[idx]).alternates);
	}
	free(array);
}/*free_marker_array*/

inline static int compare_markers(const void *markeri1, const void *markeri2) {
 if(((marker_info *)markeri1)->sep > ((marker_info *)markeri2)->sep) return -1;
   else if(((marker_info *)markeri1)->sep < ((marker_info *)markeri2)->sep) return 1;
   else return 0;
}

inline static int compare_markers_asc(const void *markeri1, const void *markeri2) {
 if(((marker_info *)markeri1)->sep < ((marker_info *)markeri2)->sep) return -1;
   else if(((marker_info *)markeri1)->sep > ((marker_info *)markeri2)->sep) return 1;
   else return 0;
}
