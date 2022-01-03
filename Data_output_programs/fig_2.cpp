//sib_recessivedrive.c to run recurrence equations of sib mating with recessive partially lethal drive  28 xi 21
// Haploids with 2 loci

#define _CRT_SECURE_NO_WARNINGS

#include <stdio.h>
#include <math.h>
#include <stdbool.h>

typedef struct equ_data {
	float d;
	float wb;
	float sm_feq;
	bool equal;
}Equ_data;


float calc_d_male(float y[2][2][2][2]) {
	float x[2][2];

	x[0][0] = y[0][0][0][0] + y[0][0][1][0] + y[0][0][0][1] + y[0][0][1][1];
	x[0][1] = y[0][1][0][0] + y[0][1][1][0] + y[0][1][0][1] + y[0][1][1][1];
	x[1][0] = y[1][0][0][0] + y[1][0][1][0] + y[1][0][0][1] + y[1][0][1][1];
	x[1][1] = y[1][1][0][0] + y[1][1][1][0] + y[1][1][0][1] + y[1][1][1][1];



	return x[0][0] * x[1][1] - x[0][1] * x[1][0];
}

float calc_d_female(float y[2][2][2][2]) {
	float x[2][2];

	x[0][0] = y[0][0][0][0] + y[1][0][0][0] + y[0][1][0][0] + y[1][1][0][0];
	x[0][1] = y[0][0][0][1] + y[1][0][0][1] + y[0][1][0][1] + y[1][1][0][1];
	x[1][0] = y[0][0][1][0] + y[1][0][1][0] + y[0][1][1][0] + y[1][1][1][0];
	x[1][1] = y[0][0][1][1] + y[1][0][1][1] + y[0][1][1][1] + y[1][1][1][1];

	return x[0][0] * x[1][1] - x[0][1] * x[1][0];
}

void run_case(float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


//	printf("\n\n\tINPUT \n\n\
//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
//sA = sibmating rate of S (wildtype is outcrossed)\n\
//sigma = fitness of offspring from sib-mated parents\n\
//K = discounting of males who sib mat (typically 0 or 1 --\n\
//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
//DD[2] = fitness of drive homozygotes \n");


	//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	printf(" number of cycles between outputting results and final number of generations for the run\n");
//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 200;
	printf("sib mating only for S\n");
	printf("gen  \tsdsd \tsdSd \tSdsd \tSdSd \tsDsD \tsDSD \tSDsD \tSDSD \tsA (%0.2f) \tWbar \tD \n", sA);



	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {


			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop

		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];

		

		printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n", gen,
			allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0],
			allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1], Wbar, Df);

		printf("sibm: \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",
			sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0],
			sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);



	}
	printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);

}


//TODO: do actual modification.
void run_case_precentages(float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


	//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");


		//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	printf(" number of cycles between outputting results and final number of generations for the run\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 200;
	printf("sib mating only for S\n");
	printf("gen  \tsdsd \tsdSd \tSdsd \tSdSd \tsDsD \tsDSD \tSDsD \tSDSD \tsA (%0.2f) \tWbar \tD \n", sA);



	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {


			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop

		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];



		printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n", gen,
			allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0],
			allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1], Wbar, Df);

		printf("sibm: \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",
			sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0],
			sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);



	}
	printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);

}






bool run_case_fixation(FILE* of, float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


	//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");


		//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	//printf(" number of cycles between outputting results and final number of generations for the run\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 1000;
	//printf("sib mating only for S\n");
	if (of != NULL) {
		fprintf(of, "gen  \tsdsd \tsdSd \tSdsd \tSdSd \tsDsD \tsDSD \tSDsD \tSDSD \tsA (%0.2f) \tWbar \tD\n", sA);
	}


	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {


			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop
		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];
		if (of != NULL) {
			fprintf(of,"%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t\t\t%1.3f \t%1.3f\n", gen,
				allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0],
				allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1], Wbar, Df);

			fprintf(of,"sibm:  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",
				sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0],
				sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);
		}


		if (Df > 0.9995) {
			return 1;
		}
	}
	//printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);

	return 0;

}

bool last_ten_equal(float lt[10]) {
	for (int i = 0; i < 9; i++)
	{
		if (roundf(lt[i]*10000)/10000 == roundf(lt[i + 1]*10000)/10000) {

		}
		else {
			return false;
		}
	}
	return true;
}

//origonal explantation 
//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");

//params
// of -- out file
// sib_mating_freq -- fequency of sibling mating
// drive_freq -- the fequency of the drive
// dib_mat_percent -- sA the sibmating rate of S
// sib_mat_fit -- sigma the fitness of the offspring from sib-mated parents

Equ_data run_case_equ(FILE* of, float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float last_ten[10] = { 0 };


	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


	//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");


		//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	//printf(" number of cycles between outputting results and final number of generations for the run\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 100000;
	//printf("sib mating only for S\n");
	if (of != NULL) {
		fprintf(of, "gen  \tsdsd \tsdSd \tSdsd \tSdSd \tsDsD \tsDSD \tSDsD \tSDSD \tsA (%0.2f) \tWbar \tD\n", sA);
	}


	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {
			

			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop
		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];
		float total_sib_mated_feq = 0;
		float total_sum = 0;
		for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
			total_sib_mated_feq += allmated[1][i2][1][i4] +

				0.5 * allmated[1][i2][0][i4] +

				0.5 * allmated[0][i2][1][i4];
			total_sum += allmated[i1][i2][i3][i4];
		}
		//printf("%.3f total sum\n", total_sum);
		total_sib_mated_feq /= 4;

		if (gen < 10) {
			last_ten[gen] = Df;
		}
		else {
			for (int i = 9; i > 0; --i)
			{
				last_ten[i] = last_ten[i - 1];
			}
			last_ten[0] = Df;
		}
		if (of != NULL) {
			fprintf(of, "%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t\t\t%1.3f \t%1.3f\n", gen,
				allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0],
				allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1], Wbar, Df);

			fprintf(of, "sibm:  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",
				sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0],
				sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);
		}


		if (gen == totcycles-1) {
			Equ_data d;
			d.equal = 1;
			d.wb = Wbar;
			d.d = Df;
			d.sm_feq = total_sib_mated_feq;
			return d;
		}
	}
	//printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);
	Equ_data rv;
	rv.equal = 0;
	return rv;

}




//the percent of femals that have childered who sibling mate has very little effect
void max_strength_vs_sib_mat_f_feq(FILE * of) {
	fprintf(of,"SMfreq, max (s-1)\n");
	for (float sib_mat_num = 0.01; sib_mat_num < 1; sib_mat_num+=0.01)
	{
		for (float drive_fit = 1; drive_fit > 0; drive_fit -= 0.01)
		{
			if (!run_case_fixation(NULL, sib_mat_num, 0.01, drive_fit, 0.01, 1)) {
				fprintf(of, "%.3f,%.3f\n",sib_mat_num, drive_fit);
				break;
			}
			
		}
	}
}

void max_strength_over_sibling_mating_precent(FILE * of) {
	fprintf(of, "Sib_mating_Percent(sA),max_drive_fit(s-1)\n");
	for (float sib_mat_per = 0.01; sib_mat_per < 1; sib_mat_per += 0.01)
	{
		for (float drive_fit = 1; drive_fit > 0; drive_fit -= 0.01)
		{
			if (!run_case_fixation(NULL, 0.01, 0.01, drive_fit, sib_mat_per, 1)) {
				fprintf(of, "%.3f,%.3f\n", sib_mat_per, drive_fit+0.01);
				break;
			}

		}
	}
}


void print_equ(FILE * of, float m) {
	fprintf(of, "drv stregth (s), final feq drive, final wbar,sm_feq\n");
	for (float i = 1; i > 0.01; i-= 0.01)
	{
		Equ_data d = run_case_equ(NULL, 0.01, 0.01, i, m, 1);
		if (d.equal) {
			fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", 1-i, d.d, d.wb,d.sm_feq);
		}
		else {
			fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", 1-i, -1.0, -1.0,d.sm_feq);
			
		}
	}
}

void print_equ_polymorphic(FILE* of, float m) {
	fprintf(of, "drv stregth (s), final feq drive, final wbar,sm_feq\n");
	for (float i = .6; i > 0.4; i -= 0.001)
	{
		Equ_data d = run_case_equ(NULL, 0.01, 0.01, i, m, 1);
		if (d.equal) {
			if (d.d != 1 || d.d != 0) {
				fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", 1 - i, d.d, d.wb, d.sm_feq);
			}
		}
		else {
			fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", 1 - i, -1.0, -1.0, d.sm_feq);

		}
	}
}



void run_case_Ds(FILE* of, float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


	//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");


		//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	//printf(" number of cycles between outputting results and final number of generations for the run\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 100;
	//printf("sib mating only for S\n");
	if (of != NULL) {
		//fprintf(of, "gen  \tsdsd \tsdSd \tSdsd \tSdSd \tsDsD \tsDSD \tSDsD \tSDSD \tsA (%0.2f) \tWbar \tD\n", sA);
		fprintf(of, "gen,mD,fD,smmD,smfD,DrvFeq,df,smr\n");
		
	}


	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {


			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop
		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];
		if (of != NULL) {
			float mD = calc_d_male(allmated);
			float fD = calc_d_female(allmated);

			float sib_sum = sibmated[0][0][0][0] + sibmated[0][0][1][0] + sibmated[1][0][0][0] + sibmated[1][0][1][0] +
				sibmated[0][1][0][1] + sibmated[0][1][1][1] + sibmated[1][1][0][1] + sibmated[1][1][1][1];
			
			float sibnor[2][2][2][2];
			


			for (int j = 0; j < 2; j++)
			{
				for (int l = 0; l < 2; l++)
				{
					for (int m = 0; m < 2; m++)
					{
						for (int n = 0; n< 2; n++)
						{
							sibnor[j][l][m][n] = sibmated[j][l][m][n] / sib_sum;

						}

					}

				}
			

			}
			float smD = calc_d_male(sibnor);
			float sfD = calc_d_female(sibnor);
			if (gen == 1) {
				fprintf(of, "%d,%.6f,%.6f,%.6f,%.6f, %.4f,", gen, mD, fD, smD, sfD, Df);
				fprintf(of, "%.2f,%.2f\n", drive_fit, sib_mate_percent);
			}
			else {
				fprintf(of, "%d,%.6f,%.6f,%.6f,%.6f, %.4f\n", gen, mD, fD, smD, sfD, Df);
			}
			/*fprintf(of, "%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t\t\t%1.3f \t%1.3f\n", gen,
				allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0],
				allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1], Wbar, Df);

			fprintf(of, "sibm:  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",
				sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0],
				sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);*/
		}


		
	}
	//printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);

	

}

void run_case_dif(FILE * of,float sib_mating_freq, float drive_freq, float drive_fit, float sib_mate_percent,
	float sib_mat_fit) {

	float transm[2][3] = { 1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
					  0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles



	int dummy, dummy1, i, j, k, totcycles, printcycles, a, b, c, d, i1, i2, i3, i4, gen;


	float
		sA, sa = 0.0, DD[3], sibprob[2], sigma, K, sum, Wbar = 1,
		allmated[2][2][2][2],   // ordered genotype for all mated females 
		outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
		sibmated[2][2][2][2],   // ordered genotype for all  females whose sires are sibs
		outmale[2][2],    // haploid males that will do the outcrossing
		outfemale[2][2];    // haploid females that will do the outcrossing

		// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
		// thus mated females have 4 alleles if we count their (haploid) sires
		// first locus is sibmating control (s, S), second locus is drive (d,D)
		// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
		// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
		// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



		// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
		allmated[i1][i2][i3][i4] = 0.0;


	//	printf("\n\n\tINPUT \n\n\
	//allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
	//sA = sibmating rate of S (wildtype is outcrossed)\n\
	//sigma = fitness of offspring from sib-mated parents\n\
	//K = discounting of males who sib mat (typically 0 or 1 --\n\
	//\t\t K=0 means sib-mated sons do not enter the random pool, 1 means they all do)\n\
	//\t\t but K=1 gives an unfair advantage to selfing that we probably don't want\n\
	//DD[2] = fitness of drive homozygotes \n");


		//scanf("%f %f %f %f %f %f", &allmated[1][0][1][0], &allmated[0][1][0][1], &sA, &sigma, &K, &DD[2]);


	allmated[1][0][1][0] = sib_mating_freq;
	allmated[0][1][0][1] = drive_freq;
	DD[2] = drive_fit;
	sA = sib_mate_percent;
	sigma = sib_mat_fit;
	K = 0;
	allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

	sibprob[0] = 0.0;
	sibprob[1] = sA;
	DD[0] = DD[1] = 1.0;

	gen = 0;


	//printf(" number of cycles between outputting results and final number of generations for the run\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = 100;
	//printf("sib mating only for S\n");
	fprintf(of, "gen,out_sdsd,out_sdSd,out_Sdsd,out_SdSd,out_sDsD,out_sDSD,out_SDsD,out_SDSD,Wbar,D,");
	fprintf(of,"sib_sdsd,sib_sdSd,sib_Sdsd,sib_SdSd,sib_sDsD,sib_sDSD,sib_SDsD,sib_SDSD\n");




	for (; gen < totcycles;) {
		for (j = 1; j <= printcycles; j++) {


			// Zero the matrix before use
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				sibmated[i1][i2][i3][i4] = 0.0;

			// sib mating from broods of all mated females -- only if she has allele 1 at first locus
			// this is set with sibprob[], which may allow us to generalize the pgm
			for (a = 0; a < 2; a++) for (b = 0; b < 2; b++) for (c = 0; c < 2; c++) for (d = 0; d < 2; d++)
				for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
					sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] * sibprob[i1] * transm[a][i1 + i3] * transm[c][i1 + i3] * transm[b][i2 + i4] * transm[d][i2 + i4] * sigma;

			// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
			// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
			// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
			// inbreeding depression taken care of by sigma fitness effect
			// this is the crux of inheritance, so we need to triple check that this does what we want
			// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
			// had the allele
			// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
			// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
			// This works if we are careful about it, but I'd like to think about it more.  
			// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
			// are themselves reduced

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) outmale[i1][i2] = outfemale[i1][i2] = 0.0;



			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {

				outmale[a][c] += allmated[i1][i2][i3][i4] * (1.0 + (K - 1.0) * sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

				outfemale[a][c] += allmated[i1][i2][i3][i4] * (1.0 - sibprob[i1]) * transm[a][i1 + i3] * transm[c][i2 + i4];

			}


			for (a = 0, sum = 0.0; a < 2; a++)  for (c = 0; c < 2; c++) sum += outmale[a][c];


			for (a = 0; a < 2; a++)  for (c = 0; c < 2; c++) outmale[a][c] /= sum;   // males normalized



			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				outmated[i1][i2][i3][i4] = outfemale[i1][i2] * outmale[i3][i4];

			// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
			// the diploid phase will be DD, thus siffering a fitness loss 
			// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 + i4 == 2) {
					outmated[i1][i2][i3][i4] *= DD[2];
					sibmated[i1][i2][i3][i4] *= DD[2];
				}
			}


			// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
			// so all progeny from Dd or dD mated females will be D
			//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++) {
				if (i2 != i4) {			// for i2 != i4, one must be D
					outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
					outmated[i1][i2][i3][i4] = 0.0;

					sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4];
					sibmated[i1][i2][i3][i4] = 0.0;

				}
			}

			// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
			// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

			// Don't need to zero allmated
			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

			for (i1 = 0, Wbar = 0.0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				Wbar += allmated[i1][i2][i3][i4];

			for (i1 = 0; i1 < 2; i1++) for (i2 = 0; i2 < 2; i2++) for (i3 = 0; i3 < 2; i3++) for (i4 = 0; i4 < 2; i4++)
				allmated[i1][i2][i3][i4] /= Wbar;


			gen++;
		}  // end of inner print loop

		float Df = allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];

		float out_sum = outmated[0][0][0][0] + outmated[0][0][1][0] + outmated[1][0][0][0] + outmated[1][0][1][0]+
			outmated[0][1][0][1] + outmated[0][1][1][1] + outmated[1][1][0][1] + outmated[1][1][1][1];

		float sib_sum = sibmated[0][0][0][0] + sibmated[0][0][1][0] + sibmated[1][0][0][0] + sibmated[1][0][1][0] +
			sibmated[0][1][0][1] + sibmated[0][1][1][1] + sibmated[1][1][0][1] + sibmated[1][1][1][1];

		float all_sum = allmated[0][0][0][0] + allmated[0][0][1][0] + allmated[1][0][0][0] + allmated[1][0][1][0] +
			allmated[0][1][0][1] + allmated[0][1][1][1] + allmated[1][1][0][1] + allmated[1][1][1][1];
		//printf("All %.3f, ", all_sum);
		fprintf(of,"%d,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,", gen,
			outmated[0][0][0][0]/out_sum, outmated[0][0][1][0]/out_sum, outmated[1][0][0][0]/out_sum, outmated[1][0][1][0]/out_sum,
			outmated[0][1][0][1]/out_sum, outmated[0][1][1][1]/out_sum, outmated[1][1][0][1]/out_sum, outmated[1][1][1][1]/out_sum, Wbar, Df);

		fprintf(of,"%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f,%1.3f\n",
			sibmated[0][0][0][0]/sib_sum, sibmated[0][0][1][0] / sib_sum, sibmated[1][0][0][0] / sib_sum, sibmated[1][0][1][0] / sib_sum,
			sibmated[0][1][0][1] / sib_sum, sibmated[0][1][1][1] / sib_sum, sibmated[1][1][0][1] / sib_sum, sibmated[1][1][1][1] / sib_sum);



	}
	printf("sA = %1.2f  sigma = %1.2f DD[2] = %1.4f  \n\n", sA, sigma, DD[2]);

}



int main()
{
	//system("dir");

	


	//run_case_fixation(stdout,0.9, 0.01, 0.5, 0.6, 1);

	//run_case_fixation(stdout, 0.9, 0.01, 0.5, 0.5, 1);

	//run_case_fixation(stdout, 0.05, 0.01, 0.6, 0.6, 1);

	//run_case_fixation(stdout, 0.3, 0.01, 0.7, 0.6, 1);
	//FILE* of = fopen("Data6.csv", "w");
	////max_strength_vs_sib_mat_f_feq(of);
	//max_strength_over_sibling_mating_precent(of);
	//run_case_fixation(stdout, 0.9, 0.01, 0.4, 0.39, 1);
	//run_case( 0.01, 0.01, 0.4, 0.5, 1);
	//run_case_equ(stdout, 0.01, 0.01, 0.39, 0.4, 1);
	/*FILE* of = fopen("Dif24.csv", "w");
	run_case_dif(of, 0.01, 0.01, 0.5, 0.5, 1);
	fclose(of);*/
	/*FILE* of = fopen("Data_equ_polymorphic2.csv", "w");

	print_equ_polymorphic(of, 0.95);

	fclose(of);*/

	//code for figure 2
	FILE* of = fopen("Datasb1.csv", "w");
	print_equ(of, 0.2);
	fclose(of);

	of = fopen("Datasb3.csv", "w");
	print_equ(of, 0.5);
	fclose(of);
	of = fopen("Datasb2.csv", "w");
	print_equ(of, 0.95);
	fclose(of);
	/*FILE* of = fopen("D1.csv", "w");
	run_case_Ds(of, 0.01, 0.01, 0.9, 0.8, 1);

	fclose(of);

	of = fopen("D2.csv", "w");
	run_case_Ds(of,0.01, 0.01, 0.7, 0.8, 1);
	
	


	of = fopen("D3.csv", "w");
	run_case_Ds(of, 0.01, 0.01, 0.1, 0.4, 1);

	fclose(of);

	of = fopen("D4.csv", "w");
	run_case_Ds(of, 0.01, 0.01, 0.35, 0.4, 1);

	fclose(of);

	of = fopen("D5.csv", "w");
	run_case_Ds(of, 0.01, 0.01, 0.7, 0.8, 1);

	fclose(of);*/


	//system("dir");
	//fclose(of);
	return 0;
	

}
