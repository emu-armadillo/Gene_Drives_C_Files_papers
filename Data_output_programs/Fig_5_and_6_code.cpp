//rK_drive_DIsurvival.c to run recurrence equations of 1 drive (males) with a resistance locus where DD homozygotes may affect #births (female only, density indep), survival (both sexes, density indep) or K (density dep).  Differs from rK_drive in that this allows you to set an all-birth suppression (intervention) and to change it after any number of generations
// 15 VII 2021  
// 
// base code modified to create data for figures




#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <math.h>


// mom, dad, progeny geotypes with 1 locus, 2 alleles
// matrix gives kid genotypes per parents

float inherit[3][3][3] =
{ 1.0,0.0,0.0,           // mom aa, dad aa kids aa Aa AA
 0.5, 0.5,0.0,		// mom aa,     Aa,     aa Aa AA
 0.0, 1.0,0.0 ,		// mom aa,     AA,     aa Aa AA

 0.5,0.5,0.0,		// mom Aa, dad aa kids aa Aa AA
 0.25,0.5,0.25,		// mom Aa, dad Aa kids aa Aa AA
 0.0,0.5,0.5,		// mom Aa, dad AA kids aa Aa AA

 0.0,1.0,0.0,		// mom AA, dad aa kids aa Aa AA
 0.0,0.5,0.5,		// mom AA, dad Aa kids aa Aa AA
 0.0,0.0,1.0 },		// mom AA, dad AA kids aa Aa AA

	inheritd[3][3][3] = 	     // for male drive only, so only matters when dad is state 1 
{ 1.0,0.0,0.0,           // mom aa, dad aa kids aa Aa AA
 0.0, 1.0,0.0,		// mom aa,     Aa,     aa Aa AA
 0.0, 1.0,0.0 ,		// mom aa,     AA,     aa Aa AA

 0.5,0.5,0.0,		// mom Aa, dad aa kids aa Aa AA
 0.0,0.5,0.5,		// mom Aa, dad Aa kids aa Aa AA
 0.0,0.5,0.5,		// mom Aa, dad AA kids aa Aa AA

 0.0,1.0,0.0,		// mom AA, dad aa kids aa Aa AA
 0.0,0.0,1.0,		// mom AA, dad Aa kids aa Aa AA
 0.0,0.0,1.0 };		// mom AA, dad AA kids aa Aa AA

// inheritd is for drive locus, which operates only in males -- second entry of matrix




//explanation of the original code
/*printf("input:  \n\
births of d- moms, birth decrement to DD moms (s_f),
viability decrement to DD progeny (s_v), -- happens before the effect of K
carrying capacity of d-, decrement to carrying capacity of DD (s_K), and
inner loop and outer loop");
printf("births of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)");
printf("births are of males and females, so females are births/2");


*/

//params: 

// br - birth rate
// S - s_f birth rate decrement to DD moms
// viablity_in - s_v viability decrment to DD progeny.
// d_cc - K of d- indivials 
// DDcc - K of the drive homozygotes
// r_f - reistanse frequancy
// total_pop - total populaiton

void run_case(float br, float S,float viability_in,float d_cc,float DDcc, int d_f, int r_f, int total_pop) {

	int dummy, dummy1, i, j, k, totcycles, inloop, outloop, gen,
		p1, p2, p3, m1, m2, m3, f1, f2, f3, patch,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	extern float inherit[3][3][3], inheritd[3][3][3];
	double exp();
	// order in matedf is drive locus, R locus,  then same order for males

	float
		s_v, viability[3], matedf[3][3][3][3], progeny[3][3], m[3][3], f[3][3], surv[3], sumprogeny, intervention,
		fsurv1[3], fsurv2[3], fsurv3[3], totalmales, initdrive, initR, sum, sumdrive, sumR, survival,
		birth[3], K[3], totalprogeny, s_f, s_K;

	// Drive is first locus, resistance second


	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything
		m[i][j] = f[i][j] = progeny[i][j] = 0.0;
	}

	for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
	{
		matedf[f1][f2][m1][m2] = matedf[f1][f2][m1][m2] = 0.0;
	}



	// msurv1[0] = msurv3[0] = 1.0;  // all wt homozygote fit =  1.0
	// fsurv1[0]  = fsurv3[0] = 1.0;

	// msurv1[1] = fsurv1[1]= 1.0;  // drive hets normal

	printf("input:  \n\
\tbirths of d- moms, birth decrement to DD moms (s_f), \n\
\tviability decrement to DD progeny (s_v), -- happens before the effect of K\n\
\tcarrying capacity of d-, decrement to carrying capacity of DD (s_K), and \n\
\tinner loop and outer loop\n");
	printf("\nbirths of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)\n");
	printf("births are of males and females, so females are births/2\n");
	//scanf("%f %f %f %f %f %d %d", &, &, &, &, &, &, &);
	inloop = 1;
	outloop = 100;
	birth[0] = br;

	s_f = S;
	s_v = viability_in;
	K[0] = d_cc;
	s_K = DDcc;
	birth[2] = (1.0 - s_f) * birth[0];

	K[2] = (1.0 - s_K) * K[0];

	viability[0] = viability[1] = 1.0;
	viability[2] = 1.0 - s_v;
	// drive only changes things for DD genotypes

	birth[1] = birth[0];
	K[1] = K[0];

	gen = 0;

	// initial numbers of mated females; we may want to allow these to be specified by input at run time

	/*
	matedf[2][0][0][0] = 1000000.0;  // drive
	matedf[0][2][0][0] = 200000.0;   // resistance
	matedf[0][0][0][0] = 98800000.0;
	*/

	matedf[2][0][0][0] = d_f; //1000000.0;  // drive at 2%
	matedf[0][2][0][0] = r_f; //1000000.0;   // resistance at 2%
	matedf[0][0][0][0] = total_pop - d_f - r_f;//98000000.0;

//target:;
//
//	printf("input 'intervention' factor that will multiply all birth rates by the same amount (0< <1)\n");
//	scanf("%f", &);
	intervention = 1;
	printf("wt births (%2.1f)\t wt carrying cap (%1.2E) \tts_f (%0.3f) \ts_v (%0.3f) \ts_K (%0.3f)\n", birth[0], K[0], s_f, s_v, s_K);
	printf("gen \tdrive \tR \tinitial-prog \tDD-progeny \trel fitness DD\n");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (I = 1; I <= outloop; I++) {
		for (J = 1; J <= inloop; J++) {


			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything except matedf
				m[i][j] = f[i][j] = progeny[i][j] = 0.0;
			}

			// make progeny in the patch -- needs 9 loops
			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)
				for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)
					for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {

						// if m2=0, then male lacks resistance gene and drive can operate, but only if male is het; Resistance is dominant

						if (m2 == 0) progeny[p1][p2] += intervention * birth[f1] * inheritd[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  // drive operates at locus 1
						else if (m2 > 0) progeny[p1][p2] += intervention * birth[f1] * inherit[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  //full resistance to drive
							// no sex diffeences in zygotes; so far numbers, not frequencies
							// births are total offspring, half male, half female
							// if drive affects births, it is a female-only effect
					}



			// scale numbers back to carrying capacity by simple truncation

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= viability[p1];

			for (totalprogeny = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalprogeny += progeny[p1][p2];

			//surv[0] = exp(0.5*(1.0 - totalprogeny/K[0]));  -- old form of density dependence, not used now

			if (K[0] >= totalprogeny) surv[0] = 1.0;
			else surv[0] = (K[0] / totalprogeny);

			if (K[1] >= totalprogeny) surv[1] = 1.0;
			else surv[1] = (K[1] / totalprogeny);

			if (K[2] >= totalprogeny) surv[2] = 1.0;
			else surv[2] = (K[2] / totalprogeny);

			// you die according to population numbers and your specific K; for now, different carrying capacity only for DD homoz,
			// as if drive effect on carrying capacity is recessive; code is general if we wish to change it



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= surv[p1];


			for (sumprogeny = 0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    sumprogeny += progeny[p1][p2];
			//printf("surv[0] = %f surv[1] = %f surv[2] = %f total = %1.4E\n",surv[0],surv[1],surv[2],sumprogeny);



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] = f[p1][p2] = progeny[p1][p2] / 2.0;
			// now split into male and female progeny  since males will be converted for freqs, will keep females as numbers

			for (totalmales = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalmales += m[p1][p2];

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] /= totalmales;


			//printf("total males = %1.2E\n\n",totalmales);   // check 

			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {
				matedf[i][j][m1][m2] = 0.0;
			}


			for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
				matedf[f1][f2][m1][m2] = f[f1][f2] * m[m1][m2];

			// don't think I need += since there none of the matedf combinations will be revisited

			/////////////////////////////////////////////
			////////////////////////////////////////////
			///////////////////////////////////////////

			gen++;
		}
		/////////////////////////////////////////////
		////////////////////////////////////////////
		///////////////////////////////////////////
									// printing goes here
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		for (sumdrive = sumR = sum = 0.0, f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)

		{
			sum += matedf[f1][f2][m1][m2];
			sumdrive += ((f1 + m1) * (matedf[f1][f2][m1][m2]) / 4.0);
			sumR += ((f2 + m2) * (matedf[f1][f2][m1][m2]) / 4.0);
		}

		printf("%d \t%1.3f \t%1.4f\t%1.3E \t%1.3E \t%1.4f\n", gen, sumdrive / sum, sumR / sum, totalprogeny, sumprogeny, birth[2] * surv[2] / (birth[0] * surv[0]));
	}  // end of I outloop 


}

//explanation of the original code
/*printf("input:  \n\
births of d- moms, birth decrement to DD moms (s_f),
viability decrement to DD progeny (s_v), -- happens before the effect of K
carrying capacity of d-, decrement to carrying capacity of DD (s_K), and
inner loop and outer loop");
printf("births of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)");
printf("births are of males and females, so females are births/2");
*/

//params: 

// br - birth rate
// S - s_f birth rate decrement to DD moms
// viablity_in - s_v viability decrment to DD progeny.
// d_cc - K of d- indivials 
// DDcc - K of the drive homozygotes
// r_f - reistanse frequancy
// total_pop - total populaiton

float run_case_csv(float br, float S, float viability_in, float d_cc, float DDcc, int d_f, int r_f, int total_pop) {

	int dummy, dummy1, i, j, k, totcycles, inloop, outloop, gen,
		p1, p2, p3, m1, m2, m3, f1, f2, f3, patch,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	extern float inherit[3][3][3], inheritd[3][3][3];
	double exp();
	// order in matedf is drive locus, R locus,  then same order for males

	float
		s_v, viability[3], matedf[3][3][3][3], progeny[3][3], m[3][3], f[3][3], surv[3], sumprogeny, intervention,
		fsurv1[3], fsurv2[3], fsurv3[3], totalmales, initdrive, initR, sum, sumdrive, sumR, survival,
		birth[3], K[3], totalprogeny, s_f, s_K;

	// Drive is first locus, resistance second


	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything
		m[i][j] = f[i][j] = progeny[i][j] = 0.0;
	}

	for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
	{
		matedf[f1][f2][m1][m2] = matedf[f1][f2][m1][m2] = 0.0;
	}



	// msurv1[0] = msurv3[0] = 1.0;  // all wt homozygote fit =  1.0
	// fsurv1[0]  = fsurv3[0] = 1.0;

	// msurv1[1] = fsurv1[1]= 1.0;  // drive hets normal

	
	//scanf("%f %f %f %f %f %d %d", &, &, &, &, &, &, &);
	inloop = 1;
	outloop = 1000;
	birth[0] = br;

	s_f = S;
	s_v = viability_in;
	K[0] = d_cc;
	s_K = DDcc;
	birth[2] = (1.0 - s_f) * birth[0];

	K[2] = (1.0 - s_K) * K[0];

	viability[0] = viability[1] = 1.0;
	viability[2] = 1.0 - s_v;
	// drive only changes things for DD genotypes

	birth[1] = birth[0];
	K[1] = K[0];

	gen = 0;

	// initial numbers of mated females; we may want to allow these to be specified by input at run time

	/*
	matedf[2][0][0][0] = 1000000.0;  // drive
	matedf[0][2][0][0] = 200000.0;   // resistance
	matedf[0][0][0][0] = 98800000.0;
	*/

	matedf[2][0][0][0] = d_f; //1000000.0;  // drive at 2%
	matedf[0][2][0][0] = r_f; //1000000.0;   // resistance at 2%
	matedf[0][0][0][0] = total_pop - d_f - r_f;//98000000.0;

//target:;
//
//	printf("input 'intervention' factor that will multiply all birth rates by the same amount (0< <1)\n");
//	scanf("%f", &);
	intervention = 1;
	///*printf("wt births (%2.1f)\t wt carrying cap (%1.2E) \tts_f (%0.3f) \ts_v (%0.3f) \ts_K (%0.3f)\n", birth[0], K[0], s_f, s_v, s_K);
	//prin*/tf("gen \tdrive \tR \tinitial-prog \tDD-progeny \trel fitness DD\n");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (I = 1; I <= outloop; I++) {
		for (J = 1; J <= inloop; J++) {


			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything except matedf
				m[i][j] = f[i][j] = progeny[i][j] = 0.0;
			}

			// make progeny in the patch -- needs 9 loops
			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)
				for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)
					for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {

						// if m2=0, then male lacks resistance gene and drive can operate, but only if male is het; Resistance is dominant

						if (m2 == 0) progeny[p1][p2] += intervention * birth[f1] * inheritd[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  // drive operates at locus 1
						else if (m2 > 0) progeny[p1][p2] += intervention * birth[f1] * inherit[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  //full resistance to drive
							// no sex diffeences in zygotes; so far numbers, not frequencies
							// births are total offspring, half male, half female
							// if drive affects births, it is a female-only effect
					}



			// scale numbers back to carrying capacity by simple truncation

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= viability[p1];

			for (totalprogeny = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalprogeny += progeny[p1][p2];

			//surv[0] = exp(0.5*(1.0 - totalprogeny/K[0]));  -- old form of density dependence, not used now

			if (K[0] >= totalprogeny) surv[0] = 1.0;
			else surv[0] = (K[0] / totalprogeny);

			if (K[1] >= totalprogeny) surv[1] = 1.0;
			else surv[1] = (K[1] / totalprogeny);

			if (K[2] >= totalprogeny) surv[2] = 1.0;
			else surv[2] = (K[2] / totalprogeny);

			// you die according to population numbers and your specific K; for now, different carrying capacity only for DD homoz,
			// as if drive effect on carrying capacity is recessive; code is general if we wish to change it



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= surv[p1];


			for (sumprogeny = 0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    sumprogeny += progeny[p1][p2];
			//printf("surv[0] = %f surv[1] = %f surv[2] = %f total = %1.4E\n",surv[0],surv[1],surv[2],sumprogeny);



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] = f[p1][p2] = progeny[p1][p2] / 2.0;
			// now split into male and female progeny  since males will be converted for freqs, will keep females as numbers

			for (totalmales = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalmales += m[p1][p2];

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] /= totalmales;


			//printf("total males = %1.2E\n\n",totalmales);   // check 

			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {
				matedf[i][j][m1][m2] = 0.0;
			}


			for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
				matedf[f1][f2][m1][m2] = f[f1][f2] * m[m1][m2];

			// don't think I need += since there none of the matedf combinations will be revisited

			/////////////////////////////////////////////
			////////////////////////////////////////////
			///////////////////////////////////////////

			gen++;
		}
		/////////////////////////////////////////////
		////////////////////////////////////////////
		///////////////////////////////////////////
									// printing goes here
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		for (sumdrive = sumR = sum = 0.0, f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)

		{
			sum += matedf[f1][f2][m1][m2];
			sumdrive += ((f1 + m1) * (matedf[f1][f2][m1][m2]) / 4.0);
			sumR += ((f2 + m2) * (matedf[f1][f2][m1][m2]) / 4.0);
		}

	//	printf("%d \t%1.3f \t%1.4f\t%1.3E \t%1.3E \t%1.4f\n", gen, sumdrive / sum, sumR / sum, totalprogeny, sumprogeny, birth[2] * surv[2] / (birth[0] * surv[0]));
	}  // end of I outloop 
	return sumdrive / sum;

}




//explanation of the original code
/*printf("input:  \n\
births of d- moms, birth decrement to DD moms (s_f),
viability decrement to DD progeny (s_v), -- happens before the effect of K
carrying capacity of d-, decrement to carrying capacity of DD (s_K), and
inner loop and outer loop");
printf("births of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)");
printf("births are of males and females, so females are births/2");


*/

//params: 

// br - birth rate
// S - s_f birth rate decrement to DD moms
// viablity_in - s_v viability decrment to DD progeny.
// d_cc - K of d- indivials 
// DDcc - K of the drive homozygotes
// r_f - reistanse frequancy
// total_pop - total populaiton

float run_case_csv_knock_down(float br, float S, float viability_in, float d_cc, float DDcc, int d_f, int r_f, int total_pop) {

	int dummy, dummy1, i, j, k, totcycles, inloop, outloop, gen,
		p1, p2, p3, m1, m2, m3, f1, f2, f3, patch,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	extern float inherit[3][3][3], inheritd[3][3][3];
	double exp();
	// order in matedf is drive locus, R locus,  then same order for males

	float
		s_v, viability[3], matedf[3][3][3][3], progeny[3][3], m[3][3], f[3][3], surv[3], sumprogeny, intervention,
		fsurv1[3], fsurv2[3], fsurv3[3], totalmales, initdrive, initR, sum, sumdrive, sumR, survival,
		birth[3], K[3], totalprogeny, s_f, s_K;

	// Drive is first locus, resistance second


	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything
		m[i][j] = f[i][j] = progeny[i][j] = 0.0;
	}

	for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
	{
		matedf[f1][f2][m1][m2] = matedf[f1][f2][m1][m2] = 0.0;
	}



	// msurv1[0] = msurv3[0] = 1.0;  // all wt homozygote fit =  1.0
	// fsurv1[0]  = fsurv3[0] = 1.0;

	// msurv1[1] = fsurv1[1]= 1.0;  // drive hets normal


	//scanf("%f %f %f %f %f %d %d", &, &, &, &, &, &, &);
	inloop = 1;
	outloop = 1000;
	birth[0] = br;

	s_f = S;
	s_v = viability_in;
	K[0] = d_cc;
	s_K = DDcc;
	birth[2] = (1.0 - s_f) * birth[0];

	K[2] = (1.0 - s_K) * K[0];

	viability[0] = viability[1] = 1.0;
	viability[2] = 1.0 - s_v;
	// drive only changes things for DD genotypes

	birth[1] = birth[0];
	K[1] = K[0];

	gen = 0;

	// initial numbers of mated females; we may want to allow these to be specified by input at run time

	/*
	matedf[2][0][0][0] = 1000000.0;  // drive
	matedf[0][2][0][0] = 200000.0;   // resistance
	matedf[0][0][0][0] = 98800000.0;
	*/

	matedf[2][0][0][0] = d_f; //1000000.0;  // drive at 2%
	matedf[0][2][0][0] = r_f; //1000000.0;   // resistance at 2%
	matedf[0][0][0][0] = total_pop - d_f - r_f;//98000000.0;

//target:;
//
//	printf("input 'intervention' factor that will multiply all birth rates by the same amount (0< <1)\n");
//	scanf("%f", &);
	intervention = 1;
	///*printf("wt births (%2.1f)\t wt carrying cap (%1.2E) \tts_f (%0.3f) \ts_v (%0.3f) \ts_K (%0.3f)\n", birth[0], K[0], s_f, s_v, s_K);
	//prin*/tf("gen \tdrive \tR \tinitial-prog \tDD-progeny \trel fitness DD\n");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (I = 1; I <= outloop; I++) {
		for (J = 1; J <= inloop; J++) {


			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything except matedf
				m[i][j] = f[i][j] = progeny[i][j] = 0.0;
			}

			// make progeny in the patch -- needs 9 loops
			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)
				for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)
					for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {

						// if m2=0, then male lacks resistance gene and drive can operate, but only if male is het; Resistance is dominant

						if (m2 == 0) progeny[p1][p2] += intervention * birth[f1] * inheritd[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  // drive operates at locus 1
						else if (m2 > 0) progeny[p1][p2] += intervention * birth[f1] * inherit[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  //full resistance to drive
							// no sex diffeences in zygotes; so far numbers, not frequencies
							// births are total offspring, half male, half female
							// if drive affects births, it is a female-only effect
					}



			// scale numbers back to carrying capacity by simple truncation

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= viability[p1];

			for (totalprogeny = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalprogeny += progeny[p1][p2];

			//surv[0] = exp(0.5*(1.0 - totalprogeny/K[0]));  -- old form of density dependence, not used now

			if (K[0] >= totalprogeny) surv[0] = 1.0;
			else surv[0] = (K[0] / totalprogeny);

			if (K[1] >= totalprogeny) surv[1] = 1.0;
			else surv[1] = (K[1] / totalprogeny);

			if (K[2] >= totalprogeny) surv[2] = 1.0;
			else surv[2] = (K[2] / totalprogeny);

			// you die according to population numbers and your specific K; for now, different carrying capacity only for DD homoz,
			// as if drive effect on carrying capacity is recessive; code is general if we wish to change it



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= surv[p1];

			if (gen % 5 == 0) {

				for (i = 0; i < 3; i++)  for (j = 0; j < 3; j++)

					progeny[i][j] /= 4.0;

			}


			for (sumprogeny = 0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    sumprogeny += progeny[p1][p2];
			//printf("surv[0] = %f surv[1] = %f surv[2] = %f total = %1.4E\n",surv[0],surv[1],surv[2],sumprogeny);



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] = f[p1][p2] = progeny[p1][p2] / 2.0;
			// now split into male and female progeny  since males will be converted for freqs, will keep females as numbers

			for (totalmales = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalmales += m[p1][p2];

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] /= totalmales;


			//printf("total males = %1.2E\n\n",totalmales);   // check 

			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {
				matedf[i][j][m1][m2] = 0.0;
			}


			for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
				matedf[f1][f2][m1][m2] = f[f1][f2] * m[m1][m2];

			// don't think I need += since there none of the matedf combinations will be revisited

			/////////////////////////////////////////////
			////////////////////////////////////////////
			///////////////////////////////////////////

			gen++;
		}
		/////////////////////////////////////////////
		////////////////////////////////////////////
		///////////////////////////////////////////
									// printing goes here
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		for (sumdrive = sumR = sum = 0.0, f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)

		{
			sum += matedf[f1][f2][m1][m2];
			sumdrive += ((f1 + m1) * (matedf[f1][f2][m1][m2]) / 4.0);
			sumR += ((f2 + m2) * (matedf[f1][f2][m1][m2]) / 4.0);
		}

		//	printf("%d \t%1.3f \t%1.4f\t%1.3E \t%1.3E \t%1.4f\n", gen, sumdrive / sum, sumR / sum, totalprogeny, sumprogeny, birth[2] * surv[2] / (birth[0] * surv[0]));
	}  // end of I outloop 
	return sumdrive / sum;

}

//explanation of the original code
/*printf("input:  \n\
births of d- moms, birth decrement to DD moms (s_f), 
viability decrement to DD progeny (s_v), -- happens before the effect of K
carrying capacity of d-, decrement to carrying capacity of DD (s_K), and 
inner loop and outer loop");
printf("births of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)");
printf("births are of males and females, so females are births/2");


*/

//params: 
// of - outfile
// br - birth rate
// S - s_f birth rate decrement to DD moms
// viablity_in - s_v viability decrment to DD progeny.
// d_cc - K of d- indivials 
// DDcc - K of the drive homozygotes
// r_f - reistanse frequancy
// total_pop - total populaiton
void run_case_csv_gen_data(FILE * of,float br, float S, float viability_in, float d_cc, float DDcc, int d_f, int r_f, int total_pop) {

	int dummy, dummy1, i, j, k, totcycles, inloop, outloop, gen,
		p1, p2, p3, m1, m2, m3, f1, f2, f3, patch,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	extern float inherit[3][3][3], inheritd[3][3][3];
	double exp();
	// order in matedf is drive locus, R locus,  then same order for males

	float
		s_v, viability[3], matedf[3][3][3][3], progeny[3][3], m[3][3], f[3][3], surv[3], sumprogeny, intervention,
		fsurv1[3], fsurv2[3], fsurv3[3], totalmales, initdrive, initR, sum, sumdrive, sumR, survival,
		birth[3], K[3], totalprogeny, s_f, s_K;

	// Drive is first locus, resistance second


	for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything
		m[i][j] = f[i][j] = progeny[i][j] = 0.0;
	}

	for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
	{
		matedf[f1][f2][m1][m2] = matedf[f1][f2][m1][m2] = 0.0;
	}



	// msurv1[0] = msurv3[0] = 1.0;  // all wt homozygote fit =  1.0
	// fsurv1[0]  = fsurv3[0] = 1.0;

	// msurv1[1] = fsurv1[1]= 1.0;  // drive hets normal


	//scanf("%f %f %f %f %f %d %d", &, &, &, &, &, &, &);
	inloop = 1;
	outloop = 100;
	birth[0] = br;

	s_f = S;
	s_v = viability_in;
	K[0] = d_cc;
	s_K = DDcc;
	birth[2] = (1.0 - s_f) * birth[0];

	K[2] = (1.0 - s_K) * K[0];

	viability[0] = viability[1] = 1.0;
	viability[2] = 1.0 - s_v;
	// drive only changes things for DD genotypes

	birth[1] = birth[0];
	K[1] = K[0];

	gen = 0;

	fprintf(of, "gen,drive,R,pupoluation\n");

	// initial numbers of mated females; we may want to allow these to be specified by input at run time

	/*
	matedf[2][0][0][0] = 1000000.0;  // drive
	matedf[0][2][0][0] = 200000.0;   // resistance
	matedf[0][0][0][0] = 98800000.0;
	*/

	matedf[2][0][0][0] = d_f; //1000000.0;  // drive at 2%
	matedf[0][2][0][0] = r_f; //1000000.0;   // resistance at 2%
	matedf[0][0][0][0] = total_pop - d_f - r_f;//98000000.0;

//target:;
//
//	printf("input 'intervention' factor that will multiply all birth rates by the same amount (0< <1)\n");
//	scanf("%f", &);
	intervention = 1;
	///*printf("wt births (%2.1f)\t wt carrying cap (%1.2E) \tts_f (%0.3f) \ts_v (%0.3f) \ts_K (%0.3f)\n", birth[0], K[0], s_f, s_v, s_K);
	//prin*/tf("gen \tdrive \tR \tinitial-prog \tDD-progeny \trel fitness DD\n");

	//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	for (I = 1; I <= outloop; I++) {
		for (J = 1; J <= inloop; J++) {


			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++) {     // zero everything except matedf
				m[i][j] = f[i][j] = progeny[i][j] = 0.0;
			}

			// make progeny in the patch -- needs 9 loops
			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)
				for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)
					for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {

						// if m2=0, then male lacks resistance gene and drive can operate, but only if male is het; Resistance is dominant

						if (m2 == 0) progeny[p1][p2] += intervention * birth[f1] * inheritd[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  // drive operates at locus 1
						else if (m2 > 0) progeny[p1][p2] += intervention * birth[f1] * inherit[f1][m1][p1] * inherit[f2][m2][p2] * matedf[f1][f2][m1][m2];  //full resistance to drive
							// no sex diffeences in zygotes; so far numbers, not frequencies
							// births are total offspring, half male, half female
							// if drive affects births, it is a female-only effect
					}



			// scale numbers back to carrying capacity by simple truncation

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= viability[p1];

			for (totalprogeny = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalprogeny += progeny[p1][p2];

			//surv[0] = exp(0.5*(1.0 - totalprogeny/K[0]));  -- old form of density dependence, not used now

			if (K[0] >= totalprogeny) surv[0] = 1.0;
			else surv[0] = (K[0] / totalprogeny);

			if (K[1] >= totalprogeny) surv[1] = 1.0;
			else surv[1] = (K[1] / totalprogeny);

			if (K[2] >= totalprogeny) surv[2] = 1.0;
			else surv[2] = (K[2] / totalprogeny);

			// you die according to population numbers and your specific K; for now, different carrying capacity only for DD homoz,
			// as if drive effect on carrying capacity is recessive; code is general if we wish to change it



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    progeny[p1][p2] *= surv[p1];


			for (sumprogeny = 0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)    sumprogeny += progeny[p1][p2];
			//printf("surv[0] = %f surv[1] = %f surv[2] = %f total = %1.4E\n",surv[0],surv[1],surv[2],sumprogeny);



			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] = f[p1][p2] = progeny[p1][p2] / 2.0;
			// now split into male and female progeny  since males will be converted for freqs, will keep females as numbers

			for (totalmales = 0.0, p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   totalmales += m[p1][p2];

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++)   m[p1][p2] /= totalmales;


			//printf("total males = %1.2E\n\n",totalmales);   // check 

			for (i = 0; i < 3; i++) for (j = 0; j < 3; j++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++) {
				matedf[i][j][m1][m2] = 0.0;
			}


			for (f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)
				matedf[f1][f2][m1][m2] = f[f1][f2] * m[m1][m2];

			// don't think I need += since there none of the matedf combinations will be revisited

			/////////////////////////////////////////////
			////////////////////////////////////////////
			///////////////////////////////////////////

			gen++;
		}
		/////////////////////////////////////////////
		////////////////////////////////////////////
		///////////////////////////////////////////
									// printing goes here
		// %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


		for (sumdrive = sumR = sum = 0.0, f1 = 0; f1 < 3; f1++) for (f2 = 0; f2 < 3; f2++)   for (m1 = 0; m1 < 3; m1++) for (m2 = 0; m2 < 3; m2++)

		{
			sum += matedf[f1][f2][m1][m2];
			sumdrive += ((f1 + m1) * (matedf[f1][f2][m1][m2]) / 4.0);
			sumR += ((f2 + m2) * (matedf[f1][f2][m1][m2]) / 4.0);
		}
		fprintf(of,"%d, %.3f, %.3f, %d\n", gen, sumdrive / sum, sumR / sum, ((int)sumprogeny));
			printf("%d \t%1.3f \t%1.4f\t%1.3E \t%1.3E \t%1.4f\n", gen, sumdrive / sum, sumR / sum, totalprogeny, sumprogeny, birth[2] * surv[2] / (birth[0] * surv[0]));
	}  // end of I outloop 
	

}


//N is the birth rate during which the simulation is run.
void knockdown_n(FILE* of, int n) {

	for (float i = 0; i <= 1; i += 0.01)
	{

		float temp = run_case_csv_knock_down(n, 0, 0, 1000000000, i, 1000000.0, 1000000.0, 1000000000);
		if (temp < 0.9995) {
			//	float temp2 = run_case_csv(10, 0, i, 1000000000, 0, 1000000.0, 1000000.0, 1000000000);
			fprintf(of, "%s%d,%.2f, ", "knockdown", n, i);
			break;
		}
	}
	for (float i = 0; i <= 1; i += 0.01)
	{

		float temp2 = run_case_csv_knock_down(n, 0, i, 1000000000, 0, 1000000.0, 1000000.0, 1000000000);
		if (temp2 < 0.9995) {
			fprintf(of, " %.2f\n", i);
			break;
		}

	}
}



//explanation of the original code
/*printf("input:  \n\
births of d- moms, birth decrement to DD moms (s_f),
viability decrement to DD progeny (s_v), -- happens before the effect of K
carrying capacity of d-, decrement to carrying capacity of DD (s_K), and
inner loop and outer loop");
printf("births of DD will be (1-s_f), viability of DD (1-s_v), carrying capacity of DD will be (1-s_K)");
printf("births are of males and females, so females are births/2");


*/

//params: 
// br - birth rate
// S - s_f birth rate decrement to DD moms
// viablity_in - s_v viability decrment to DD progeny.
// d_cc - K of d- indivials 
// DDcc - K of the drive homozygotes
// r_f - reistanse frequancy
// total_pop - total populaiton

void no_knockdown_n(FILE* of, int n) {

	for (float i = 0; i <= 1; i += 0.01)
	{

		float temp = run_case_csv(n, 0, 0, 1000000000, i, 1000000.0, 1000000.0, 1000000000);
		if (temp < 0.9995) {
			//	float temp2 = run_case_csv(10, 0, i, 1000000000, 0, 1000000.0, 1000000.0, 1000000000);
			fprintf(of, "%d,%.2f, ", n, i);
			break;
		}
	}
	for (float i = 0; i <= 1; i += 0.01)
	{

		float temp2 = run_case_csv(n, 0, i, 1000000000, 0, 1000000.0, 1000000.0, 1000000000);
		if (temp2 < 0.9995) {
			fprintf(of, " %.2f\n", i);
			break;
		}

	}
}


//
// the permater s_v is used in the program
// s_f is not used however it is included for the case where the fitness effect is
// in females only.
int main()
{
	//This code was used to generate the data for figure 6.
	//uncomment and run to get data files
	//
	//**************************************************
	//system("dir");
	//FILE* of = fopen("DataGen_1.csv", "w");
	//run_case_csv_gen_data(of, 3, 0, 0, 1000000000, 0.4, 1000000.0, 1000000.0, 1000000000);
	//
	//fclose(of);

	//of = fopen("DataGen_2.csv", "w");
	//run_case_csv_gen_data(of, 3, 0, 0.4, 1000000000, 0, 1000000.0, 1000000.0, 1000000000);
	//fclose(of);
	//***********************************************************


	//
	//This produces two outpute files one 
	//One has data file without knockdown, the other data file with knockdown
	//
	//
	FILE* of = fopen("Data3.csv", "w");
	fprintf(of, "status, k, sv \n");

	for (int j = 2; j < 50; j++)
	{


		knockdown_n(of, j);
		printf("%d\n", j);
	}

	fclose(of);

	of = fopen("no_knockdown_data.csv", "w");

//	fprintf(of, "status, k, sv \n");
	for (int j = 2; j < 50; j++)
	{


		no_knockdown_n(of, j);
	}

	fclose(of);
	return 0;
}
