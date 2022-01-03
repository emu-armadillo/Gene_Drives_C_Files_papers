// ToUploadGraph_1.cpp : This file contains the 'main' function. Program execution begins and ends there.
//
#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#define MAX_GEN 300

double t[3][2] = { {1.0,0.0},
	{0.5, 0.5},
	{0.0, 1.0} };




typedef struct gen_data {
	bool fixed;
	double d1;
	double d2;
	double r;
	double wb;
	double data110;
	int gen;

} gen_data_t;

void print_data(const gen_data_t* const print, FILE* f) {
	/*if (print->fixed) {
		fprintf(f,"fixation: True\n");
	}
	else {
		fprintf(f,"fixation: False\n");
	}*/
	fprintf(f, "gen \tD1 \tD2 \tR \tWb\n");
	fprintf(f, "%d  \t%.3f\t%.3f\t%.3f\t%.3f\n", print->gen, print->d1, print->d2, print->r, print->wb);
}

void init_tans_maxtrics(double tm1[3][3][2], double tm2[3][3][2]) {

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 2; j++) {   // set transmission matrices
			tm1[0][i][j] = tm2[0][i][j] = t[i][j];
			tm1[1][i][j] = tm2[1][i][j] = t[i][j];
			tm1[2][i][j] = tm2[2][i][j] = t[i][j];
		}
	}
	/* above; tm1 is a transmission matrix for drive locus 1, tm2 for drive locus 2.  both in males,
	tm1[i][j][k]:  The first entry 0,1,2 gives the number of resistance alleles at locus 3 in the parental genotype.
	It will be envforced below that R is dominant, so drive only operates for 0
	i is the number of D alleles at the parental locus (0, 1, 2); 1 means a heterozygote, which is subject to drive
	j is the number of D alleles in the gamete

	They are all set here to follow Mendel; m1[0][1][] and m2[0][1][] will be reset below.

	*/

	tm1[0][1][1] = 1.0;
	tm1[0][1][0] = 0.0;
	tm2[0][1][1] = 1.0;
	tm2[0][1][0] = 0.0;	// Assume perfect drive

	// above:  m[i][j][k] the transmission matrices for both drive loci are modified so that R resistance allele is dominant,
	//and drive operates in the heterogygotes (when j=1); dominance is due to the fact that distortion only operates when i=0


}

void zero_gamete(double fg[2][2][2], double mg[2][2][2]) {
	for (int i = 0; i < 2; i++) {
		for (int j = 0; j < 2; j++) {
			for (int k = 0; k < 2; k++) {
				fg[i][j][k] = mg[i][j][k] = 0.0; // zero matrices
			}
		}
	}
}

void set_initial_diploid_vals(double dipmales[3][3][3], double dipfemales[3][3][3],
	double dr_1_feq, double dr_2_feq, double r_feq) {

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			for (int k = 0; k < 3; k++)
			{
				dipfemales[i][j][k] = 0.0;
				dipmales[i][j][k] = 0.0;

			}
		}
	}

	//this counteracts the halving of resisance frequency that occures as a result of making it heterozygous.
	r_feq *= 2;


	dipmales[1][0][0] = dipfemales[1][0][0] = dr_1_feq;
	dipmales[0][1][0] = dipfemales[0][1][0] = dr_2_feq;
	//making resistance homozygous necesary to prevent it from being halved in frequency
	dipmales[0][0][1] = dipfemales[0][0][1] =( r_feq);

	dipmales[0][0][0] = dipfemales[0][0][0] = 1.0 - dr_1_feq - dr_2_feq - (r_feq);
}

void set_homoz_fitness_values(double mfec1[3], double mfec2[3], double ffec1[3], double ffec2[3],
	double dr_1_m_fit, double dr_2_m_fit, double dr_1_f_fit, double dr_2_f_fit) {
	
	mfec1[0] = mfec2[0] = mfec1[1] = mfec2[1] = 1.0;
	
	mfec1[2] = dr_1_m_fit;
	mfec2[2] = dr_2_m_fit;


	ffec1[0] = ffec2[0] = ffec1[1] = ffec2[1] = 1.0;
	
	ffec1[2] = dr_1_f_fit;
	ffec2[2] = dr_2_f_fit;

}


gen_data_t run_case_untill_fixation_sequential(FILE* of, double d1_feq, double d2_feq, double r_feq, double dr_1_fit_m, double dr_2_fit_m, double dr_1_fit_f,
	double dr_2_fit_f)
{
	

	int dummy, dummy1, i, j, k, totcycles, printcycles, gen,
		p1, p2, p3,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	double
		tmpmg[2][2][2], mg[2][2][2], tm1[3][3][2], tm2[3][3][2], tmp,
		tmpfg[2][2][2], fg[2][2][2], dipmales[3][3][3], dipfemales[3][3][3],
		seg1, seg2, X, Y, tot,
		ffec1[3], ffec2[3], ffec3[3], ffectot,
		mfec1[3], mfec2[3], mfec3[3], mfectot,
		sum, x1, x2, x3, malpha, falpha, z, D1, D2, R, Wbar = 0;

	bool d2_introduced = false;

	for (i = 0; i < 3; i++)
	{
		//ffec1[i] = 0;
		//ffec2[i] = 0;
		//ffec3[i] = 0;
		//mfec1[i] = 0;
		//mfec2[i] = 0;
		//mfec3[i] = 0;
	}



	init_tans_maxtrics(tm1, tm2);


	set_initial_diploid_vals(dipmales, dipfemales, d1_feq, 0, r_feq);


	// above:  input and set initial conditions

	set_homoz_fitness_values(mfec1, mfec2, ffec1, ffec2, dr_1_fit_m, dr_2_fit_m, dr_1_fit_f, dr_2_fit_f);

	if (of != NULL) {

		fprintf(of, "mfec1[1] = %1.2f ,mfec1[2] = %1.2f ,mfec2[1] = %1.2f ,mfec2[2] = %1.2f \n",
			mfec1[1], mfec1[2], mfec2[1], mfec2[2]);

		fprintf(of, "ffec1[1] = %1.2f ,ffec1[2] = %1.2f ,ffec2[1] = %1.2f ,ffec2[2] = %1.2f \n",
			ffec1[1], ffec1[2], ffec2[1], ffec2[2]);
	}

	//target1:;
		//printf(" number of cycles before printing and number of times to print\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = MAX_GEN * 2; //MAX_GEN
	if (of != NULL) {
		fprintf(of, "\ninitial vals\n");

		fprintf(of, "D1 feq\t D2 feq \t R feq \t\tD1 mfit \t D1 ffit \t D2 mfit \t D2 ffit\n");
		fprintf(of, "%1.3f \t %1.3f  \t %1.3f \t\t%1.3f \t\t %1.3f  \t %1.3f   \t %1.3f \n",
			d1_feq, d2_feq, r_feq, dr_1_fit_m, dr_1_fit_f, dr_2_fit_m, dr_2_fit_f);
	}


	for (I = 1, gen = 0; I <= totcycles; I++) {
		for (J = 1; J <= printcycles; J++) {

			zero_gamete(fg, mg);
			Wbar = 0.0;

			// Guts start here  get gametes from diploids; fecundity has its effect in diploid to gametes


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {  // diploid genotypes to gametes

				X = dipfemales[p1][p2][p3];

				z = ffec1[p1] * ffec2[p2];   // female fitness step 1


				Wbar += z * X;

				// Now generate haploid eggs

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     fg[o1][o2][o3] += t[p1][o1] * t[p2][o2] * t[p3][o3] * z * X;


				// Now males, which have SD

				Y = dipmales[p1][p2][p3];

				z = mfec1[p1] * mfec2[p2];   // male fitness step 1


				// Now generate haploid sperm
				RR = 1 * (p3 > 0);  // RR = 0 means resistanc allele absent

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     mg[o1][o2][o3] +=
					tm1[RR][p1][o1] * tm2[RR][p2][o2] * t[p3][o3] * z * Y;



			}  // gametes done, can normalize in diploids

// make diploids

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++)     dipmales[p1][p2][p3] = dipfemales[p1][p2][p3] = 0.0;

			for (p1 = 0; p1 < 2; p1++) for (p2 = 0; p2 < 2; p2++) for (p3 = 0; p3 < 2; p3++)     for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     dipfemales[p1 + o1][p2 + o2][p3 + o3] += mg[p1][p2][p3] * fg[o1][o2][o3];


			for (p1 = 0, tot = 0.0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) 	 tot += dipfemales[p1][p2][p3];


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {
				dipfemales[p1][p2][p3] /= tot;
				dipmales[p1][p2][p3] = dipfemales[p1][p2][p3];     // normalizes
			}

			gen++;
		}  // end of inner printcycles

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += fg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) fg[i][j][k] /= tmp;

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += mg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) mg[i][j][k] /= tmp;

		//printf("tot = %1.3f\n",tot);

		//if (I == 1) fprintf(of, "gen \t000 \t001 \t010 \t100 \t011 \t101 \t110 \t111    \tD1 \t\tD2 \t\tR \t\tWb\n");
		if (I == 1 && of != NULL) {
			fprintf(of, "\n");
		}
		D1 = 0.5 * dipmales[1][0][0] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][0][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[1][2][0] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[1][2][2] + dipmales[2][0][0] + dipmales[2][0][1] + dipmales[2][0][2] + dipmales[2][1][0] + dipmales[2][1][1] + dipmales[2][1][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		D2 = 0.5 * dipmales[0][1][0] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][1][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[2][1][0] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][1][2] + dipmales[0][2][0] + dipmales[0][2][1] + dipmales[0][2][2] + dipmales[1][2][0] + dipmales[1][2][1] + dipmales[1][2][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		R = 0.5 * dipmales[0][0][1] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][2][1] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[2][0][1] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][2][1] + dipmales[0][0][2] + dipmales[0][1][2] + dipmales[0][2][2] + dipmales[1][0][2] + dipmales[1][1][2] + dipmales[1][2][2] + dipmales[2][0][2] + dipmales[2][1][2] + dipmales[2][2][2];

		if (!d2_introduced && D1 > 0.9995)
		{
			//TEMP_CHANGE
			dipmales[0][1][0] = d2_feq;
			dipfemales[0][1][0] = d2_feq;
			for (p1 = 0, tot = 0.0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) 	 tot += dipfemales[p1][p2][p3];


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {
				dipfemales[p1][p2][p3] /= tot;
				dipmales[p1][p2][p3] = dipfemales[p1][p2][p3];     // normalizes
			}
			d2_introduced = true;
		}

		//if fixation has occured print everything and exit
		//prints and exits if on the last loop and fixation has not occured.
		if (D1 > 0.9995 && D2 > 0.9995 || I == totcycles - 1) {
			

			gen_data_t to_r;

			to_r.fixed = (D1 > 0.9995 && D2 > 0.9995);
			to_r.d1 = D1;
			to_r.d2 = D2;
			to_r.r = R;
			to_r.wb = Wbar;
			to_r.gen = I * printcycles;


			if (of != NULL) {
				print_data(&to_r, of);
			}

			return to_r;
		}



	}  // end of outer print loops

}


//params
//printfile, drive 1 feq, drive 2 feq, resistance feq, drive 1 male fit, drive 2 male fit, drive 1 female fit, 
// drive 2 female fit
//the print file must already be open and must be closed outside the fuction
gen_data_t run_case_untill_fixation(FILE* of, double d1_feq, double d2_feq, double r_feq, double dr_1_fit_m, double dr_2_fit_m, double dr_1_fit_f,
	double dr_2_fit_f)
{
	//int dummy = 0, dummy1 = 0, i = 0, j = 0, 
	//	k = 0, totcycles = 0, printcycles = 0, gen = 0,
	//	p1 = 0, p2 = 0, p3 = 0,
	//	o1 = 0, o2 = 0, o3 = 0, I = 0, J = 0,
	//	i1 = 0, i2 = 0, i3 = 0, j1j2 = 0, j3 = 0, k1 = 0, k2 = 0, k3 = 0, RR = 0;

	//float
	//	tmpmg[2][2][2], mg[2][2][2], tm1[3][3][2], tm2[3][3][2], tmp,
	//	tmpfg[2][2][2], fg[2][2][2], dipmales[3][3][3], dipfemales[3][3][3],
	//	seg1 = 0, seg2 = 0, X = 0, Y = 0, tot = 0,
	//	ffec1[3], ffec2[3], ffec3[3], ffectot = 0,
	//	mfec1[3], mfec2[3], mfec3[3], mfectot = 0,
	//	sum = 0, x1 = 0, x2 = 0, x3 = 0, malpha = 0, 
	//	falpha = 0, z = 0, D1 = 0, D2 = 0, R = 0, Wbar = 0;

	int dummy, dummy1, i, j, k, totcycles, printcycles, gen,
		p1, p2, p3,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	double
		tmpmg[2][2][2], mg[2][2][2], tm1[3][3][2], tm2[3][3][2], tmp,
		tmpfg[2][2][2], fg[2][2][2], dipmales[3][3][3], dipfemales[3][3][3],
		seg1, seg2, X, Y, tot,
		ffec1[3], ffec2[3], ffec3[3], ffectot,
		mfec1[3], mfec2[3], mfec3[3], mfectot,
		sum, x1, x2, x3, malpha, falpha, z, D1, D2, R, Wbar = 0;


	for (i = 0; i < 3; i++)
	{
		//ffec1[i] = 0;
		//ffec2[i] = 0;
		//ffec3[i] = 0;
		//mfec1[i] = 0;
		//mfec2[i] = 0;
		//mfec3[i] = 0;
	}


	init_tans_maxtrics(tm1, tm2);


	set_initial_diploid_vals(dipmales, dipfemales, d1_feq, d2_feq, r_feq);


	// above:  input and set initial conditions

	set_homoz_fitness_values(mfec1, mfec2, ffec1, ffec2, dr_1_fit_m, dr_2_fit_m, dr_1_fit_f, dr_2_fit_f);

	if (of != NULL) {

		fprintf(of, "mfec1[1] = %1.2f ,mfec1[2] = %1.2f ,mfec2[1] = %1.2f ,mfec2[2] = %1.2f \n",
			mfec1[1], mfec1[2], mfec2[1], mfec2[2]);

		fprintf(of, "ffec1[1] = %1.2f ,ffec1[2] = %1.2f ,ffec2[1] = %1.2f ,ffec2[2] = %1.2f \n",
			ffec1[1], ffec1[2], ffec2[1], ffec2[2]);
	}

	//target1:;
		//printf(" number of cycles before printing and number of times to print\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = MAX_GEN; //MAX_GEN
	if (of != NULL) {
		fprintf(of, "\ninitial vals\n");

		fprintf(of, "D1 feq\t D2 feq \t R feq \t\tD1 mfit \t D1 ffit \t D2 mfit \t D2 ffit\n");
		fprintf(of, "%1.3f \t %1.3f  \t %1.3f \t\t%1.3f \t\t %1.3f  \t %1.3f   \t %1.3f \n",
			d1_feq, d2_feq, r_feq, dr_1_fit_m, dr_1_fit_f, dr_2_fit_m, dr_2_fit_f);
	}


	for (I = 1, gen = 0; I <= totcycles; I++) {
		for (J = 1; J <= printcycles; J++) {

			zero_gamete(fg, mg);
			Wbar = 0.0;

			// Guts start here  get gametes from diploids; fecundity has its effect in diploid to gametes


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {  // diploid genotypes to gametes

				X = dipfemales[p1][p2][p3];

				z = ffec1[p1] * ffec2[p2];   // female fitness step 1


				Wbar += z * X;

				// Now generate haploid eggs

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     fg[o1][o2][o3] += t[p1][o1] * t[p2][o2] * t[p3][o3] * z * X;


				// Now males, which have SD

				Y = dipmales[p1][p2][p3];

				z = mfec1[p1] * mfec2[p2];   // male fitness step 1


				// Now generate haploid sperm
				RR = 1 * (p3 > 0);  // RR = 0 means resistanc allele absent

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     mg[o1][o2][o3] +=
					tm1[RR][p1][o1] * tm2[RR][p2][o2] * t[p3][o3] * z * Y;



			}  // gametes done, can normalize in diploids

// make diploids

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++)     dipmales[p1][p2][p3] = dipfemales[p1][p2][p3] = 0.0;

			for (p1 = 0; p1 < 2; p1++) for (p2 = 0; p2 < 2; p2++) for (p3 = 0; p3 < 2; p3++)     for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     dipfemales[p1 + o1][p2 + o2][p3 + o3] += mg[p1][p2][p3] * fg[o1][o2][o3];


			for (p1 = 0, tot = 0.0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) 	 tot += dipfemales[p1][p2][p3];


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {
				dipfemales[p1][p2][p3] /= tot;
				dipmales[p1][p2][p3] = dipfemales[p1][p2][p3];     // normalizes
			}

			gen++;
		}  // end of inner printcycles

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += fg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) fg[i][j][k] /= tmp;

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += mg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) mg[i][j][k] /= tmp;

		//printf("tot = %1.3f\n",tot);

		//if (I == 1) fprintf(of, "gen \t000 \t001 \t010 \t100 \t011 \t101 \t110 \t111    \tD1 \t\tD2 \t\tR \t\tWb\n");
		if (I == 1 && of != NULL) {
			fprintf(of, "\n");
		}
		D1 = 0.5 * dipmales[1][0][0] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][0][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[1][2][0] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[1][2][2] + dipmales[2][0][0] + dipmales[2][0][1] + dipmales[2][0][2] + dipmales[2][1][0] + dipmales[2][1][1] + dipmales[2][1][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		D2 = 0.5 * dipmales[0][1][0] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][1][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[2][1][0] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][1][2] + dipmales[0][2][0] + dipmales[0][2][1] + dipmales[0][2][2] + dipmales[1][2][0] + dipmales[1][2][1] + dipmales[1][2][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		R = 0.5 * dipmales[0][0][1] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][2][1] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[2][0][1] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][2][1] + dipmales[0][0][2] + dipmales[0][1][2] + dipmales[0][2][2] + dipmales[1][0][2] + dipmales[1][1][2] + dipmales[1][2][2] + dipmales[2][0][2] + dipmales[2][1][2] + dipmales[2][2][2];


		//if fixation has occured print everything and exit
		//prints and exits if on the last loop and fixation has not occured.
		if (D1 > 0.9995 && D2 > 0.9995 || I == totcycles - 1) {
			/*fprintf(of, "%d\nfemale  %1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ", I * printcycles,
				fg[0][0][0], fg[0][0][1], fg[0][1][0], fg[1][0][0], fg[0][1][1], fg[1][0][1], fg[1][1][0], fg[1][1][1]);
			fprintf(of, "\nmale  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ",
				mg[0][0][0], mg[0][0][1], mg[0][1][0], mg[1][0][0], mg[0][1][1], mg[1][0][1], mg[1][1][0], mg[1][1][1]);*/



				/* for(i=0,tmp=0.0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++) tmp += mg[i][j][k];
				 for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++) mg[i][j][k] /= tmp;

				if(I==1) printf("gen \t000 \t001 \t010 \t100 \t011 \t101 \t110 \t111 \tD1 \tD2 \tR \tWb (male)\n");

				printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ",I*printcycles,
				mg[0][0][0],mg[0][0][1],mg[0][1][0],mg[1][0][0],mg[0][1][1],mg[1][0][1],mg[1][1][0],mg[1][1][1]);
				*/





				//fprintf(of, "\t%1.3f \t%1.3f \t%1.3f \t%1.3f\n", D1, D2, R, Wbar);

			gen_data_t to_r;

			to_r.fixed = (D1 > 0.9995 && D2 > 0.9995);
			to_r.d1 = D1;
			to_r.d2 = D2;
			to_r.r = R;
			to_r.wb = Wbar;
			to_r.gen = I * printcycles;


			if (of != NULL) {
				print_data(&to_r, of);
			}

			return to_r;
		}



	}  // end of outer print loops

}

//same as run_case_untill fixation accept for the flags
//flags are boolian values that get ored with the drives fixation state.
//for example, if drive 2 has a frequency of zero and you just care about drive 1 fixing you would 
//set d2_flag = true, and d1_flag = false.
gen_data_t run_case_untill_fixation_flags(FILE* of, double d1_feq, double d2_feq, double r_feq, double dr_1_fit_m, double dr_2_fit_m, double dr_1_fit_f,
	double dr_2_fit_f, bool d1_flag, bool d2_flag)
{
	//int dummy = 0, dummy1 = 0, i = 0, j = 0, 
	//	k = 0, totcycles = 0, printcycles = 0, gen = 0,
	//	p1 = 0, p2 = 0, p3 = 0,
	//	o1 = 0, o2 = 0, o3 = 0, I = 0, J = 0,
	//	i1 = 0, i2 = 0, i3 = 0, j1j2 = 0, j3 = 0, k1 = 0, k2 = 0, k3 = 0, RR = 0;

	//float
	//	tmpmg[2][2][2], mg[2][2][2], tm1[3][3][2], tm2[3][3][2], tmp,
	//	tmpfg[2][2][2], fg[2][2][2], dipmales[3][3][3], dipfemales[3][3][3],
	//	seg1 = 0, seg2 = 0, X = 0, Y = 0, tot = 0,
	//	ffec1[3], ffec2[3], ffec3[3], ffectot = 0,
	//	mfec1[3], mfec2[3], mfec3[3], mfectot = 0,
	//	sum = 0, x1 = 0, x2 = 0, x3 = 0, malpha = 0, 
	//	falpha = 0, z = 0, D1 = 0, D2 = 0, R = 0, Wbar = 0;

	int dummy, dummy1, i, j, k, totcycles, printcycles, gen,
		p1, p2, p3,
		o1, o2, o3, I, J,
		i1, i2, i3, j1j2, j3, k1, k2, k3, RR;

	double
		tmpmg[2][2][2], mg[2][2][2], tm1[3][3][2], tm2[3][3][2], tmp,
		tmpfg[2][2][2], fg[2][2][2], dipmales[3][3][3], dipfemales[3][3][3],
		seg1, seg2, X, Y, tot,
		ffec1[3], ffec2[3], ffec3[3], ffectot,
		mfec1[3], mfec2[3], mfec3[3], mfectot,
		sum, x1, x2, x3, malpha, falpha, z, D1, D2, R, Wbar = 0;


	for (i = 0; i < 3; i++)
	{
		//ffec1[i] = 0;
		//ffec2[i] = 0;
		//ffec3[i] = 0;
		//mfec1[i] = 0;
		//mfec2[i] = 0;
		//mfec3[i] = 0;
	}


	init_tans_maxtrics(tm1, tm2);


	set_initial_diploid_vals(dipmales, dipfemales, d1_feq, d2_feq, r_feq);


	// above:  input and set initial conditions

	set_homoz_fitness_values(mfec1, mfec2, ffec1, ffec2, dr_1_fit_m, dr_2_fit_m, dr_1_fit_f, dr_2_fit_f);


	if (of != NULL) {
		fprintf(of, "mfec1[1] = %1.2f ,mfec1[2] = %1.2f ,mfec2[1] = %1.2f ,mfec2[2] = %1.2f \n",
			mfec1[1], mfec1[2], mfec2[1], mfec2[2]);

		fprintf(of, "ffec1[1] = %1.2f ,ffec1[2] = %1.2f ,ffec2[1] = %1.2f ,ffec2[2] = %1.2f \n",
			ffec1[1], ffec1[2], ffec2[1], ffec2[2]);
	}


	//target1:;
		//printf(" number of cycles before printing and number of times to print\n");
	//	scanf("%d %d", &printcycles, &totcycles);
	printcycles = 1;
	totcycles = MAX_GEN; //MAX_GEN
	if (of != NULL) {
		fprintf(of, "\ninitial vals\n");

		fprintf(of, "D1 feq\t D2 feq \t R feq \t\tD1 mfit \t D1 ffit \t D2 mfit \t D2 ffit\n");
		fprintf(of, "%1.3f \t %1.3f  \t %1.3f \t\t%1.3f \t\t %1.3f  \t %1.3f   \t %1.3f \n",
			d1_feq, d2_feq, r_feq, dr_1_fit_m, dr_1_fit_f, dr_2_fit_m, dr_2_fit_f);
	}

	for (I = 1, gen = 0; I <= totcycles; I++) {
		for (J = 1; J <= printcycles; J++) {

			zero_gamete(fg, mg);
			Wbar = 0.0;

			// Guts start here  get gametes from diploids; fecundity has its effect in diploid to gametes


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {  // diploid genotypes to gametes

				X = dipfemales[p1][p2][p3];

				z = ffec1[p1] * ffec2[p2];   // female fitness step 1


				Wbar += z * X;

				// Now generate haploid eggs

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     fg[o1][o2][o3] += t[p1][o1] * t[p2][o2] * t[p3][o3] * z * X;


				// Now males, which have SD

				Y = dipmales[p1][p2][p3];

				z = mfec1[p1] * mfec2[p2];   // male fitness step 1


				// Now generate haploid sperm
				RR = 1 * (p3 > 0);  // RR = 0 means resistanc allele absent

				for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     mg[o1][o2][o3] +=
					tm1[RR][p1][o1] * tm2[RR][p2][o2] * t[p3][o3] * z * Y;



			}  // gametes done, can normalize in diploids

// make diploids

			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++)     dipmales[p1][p2][p3] = dipfemales[p1][p2][p3] = 0.0;

			for (p1 = 0; p1 < 2; p1++) for (p2 = 0; p2 < 2; p2++) for (p3 = 0; p3 < 2; p3++)     for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     dipfemales[p1 + o1][p2 + o2][p3 + o3] += mg[p1][p2][p3] * fg[o1][o2][o3];


			for (p1 = 0, tot = 0.0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) 	 tot += dipfemales[p1][p2][p3];


			for (p1 = 0; p1 < 3; p1++) for (p2 = 0; p2 < 3; p2++) for (p3 = 0; p3 < 3; p3++) {
				dipfemales[p1][p2][p3] /= tot;
				dipmales[p1][p2][p3] = dipfemales[p1][p2][p3];     // normalizes
			}

			gen++;
		}  // end of inner printcycles

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += fg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) fg[i][j][k] /= tmp;

		for (i = 0, tmp = 0.0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) tmp += mg[i][j][k];
		for (i = 0; i < 2; i++) for (j = 0; j < 2; j++) for (k = 0; k < 2; k++) mg[i][j][k] /= tmp;

		//printf("tot = %1.3f\n",tot);

		//if (I == 1) fprintf(of, "gen \t000 \t001 \t010 \t100 \t011 \t101 \t110 \t111    \tD1 \t\tD2 \t\tR \t\tWb\n");
		if (I == 1 && of != NULL) {
			fprintf(of, "\n");
		}
		D1 = 0.5 * dipmales[1][0][0] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][0][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[1][2][0] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[1][2][2] + dipmales[2][0][0] + dipmales[2][0][1] + dipmales[2][0][2] + dipmales[2][1][0] + dipmales[2][1][1] + dipmales[2][1][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		D2 = 0.5 * dipmales[0][1][0] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][1][2] + 0.5 * dipmales[1][1][0] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][1][2] + 0.5 * dipmales[2][1][0] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][1][2] + dipmales[0][2][0] + dipmales[0][2][1] + dipmales[0][2][2] + dipmales[1][2][0] + dipmales[1][2][1] + dipmales[1][2][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

		R = 0.5 * dipmales[0][0][1] + 0.5 * dipmales[0][1][1] + 0.5 * dipmales[0][2][1] + 0.5 * dipmales[1][0][1] + 0.5 * dipmales[1][1][1] + 0.5 * dipmales[1][2][1] + 0.5 * dipmales[2][0][1] + 0.5 * dipmales[2][1][1] + 0.5 * dipmales[2][2][1] + dipmales[0][0][2] + dipmales[0][1][2] + dipmales[0][2][2] + dipmales[1][0][2] + dipmales[1][1][2] + dipmales[1][2][2] + dipmales[2][0][2] + dipmales[2][1][2] + dipmales[2][2][2];


		//if fixation has occured print everything and exit
		//prints and exits if on the last loop and fixation has not occured.
		if ((D1 > 0.9995 || d1_flag) && (D2 > 0.9995 || d2_flag) || I == totcycles - 1) {
			/*fprintf(of, "%d\nfemale  %1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ", I * printcycles,
				fg[0][0][0], fg[0][0][1], fg[0][1][0], fg[1][0][0], fg[0][1][1], fg[1][0][1], fg[1][1][0], fg[1][1][1]);
			fprintf(of, "\nmale  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ",
				mg[0][0][0], mg[0][0][1], mg[0][1][0], mg[1][0][0], mg[0][1][1], mg[1][0][1], mg[1][1][0], mg[1][1][1]);*/



				/* for(i=0,tmp=0.0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++) tmp += mg[i][j][k];
				 for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++) mg[i][j][k] /= tmp;

				if(I==1) printf("gen \t000 \t001 \t010 \t100 \t011 \t101 \t110 \t111 \tD1 \tD2 \tR \tWb (male)\n");

				printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f ",I*printcycles,
				mg[0][0][0],mg[0][0][1],mg[0][1][0],mg[1][0][0],mg[0][1][1],mg[1][0][1],mg[1][1][0],mg[1][1][1]);
				*/





				//fprintf(of, "\t%1.3f \t%1.3f \t%1.3f \t%1.3f\n", D1, D2, R, Wbar);

			gen_data_t to_r;

			to_r.fixed = (D1 > 0.9995 && D2 > 0.9995);
			to_r.d1 = D1;
			to_r.d2 = D2;
			to_r.r = R;
			to_r.wb = Wbar;
			to_r.gen = I * printcycles;
			to_r.data110 = fg[1][1][0];
			if (of != NULL) {
				print_data(&to_r, of);
			}

			return to_r;
		}



	}  // end of outer print loops

}


void max_reduction_three_drives_graph_data(FILE* of) {
	fprintf(of, "R0,con,seq,one\n");
	for (float R0 = 0.005; R0 < 0.4; R0 += 0.01)
	{


		gen_data_t seq_max;
		seq_max.wb = 111;
		seq_max.r = -1;
		gen_data_t con_max;
		con_max.wb = 111;
		con_max.r = -1;
		gen_data_t oned_max;
		oned_max.wb = 111;
		oned_max.r = -1;
		float D1_max_s_con = 10;
		float D2_max_s_con = 10;
		float D1_max_s_seq = 10;
		float D2_max_s_seq = 10;
		for (float D1s = 1; D1s > 0.5; D1s -= 0.01)
		{
			for (float D2s = 1; D2s > 0.5; D2s -= 0.01)
			{
				gen_data_t seq_curr = run_case_untill_fixation_sequential(NULL, 0.01, 0.01, R0, D1s, D2s, D1s, D2s);
				gen_data_t con_curr = run_case_untill_fixation(NULL, 0.01, 0.01, R0, D1s, D2s, D1s, D2s);
				if (seq_curr.wb < seq_max.wb && (seq_curr.d1 > 0.9995 && seq_curr.d2 > 0.9995)) {
					seq_max = seq_curr;
					D2_max_s_seq = D2s;
					D1_max_s_seq = D1s;
				}
				if (con_curr.wb < con_max.wb && (con_curr.d1 > 0.9995 && con_curr.d2 > 0.9995)) {
					con_max = con_curr;
					D1_max_s_con = D1s;
					D2_max_s_con = D2s;
				}
			}
			gen_data_t one_d_curr = run_case_untill_fixation_flags(NULL, 0.02, 0, R0, D1s, 1, D1s, 1, 0, 1);
			if (one_d_curr.wb < oned_max.wb && (one_d_curr.d1 > 0.9995)) {
				oned_max = one_d_curr;
			}
		}
		//r,con,seq,one
		fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", R0, con_max.wb, seq_max.wb, oned_max.wb);
		printf("%.3f\n", R0);
	}
}

void max_reduction_three_drives_graph_data_f_only(FILE* of) {
	fprintf(of, "R0,con,seq,one\n");
	for (float R0 = 0.005; R0 < 0.4; R0 += 0.01)
	{


		gen_data_t seq_max;
		seq_max.wb = 111;
		seq_max.r = -1;
		gen_data_t con_max;
		con_max.wb = 111;
		con_max.r = -1;
		gen_data_t oned_max;
		oned_max.wb = 111;
		oned_max.r = -1;
		float D1_max_s_con = 10;
		float D2_max_s_con = 10;
		float D1_max_s_seq = 10;
		float D2_max_s_seq = 10;
		for (float D1s = 1; D1s > 0.4; D1s -= 0.01)
		{
			for (float D2s = 1; D2s > 0.4; D2s -= 0.01)
			{
				gen_data_t seq_curr = run_case_untill_fixation_sequential(NULL, 0.01, 0.01, R0, 1, 1, D1s, D2s);
				gen_data_t con_curr = run_case_untill_fixation(NULL, 0.01, 0.01, R0, 1, 1, D1s, D2s);
				if (seq_curr.wb < seq_max.wb && (seq_curr.d1 > 0.9995 && seq_curr.d2 > 0.9995)) {
					seq_max = seq_curr;
					D2_max_s_seq = D2s;
					D1_max_s_seq = D1s;
				}
				if (con_curr.wb < con_max.wb && (con_curr.d1 > 0.9995 && con_curr.d2 > 0.9995)) {
					con_max = con_curr;
					D1_max_s_con = D1s;
					D2_max_s_con = D2s;
				}
			}
			gen_data_t one_d_curr = run_case_untill_fixation_flags(NULL, 0.02, 0, R0, 1, 1, D1s, 1, 0, 1);
			if (one_d_curr.wb < oned_max.wb && (one_d_curr.d1 > 0.9995)) {
				oned_max = one_d_curr;
			}
		}
		//r,con,seq,one
		fprintf(of, "%.3f,%.3f,%.3f,%.3f\n", R0, con_max.wb, seq_max.wb, oned_max.wb);
		printf("%.3f\n", R0);
	}
}


int main(void)
{
	//system("dir");
	FILE* of = fopen("g1_data_both2.csv", "w");
	max_reduction_three_drives_graph_data(of);
	fclose(of);
	of = fopen("g1_data_f_only2.csv", "w");
	max_reduction_three_drives_graph_data_f_only(of);
	fclose(of);



	return 0;
}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
