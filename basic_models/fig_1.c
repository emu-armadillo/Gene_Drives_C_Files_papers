/*
fc2drives_het_loop_init.c to run recurrence equations of 2 independent drives with a resistance locus
3 loci with males and females; drive in male only  3V21
assigns drive hets a fitness decrement and loops through different levels of drive imperfection
loops through to find maximum s
this program varies initial drive freq (x3)
both drives are introduced simultaneously at equal frequencies

This program generates outputs that cover more conditions than given in Fig. 1 (delta >0 is imperfect drive, but the figure assumed
delta =0).  Likewise, the program allows specifying a heterozygote fitness effect of the drive, but the figure assumed no heterozygous 
fitness effect.  The output does not span the density of initial R values used in Fig. 1, but that is trivialy changed.  
To generate results for a single drive, set x2 or x1 = 0 (search on **** in the code);
to generate results for female-only fitness effects, set mfec[2] and mfec[1] = 1.0 (search on **** in the code);
*/

#include <stdio.h>
#include <math.h>

// transmission matrix m[i,j] for a diploid parental genotype i (0,1,2 copies of the A allele ) and offspring gamete type j (0,1 copies of the A allele) ; matrix entries give the fraction of gametes for each parent-gamete combination.  m[0][0] is 1.0 because the aa parent produces only a gametes

float t[3][2] = {1.0,0.0,
		 0.5, 0.5,
		0.0, 1.0 };

int main()
{

int dummy,dummy1, i,j,k,totcycles,printcycles,gen,
p1,p2,p3,recessive,
o1, o2, o3,I,J,
i1,i2,i3,j1j2,j3,k1,k2,k3,RR;

extern float t[3][2];

float 
tmpmg[2][2][2],mg[2][2][2],  tm1[3][3][2], tm2[3][3][2], tmp,
tmpfg[2][2][2],fg[2][2][2],  dipmales[3][3][3],dipfemales[3][3][3],
seg1,seg2,X,Y,tot,delta,s,h,
ffec1[3],ffec2[3],ffec3[3],ffectot,
mfec1[3],mfec2[3],mfec3[3],mfectot,
sum, x1,x2,x3,malpha,falpha,z,D1,D2,R,Wbar;

double pow();


for(i=0;i<3;i++) for(j=0;j<2;j++) {   // set transmission matrices
tm1[0][i][j] = tm2[0][i][j] = t[i][j];
tm1[1][i][j] = tm2[1][i][j] = t[i][j];
tm1[2][i][j] = tm2[2][i][j] = t[i][j];
				}
/* above; tm1 is a transmission matrix for drive locus 1, tm2 for drive locus 2.  both in males,
tm1[i][j][k]:  The first entry 0,1,2 gives the number of resistance alleles at locus 3 in the parental genotype.  
It will be envforced below that R is dominant, so drive only operates for 0
i is the number of D alleles at the parental locus (0, 1, 2); 1 means a heterozygote, which is subject to drive
j is the number of D alleles in the gamete

They are all set here to follow Mendel; m1[0][1][] and m2[0][1][] will be reset below.

// above:  tm[i][j][k] the transmission matrices for both drive loci are modified so that R resistance allele is dominant, and drive operates in the heterogygotes (when j=1); dominance is due to the fact that distortion only operates when i=0

*/
mfec1[0] = mfec2[0] = mfec1[1] = mfec2[1]  = 1.0;
ffec1[0] = ffec2[0]  = ffec1[1] = ffec2[1]  = 1.0;

x1 =  0.001;
x2  = 0.001;  // set either of these to 0.0 to output for a single drive   ****

recessive = 0;
printf("input fitness decrement of drive heterozygotes\n");
scanf("%f",&h);
printf("fc2drives_het_loop.c   h = %0.3f\n",h);

printf("\n\t In the printout below, s= indicates the smallest s that blocked drive evolution; \n\
the number in parens is 1-(s-0.005)(s-0.005), hence the largest combined suppression for which drive evolved. \n\
The printout must be changed to s-0.005 when a single drive is modeled\n\n");

printf("Initial conditions (R0, delta) \t\t Final states (D1, D2, R) serve as a check that resistance evolved and drive did not \n\t\
 \t\t\t\t\t(I=2991 also indicates drive did not evolve) \n\n");

for(x3= 0.03;x3<0.80;x3+= 0.05)	{     // #1


	for (delta= 0.0; delta<0.16; delta += 0.05)	{  // #2

tm1[0][1][1] = 1.0 - delta;
tm1[0][1][0] = 0.0 + delta;
tm2[0][1][1] = 1.0 - delta;
tm2[0][1][0] = 0.0 + delta;	// perfect drive if delta = 0

// above:  tm[i][j][k] the transmission matrices for both drive loci are modified so that R resistance allele is dominant, and drive operates in the heterogygotes (when j=1); dominance is due to the fact that distortion only operates when i=0


//printf("\t%1.3f \t%1.3f \t%1.3f \tdelta = %0.4f\n",D1,D2,R,delta);


// above:  input and set initial conditions

		for(s=0.0;s<0.5;s+= 0.005)		{   // #3  s is the fitness effect in drive homozygotes of a single drive

mfec1[2] = mfec2[2] = 1.0 -s;  // set to 1.0 for a drive with effects just in females  ****
mfec1[1] = mfec2[1] = 1.0 -h;  // set to 1.0 for a drive with effects just in females  ****

ffec1[2] = ffec2[2] = 1.0 -s;
ffec1[1] = ffec2[1] = 1.0 -h;


for(i=0;i<3;i++) for(j=0;j<3;j++) for(k=0;k<3;k++) dipmales[i][j][k] = dipfemales[i][j][k] = 0.0;

dipmales[1][0][0] = dipfemales[1][0][0] = x1;  // initial frequency of drive 1 is x1/2 as these are heterozygotes
dipmales[0][1][0] = dipfemales[0][1][0] = x2;  // initial frequency of drive 2 is x2/2 as these are heterozygotes
dipmales[0][0][1] = dipfemales[0][0][1] = x3;  // initial frequency of resistance is x3/2 as these are heterozygotes
dipmales[0][0][0] = dipfemales[0][0][0] = 1.0 -x1 -x2 -x3;

gen = 0;

target1:;

for(I=1;I<=2990;I++)	{

	for(i=0;i<2;i++) for(j=0;j<2;j++) for(k=0;k<2;k++) fg[i][j][k] = mg[i][j][k] = 0.0; // zero matrices
Wbar = 0.0;

// Guts start here  get gametes from diploids; fecundity has its effect in diploid to gametes


for(p1=0;p1<3;p1++) for(p2=0;p2<3;p2++) for(p3=0;p3<3;p3++) 	{  // diploid genotypes to gametes

X  = dipfemales[p1][p2][p3];

// if (recessive ==0)  

z = ffec1[p1]*ffec2[p2];   // female fitness step 1 -- in AA-- --BB or AABB

/*
if (recessive ==1) 	{
z = 1.0;
if (p1 == 2 && p2 ==2) z = ffec1[p1]*ffec2[p2];   // female fitness step 1 -- only in AA BB
			}
*/


Wbar += z*X;

// Now generate haploid eggs

for(o1=0;o1<2;o1++) for(o2=0;o2<2;o2++) for(o3=0;o3<2;o3++)     fg[o1][o2][o3] += t[p1][o1] * t[p2][o2] * t[p3][o3]*z*X;


// Now males, which have SD

Y  = dipmales[p1][p2][p3];

// if (recessive ==0)   
 z = mfec1[p1]*mfec2[p2];   // male fitness step 1

/*
if (recessive ==1)	{
z = 1.0;
if (p1 ==2 && p2 ==2) z = mfec1[p1]*mfec2[p2];
			}
*/

// Now generate haploid sperm
RR = 1*(p3>0);  // RR = 0 means resistanc allele absent

for(o1=0;o1<2;o1++) for(o2=0;o2<2;o2++) for(o3=0;o3<2;o3++)     mg[o1][o2][o3] += 
			tm1[RR][p1][o1]*tm2[RR][p2][o2]* t[p3][o3] * z*Y;



								}  // gametes done, can normalize in diploids

// make diploids

for(p1=0;p1<3;p1++) for(p2=0;p2<3;p2++) for(p3=0;p3<3;p3++)     dipmales[p1][p2][p3] = dipfemales[p1][p2][p3] = 0.0;

for(p1=0;p1<2;p1++) for(p2=0;p2<2;p2++) for(p3=0;p3<2;p3++)     for(o1=0;o1<2;o1++) for(o2=0;o2<2;o2++) for(o3=0;o3<2;o3++)     dipfemales[p1+o1][p2+o2][p3+o3] += mg[p1][p2][p3]*fg[o1][o2][o3];


for(p1=0,tot=0.0;p1<3;p1++) for(p2=0;p2<3;p2++) for(p3=0;p3<3;p3++) 	 tot += dipfemales[p1][p2][p3];
								

for(p1=0;p1<3;p1++) for(p2=0;p2<3;p2++) for(p3=0;p3<3;p3++) 	{
dipfemales[p1][p2][p3] /= tot;
dipmales[p1][p2][p3] = dipfemales[p1][p2][p3];     // normalizes
								}
									
gen++;

D1 = 0.5*dipmales[1][0][0] + 0.5*dipmales[1][0][1] + 0.5*dipmales[1][0][2] + 0.5*dipmales[1][1][0] + 0.5*dipmales[1][1][1] + 0.5*dipmales[1][1][2] + 0.5*dipmales[1][2][0] + 0.5*dipmales[1][2][1] + 0.5*dipmales[1][2][2] + dipmales[2][0][0] + dipmales[2][0][1] + dipmales[2][0][2] + dipmales[2][1][0] + dipmales[2][1][1] + dipmales[2][1][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

D2 = 0.5*dipmales[0][1][0] + 0.5*dipmales[0][1][1] + 0.5*dipmales[0][1][2] + 0.5*dipmales[1][1][0] + 0.5*dipmales[1][1][1] + 0.5*dipmales[1][1][2] + 0.5*dipmales[2][1][0] + 0.5*dipmales[2][1][1] + 0.5*dipmales[2][1][2] + dipmales[0][2][0] + dipmales[0][2][1] + dipmales[0][2][2] + dipmales[1][2][0] + dipmales[1][2][1] + dipmales[1][2][2] + dipmales[2][2][0] + dipmales[2][2][1] + dipmales[2][2][2];

if (D1 > 0.999 || D2 > 0.999) I = 3000;  // one or both drives evolved; move to next s
// want first s for which resistance evolved, meaning I=2991
									}  // end of I loop

// if get to here and I = 2991, we know that drives got shut down

if (I==2991)	{


R= 0.5*dipmales[0][0][1] + 0.5*dipmales[0][1][1] + 0.5*dipmales[0][2][1] + 0.5*dipmales[1][0][1] + 0.5*dipmales[1][1][1] + 0.5*dipmales[1][2][1] + 0.5*dipmales[2][0][1] + 0.5*dipmales[2][1][1] + 0.5*dipmales[2][2][1] + dipmales[0][0][2] + dipmales[0][1][2] + dipmales[0][2][2] + dipmales[1][0][2] + dipmales[1][1][2] + dipmales[1][2][2] + dipmales[2][0][2] + dipmales[2][1][2] + dipmales[2][2][2];


//printf("I = %d \tR0=%0.3f \tdelta=%0.3f \tD1=%1.3f \tD2=%1.3f \tR=%1.3f \ts=%0.3f (%0.3f)\n",I,x3/2.0,delta,D1,D2,R,s,(1.0-(1.0-s + 0.005)*(1.0-s + 0.005)));

printf("R0=%0.3f \tdelta=%0.3f \ts=%0.3f (%0.3f)   \tD1=%1.3f \tD2=%1.3f \tR=%1.3f \tI=%d\n",x3/2.0,delta,s,(1.0-(1.0-s + 0.005)*(1.0-s + 0.005)),D1,D2,R,I);

goto target2;
		}

				}  // end of s loop

target2:;
					} // end of delta loop
printf("\n");

}  // end of x3 (R) loop

}
