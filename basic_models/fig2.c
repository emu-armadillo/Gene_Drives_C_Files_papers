/*
formerly sib_recessivedrive_partdrive.c 
See associated README file to understand how this code generates recursion
to run recurrence equations of sib mating with recessive partially lethal drive  28 xi 21
 Haploids with 2 loci
this program automatically includes perfect and imperfect drive (alpha = 0 is perfect drive)
and includes a fitness cost or no fitness cost assigned to drive heterozygotes (h = 0 is no cost)
The output is not quite as dense as used in Fig. 2
Output is sent to a file with the name 'out'; what gets printed to the screen is the progress
*/

#include <stdio.h>
#include <math.h>
FILE *fp, *fopen();

float transm[2][3] = {1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
                      0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles

int main()
{

int dummy,dummy1, i,j,k,totcycles,printcycles,a,b,c,d,i1,i2,i3,i4,gen;

char out[] = "out";

extern float  transm[2][3];

float 
h,s,sA,sa=0.0,DD[3],sibprob[2],sigma,K,sum,Wbar,Outbar,Sibbar,alpha,sumS,sumD,store[101][2][2][4],
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
for(i1=0;i1<100;i1++)  for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<3;i4++) store[i1][i2][i3][i4]=0.0;	

printf("alpha has half the effect as delta -- alpha = 1 is no distortion\n\
so to have 90%% distortion set alpha = 0.2\n");

//scanf("%f",&alpha);

//printf("input fitness cost to drive heterozygote\n");
//scanf("%f",&h);
K = 0.0;
sigma = 1.0;
sibprob[0] = 0.0;

// uncomment the appropriate line of the next 3 to set 'm' in Fig. 2, m = sibprob[1]

//sibprob[1] = 0.5;
sibprob[1] = 0.95;
//sibprob[1] = 0.2;

printf("m = %0.2f\n",sibprob[1]);

//printf("s \tsumD \tsumS \tWbar (m = %0.3f  alpha = %0.3f, h = %0.3f)\n",sibprob[1],alpha,h);

for(alpha = 0.0;alpha < 0.21;alpha+= 0.2)	{  // alpha will be either 0.0 or 0.2, if =.2, then drive is 90%
for(h=0.0; h<0.11;h+= 0.1)			{  // h will be either 0.0 or 0.1
for (s=0.0;s<0.9;s+= .01)	{

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)
allmated[i1][i2][i3][i4] = 0.0;

allmated[1][0][1][0] = 0.01;  // SdSd mated females
allmated[0][1][0][1] = 0.01;  // sDsD mated females
allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];


//sibprob[1] = sA;
DD[0] = DD[1] = 1.0;
DD[2] = 1.0 -s;
DD[1] = 1.0 -h;


for(gen=0;gen<100000;)	{


// Zero the matrix before use
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)
sibmated[i1][i2][i3][i4] = 0.0;

// sib mating from broods of all mated females -- only if she has allele 1 at first locus
// this is set with sibprob[], which may allow us to generalize the pgm
for(a=0;a<2;a++) for(b=0;b<2;b++) for(c=0;c<2;c++) for(d=0;d<2;d++)
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)
sibmated[a][b][c][d] += allmated[i1][i2][i3][i4] *sibprob[i1]*transm[a][i1+i3]*transm[c][i1+i3]*transm[b][i2+i4]*transm[d][i2+i4]*sigma;

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

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) outmale[i1][i2] = outfemale[i1][i2] =0.0;



for(a=0;a<2;a++)  for(c=0;c<2;c++) for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{

outmale[a][c] += allmated[i1][i2][i3][i4] *(1.0 + (K-1.0)*sibprob[i1]) * transm[a][i1+i3]*transm[c][i2+i4]; 

outfemale[a][c] += allmated[i1][i2][i3][i4] *(1.0 -sibprob[i1]) * transm[a][i1+i3]*transm[c][i2+i4]; 

										}


for(a=0,sum=0.0;a<2;a++)  for(c=0;c<2;c++) sum += outmale[a][c];


for(a=0;a<2;a++)  for(c=0;c<2;c++) outmale[a][c] /= sum;   // males normalized



for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	
outmated[i1][i2][i3][i4] = outfemale[i1][i2]*outmale[i3][i4];

// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
// Also reduce fitness of Dd drive 'heterozygotes'  -- D females mated with d sires or d x D
// the diploid phase will be DD, (or Dd or dD) thus suffering a fitness loss 
// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
if(i2+i4 == 2)	{
outmated[i1][i2][i3][i4] *= DD[2];
sibmated[i1][i2][i3][i4] *= DD[2];
		}
if(i2+i4 == 1)	{
outmated[i1][i2][i3][i4] *= DD[1];
sibmated[i1][i2][i3][i4] *= DD[1];
		}
										}


// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;
// so all progeny from Dd or dD mated females will be D
//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
if (i2 != i4) 	{			// for i2 != i4, one must be D
	outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4]*(1.0 - alpha);
	 outmated[i1][i2][i3][i4] = alpha*outmated[i1][i2][i3][i4];

	sibmated[i1][1][i3][1] += sibmated[i1][i2][i3][i4]*(1.0 - alpha);
	sibmated[i1][i2][i3][i4] = alpha*sibmated[i1][i2][i3][i4];

// Note that this zero cannot overwrite the previous cell because i2 and i4 are not both 1

		}								}
	
// All fitness effects should be accounted for above this line -- DD fitness loss and inbreeding depression
// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

// Don't need to zero allmated
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	
allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4] + sibmated[i1][i2][i3][i4];

for(Wbar = 0.0,i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
Wbar += allmated[i1][i2][i3][i4];
												}
gen++;

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)		{
 allmated[i1][i2][i3][i4] /= Wbar;
 sibmated[i1][i2][i3][i4] /= Wbar;
											}

					}  // end of gen print loop
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for(sumD=0.0,sumS=0.0,i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)		{
sumD +=  allmated[i1][i2][i3][i4] *((float)i2 + (float)i4)/2.0;
sumS +=  allmated[i1][i2][i3][i4] *((float)i1 + (float)i3)/2.0;
											}
i1 = (int)(s*100.1);
i2 = (int)(alpha>0.05);
i3 = (int)(h > 0.05);
store[i1][i2][i3][0] = sumD;
store[i1][i2][i3][1] = sumS;
store[i1][i2][i3][2] = Wbar;
printf("%d %d %d \n",i1,i2,i3);




				}  // end of s loop
}  // alpha loop
}  // h loop
// Want printout to be sumD, sumD sumD sumD, sumS sumS, sumS, sumS, wbar

fp = fopen(out,"a+");
fprintf(fp,"m = %0.2f\n",sibprob[1]);
fprintf(fp,"sumD is final frequency of drive allele, sumS is final frequency of sib mating allele, Wbar is final fitness\n");
fprintf(fp,"s    sumD (00  10  01  11)  sumS (00  10  01 11)  Wbar (00  10  01  11)\n\
where 00 is alpha =0.0, h = 0.0, ...,   11 is alpha = 0.2, h = 0.1\n");

for(i1=0;i1<100;i1++) {
	   fprintf(fp,"%0.3f  %0.3f %0.3f %0.3f %0.3f  %0.3f %0.3f %0.3f %0.3f  %0.3f %0.3f %0.3f %0.3f\n",
(float)i1/100.0,
store[i1][0][0][0],store[i1][1][0][0],store[i1][0][1][0],store[i1][1][1][0],
store[i1][0][0][1],store[i1][1][0][1],store[i1][0][1][1],store[i1][1][1][1],	
store[i1][0][0][2],store[i1][1][0][2],store[i1][0][1][2],store[i1][1][1][2]	);
		}
fprintf(fp,"\n");

fclose(fp);
}

