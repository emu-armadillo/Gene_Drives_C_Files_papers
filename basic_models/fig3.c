/*fig3.c, formerly sib_clvR.c to run recurrence equations of sib mating with ClvR drive that has fitness effects 14 vii 21
 Haploids with 2 loci
In contrast to the data presented in Fig. 3, this program runs single conditions specified by input and generates the evolution over time
at intervals specified on input
It is easily modified to loop through different conditions and to create output at a specified number of generations or after conditions have been
met.


*/



#include <stdio.h>
#include <math.h>

float transm[2][3] = {1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
                      0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles

int main()
{

int dummy,dummy1, i,j,k,totcycles,printcycles,a,b,c,d,e,f,i1,i2,i3,i4,i5,i6,gen;

extern float  transm[2][3];

float 
sA,sa=0.0,CC[3],sibprob[2],sigma,K,sum,Wbar,Outbar,Sibbar,total[2][2][2],totals[2][2][2],
allmated[2][2][2][2][2][2],   // ordered genotype for all mated females 
outmated[2][2][2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
sibmated[2][2][2][2][2][2],   // ordered genotype for all  females whose sires are sibs
outmale[2][2][2],    // haploid males that will do the outcrossing
outfemale[2][2][2];    // haploid females that will do the outcrossing

// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
// thus mated females have 6 alleles if we count their (haploid) sires
// first locus is sibmating control (s, S), second locus is ClvR (c,C), third is target (g,G)
// GG dies w/o at least one C; Cc, cC, and CC conert all g to G in diploid
// if mom has 's' then all kids go to outcross pool; if she has 'S' then sA are reserved for sibmating
// mated female indicies are ordered:  female-sib locus, female-drive,sire-sib,sire-drive
// Note that K=1 provides an 'unfair' advantage for selfing, so our runs should have K=0



// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 
allmated[i1][i2][i3][i4][i5][i6] = 0.0;


printf("\n\n A program to model the evolution of sib mating as a possible block to a Cleave-Rescue gene drive\n\n\
\tINPUT 6 numbers \n\n\
allmated[1][0][0][1][0][0], allmated[0][1][0][0][1][0]  = starting frequencies of Sd-Sd, and sC-sC mated females\n\
\t S is sib mating allele, C is ClvR allele\n\
sigma = fitness of offspring from sib-mated parents (assumed to be 1 in Fig. 3)\n\
K = discounting of males who sib mat (assumed to 0 in Fig. 3, max of 1) \n\
m = sibmating rate of S (wildtype is outcrossed), formerly sA\n\
CC[2] = fitness of ClvR homozygotes, hence the level of suppression\n");


scanf("%f %f %f %f %f %f",&allmated[1][0][0][1][0][0],&allmated[0][1][0][0][1][0],&sigma,&K,&sA,&CC[2]);

allmated[0][0][0][0][0][0] = 1.0 - allmated[1][0][0][1][0][0] - allmated[0][1][0][0][1][0];

sibprob[0] = 0.0;
sibprob[1] = sA;
CC[0] = CC[1] = 1.0;

gen = 0;
//printf("%1.3f %1.3f %1.3f (all 000000, 100100, 010010)\n", allmated[0][0][0][0][0][0],allmated[1][0][0][1][0][0], allmated[0][1][0][0][1][0]);

target1:;

printf(" number of cycles between outputting results and final number of generations for the run\n");
scanf("%d %d",&printcycles,&totcycles);

printf("sib mating only for S\n");
printf("gen  \tscg \tscG \tsCg \tsCG \t \tScg \tScG \tSCg \tSCG \tWbar  \tm (%0.2f) \n",sA);


for(;gen<totcycles;)	{
	for(j=1;j<=printcycles;j++)	{


//for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) 
//for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++)   
//printf("%0.2f ",allmated[i1][i2][i3][i4][i5][i6]);
//printf("\n");

// Zero the matrix before use
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) 
for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++)   sibmated[i1][i2][i3][i4][i5][i6] = 0.0;

// sib mating from broods of all mated females -- only if she has allele 1 at first locus
// this is set with sibprob[], which may allow us to generalize the pgm
for(a=0;a<2;a++) for(b=0;b<2;b++) for(c=0;c<2;c++) for(d=0;d<2;d++) for(e=0;e<2;e++) for(f=0;f<2;f++) {
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++)  {

sibmated[a][b][c][d][e][f] += allmated[i1][i2][i3][i4][i5][i6] *sibprob[i1]*transm[a][i1+i4]*transm[d][i1+i4]*transm[b][i2+i5]*transm[e][i2+i5]*transm[c][i3+i6]*transm[f][i3+i6]*sigma;
	}		}

//printf("%1.3f %1.3f %1.3f (all 000000, 100100, 010010)\n", allmated[0][0][0][0][0][0],allmated[1][0][0][1][0][0], allmated[0][1][0][0][1][0]);
//printf("%1.3f (sib 000000)\n",sibmated[0][0][0][0][0][0]);
// diploid will have i1+i4 non-wt alleles at locus 1, i2+i5 non-wt alleles at locus 2; i3+i6 at locus 3
// violations of Mendel are not present here because we converted any G/g or g/G or G/G mated females to g/g at a different stage
// starting frequencies were of C/C, 
// inbreeding depression taken care of by sigma fitness effect
// no male-female difference except in control of sib mating, so progeny inherit same whether mom or dad 
// had the allele
// Note that we are assigning inbreeding depression early -- to a female whose sire was a sib.  This anticipates 
// that her progeny will suffer (in survival), even though it is not she who has the problem with inbreeding.
// This works if we are careful about it, but I'd like to think about it more.  
// In particular, we discount sib-mated females in advance because they will have fewer offspring, not because they 
// are themselves reduced

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++)  outmale[i1][i2][i3] = outfemale[i1][i2][i3] =0.0;



for(a=0;a<2;a++) for(b=0;b<2;b++)  for(c=0;c<2;c++) for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++)	{

outmale[a][b][c] += allmated[i1][i2][i3][i4][i5][i6] *(1.0 + (K-1.0)*sibprob[i1]) * transm[a][i1+i4]*transm[b][i2+i5]*transm[c][i3+i6]; 

outfemale[a][b][c] += allmated[i1][i2][i3][i4][i5][i6]*(1.0 -sibprob[i1]) * transm[a][i1+i4]*transm[b][i2+i5]*transm[c][i3+i6]; 

										}
//printf("%1.3f  %1.3f outfemale (000), outmale\n",outfemale[0][0][0],outmale[0][0][0]);

for(a=0,sum=0.0;a<2;a++)  for(b=0;b<2;b++) for(c=0;c<2;c++) sum += outmale[a][b][c];

for(a=0;a<2;a++)  for(b=0;b<2;b++) for(c=0;c<2;c++)  outmale[a][b][c] /= sum;   // males normalized



for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{
outmated[i1][i2][i3][i4][i5][i6] = outfemale[i1][i2][i3]*outmale[i4][i5][i6];
	}

//printf("%1.3f  outmated (000)\n",outmated[0][0][0][0][0][0]);

// Now reduce fitness of CC drive 'homozygotes'  -- C females mated with C sires
// the diploid phase will be CC, thus suffering a fitness loss 
// Note that we do this BEFORE we convert g/g, G/g and g/G mated females to G/G
// We also kill G/G here, since G/G die only if C is absent, so there is no worry of precedence

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{
if((i2==1) && (i5 == 1))	{
outmated[i1][i2][i3][i4][i5][i6] *= CC[2];
sibmated[i1][i2][i3][i4][i5][i6] *= CC[2];
		}

if((i3+i6 ==2) && (i2+i5 == 0))	{   // if all G and no C, we will kill them off; no overlap with CC lines above
outmated[i1][i2][i3][i4][i5][i6]  = 0.0;
sibmated[i1][i2][i3][i4][i5][i6]  = 0.0;
				}
										}


// A trick:  we Now convert all g/g G/g and g/G mated females to G/G if C present -- because of ClvR
// all progeny from Cc or cC mated females will be G once the diploid is formed
//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on C/C and killing G/G; note that we are killing G/G zygotes in the absence of C, but would not kill G/G conversions; in any case, C is present for these conversions


for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{

if (((i2 + i5)>0) && ((i3+i6) <2)) 	{			// i3+16<2 means that at least one allele is wt g
			// i2 +i5 >0 means at least one is ClvR 
			// Otherwise the zeroing would erase the coverted genotype
	outmated[i1][i2][1][i4][i5][1] += outmated[i1][i2][i3][i4][i5][i6];  // converts any g to G
	 outmated[i1][i2][i3][i4][i5][i6] = 0.0;

	sibmated[i1][i2][1][i4][i5][1] += sibmated[i1][i2][i3][i4][i5][i6];
	sibmated[i1][i2][i3][i4][i5][i6] = 0.0;

		}								}
	
// All fitness effects should be accounted for above this line -- CC fitness loss and inbreeding depression
// Thus, we can now combine mated females from outcrosses and sibmating;  best to think about this

// Don't need to zero allmated
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{
allmated[i1][i2][i3][i4][i5][i6] = outmated[i1][i2][i3][i4][i5][i6] + sibmated[i1][i2][i3][i4][i5][i6];
		}
//printf("allmated 000000 = %0.3f\n",allmated[0][0][0][0][0][0]);

for(i1=0,Wbar=0.0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++) for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{
Wbar += allmated[i1][i2][i3][i4][i5][i6];
			}
//printf("wbar %0.3f\n",Wbar);

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	 for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{
 allmated[i1][i2][i3][i4][i5][i6] /= Wbar;
 sibmated[i1][i2][i3][i4][i5][i6] /= Wbar;
											}

gen++;
					}  // end of inner print loop

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) total[i1][i2][i3] = 0.0;

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	 for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{ 
total[i1][i2][i3] += allmated[i1][i2][i3][i4][i5][i6]/2.0 ;
total[i4][i5][i6] += allmated[i1][i2][i3][i4][i5][i6]/2.0 ;
		}

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) totals[i1][i2][i3] = 0.0;

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	 for(i5=0;i5<2;i5++) for(i6=0;i6<2;i6++) 	{ 
totals[i1][i2][i3] += sibmated[i1][i2][i3][i4][i5][i6]/2.0 ;
totals[i4][i5][i6] += sibmated[i1][i2][i3][i4][i5][i6]/2.0 ;
		}


printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t(s/S) \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",gen,
total[0][0][0], total[0][0][1], total[0][1][0], total[0][1][1], 
total[1][0][0], total[1][0][1], total[1][1][0], total[1][1][1], Wbar);

printf("sibm:  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n\n",
totals[0][0][0], totals[0][0][1], totals[0][1][0], totals[0][1][1], 
totals[1][0][0], totals[1][0][1], totals[1][1][0], totals[1][1][1]);

/*
// Now calculate freqs w/in groups
for(i1=0,Sibbar = Outbar=0.0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
Outbar += outmated[i1][i2][i3][i4];
Sibbar += sibmated[i1][i2][i3][i4];
													}

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
outmated[i1][i2][i3][i4] /=Outbar;
sibmated[i1][i2][i3][i4] /=Sibbar;
													}

printf("Out  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \tdd/CC \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t\t freqs per Group\n",
outmated[0][0][0][0], outmated[0][0][1][0], outmated[1][0][0][0], outmated[1][0][1][0], 
outmated[0][1][0][1], outmated[0][1][1][1], outmated[1][1][0][1], outmated[1][1][1][1],Wbar);

printf("sibm:  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \tdd/CC \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t\t freqs per Group\n\n",
sibmated[0][0][0][0], sibmated[0][0][1][0], sibmated[1][0][0][0], sibmated[1][0][1][0], 
sibmated[0][1][0][1], sibmated[0][1][1][1], sibmated[1][1][0][1], sibmated[1][1][1][1]);
*/


				}

printf("eqm C = %1.3f   m = %1.2f  sigma = %1.2f CC[2] = %1.4f  \n\n\n",1.0/(2.0-CC[2]),sA,sigma,CC[2]);
 goto target1;

}

