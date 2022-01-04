//type1.c to run recurrence equations of recessive partially lethal drive  with dominant resistance 27 vii 21
// Haploids with 2 loci

//This program can be used to manually find the bound of type1 resistance (grey box)

#include <stdio.h>
#include <math.h>

float transm[2][3] = {1.0,0.5,0.0,    //  for [0][0], [0][1], [0][2]
                      0.0, 0.5, 1.0 };   // [1][0], [1][1], [1][2]
// first index is what goes to the gamete/kid, second index is parent/diploid type
// diploid genotype is not ordered, we just count number of alleles

int main()
{

int dummy,dummy1, i,j,k,totcycles,printcycles,a,b,c,d,i1,i2,i3,i4,gen;

extern float  transm[2][3];

float 
sA,sa=0.0,DD[3],sibprob[2],sigma,K,sum,Wbar,Outbar,Sibbar,
allmated[2][2][2][2],   // ordered genotype for all mated females 
outmated[2][2][2][2],   // ordered genotype for all  females whose sires are not sibs
outmale[2][2],    // haploid males that will do the outcrossing
outfemale[2][2];    // haploid females that will do the outcrossing

// Note that individuals are haploid, so they have only allele 0 or 1 at each locus
// thus mated females have 4 alleles if we count their (haploid) sires
// first locus is resistance control (s, S), second locus is drive (d,D)
// if either parent has S, then drive blocked
// mated female indicies are ordered:  female-resistance locus, female-drive,sire-resistance,sire-drive



// Set to zero  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)
allmated[i1][i2][i3][i4] = 0.0;


printf("\n\n\tINPUT 3 numbers \n\n\
allmated[1][0][1][0], allmated[0][1][0][1]  = starting frequencies of Sd-Sd, and sD-sD mated females\n\
DD[2] = fitness of drive homozygotes \n");


scanf("%f %f %f",&allmated[1][0][1][0],&allmated[0][1][0][1],&DD[2]);

allmated[0][0][0][0] = 1.0 - allmated[1][0][1][0] - allmated[0][1][0][1];

DD[0] = DD[1] = 1.0;

gen = 0;

target1:;

printf(" number of cycles between outputting results and final number of generations for the run\n");
scanf("%d %d",&printcycles,&totcycles);

printf("resistance operates only if S is absent from both parents\n");
printf("gen  \tsdsd \tsdSd \tSdsd \tSdSd \t \tsDsD \tsDSD \tSDsD \tSDSD \tWbar   \n");


for(;gen<totcycles;)	{
	for(j=1;j<=printcycles;j++)	{



// diploid will have i1+i3 non-wt alleles at locus 1, i2+i4 non-wt alleles at locus 2
// drive is not present here because we converted any D/d or d/D mated females to D/D at a different stage
// starting frequencies were of D/D, but we can think of that as D/d then converted by drive
// this is the crux of inheritance, so we need to triple check that this does what we want
// no male-female difference,  so progeny inherit same whether mom or dad 
// had the allele

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) outmale[i1][i2] = outfemale[i1][i2] =0.0;



for(a=0;a<2;a++)  for(c=0;c<2;c++) for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{

outmale[a][c] += allmated[i1][i2][i3][i4] *transm[a][i1+i3]*transm[c][i2+i4]; 

outfemale[a][c] += allmated[i1][i2][i3][i4] * transm[a][i1+i3]*transm[c][i2+i4]; 

										}


for(a=0,sum=0.0;a<2;a++)  for(c=0;c<2;c++) sum += outmale[a][c];


for(a=0;a<2;a++)  for(c=0;c<2;c++) outmale[a][c] /= sum;   // males normalized



for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	
outmated[i1][i2][i3][i4] = outfemale[i1][i2]*outmale[i3][i4];

// Now reduce fitness of DD drive 'homozygotes'  -- D females mated with D sires
// the diploid phase will be DD, thus siffering a fitness loss 
// Note that we do this BEFORE we convert D/d and d/D mated females to D/D

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
if(i2+i4 == 2)	{
outmated[i1][i2][i3][i4] *= DD[2];
		}
										}


// A trick:  we Now convert all D/d and d/D mated females to D/D -- because 100% drive will operate in all diploids;  but only if S absent
// so all progeny from Dd or dD mated females will be D
//We do it now to avoid having to keep track of segregation, and do it after we imposed a fitness effect on D/D 


for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	{
if ((i1+i3==0) &&(i2 != i4)) 	{			
// for i2 != i4, one must be D; S absent only if i1+i3 =0
	outmated[i1][1][i3][1] += outmated[i1][i2][i3][i4];
	 outmated[i1][i2][i3][i4] = 0.0;

// Note that this zero cannot overwrite the previous cell because i2 and i4 are not both 1

		}								}
	
// All fitness effects should be accounted for above this line -- DD fitness loss 

// Don't need to zero allmated
for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	
allmated[i1][i2][i3][i4] = outmated[i1][i2][i3][i4];

for(i1=0,Wbar=0.0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)	
Wbar += allmated[i1][i2][i3][i4];

for(i1=0;i1<2;i1++) for(i2=0;i2<2;i2++) for(i3=0;i3<2;i3++) for(i4=0;i4<2;i4++)		{
 allmated[i1][i2][i3][i4] /= Wbar;
											}

gen++;
					}  // end of inner print loop

printf("%d  \t%1.3f \t%1.3f \t%1.3f \t%1.3f \tdd/DD \t%1.3f \t%1.3f \t%1.3f \t%1.3f \t%1.3f\n",gen,
allmated[0][0][0][0], allmated[0][0][1][0], allmated[1][0][0][0], allmated[1][0][1][0], 
allmated[0][1][0][1], allmated[0][1][1][1], allmated[1][1][0][1], allmated[1][1][1][1],Wbar);


				}
printf(" DD[2] = %1.4f  \n\n\n",DD[2]);
 goto target1;

}
