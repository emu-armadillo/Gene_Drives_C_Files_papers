I 2022

This file is a guide to the code underlying figures and results for Cook, Bull, Gomulkiewicz:  Gene drive escape from resistance depends on mechanism and ecology.  Several files of code are provided, each corresponding to a different figure. The files included contain all the code necessary to generate the data for all figures, however, some simple mathematics and/or interpretation of output data may be necessary in order to recreate the figures. The repository contains two directories one with the code used to generate the figures excluding figure two, which takes data from the figure 3 file and a trival modification to the figure 1 file. the second contains the underling models for the first

The gene frequency recursions in these programs evolution calculate offspring genotype numbers by exhaustive enumeration of all possible parental matings, multiplied by the fraction of gamete or offspring types produced by the parental genotypes.

In the code for 2 concurrent drives (Fig. 1), the key recursion for transmission from (haploid) fathers to progeny is in the 3-level, nested loop:

for (o1 = 0; o1 < 2; o1++) for (o2 = 0; o2 < 2; o2++) for (o3 = 0; o3 < 2; o3++)     
	mg[o1][o2][o3] += tm1[RR][p1][o1] * tm2[RR][p2][o2] * ttt[p3][o3] * z * Y;

where mg is the progeny gamete genotype at 3 haploid loci, Y is the paternal genotype frequency, z is its fitness, with tm1 and tm2 as transmission vectors that depend on the resistance genotype of the father (RR), the genotype at drive locus 1 or 2 in the father (p1, p2)) and the allele inherited by the progeny at that locus (o1 or o2).  The transmission vectors are trivially set and ensure that progeny inherit only the parental alleles at Mendelian or drive frequencies -- tm1 and tm2 are specifically for drive.


In contrast to the 2drive program, code for sib mating evolutino uses counts of mated (haploid) females. Under haploidy, her state represented by 4 indices  The first index is her allele at the sib-mating locus, the second index is her allele at the drive locus, and the third and fourth indices are the same, repectively, for her sire.   With 2 alleles per locus, each element in the vector can take two states, 0 or 1.  Drive can be manifested simply as converting a (D,d) female mated into a (D,D) mated female, and likewise for a (d,D) mated female, where D is the drive allele and d its wild-type allele.  This shortcut is based on the mated female producing a brief diploid progeny intermediate that then becomes haploid according to the rules of segregation.  Fitness is easily assigned to the mated female (accruing to her progeny) based on what progeny genoytpe she will produce, but fitness must be assigned before the code converts the drive locus (so that a Dd diploid is assigned the fitness for a heterozygote and not a DD homozygote).


Other programs follow similar rules.  

They compile on Mac OSX 10.14 with 'gcc <filename> -lm'

The programs and their associated figures are:

Fig. 1 (2 homing drives and type-M resistance)

Fig. 2 (1 homing drive and sib mating)

Fig. 3 (ClvR and sib mating)

Fig. 5, 6 (ecology)
	
	

