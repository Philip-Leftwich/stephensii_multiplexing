# stephensii_multiplexing

- 1590 = Cas9 (B element)

- A = gRNA elements

Classing multiplexing strategy

•	Cutting assay of the singleplex at target site A in Cd (Crossing data summary A1759 B_1590): you might have seen this data before for the PolII promoters. However I think it is important to say that I changed the order of the crosses to female:male to match the rest of the crosses. 
•	Cutting assay of the singleplex deletion at target site A in Cd ( Crossing data summary A_2072 B_1590)
•	Cutting assay of the multiplex deletion at target site A in Cd with 3 additional sgRNAs ( Crossing data summary A_2301 B_1590). 
•	Cutting assay of the singleplex at target site A in Cd in background of a homozygous resistance allele at target site A (Crossing data summary A_1759 B_1590xQA383P): resistance is at the sgRNA 1 target site so no cutting should occur. 
•	Cutting assay of the multiplex deletion at target site A in Cd with 3 additional sgRNAs in background of a homozygous resistance allele at target site A (Crossing data summary A_2301 B_1590xQA383P):  resistance is at the sgRNA 1 target site so cutting should be due to the 3 additional gRNAs (which were not individually tested). 

Separate multiplexing strategy: 
•	Cutting assay of the singleplex at target site B in cd ( Crossing data summary A_2273 B_1590): I am not sure if it is relevant but this singleplex line was a rescue in comparison to the others. 
•	Cutting assay of the singleplex at target site B in cd, in the background of a homozygous resistance allele at target site A ( Crossing data summary A_2273 B_1590xQA383P).
•	Cutting assay of singleplex at target site A in cd ( Crossing data summary A_1759 B_1590): is the same cross used for the classign multiplexing strategy. 
•	Cutting assay of the singleplex at target site A in cd, in the background of a homozygous resistance allele at target site 2 ( Crossing data summary A_1759 B_1590xD251). 

Note that I cleaned all the data following the same template and all the crosses are set female:male order. Crosses that only have one direction are due to fitness cost from the Cd lines. The common direction in all the crosses is: grandmother Cas9 and father trans-het since the fitness cost mainly affects females. 


## Indel data


3m is a variant of sgRNA 3 with a SNP which is common in the wild type (SDA-500) (so samples lacking sgRNA3 would be those lacking this variant and exact sgRNA 3 copy) 
4m is a variant of sgRNA 4 with a SNP which is common in the wild type (SDA-500), (so samples lacking sgRNA4 would be those lacking this variant and exact sgRNA 4 copy). 

I have included the raw file used and made explanations in the file.

We considered only modified reads and alleles having more than 10 reads.


