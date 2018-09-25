#!/usr/bin/env bash


#Homer Path
PATH=$PATH:/Users/M/Google_Drive/Scripts/motiff/homer/.//bin/

#homer2 denovo -i HS_combined_ADOG_Annotated_sequences_500.fa -b Random_Intergenic_500.fa > HS_homer_500.txt
#findMotifs.pl HS_combined_ADOG_Annotated_sequences_500.fa fasta HS_homer_500 -fasta <background.fa> [options]



############ HS MOTIFS
#homer2 denovo -i HS_combined_ADOG_Annotated_sequences_BO_500.fa -b Random_Intergenic_500.fa > OUT/HS_combined_ADOG_Annotated_sequences_BO_500_DENOVO.txt
#findMotifs.pl HS_combined_ADOG_Annotated_sequences_BO_500.fa worm OUT/HS_combined_ADOG_Annotated_sequences_BO_500_findmotif_OUT -fasta Random_Intergenic_500.fa
#
#homer2 denovo -i HS_combined_ADOG_Annotated_sequences_HS_500.fa -b Random_Intergenic_500.fa > OUT/HS_combined_ADOG_Annotated_sequences_HS_500_DENOVO.txt
#findMotifs.pl HS_combined_ADOG_Annotated_sequences_HS_500.fa worm OUT/HS_combined_ADOG_Annotated_sequences_HS_500_findmotif_OUT -fasta Random_Intergenic_500.fa
#
#homer2 denovo -i HS_combined_ADOG_Annotated_sequences_WT_500.fa -b Random_Intergenic_500.fa > OUT/HS_combined_ADOG_Annotated_sequences_WT_500_DENOVO.txt
#findMotifs.pl HS_combined_ADOG_Annotated_sequences_WT_500.fa worm OUT/HS_combined_ADOG_Annotated_sequences_WT_500_findmotif_OUT -fasta Random_Intergenic_500.fa
