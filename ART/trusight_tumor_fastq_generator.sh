#!/bin/bash

##MODIFY AS NECESSARY FOR YOUR OWN FOLDER
#illumina test examples
#art=/home/tom/Documents/bioinformatics/ART/bin/art_illumina
art=/home/tom/Documents/bioinformatics/art_bin_ChocolateCherryCake/art_illumina
reference3=/home/tom/Documents/bioinformatics/ART/Testing/FPA_Complex_Mutations_1_thru_11.fa
reference4=/home/tom/Documents/bioinformatics/ART/Testing/FPB_Complex_Mutations_1_thru_11.fa
#################################


# 7) amplicaton read simulation: generate one 121bp paired-end reads from both ends for each amplicon reference


#####
###DEPTH OF 10 TO MAKE REALLY SMALL FILE FOR PURPOSES OF DEBUGGING PIPELINE AS IT IS VERY
###FAST TO PROCESS
#$art -i $reference3  -amp  -o ./Complex_Mutations_1_thru_11_small_S0_R -p -l 121 -f 10 --seqSys MS 
#######
$art -i $reference3  -amp  -o ./Complex_Mutations_1_thru_11_large_S0_R -p -l 121 -f 500 --seqSys MS 

rm ./Complex_Mutations_1_thru_11_large_S0_R1.aln 
rm ./Complex_Mutations_1_thru_11_large_S0_R2.aln


#######
###DEPTH OF 10 TO MAKE REALLY SMALL FILE FOR PURPOSES OF DEBUGGING PIPELINE AS IT IS VERY
###FAST TO PROCESS
#$art -i $reference4  -amp  -o ./Complex_Mutations_1_thru_11_small_S1_R -p -l 121 -f 10 --seqSys MS 
#######

$art -i $reference4  -amp  -o ./Complex_Mutations_1_thru_11_large_S1_R -p -l 121 -f 500 --seqSys MS  

#rm ./FPB_EGFR_deletion_BRAF_snp_S0.sam 
rm ./Complex_Mutations_1_thru_11_large_S1_R2.aln 
rm ./Complex_Mutations_1_thru_11_large_S1_R1.aln




