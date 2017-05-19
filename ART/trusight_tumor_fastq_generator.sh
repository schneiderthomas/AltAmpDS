#!/bin/bash

##MODIFY AS NECESSARY FOR YOUR OWN FOLDER
#illumina test examples
#art=/home/tom/Documents/bioinformatics/ART/bin/art_illumina
basedir=$(cd ../ && pwd)
art=$basedir/art_bin_ChocolateCherryCake/art_illumina
#reference1=$PWD/FPA_normal_with_mutations.fa
#reference2=$PWD/FPB_normal_with_mutations.fa
reference1=$PWD/FPA_normal.fa
reference2=$PWD/FPB_normal.fa
#################################


# 7) amplicaton read simulation: generate one 121bp paired-end reads from both ends for each amplicon reference


#####
###DEPTH OF 10 TO MAKE REALLY SMALL FILE FOR PURPOSES OF DEBUGGING PIPELINE AS IT IS VERY
###FAST TO PROCESS
#$art -i $reference1  -amp  -o ./normal_with_mutations_S0_R -p -l 121 -f 10 --seqSys MS 
#######
#$art -i $reference1  -amp  -o ./normal_with_mutations_S0_R -p -l 121 -f 10 --seqSys MS 

#rm ./normal_with_mutations_S0_R1.aln 
#rm ./normal_with_mutations_S0_R2.aln


$art -i $reference1  -amp  -o ./normal_S0_R -p -l 121 -f 10 --seqSys MS 

rm ./normal_S0_R1.aln 
rm ./normal_S0_R2.aln

mv normal_S0_R1.fq normal_S0_R1.fastq
mv normal_S0_R2.fq normal_S0_R2.fastq


#######
###DEPTH OF 10 TO MAKE REALLY SMALL FILE FOR PURPOSES OF DEBUGGING PIPELINE AS IT IS VERY
###FAST TO PROCESS
#$art -i $reference2  -amp  -o ./normal_with_mutations_S1_R -p -l 121 -f 10 --seqSys MS 
#######

#$art -i $reference2  -amp  -o ./normal_S1_R -p -l 121 -f  --seqSys MS  


#rm ./normal_with_mutations_S1_R1.aln 
#rm ./normal_with_mutations_S1_R2.aln


$art -i $reference2  -amp  -o ./normal_S1_R -p -l 121 -f 10 --seqSys MS  


rm ./normal_S1_R1.aln 
rm ./normal_S1_R2.aln

mv normal_S1_R1.fq normal_S1_R1.fastq
mv normal_S1_R2.fq normal_S1_R2.fastq




