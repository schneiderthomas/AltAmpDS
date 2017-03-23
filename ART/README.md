# ART Folder

example bash script to generate artifical Trusight Tumor 26 fastq files


## FPA_normal.fa and FPB_normal.fa
these are fasta files created from the Trusight manifest files that include 
both the probes and the target region sequence together, you can use this as
your reference when using the trusight_tumor_fastq_generator.sh


## FPA_normal_with_mutations.fa and FPB_normal_with_mutations.fa
exactly the same as FPA_normal.fa and FPB_normal.fa however the last
three fasta sequences modify EGFR, KIT, or PTEN amplicons to create
the mutations described in the fasta description line


## trusight_tumor_fastq_generator.sh
a simple script file calling the ART program that makes an artifical fastq file
similar to trusight tumor 26, users will need to modify the top portion of the 
script to call the appropriate fasta files and ART program
