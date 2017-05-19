#!/bin/sh



####
#CONSTANT VARIABLES
#
####
####
#change when moving to different folder system
####
HOME_DIR=$PIPELINE_DIR
THREADS=25
MEM=16
###
export PATH=$PATH:$HOME_DIR/bwa-0.7.10
BEDTOOLS_DIR=$HOME_DIR/bedtools2/bin/
TRIMMOMATIC_JAR=$HOME_DIR/Trimmomatic-0.33/trimmomatic-0.33.jar
FASTQC_DIR=$HOME_DIR/FastQC
PICARD_DIR=$HOME_DIR/picard
SAMTOOLS_DIR=$HOME_DIR/samtools-1.2
IGVTOOLS_JAR=$HOME_DIR/IGVTools/igvtools.jar
VARSCAN_JAR=$HOME_DIR/VarScan/VarScan.v2.4.0.jar
VCFLIB_HOME=$HOME_DIR/vcflib/bin
FREEBAYES_HOME=$HOME_DIR/freebayes/bin
BCFTOOLS_DIR=$HOME_DIR/bcftools-1.2/
export PATH=$PATH:$HOME_DIR/htslib-1.2.1
SNPEFF_DIR=$HOME_DIR/snpEff
ANNOVAR_HOME=$HOME_DIR/annovar
AMPLICONS_BED_VARIANTS=$HOME_DIR/amplicons_variants_tmp_from_manifest_intervals.bed
AMPLICONS_BED_SAMTOOLS=$HOME_DIR/amplicons_samtools_tmp_from_manifest_intervals_merged.bed
HG19=$HOME_DIR/hg19/WholeGenomeFASTA/genome.fa
COVERAGEQC_DIR=$HOME_DIR/Alternative_QC_Reporter/CoverageQC
GATK_DIR=$HOME_DIR/GenomeAnalysisTK-3.5
ANNOVAR_DB_HOME=$HOME_DIR/annovar/humandb
manifest_libA=$HOME_DIR/Manifest_Folder/TruSightTumor-FPA-Manifest_RevB.txt
manifest_libB=$HOME_DIR/Manifest_Folder/TruSightTumor-FPB-Manifest_RevB.txt
AMPLICONS_BED_VARIANTS_NAMED=$HOME_DIR/amplicons_variants_named.bed
consensus_RefSeq_CMP26=$HOME_DIR/consensus_RefSeq_CMP26_Original.csv
validation=false
debugging=false

usage ()
{
  echo 'Usage : Script  -fqrdir <directory where fastq files are located> -fqr1 <library 1 read 1 fastq> -fqr2 <library 1 read 2 fastq> -fqr1 <library 2 read 1 fastq> -fqr2 <library 2 read 2 fastq> -s <save directory> -sname <save name> -vcforig <trusight vcf file> -sequencer <sequencer using, either MiSeq or NextSeq> -manifest_libA <the illumina provided manifest file for library A, trusight tumor if not provided> -manifest_libB <the illumina provided manifest file for library B, trusight tumor default if not provided> -validation <if validating> -debug <if debugging>'
  exit
}




while [ "$1" != "" ]; do
case $1 in
        -fqrdir )   shift
                    FQRDIR=$1/
                    ;;
        -fqr1 )           shift
                    tempFQR1=$1
                       ;;
        -fqr2 )        shift
                    tempFQR2=$1
                       ;;
        -fqr3 )           shift
                       tempFQR3=$1
                       ;;
        -fqr4 )        shift
                       tempFQR4=$1
                       ;;         
        -s )        shift
                       DIR=$1
                       ;;
        -sname )        shift
                    OUTNAME=$1
                       ;; 
        -vcforig )    shift
                    VCFORIG=$1
                       ;;
        -home_dir )    shift
                    HOME_DIR=$1
                       ;;
        -sequencer )   shift
		    SEQUENCER=$1
		    ;;
        -manifest_libA ) shift
		    manifest_libA=$1
		    ;;
        -manifest_libB ) shift
		    manifest_libB=$1
		    ;;
        -validation ) shift
		    validation=$1
		    ;;
        -debugging ) shift
		    debugging=$1
		    ;;
        -threads ) shift
		   THREADS=$1
		    ;;
	-memory ) shift
		   MEM=$1
		    ;;	    
    esac
    shift
done

if [ "$FQRDIR" = "" ]
then
    echo "line 100"
    usage
fi
if [ "$tempFQR1" = "" ]
then
    echo "line 105"
    usage
fi
if [ "$tempFQR2" = "" ]
then
    echo "line 110"
    usage
fi
if [ "$tempFQR3" = "" ]
then
    echo "line 115"
    usage
fi
if [ "$tempFQR4" = "" ]
then
    echo "line 120"
    usage
fi
if [ "$DIR" = "" ]
then
    echo "line 125"
    usage
fi
if [ "$OUTNAME" = "" ]
then
    echo "line 130"
    usage
fi

if [ "$SEQUENCER" = "" ] && [ "$SEQUENCER" != MiSeq ] && [ "$SEQUENCER" != NextSeq ]
then
    echo "line 136"
    usage
fi
#NOT DOING BLANK CHECK FOR VCFORIG SINCE IT IS GOING TO BE AN OPTIONAL INPUT
#


echo "ALL IS WELL.  Library 1 Read 1 FASTQ =$tempFQR1, library 1 Read 2 FASTQ =$tempFQR2, Library 2 Read 1 FASTQ =$tempFQR3, library 2 Read 2 FASTQ =$tempFQR4 OUTNAME= $OUTNAME DIRECTORY = $DIR FASTQ DIRECTORY =$FQRDIR" 







echo "`date`: start"




#unzipping if necessary

if   [ -f "$DIR/$tempFQR1" ]  &&  [ -s "$DIR/$tempFQR1" ]  &&  [ -f "$DIR/$tempFQR2" ]  &&   [ -s "$DIR/$tempFQR2" ]  &&  [ -f "$DIR/$tempFQR3" ]  &&  [ -s "$DIR/$tempFQR3" ]  &&  [ -f "$DIR/$tempFQR4" ]  &&  [ -s "$DIR/$tempFQR4" ] 
#&& [ ! [[-f "$DIR/$tempFQR2"] && [-s "$DIR/$tempFQR2"]] && ! [[-f "$DIR/$tempFQR3"] && [-s "$DIR/$tempFQR3"]] && ! [[-f "$DIR/$tempFQR4"] && [-s "$DIR/$tempFQR4"]] ]
then
echo "already moved fasqs and unzipped"
elif [ -f "$FQRDIR$tempFQR1" ]  &&  [ -s "$FQRDIR$tempFQR1" ]  &&  [ -f "$FQRDIR$tempFQR2" ]  &&   [ -s "$FQRDIR$tempFQR2" ]  &&  [ -f "$FQRDIR$tempFQR3" ]  &&  [ -s "$FQRDIR$tempFQR3" ]  &&  [ -f "$FQRDIR$tempFQR4" ]  &&  [ -s "$FQRDIR$tempFQR4" ] 
then
echo "already unzipped, moving files"
cp $FQRDIR$tempFQR1 $DIR/$tempFQR1 &
cp $FQRDIR$tempFQR2 $DIR/$tempFQR2 &
cp $FQRDIR$tempFQR3 $DIR/$tempFQR3 &
cp $FQRDIR$tempFQR4 $DIR/$tempFQR4 &
wait

elif  [ -f "$FQRDIR$tempFQR1.gz" ]  &&  [ -s "$FQRDIR$tempFQR1.gz" ]  &&  [ -f "$FQRDIR$tempFQR2.gz" ]  &&   [ -s "$FQRDIR$tempFQR2.gz" ]  &&  [ -f "$FQRDIR$tempFQR3.gz" ]  &&  [ -s "$FQRDIR$tempFQR3.gz" ]  &&  [ -f "$FQRDIR$tempFQR4.gz" ]  &&  [ -s "$FQRDIR$tempFQR4.gz" ] 
then
gunzip -c $FQRDIR$tempFQR1.gz > $DIR/$tempFQR1 &
gunzip -c $FQRDIR$tempFQR2.gz > $DIR/$tempFQR2 &
gunzip -c $FQRDIR$tempFQR3.gz > $DIR/$tempFQR3 &
gunzip -c $FQRDIR$tempFQR4.gz > $DIR/$tempFQR4 &
echo "UNZIPPINGS"
wait

else
echo "There are no fastqs"
echo $FQRDIR$tempFQR1
echo $FQRDIR$tempFQR2
echo $FQRDIR$tempFQR3
echo $FQRDIR$tempFQR4
exit

 
fi


wait
FQR1=$DIR/$tempFQR1
FQR2=$DIR/$tempFQR2
FQR3=$DIR/$tempFQR3
FQR4=$DIR/$tempFQR4

currentR1number=""
currentR3number=""

###in case order is reversed because of incorrect sorting of s9 and s10
if [[ $FQR1 =~([[:digit:]]+)_L ]]
then
	currentR1number=${BASH_REMATCH[1]}
elif [[ $FQR1 =~([[:digit:]]+)_R ]]
then 
        currentR1number=${BASH_REMATCH[1]}
else 
	echo "can't double check sorted correctly, check the name of the fastq file, follow the usual illumina convention ie samplehnnum_S9_L001_R1_001.fastq or samplehnnum_S9_R1_001.fastq"
	echo "the fastq file name is below"
	echo $FQR1
	exit
fi
if [[ $FQR3 =~([[:digit:]]+)_L ]]
then
	currentR3number=${BASH_REMATCH[1]}
elif [[ $FQR3 =~([[:digit:]]+)_R ]]
then 
        currentR3number=${BASH_REMATCH[1]}
else 
	echo "can't double check sorted correctly, check the name of the fastq file, follow the usual illumina convention ie samplehnnum_S9_L001_R1_001.fastq or samplehnnum_S9_R1_001.fastq"
	echo "the fastq file name is below"
	echo $FQR3
	exit
fi
echo "The currentR1number is: $currentR1number"
echo "The currentR3number is: $currentR3number"
LIBA_SNUMBER=$currentR1number
LIBB_SNUMBER=$currentR3number
# -gt means greater than, -lt means less than
if [[ $currentR1number -gt $currentR3number ]]
then
#need to switch the order of the fastqs
echo $currentR1number
echo $currentR3number
LIBA_SNUMBER=$currentR3number
LIBB_SNUMBER=$currentR1number
FQR1=$DIR/$tempFQR3
FQR2=$DIR/$tempFQR4
FQR3=$DIR/$tempFQR1
FQR4=$DIR/$tempFQR2
echo $FQR1
echo $FQR2
echo $FQR3
echo $FQR4
echo "switched"
fi





python $HOME_DIR/create_bed_from_manifest.py -i $manifest_libA -o $HOME_DIR/amplicons_samtools_tmp_from_manifest_intervals.bed -e 1
python $HOME_DIR/create_bed_from_manifest.py -i $manifest_libA -o $AMPLICONS_BED_VARIANTS -e 0
#$BEDTOOLS_DIR/bedtools merge -i $HOME_DIR/amplicons_samtools_tmp_from_manifest_intervals.bed > $AMPLICONS_BED_SAMTOOLS
$BEDTOOLS_DIR/mergeBed -i $HOME_DIR/amplicons_samtools_tmp_from_manifest_intervals.bed > $AMPLICONS_BED_SAMTOOLS

perl $FASTQC_DIR/fastqc -o $DIR $FQR1 &
perl $FASTQC_DIR/fastqc -o $DIR $FQR2 &
perl $FASTQC_DIR/fastqc -o $DIR $FQR3 &
perl $FASTQC_DIR/fastqc -o $DIR $FQR4 &
wait

if [[ $SEQUENCER = "NextSeq" ]]
then

java -Xmx${MEM}g -jar ${TRIMMOMATIC_JAR} PE -phred33 -threads ${THREADS} $FQR1 $FQR2 $FQR1.cut.fastq $FQR1.unpaired.fastq $FQR2.cut.fastq $FQR2.unpaired.fastq SLIDINGWINDOW:5:15 MINLEN:100 &
java -Xmx${MEM}g -jar ${TRIMMOMATIC_JAR} PE -phred33 -threads ${THREADS} $FQR3 $FQR4 $FQR3.cut.fastq $FQR3.unpaired.fastq $FQR4.cut.fastq $FQR4.unpaired.fastq SLIDINGWINDOW:5:15 MINLEN:100 &
wait


FQR1BWAMEM=$FQR1.cut.fastq
FQR2BWAMEM=$FQR2.cut.fastq
FQR3BWAMEM=$FQR3.cut.fastq
FQR4BWAMEM=$FQR4.cut.fastq
tID=${SEQUENCER}.AllLanes
else

FQR1BWAMEM=$FQR1
FQR2BWAMEM=$FQR2
FQR3BWAMEM=$FQR3
FQR4BWAMEM=$FQR4
tID=${SEQUENCER}.Lane1

fi

if [ "$debugging" = false ]
then
rm -f $FQR1.unpaired.fastq
rm -f $FQR2.unpaired.fastq 
rm -f $FQR3.unpaired.fastq
rm -f $FQR4.unpaired.fastq 
fi
##ONLY NEED TO DO THIS ONCE
##bwa index $HG19


#note the -M tag is important if you are to use picard
#-L 50, increases the clipping penalty -> discourages clipping by BWA-MEM, likely will increase indel sensitivity (especially insertions as already proven in a case) 
# while hurting SNP accuracy
#giving exact same read group ID for baserecalibration so can merge the bam files (note: you should only do this if the samples are on the same lane which is only true for miseq)
#baserecalibration has to be done on the same lane
#see https://toolshed.g2.bx.psu.edu/repository/display_tool?repository_id=c45d6c51a4fcfc6c&tool_config=%2Fsrv%2Ftoolshed%2Fmain%2Fvar%2Fdata%2Frepos%2F000%2Frepo_259%2Fpicard_AddOrReplaceReadGroups.xml&changeset_revision=ab1f60c26526
#for explanation on read groups
bwa mem -M -R "@RG\tID:${tID}\tLB:TrusightTumor\tPL:illumina\tPU:1\tSM:$OUTNAME$LIBA_SNUMBER" -L 50 -t $THREADS $HG19 $FQR1BWAMEM $FQR2BWAMEM > $DIR/${OUTNAME}.A.bwa-mem.sam &
bwa mem -M -R "@RG\tID:${tID}\tLB:TrusightTumor\tPL:illumina\tPU:1\tSM:$OUTNAME$LIBB_SNUMBER" -L 50 -t $THREADS $HG19 $FQR3BWAMEM $FQR4BWAMEM > $DIR/${OUTNAME}.B.bwa-mem.sam &
wait 




#for ALIGNER in bwa-mem bowtie2

for ALIGNER in bwa-mem
do



################WITHOUT BASERECALIBRATOR (FOUND IT NOT TO WORK WELL WITH TruSightTumor samples)
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SamFormatConverter INPUT=$DIR/${OUTNAME}.A.$ALIGNER.sam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.bam &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SamFormatConverter INPUT=$DIR/${OUTNAME}.B.$ALIGNER.sam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.bam &
wait
########################



if [ "$debugging" = false ]
then
####REMOVE SAM FILE TO SAVE SPACE
rm -f $DIR/${OUTNAME}.A.$ALIGNER.sam 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.sam 
############
fi

if [ "$validation" = true ]
then
#########################FOR DEBUGGING#################################
#getting bam index of the raw bam for debugging purpses
#NOW SORTING BECAUSE BUILD BAM INDEX NEEDS THIS
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.A.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.A_sorted.$ALIGNER.bam SORT_ORDER=coordinate &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.B.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.B_sorted.$ALIGNER.bam SORT_ORDER=coordinate &
wait
#NOW BUILDING BAM INDEX BECAUSE VARIANT CALLERS NEED THIS 
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.A_sorted.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.A_sorted.$ALIGNER.bai &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.B_sorted.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.B_sorted.$ALIGNER.bai &
wait
######################################################################
fi




##NOW I AM FILTERING THE READS USING SAMTOOLS
##https://www.biostars.org/p/56246/
## see website for explanation, only keeping proper pair and mapped read, because extension ligation method used this is the best, doesn't make sense if didn't map, then probe offtarget getting junk
#-L option limits it to amplicons that are only mapped to the region specified
#-q 1 ALSO GETS RID OF READS IN WHICH THE MATE IS MATCHED TO DIFFERENT CHROMOSOME, DON'T TRUST A READ LIKE THAT BECAUSE NOT HOW PANEL DESIGNED
#-F INT Do not output alignments with any bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].

#-f INT Only output alignments with all bits set in INT present in the FLAG field. INT can be specified in hex by beginning with `0x' (i.e. /^0x[0-9A-F]+/) or in octal by beginning with `0' (i.e. /^0[0-7]+/) [0].
#-f 2 means it is a proper pair - getting rid of bwa-mem can sometime not properly assign this value
#-F 8 the mate is unmapped, -F 4 means the first read is unmapped
$SAMTOOLS_DIR/samtools view -b -q 1 -F 4 -F 8 -L $AMPLICONS_BED_SAMTOOLS $DIR/${OUTNAME}.A.$ALIGNER.bam > $DIR/${OUTNAME}.A_filtered.$ALIGNER.bam &
$SAMTOOLS_DIR/samtools view -b -q 1 -F 4 -F 8 -L $AMPLICONS_BED_SAMTOOLS $DIR/${OUTNAME}.B.$ALIGNER.bam > $DIR/${OUTNAME}.B_filtered.$ALIGNER.bam &

wait

if [ "$debugging" = false ]
then
###TO FREE SPACE
rm -f $DIR/${OUTNAME}.A.$ALIGNER.bam 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.bam
rm -f $DIR/${OUTNAME}.A.$ALIGNER.pre.bam 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.pre.bam 
####
fi


###############
##CODE FOR NO CLIPPED READS
python $HOME_DIR/clipPrimers_both_ends.py -i $DIR/${OUTNAME}.A_filtered.$ALIGNER.bam -o $DIR/${OUTNAME}.A_filtered_clipped.$ALIGNER.bam -t "A" -m $manifest_libA &
python $HOME_DIR/clipPrimers_both_ends.py -i $DIR/${OUTNAME}.B_filtered.$ALIGNER.bam -o $DIR/${OUTNAME}.B_filtered_clipped.$ALIGNER.bam -t "B" -m $manifest_libB &
wait

#python $HOME_DIR/clipPrimers_manifest_speed.py -i $DIR/${OUTNAME}.A_filtered.$ALIGNER.bam -o $DIR/${OUTNAME}.A_filtered_clipped.$ALIGNER.bam -t "A" -m $manifest_libA &
#python $HOME_DIR/clipPrimers_manifest_speed.py -i $DIR/${OUTNAME}.B_filtered.$ALIGNER.bam -o $DIR/${OUTNAME}.B_filtered_clipped.$ALIGNER.bam -t "B" -m $manifest_libB &
#wait
######################

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar CleanSam INPUT=$DIR/${OUTNAME}.A_filtered_clipped.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.A_filtered_clipped_cleaned.$ALIGNER.bam &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar CleanSam INPUT=$DIR/${OUTNAME}.B_filtered_clipped.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.B_filtered_clipped_cleaned.$ALIGNER.bam &
wait

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar ValidateSamFile INPUT=$DIR/${OUTNAME}.A_filtered_clipped_cleaned.$ALIGNER.bam  &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar ValidateSamFile INPUT=$DIR/${OUTNAME}.B_filtered_clipped_cleaned.$ALIGNER.bam  &
wait

if [ "$debugging" = false ]
then
###TO FREE SPACE
rm -f $DIR/${OUTNAME}.A_filtered.$ALIGNER.bam 
rm -f $DIR/${OUTNAME}.B_filtered.$ALIGNER.bam 
rm -f $DIR/${OUTNAME}.A_filtered_clipped.$ALIGNER.bam 
rm -f $DIR/${OUTNAME}.B_filtered_clipped.$ALIGNER.bam 
####
fi
#java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${HG19} -I $DIR/${OUTNAME}.A.$ALIGNER.final.bam -o $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -dcov 100000 -minIndelFrac 0.03 -glm BOTH -L $AMPLICONS_BED_VARIANTS -dt NONE &
#wait

#NOW SORTING BECAUSE BUILD BAM INDEX NEEDS THIS
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.A_filtered_clipped_cleaned.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.bam SORT_ORDER=coordinate &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.B_filtered_clipped_cleaned.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.bam SORT_ORDER=coordinate &
wait

#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.A_filtered_clipped_cleaned.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.pre.bam SORT_ORDER=coordinate &
#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SortSam INPUT=$DIR/${OUTNAME}.B_filtered_clipped_cleaned.$ALIGNER.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.pre.bam SORT_ORDER=coordinate &
#wait
#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.A.$ALIGNER.pre.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.pre.bai &
#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.B.$ALIGNER.pre.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.pre.bai &
#wait
#java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T PrintReads -R ${HG19} -I $DIR/${OUTNAME}.A.$ALIGNER.pre.bam -o $DIR/${OUTNAME}.A.$ALIGNER.final.bam -filterNoBases -rf BadCigar &
#java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T PrintReads -R ${HG19} -I $DIR/${OUTNAME}.B.$ALIGNER.pre.bam -o $DIR/${OUTNAME}.B.$ALIGNER.final.bam -filterNoBases -rf BadCigar &
#wait


if [ "$debugging" = false ]
then
###TO FREE SPACE
rm -f $DIR/${OUTNAME}.A_filtered_clipped_cleaned.$ALIGNER.bam 
rm -f $DIR/${OUTNAME}.B_filtered_clipped_cleaned.$ALIGNER.bam 
#rm -f $DIR/${OUTNAME}.A.$ALIGNER.pre.bam
#rm -f $DIR/${OUTNAME}.B.$ALIGNER.pre.bam
#rm -f $DIR/${OUTNAME}.A.$ALIGNER.pre.bai
#rm -f $DIR/${OUTNAME}.B.$ALIGNER.pre.bai
####
fi


#NOW BUILDING BAM INDEX BECAUSE VARIANT CALLERS NEED THIS 

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.bai &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.bai &

wait





###DEBUGGING ONLY
#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SamFormatConverter INPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.sam &
#java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar SamFormatConverter INPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.sam &
#wait


#################VARSCAN##################################
#suggested code for pileup by varscan, samtools mpileup -f [reference sequence] [BAM file(s)] >myData.mpileup
#NOTE: NEED THE -d100000 to make the read depth > 8000 to 100,0000 which can happen in trusight tumor and potentially in the nextseq
#I am disabiling BAQ because I missed a snp next to a deletion in a case and it was due to this aggressive filter
#-S, to keep strand bias or -t SP since -S is depreciated
#-Q, minimum base quality to be considered
#-A, --count-orphans Do not skip anomalous read pairs in variant calling (need to do this as some reads are improperly flagged as not a proper pair)

#note specifying Q 20 here to maximize sensitivity
$SAMTOOLS_DIR/samtools mpileup -d 1000000 -B -A -f $HG19 $DIR/${OUTNAME}.A.$ALIGNER.final.bam -t SP -o $DIR/${OUTNAME}.A.$ALIGNER.mpileup &
$SAMTOOLS_DIR/samtools mpileup -d 1000000 -B -A -f $HG19 $DIR/${OUTNAME}.B.$ALIGNER.final.bam -t SP -o $DIR/${OUTNAME}.B.$ALIGNER.mpileup &
wait

#########
#for CoverageQC
$SAMTOOLS_DIR/samtools mpileup -d 1000000 -B -A -Q 20 -f $HG19 $DIR/${OUTNAME}.A.$ALIGNER.final.bam $DIR/${OUTNAME}.B.$ALIGNER.final.bam -t DP,DV,DPR,DP4,SP -u -v -o $DIR/${OUTNAME}.$ALIGNER.mpileup_all.noe.genome.vcf





java -Xmx${MEM}g -jar $IGVTOOLS_JAR index $DIR/${OUTNAME}.$ALIGNER.mpileup_all.noe.genome.vcf
#USAGE: java -jar VarScan.jar readcounts [pileup file] OPTIONS
#	pileup file - The SAMtools pileup file
	
#        OPTIONS:
#        --variants-file A list of variants at which to report readcounts
#        --output-file   Output file to contain the readcounts
#        --min-coverage  Minimum read depth at a position to make a call [8]
#        --min-base-qual Minimum base quality at a position to count a read [30]
#java -jar $VARSCAN_JAR readcounts $DIR/${OUTNAME}.A.$ALIGNER.mpileup --min-base-qual 20 --output-file $DIR/${OUTNAME}.A.$ALIGNER.readcounts &
#java -jar $VARSCAN_JAR readcounts $DIR/${OUTNAME}.B.$ALIGNER.mpileup --min-base-qual 20 --output-file $DIR/${OUTNAME}.B.$ALIGNER.readcounts &
#wait

######USAGE######mpileup2cns
#USAGE: java -jar VarScan.jar mpileup2cns [pileup file] OPTIONS
#        mpileup file - The SAMtools mpileup file

#        OPTIONS:
#        --min-coverage  Minimum read depth at a position to make a call [8]
#        --min-reads2    Minimum supporting reads at a position to call variants [2]
#        --min-avg-qual  Minimum base quality at a position to count a read [15]
#        --min-var-freq  Minimum variant allele frequency threshold [0.01]
#        --min-freq-for-hom      Minimum frequency to call homozygote [0.75]
#        --p-value       Default p-value threshold for calling variants [99e-02]
#        --strand-filter Ignore variants with >90% support on one strand [1]
#        --output-vcf    If set to 1, outputs in VCF format
#        --vcf-sample-list       For VCF output, a list of sample names in order, one per line
#        --variants      Report only variant (SNP/indel) positions [0]
######        


##############REAL CODE
##the code below makes a genome vcf

#You want the p-values to be lower rather than higher (so, 0.01 not 0.99) to decrease false positive, for sensitivity want it higher
#https://www.biostars.org/p/93153/
#if set high like .99 (which supposedly it sets it to this if p-value is NOT specified BUT BUGGY SO HAVE TO SPECIFIY) it does not calculate the p-value (this will likely give you the most results) here it will based calls solely on depth, although appears somewhat buggy and it still seems to be applying a p-value
#most specific .01, .05 will increase calls at risk of increaing false positive (ie 1% chance due to noise versus 5% chance due to noise)
#The variant P-value is computed by one-tailed Fisher's exact test (FET) in a two-by-two table, comparing the total number of reference-supporting reads and the total number of variant supporting reads (normal and tumor values are combined) to the expected distribution for a nonvariant position due to sequencing error (0.01%). For example, the expected read distribution for a nonvariant position with 500× coverage in each sample would be 999 reference-supporting reads, and one variant-supporting read due to sequencing error.
#significant or due to sequencing error (usually always have high p-value given lab cut off of 5% at 500X coverage; if wanted to call more would need a higher cutoff


######################

if [ "$debugging" = false ]
then
#need to turn of strand filter because it will throw out read instead of gives an SB in the filter portion of the vcf
java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.A.$ALIGNER.mpileup --strand-filter 0 --min-coverage 200 --min-reads2 10 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.1 --variants --output-vcf 1 >  $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre1.vcf &
java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.B.$ALIGNER.mpileup --strand-filter 0 --min-coverage 200 --min-reads2 10 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.1 --variants --output-vcf 1 >  $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre1.vcf &
#wait

###below code is with the p-value specified (NOTE HAVE TO SPECIFY buggy if don't, but didn't seem to make huge difference when compared cases or compared case where indel was missed
###can use if want to adjust, but validation on above parameters with p-vale 
#java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.A.$ALIGNER.mpileup --strand-filter 0 --min-coverage 200 --min-reads2 10 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.99 --variants --output-vcf 1 >  $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre1.vcf &
#java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.B.$ALIGNER.mpileup --strand-filter 0 --min-coverage 200 --min-reads2 10 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.99 --variants --output-vcf 1 >  $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre1.vcf &
wait


else

#################DEBUGGING CODE
###so can work with really small fastq files
#java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.A.$ALIGNER.mpileup --output-vcf 1 --strand-filter 0 --min-coverage 2 --min-reads2 1 --min-var-freq 0.025 --p-value 0.01 >  $DIR/${OUTNAME}.A.varscan.$ALIGNER.genome.vcf &
#java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.B.$ALIGNER.mpileup --output-vcf 1 --strand-filter 0 --min-coverage 2 --min-reads2 1 --min-var-freq 0.025 --p-value 0.01 >  $DIR/${OUTNAME}.B.varscan.$ALIGNER.genome.vcf &
#wait
##the code below makes a vcf with both snps and indels
java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.A.$ALIGNER.mpileup --strand-filter 0 --min-coverage 2 --min-reads2 1 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.99 --variants --output-vcf 1 >  $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre1.vcf &
java -jar $VARSCAN_JAR mpileup2cns $DIR/${OUTNAME}.B.$ALIGNER.mpileup --strand-filter 0 --min-coverage 2 --min-reads2 1 --min-avg-qual 20 --min-var-freq 0.025 --min-freq-for-hom 0.75 --p-value 0.99 --variants --output-vcf 1 >  $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre1.vcf &


wait
########################


fi

#breakup multiallele sites in case made because it makes future processing difficult
$VCFLIB_HOME/vcfbreakmulti $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre1.vcf > $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf
$VCFLIB_HOME/vcfbreakmulti $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre1.vcf > $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf



bgzip -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf &
bgzip -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf &
wait
tabix -f -p vcf $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf.gz &
tabix -f -p vcf $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf.gz &
wait

#applying strand bias to reads, if 90% of the reads come from one strand then there is a strand bias
$BCFTOOLS_DIR/bcftools filter -R $AMPLICONS_BED_SAMTOOLS -s "SB" $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf.gz -e "((FORMAT/RDF)+(FORMAT/ADF) > 0 && (FORMAT/ADF)/((FORMAT/RDF)+(FORMAT/ADF)) < .1 * (FORMAT/AD)/((FORMAT/RD)+(FORMAT/AD))) || ((FORMAT/RDR)+(FORMAT/ADR)>0 && (FORMAT/ADR)/((FORMAT/RDR)+(FORMAT/ADR)) < .1 * (FORMAT/AD)/((FORMAT/RD)+(FORMAT/AD)))" > $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf &
$BCFTOOLS_DIR/bcftools filter -R $AMPLICONS_BED_SAMTOOLS -s "SB" $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf.gz -e "((FORMAT/RDF)+(FORMAT/ADF) > 0 && (FORMAT/ADF)/((FORMAT/RDF)+(FORMAT/ADF)) < .1 * (FORMAT/AD)/((FORMAT/RD)+(FORMAT/AD))) || ((FORMAT/RDR)+(FORMAT/ADR)>0 && (FORMAT/ADR)/((FORMAT/RDR)+(FORMAT/ADR)) < .1 * (FORMAT/AD)/((FORMAT/RD)+(FORMAT/AD)))" > $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf &
wait


if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf.gz
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf.gz
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre1.vcf
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.pre2.vcf.idx
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre1.vcf
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.pre2.vcf.idx
fi

bgzip -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf &
bgzip -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf &

wait

tabix -f -p vcf $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf.gz &
tabix -f -p vcf $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf.gz &

wait


#see documentation https://samtools.github.io/bcftools/bcftools.html#merge 
#-m none   - this setting means that no multiallelic sites will be used which is EXACTLY what I want
$BCFTOOLS_DIR/bcftools merge $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf.gz $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf.gz --force-samples -m none > $DIR/${OUTNAME}.varscan_all.$ALIGNER.vcf &
wait


if [ "$debugging" = false ]
then
####TO FREE UP SPACE
rm -f $DIR/${OUTNAME}.varscan_all.$ALIGNER.pre.vcf
rm -f $DIR/${OUTNAME}.A.$ALIGNER.mpileup
rm -f $DIR/${OUTNAME}.B.$ALIGNER.mpileup
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf.gz
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf.gz
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf.idx
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf.idx
rm -f $DIR/${OUTNAME}.A.varscan.$ALIGNER.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.B.varscan.$ALIGNER.genome.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.varscan_all.$ALIGNER.genome_temp.vcf
fi

######



##########################################################


####################################################FREEBAYES########################################

#trying freebayes,  freebayes appears to work better for small snps 
if [ "$debugging" = false ]
then
###############REGULAR_PARAMETERS
#note: min-base-quality default is 0, changed to 20 , ie don't count a base <20
#-F --min-alternate-fraction set to 0.025 Require at least this fraction of observations supporting an alternate allele within a single individual in the in order to evaluate the position.  default: 0.2
#rest are freebayes defaults
#-m --min-mapping-quality  Exclude alignments from analysis if they have a mapping quality less than Q.  default: 1
# -R --min-supporting-allele-qsum Q Consider any allele in which the sum of qualities of supporting observations is at least Q.  default: 0
# -e --read-indel-limit N Exclude reads with more than N separate gaps. default: ~unbounded
$FREEBAYES_HOME/freebayes --bam $DIR/${OUTNAME}.A.$ALIGNER.final.bam --fasta-reference $HG19 --vcf $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre1.vcf --targets $AMPLICONS_BED_SAMTOOLS --min-mapping-quality "1" --min-base-quality "20" --min-supporting-allele-qsum "0" --min-supporting-mapping-qsum "0" --read-indel-limit "1000" --min-alternate-fraction "0.03" --min-alternate-count "6" --min-alternate-qsum "0" --min-alternate-total "6" --min-coverage "200" &
$FREEBAYES_HOME/freebayes --bam $DIR/${OUTNAME}.B.$ALIGNER.final.bam --fasta-reference $HG19 --vcf $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre1.vcf --targets $AMPLICONS_BED_SAMTOOLS --min-mapping-quality "1" --min-base-quality "20" --min-supporting-allele-qsum "0" --min-supporting-mapping-qsum "0" --read-indel-limit "1000" --min-alternate-fraction "0.03" --min-alternate-count "6" --min-alternate-qsum "0" --min-alternate-total "6" --min-coverage "200" &
####################################
wait


else

$FREEBAYES_HOME/freebayes --bam $DIR/${OUTNAME}.A.$ALIGNER.final.bam --fasta-reference $HG19 --vcf $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre1.vcf --targets $AMPLICONS_BED_SAMTOOLS --min-mapping-quality "1" --min-base-quality "20" --min-supporting-allele-qsum "0" --min-supporting-mapping-qsum "0" --read-indel-limit "1000" --min-alternate-fraction "0.03" --min-alternate-count "6" --min-alternate-qsum "0" --min-alternate-total "1" --min-coverage "2" &
$FREEBAYES_HOME/freebayes --bam $DIR/${OUTNAME}.B.$ALIGNER.final.bam --fasta-reference $HG19 --vcf $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre1.vcf --targets $AMPLICONS_BED_SAMTOOLS --min-mapping-quality "1" --min-base-quality "20" --min-supporting-allele-qsum "0" --min-supporting-mapping-qsum "0" --read-indel-limit "1000" --min-alternate-fraction "0.03" --min-alternate-count "6" --min-alternate-qsum "0" --min-alternate-total "1" --min-coverage "2" &
####################################
wait


fi


#freebayes, makes multiple alt allles, trying to breakup because it makes future processing difficult
$VCFLIB_HOME/vcfbreakmulti $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre1.vcf > $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf
$VCFLIB_HOME/vcfbreakmulti $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre1.vcf > $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf

bgzip -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf &
bgzip -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf &
wait

tabix -f -p vcf $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf.gz &
tabix -f -p vcf $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf.gz &
wait

#applying strand bias to reads, if 90% of the reads come from one strand then there is a strand bias
$BCFTOOLS_DIR/bcftools filter -s "SB" $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf.gz -e "((INFO/SAF)+(INFO/SRF)>0 && (INFO/SAF)/((INFO/SAF)+(INFO/SRF)) < .1 * ((INFO/SAF)+(INFO/SAR))/(INFO/DP)) || ((INFO/SAR)+(INFO/SRR) > 0  && (INFO/SAR)/((INFO/SAR)+(INFO/SRR)) < .1 * ((INFO/SAF)+(INFO/SAR))/(INFO/DP))" > $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf &
$BCFTOOLS_DIR/bcftools filter -s "SB" $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf.gz -e "((INFO/SAF)+(INFO/SRF)>0 && (INFO/SAF)/((INFO/SAF)+(INFO/SRF)) < .1 * ((INFO/SAF)+(INFO/SAR))/(INFO/DP)) || ((INFO/SAR)+(INFO/SRR) > 0  && (INFO/SAR)/((INFO/SAR)+(INFO/SRR)) < .1 * ((INFO/SAF)+(INFO/SAR))/(INFO/DP))" > $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf &
wait



if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf.gz
rm -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf.gz
rm -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre2.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre2.vcf.gz.tbi
rm -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.pre1.vcf
rm -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.pre1.vcf
fi

echo "Freebayes is going"
wait



bgzip -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf &
bgzip -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf &
wait

tabix -f -p vcf $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf.gz &
tabix -f -p vcf $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf.gz &
wait




#see documentation https://samtools.github.io/bcftools/bcftools.html#merge 
#-m none   - this setting means that no multiallelic sites will be made which is EXACTLY what I want, other tools like vcftools could not do this
$BCFTOOLS_DIR/bcftools merge $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf.gz $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf.gz --force-samples -m none > $DIR/${OUTNAME}.freebayes_all.$ALIGNER.vcf &
wait



if [ "$debugging" = false ]
then
###TO FREE UP SPACE
rm -f $DIR/${OUTNAME}.freebayes_all.$ALIGNER.pre.vcf
rm -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf.gz 
rm -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf.gz
rm -f $DIR/${OUTNAME}.A.freebayes.$ALIGNER.vcf.gz.tbi 
rm -f $DIR/${OUTNAME}.B.freebayes.$ALIGNER.vcf.gz.tbi 
####
fi

#########################################################################



##################################GATK########################################


java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${HG19} -I $DIR/${OUTNAME}.A.$ALIGNER.final.bam -o $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -dcov 1000000 -minIndelFrac 0.03 -mbq 20 -glm BOTH -L $AMPLICONS_BED_SAMTOOLS -dt NONE &
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T UnifiedGenotyper -R ${HG19} -I $DIR/${OUTNAME}.B.$ALIGNER.final.bam -o $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -dcov 1000000 -minIndelFrac 0.03 -mbq 20 -glm BOTH -L $AMPLICONS_BED_SAMTOOLS -dt NONE &
wait


if [ "$debugging" = false ]
then
#############REGULAR PARAMETERS
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants -R ${HG19} -o $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --variant $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -select "DP > 200" &
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants -R ${HG19} -o $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --variant $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -select "DP > 200" &
############################

wait


else
####################################Debugging only
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants -R ${HG19} -o $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --variant $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -select "DP > 2" &
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T SelectVariants -R ${HG19} -o $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --variant $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf -select "DP > 2" &
############################
wait

fi


#adding a strand bias filter
#can't do strand filter similar to the others because the vcf doesn't give depth read forward and read reverse
##INFO=<ID=FS,Number=1,Type=Float,Description="Phred-scaled p-value using Fisher's exact test to detect strand bias">
# #where the higher FS the higher probability of strand bias
#https://www.broadinstitute.org/gatk/events/slides/1506/GATKwr8-D-4-Manual_filtering.pdf; for snps suggest > 60.0 so will set as such since never should trust strand bias
#with amplicon sequencing when it comes to indels
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T VariantFiltration -R ${HG19} -o $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf --variant $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --filterExpression "FS > 60.0" --filterName "SB" &
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T VariantFiltration -R ${HG19} -o $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf --variant $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf --filterExpression "FS > 60.0" --filterName "SB" &
wait


if [ "$debugging" = false ]
then
###TO FREE UP SPACE
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf.idx
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.pre.vcf.idx
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf.idx 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.temp.vcf.idx 
####
fi


#Not using combine variants as there is no way I specify not to make multiallelic variants and can't use vcfbreakmulti as it does not properly break one allele percent and the other in multisample vcf
#java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -R $HG19  -T CombineVariants --variant $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf --variant $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf -o $DIR/${OUTNAME}.gatk_all.temp.$ALIGNER.vcf -genotypeMergeOptions UNIQUIFY &
#wait

bgzip -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf &
bgzip -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf &
wait

tabix -f -p vcf $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz &
tabix -f -p vcf $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz &
wait

$BCFTOOLS_DIR/bcftools merge $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz --force-samples -m none > $DIR/${OUTNAME}.gatk_all.$ALIGNER.vcf &
wait


if [ "$debugging" = false ]
then
###TO FREE UP SPACE
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz.tbi 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.gz.tbi 
rm -f $DIR/${OUTNAME}.gatk_all.temp.$ALIGNER.vcf
rm -f $DIR/${OUTNAME}.gatk_all.temp.$ALIGNER.vcf.idx
rm -f $DIR/${OUTNAME}.A.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.idx 
rm -f $DIR/${OUTNAME}.B.$ALIGNER.filtered_clipped_sorted.gatk.indel.snp.vcf.idx
####
fi

##################################################################################################

#NOTE: GETTING RID OF THE -onlyTr parameter, ie -onlyTr $HOME_DIR/consensus_CMP26_for_snpeff.txt, now will annotate for all transcripts
java -Xmx${MEM}g -jar $SNPEFF_DIR/snpEff.jar hg19 -formatEff $DIR/${OUTNAME}.gatk_all.$ALIGNER.vcf  > $DIR/${OUTNAME}.gatk_all.$ALIGNER.pre.vcf &
java -Xmx${MEM}g -jar $SNPEFF_DIR/snpEff.jar hg19 -formatEff $DIR/${OUTNAME}.freebayes_all.$ALIGNER.vcf  > $DIR/${OUTNAME}.freebayes_all.$ALIGNER.pre.vcf &
java -Xmx${MEM}g -jar $SNPEFF_DIR/snpEff.jar hg19 -formatEff $DIR/${OUTNAME}.varscan_all.$ALIGNER.vcf  > $DIR/${OUTNAME}.varscan_all.$ALIGNER.pre.vcf &
wait

##need to remove contigs, as bcftools adds this to the vcf it makes, however, pyVcf which modify_annovar uses does not like the contig feature and crashes if a vcf has
##this saying for example : "One of the contig lines is malformed: %s" % contig_string)
##SyntaxError: One of the contig lines is malformed: ##contig=<ID=chr2>
##work around for right now is to just get rid of as likely not necessary as before bcftools was being used this was not even present 
perl -pe 's/^##contig.*\n//' $DIR/${OUTNAME}.gatk_all.$ALIGNER.pre.vcf > $DIR/${OUTNAME}.gatk_all.$ALIGNER.final.vcf
perl -pe 's/^##contig.*\n//' $DIR/${OUTNAME}.freebayes_all.$ALIGNER.pre.vcf > $DIR/${OUTNAME}.freebayes_all.$ALIGNER.final.vcf
perl -pe 's/^##contig.*\n//' $DIR/${OUTNAME}.varscan_all.$ALIGNER.pre.vcf > $DIR/${OUTNAME}.varscan_all.$ALIGNER.final.vcf

if [ "$debugging" = false ]
then
###TO FREE SPACE
#rm -f $DIR/${OUTNAME}.vardict_all.$ALIGNER.vcf
rm -f $DIR/${OUTNAME}.gatk_all.$ALIGNER.pre.vcf
rm -f $DIR/${OUTNAME}.freebayes_all.$ALIGNER.pre.vcf
rm -f $DIR/${OUTNAME}.varscan_all.$ALIGNER.pre.vcf
rm -f $DIR/${OUTNAME}.gatk_all.$ALIGNER.vcf 
rm -f $DIR/${OUTNAME}.freebayes_all.$ALIGNER.vcf
rm -f $DIR/${OUTNAME}.varscan_all.$ALIGNER.vcf 
####
fi


##NOW CHECKING IF WE CAN FIND THE TRUSIGHT VCF
if [ "$VCFORIG" != "" ]
then


cp $VCFORIG $DIR/${OUTNAME}.illumina.pre.vcf
java -Xmx${MEM}g -jar $SNPEFF_DIR/snpEff.jar hg19 -formatEff $DIR/${OUTNAME}.illumina.pre.vcf  > $DIR/${OUTNAME}.illumina.final.vcf

if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.illumina.pre.vcf
fi

fi	






#DOING ANNOVAR ON JUST THE GATK, FREEBAYES AND VARSCAN
perl ${ANNOVAR_HOME}/convert2annovar.pl $DIR/${OUTNAME}.gatk_all.$ALIGNER.final.vcf --includeinfo -format vcf4 -allsample -withfreq --outfile $DIR/${OUTNAME}.gatk_all.$ALIGNER.annovar &
perl ${ANNOVAR_HOME}/convert2annovar.pl $DIR/${OUTNAME}.freebayes_all.$ALIGNER.final.vcf --includeinfo -format vcf4 -allsample -withfreq --outfile $DIR/${OUTNAME}.freebayes_all.$ALIGNER.annovar &
perl ${ANNOVAR_HOME}/convert2annovar.pl $DIR/${OUTNAME}.varscan_all.$ALIGNER.final.vcf --includeinfo -format vcf4 -allsample -withfreq --outfile $DIR/${OUTNAME}.varscan_all.$ALIGNER.annovar &
wait


##NOW CHECKING IF WE CAN FIND THE TRUSIGHT VCF
if [ "$VCFORIG" != "" ]
then


perl ${ANNOVAR_HOME}/convert2annovar.pl $DIR/${OUTNAME}.illumina.final.vcf --includeinfo -format vcf4 -allsample -withfreq --outfile $DIR/${OUTNAME}.illumina.miseq.annovar &
wait


fi


perl ${ANNOVAR_HOME}/table_annovar_splicing_5.pl -buildver hg19 -out $DIR/${OUTNAME}.gatk_all.$ALIGNER -otherinfo -remove -protocol refGene,1000g2015aug_all,snp138,cosmic70,clinvar_20150629,exac03 -operation g,f,f,f,f,f $DIR/${OUTNAME}.gatk_all.$ALIGNER.annovar ${ANNOVAR_DB_HOME} &
perl ${ANNOVAR_HOME}/table_annovar_splicing_5.pl -buildver hg19 -out $DIR/${OUTNAME}.freebayes_all.$ALIGNER -otherinfo -remove -protocol refGene,1000g2015aug_all,snp138,cosmic70,clinvar_20150629,exac03 -operation g,f,f,f,f,f $DIR/${OUTNAME}.freebayes_all.$ALIGNER.annovar ${ANNOVAR_DB_HOME} &
perl ${ANNOVAR_HOME}/table_annovar_splicing_5.pl -buildver hg19 -out $DIR/${OUTNAME}.varscan_all.$ALIGNER -otherinfo -remove -protocol refGene,1000g2015aug_all,snp138,cosmic70,clinvar_20150629,exac03 -operation g,f,f,f,f,f $DIR/${OUTNAME}.varscan_all.$ALIGNER.annovar ${ANNOVAR_DB_HOME} &

wait
 

##NOW CHECKING IF WE CAN FIND THE TRUSIGHT VCF
if [ "$VCFORIG" != "" ]
then


perl ${ANNOVAR_HOME}/table_annovar_splicing_5.pl -buildver hg19 -out $DIR/${OUTNAME}.illumina.miseq -otherinfo -remove -protocol refGene,1000g2015aug_all,snp138,cosmic70,clinvar_20150629,exac03 -operation g,f,f,f,f,f $DIR/${OUTNAME}.illumina.miseq.annovar ${ANNOVAR_DB_HOME} 
python $HOME_DIR/modify_annovar.py -i $DIR/${OUTNAME}.illumina.miseq.hg19_multianno.txt -o $DIR/${OUTNAME}.illumina.miseq.final.hg19_multianno.txt -v $DIR/${OUTNAME}.illumina.final.vcf -p "MiSeq->MiSeq" -c $consensus_RefSeq_CMP26 

if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.illumina.miseq.annovar &
rm -f $DIR/${OUTNAME}.illumina.miseq.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.illumina.miseq.hg19_multianno.txt2.txt &
fi


wait


#######FOR VALIDATION PURPOSES

if [ "$validation" = true ]
then
BAM_A=$DIR/${OUTNAME}_S${LIBA_SNUMBER}
BAM_B=$DIR/${OUTNAME}_S${LIBB_SNUMBER}

echo $BAM_A
echo $BAM_B

cp ${FQRDIR}Alignment/${OUTNAME}_S${LIBA_SNUMBER}.bam $BAM_A.bam &
cp ${FQRDIR}Alignment/${OUTNAME}_S${LIBB_SNUMBER}.bam $BAM_B.bam &
wait


java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups INPUT=$BAM_A.bam RGID=1 RGLB=TruSightTumor RGPL=illumina RGPU=unknown RGSM=total OUTPUT=$BAM_A.newheader.bam VALIDATION_STRINGENCY=SILENT &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups INPUT=$BAM_B.bam RGID=1 RGLB=TruSightTumor RGPL=illumina RGPU=unknown RGSM=total OUTPUT=$BAM_B.newheader.bam VALIDATION_STRINGENCY=SILENT &
wait

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$BAM_A.newheader.bam OUTPUT=$BAM_A.newheader.bai VALIDATION_STRINGENCY=SILENT &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$BAM_B.newheader.bam OUTPUT=$BAM_B.newheader.bai VALIDATION_STRINGENCY=SILENT &
wait

if [ "$debugging" = false ]
then
rm -f $BAM_A.bam
rm -f $BAM_B.bam
fi


###NOTE: TO CALCULATE THE DEPTH ILLUMINA USES A MINIMUM BASE QUALITY OF 20 (maybe, haven't tested 18 though)
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage -R ${HG19} -o $DIR/${OUTNAME}.illumina.gatk.doc -I $BAM_A.newheader.bam -I $BAM_B.newheader.bam -L $AMPLICONS_BED_VARIANTS --nBins 4999 --stop 5000 -ct 500 -ct 1000 --minBaseQuality 20


if [ "$debugging" = false ]
then
rm -f $BAM_A.newheader.bam
rm -f $BAM_B.newheader.bam
rm -f $BAM_A.newheader.bai
rm -f $BAM_B.newheader.bai
fi

python $HOME_DIR/append_bed_name.py -b $AMPLICONS_BED_VARIANTS_NAMED -i $DIR/${OUTNAME}.illumina.gatk.doc.sample_interval_summary -o $DIR/${OUTNAME}.illumina.gatk.doc.final.interval_summary.txt

if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_cumulative_coverage_proportions
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_cumulative_coverage_counts
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_interval_summary
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_interval_statistics
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_statistics
rm -f $DIR/${OUTNAME}.illumina.gatk.doc.sample_summary
rm -f $DIR/${OUTNAME}.illumina.gatk.doc
fi

fi # if validation check

fi #if vcf original check



if [ "$validation" = true ]
then
#now calculating the custom pipeline

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups INPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.bam RGID=1 RGLB=TruSightTumor RGPL=$ALIGNER RGPU=unknown RGSM=total OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bam VALIDATION_STRINGENCY=SILENT &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar AddOrReplaceReadGroups INPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.bam RGID=1 RGLB=TruSightTumor RGPL=$ALIGNER RGPU=unknown RGSM=total OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bam VALIDATION_STRINGENCY=SILENT &
wait

java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bam OUTPUT=$DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bai VALIDATION_STRINGENCY=SILENT &
java -Xmx${MEM}g -jar $PICARD_DIR/picard.jar BuildBamIndex INPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bam OUTPUT=$DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bai VALIDATION_STRINGENCY=SILENT &
wait


###NOTE: TO CALCULATE THE DEPTH ILLUMINA USES A MINIMUM BASE QUALITY OF 20
java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage -R ${HG19} -o $DIR/${OUTNAME}.$ALIGNER.gatk.doc -I $DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bam -I $DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bam -L $AMPLICONS_BED_VARIANTS --nBins 4999 --stop 5000 -ct 500 -ct 1000 --minBaseQuality 20

java -Xmx${MEM}g -jar $GATK_DIR/GenomeAnalysisTK.jar -T DepthOfCoverage -R ${HG19} -o $DIR/${OUTNAME}.$ALIGNER.gatk.doc.genome.txt -I $DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bam -I $DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bam -L $AMPLICONS_BED_SAMTOOLS --omitPerSampleStats --omitIntervalStatistics --nBins 4999 --stop 5000 -ct 500 -ct 1000 --minBaseQuality 20


if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bam
rm -f $DIR/${OUTNAME}.B.$ALIGNER.final.newheader.bam
rm -f $DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bai
rm -f $DIR/${OUTNAME}.A.$ALIGNER.final.newheader.bai
fi


python $HOME_DIR/append_bed_name.py -b $AMPLICONS_BED_VARIANTS_NAMED -i $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_interval_summary -o $DIR/${OUTNAME}.$ALIGNER.gatk.doc.final.interval_summary.txt


fi


if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_cumulative_coverage_proportions
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_cumulative_coverage_counts
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_interval_summary
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_interval_statistics
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_statistics
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.sample_summary
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.genome.txt.sample_cumulative_coverage_counts
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc.genome.txt.sample_cumulative_coverage_proportions
rm -f $DIR/${OUTNAME}.$ALIGNER.gatk.doc
fi








if [ "$debugging" = false ]
then
rm -f $DIR/${OUTNAME}.gatk_all.$ALIGNER.annovar &
rm -f $DIR/${OUTNAME}.freebayes_all.$ALIGNER.annovar &
rm -f $DIR/${OUTNAME}.varscan_all.$ALIGNER.annovar &
fi

wait




python $HOME_DIR/modify_annovar.py -i $DIR/${OUTNAME}.gatk_all.$ALIGNER.hg19_multianno.txt -o $DIR/${OUTNAME}.gatk_all.$ALIGNER.final.hg19_multianno.txt -v $DIR/${OUTNAME}.gatk_all.$ALIGNER.final.vcf -p "$ALIGNER->gatk_unifiedgenotyper" -c $consensus_RefSeq_CMP26 &
python $HOME_DIR/modify_annovar.py -i $DIR/${OUTNAME}.freebayes_all.$ALIGNER.hg19_multianno.txt -o $DIR/${OUTNAME}.freebayes_all.$ALIGNER.final.hg19_multianno.txt -v $DIR/${OUTNAME}.freebayes_all.$ALIGNER.final.vcf -p "$ALIGNER->freebayes" -c $consensus_RefSeq_CMP26 &
python $HOME_DIR/modify_annovar.py -i $DIR/${OUTNAME}.varscan_all.$ALIGNER.hg19_multianno.txt -o $DIR/${OUTNAME}.varscan_all.$ALIGNER.final.hg19_multianno.txt -v $DIR/${OUTNAME}.varscan_all.$ALIGNER.final.vcf -p "$ALIGNER->varscan" -c $consensus_RefSeq_CMP26 &
wait



python $HOME_DIR/merge_annovar.py -i $DIR/${OUTNAME}.gatk_all.$ALIGNER.final.hg19_multianno.txt -i $DIR/${OUTNAME}.freebayes_all.$ALIGNER.final.hg19_multianno.txt -i $DIR/${OUTNAME}.varscan_all.$ALIGNER.final.hg19_multianno.txt -o $DIR/${OUTNAME}.$ALIGNER.allvariantcallers.final.hg19_multianno.txt &

wait 

if [ "$debugging" = false ]
then
#####TO-FREE UP SPACE#######
rm -f $DIR/${OUTNAME}.gatk_all.bwa-mem.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bwa-mem.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bwa-mem.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bowtie2.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bowtie2.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bowtie2.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bwa-mem.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bwa-mem.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bwa-mem.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bowtie2.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bowtie2.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bowtie2.final.hg19_multianno.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bwa-mem.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bwa-mem.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bwa-mem.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bowtie2.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bowtie2.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bowtie2.final.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bwa-mem.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bwa-mem.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bwa-mem.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.gatk_all.bowtie2.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.freebayes_all.bowtie2.hg19_multianno.txt2.txt &
rm -f $DIR/${OUTNAME}.varscan_all.bowtie2.hg19_multianno.txt2.txt &
fi

wait
############



####DEBUGGING_ONLY_FOR_VALIDATION_PURPOSES
###to prevent a weird Can't connect to X11 window server using 'localhost:10.0' as the value of the DISPLAY variable
###have to add -Djava.awt.headless=true at runtime (ie when running the java code)
if [ "$validation" = true ]
then
export DISPLAY=:0.0
java -Djava.awt.headless=true -Xmx${MEM}g -jar $COVERAGEQC_DIR/coverageQc.jar $DIR/${OUTNAME}.$ALIGNER.mpileup_all.noe.genome.vcf $COVERAGEQC_DIR/cancer_panel_26.20151130.exons.bed $COVERAGEQC_DIR/cancer_panel_26.20160216.amplicons.bed $COVERAGEQC_DIR/20150702_aligners.csv $COVERAGEQC_DIR/20150810_filestolookfor.csv $COVERAGEQC_DIR/20150705_genes_excluding.csv $COVERAGEQC_DIR/'Do not call_26.20160216.list.xlsx'
fi
done


 

if [ "$debugging" = false ]
then
##removing the temporary fastqs to save space, when debugging, may want to comment this out
rm -f $FQR1
rm -f $FQR2
rm -f $FQR3
rm -f $FQR4
rm -f $FQR1.cut.fastq
rm -f $FQR2.cut.fastq
rm -f $FQR3.cut.fastq
rm -f $FQR4.cut.fastq
fi



echo "DONE with $FQR1"


