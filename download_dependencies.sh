#!/bin/sh



DIR=$PWD

#download dependencies calls for for debian environments, calls will need to be modified if
#working with red-hat

#if not done (ubuntu version, see readme for red-hat equivalent
#sudo apt-get install git #ubuntu
#sudo apt-get install unzip #ubuntu
#sudo apt-get install openjdk-8-jdk #Ubuntu
#sudo apt-get install wget #Ubuntu
#sudo apt-get install gcc #ubuntu
#sudo yum install python-dev #Ubuntu
#sudo apt-get install zlib1g-dev #ubuntu
#sudo apt-get install g++ #ubuntu
#sudo apt-get install python-pip
#sudo pip install biopython==1.66
#sudo pip install pysam==0.8.4
#sudo pip install pyvcf==0.6.7
#sudo pip install pandas==0.16.2
#sudo pip install regex==2015.3.18
#sudo apt-get install zenity
#sudo apt-get install xterm xorg dbus #ubuntu
#sudo apt-get install alien dpkg-dev debhelper build-essential #needed for unbuntu
#must manually download ANNOVAR and GATK after accepting license aggrement

#install freebayes
#git clone --recursive git://github.com/ekg/freebayes.git
#cd $DIR/freebayes
#make

cd $DIR

#install bedtools
#wget https://github.com/arq5x/bedtools2/releases/download/v2.19.1/bedtools-2.19.1.tar.gz

#tar -zxvf bedtools-2.19.1.tar.gz
#mv bedtools2-2.19.1 bedtools2
#rm bedtools-2.19.1.tar.gz
#cd $DIR/bedtools
#make

cd $DIR


#install samtools
#wget https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2
#tar -vxjf samtools-1.2.tar.bz2
#rm samtools-1.2.tar.bz2
#cd samtools-1.2
#make

cd $DIR

#install bcftools
#wget https://github.com/samtools/bcftools/releases/download/1.2/bcftools-1.2.tar.bz2
#tar -vxjf bcftools-1.2.tar.bz2
#rm bcftools-1.2.tar.bz2
#cd bcftools-1.2
#make

cd $DIR

#install bcl2fastq
#wget 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.zip'
#unzip bcl2fastq2-v2.17.1.14-Linux-x86_64.zip
#sudo alien bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm
#sudo dpkg -i bcl2fastq2_0v2.17.1.14-2_amd64.deb
rm bcl2fastq2_0v2.17.1.14-2_amd64.deb
rm bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm
rm bcl2fastq2-v2.17.1.14-Linux-x86_64.zip

cd $DIR

#install htslib
#wget https://github.com/samtools/htslib/releases/download/1.2.1/htslib-1.2.1.tar.bz2
#tar -vxjf htslib-1.2.1.tar.bz2
#rm htslib-1.2.1.tar.bz2
#cd htslib-1.2.1
#make

cd $DIR


#install FASTQC
wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.3.zip
unzip fastqc_v0.11.3.zip
rm fastqc_v0.11.3.zip

cd $DIR

#install IGVTools 2.3.57 
wget http://data.broadinstitute.org/igv/projects/downloads/igvtools_2.3.57.zip
unzip igvtools_2.3.57.zip
rm igvtools_2.3.57.zip

#install picard
wget https://github.com/broadinstitute/picard/releases/download/2.1.0/picard-tools-2.1.0.zip
unzip picard-tools-2.1.0.zip
mv picard-tools-2.1.0 picard
rm picard-tools-2.1.0.zip

#install snpeff 4.1g 2015-05-17 
wget https://sourceforge.net/projects/snpeff/files/snpEff_v4_1g_core.zip/download
unzip download
rm download
cd snpeff
mkdir data 
java -jar snpEff.jar download hg19

cd $DIR

#install varscan v2.4.0
mkdir VarScan
cd VarScan
wget https://github.com/dkoboldt/varscan/raw/master/VarScan.v2.4.0.jar

cd $DIR

#install CoverageQC
git clone https://github.com/ghsmith/coverageQc -b agnostic_v2
mkdir Alternative_QC_Reporter
cd Alternative_QC_Reporter
mkdir CoverageQC
cd $DIR
cp coverageQc/dist/* Alternative_QC_Reporter/CoverageQC
cp -r coverageQC/lib/ Alternative_QC_Reporter/CoverageQC

#download hg19
mkdir hg19
cd hg19
mkdir WholeGenomeFASTA
cd WholeGenomeFASTA
cd ./hg19/WholeGenomeFASTA
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrM.fa.gz' -O chrM.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr1.fa.gz' -O chr1.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr2.fa.gz' -O chr2.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr3.fa.gz' -O chr3.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr4.fa.gz' -O chr4.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr5.fa.gz' -O chr5.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr6.fa.gz' -O chr6.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr7.fa.gz' -O chr7.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr8.fa.gz' -O chr8.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr9.fa.gz' -O chr9.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr10.fa.gz' -O chr10.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr11.fa.gz' -O chr11.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr12.fa.gz' -O chr12.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr13.fa.gz' -O chr13.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr14.fa.gz' -O chr14.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr15.fa.gz' -O chr15.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr16.fa.gz' -O chr16.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr17.fa.gz' -O chr17.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr18.fa.gz' -O chr18.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr19.fa.gz' -O chr19.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr20.fa.gz' -O chr20.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr21.fa.gz' -O chr21.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chr22.fa.gz' -O chr22.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrX.fa.gz' -O chrX.fa.gz
wget --timestamping 'ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/chromosomes/chrY.fa.gz' -O chrY.fa.gz
gunzip ./chr*
cat chrM.fa chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa > genome.fa
rm -f chr*

#install bwa
wget https://sourceforge.net/projects/bio-bwa/files/bwa-0.7.10.tar.bz2/download
tar -vxjf download
rm download
cd bwa-0.7.10

cd $DIR


#index hg19
./bwa-0.7.10/bwa index ./hg19/WholeGenomeFASTA/genome.fa

#make picard index
java -jar $DIR/picard/picard.jar  CreateSequenceDictionary R=$DIR/hg19/WholeGenomeFASTA/genome.fa O=$DIR/hg19/WholeGenomeFASTA/genome.dict &

#install Trimmomatic
wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip
unzip Trimmomatic-0.33.zip

#install vcflib 
#already installed with freebayes

#download ANNOVAR databases NOTE: this assumes ANNOVAR was downloaded
cd annovar
rc=$?
if [[ $rc != 0 ]]; then 
echo "annovar needs to be installed and the annovar directory placed in this folder" 
exit $rc
fi

mkdir humandb
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar refGene humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar clinvar_20150629 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar cosmic70 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar exac03 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar snp138 humandb/
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar 1000g2015aug humandb/




