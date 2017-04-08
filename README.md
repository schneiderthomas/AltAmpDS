# AltAmpDS
Alternative Bioinformatic Pipeline for AmpliconDS

This program is designed to run through a NextSeq or MiSeq run directory looking for 
fastq files located in ${current directory or specified directory}/Data/Intensities/BaseCalls/ 

The main script file to run is runAltPipeline.sh.

## example usage
```bash
bash /<AltAmpDS directory>/runAltPipeline -h #to get help and see the different parameters
bash /<AltAmpDS directory>/runAltPipeline -s /<AltAmpDS directory>/trusight_tumor_pipeline.sh > output_alt_pipeline_run.txt 2>&1&
nohup sh /<AltAmpDS directory>/runAltPipeline.sh -debugging true -validation true> output_alt_pipeline_run.txt 2>&1&
```

It is highlest suggested to make script alias to make running the pipeline easier
```bash
cd ~
vim ./.bashrc
```
in the bashrc file under the # User specific aliases and functions section (modify as appropriate for your machine)<br />

```
alias runAltPipeline='nohup bash /<AltAmpDS directory>/runAltPipeline.sh > output_alt_pipeline_run.txt 2>&1&'
alias debugRunAltPipeline='nohup bash /<AltAmpDS directory>/runAltPipeline.sh -debugging true -validation true> output_alt_pipeline_run.txt 2>&1&'
alias validationRunAltPipeline='nohup sh /<AltAmpDS directory>/runAltPipeline.sh -validation true > output_alt_pipeline_run.txt 2>&1&'

```
where runAltPipeline is the default, debugRunAltPipeline and validationRunAltPipeline do not get rid of temporary files, debugRunAltPipeline has less restrictions region depth (to use when testing pipeline with very small artifical fastqs) <br />

**Note:when running the above code the user needs to be in the top directory of a NextSeq or MiSeq folder
as a home_dir was not specified** <br />









## Parameters
PIPELINLE_DIR, this variable needs to be set in ./.bashrc file
```bash
cd ~
vim ./.bashrc
###add this under the alias section, modify to point to the main directory folder of this repository
PIPELINE_DIR=/home/ec2-user/ampDsTs;export PIPELINE_DIR
``` 
THREADS - number of threads to use when calling functions that support multi-threaded workflow (default 25) <br />
this parameter can be changed by specifying the -threads parameter when calling runAltPipeline.sh <br />
MEMORY - integer: the amount of memory to specify for the java virtual manager to use: default 16 <br />
active_case_limit - integer: number of cases to process at one time, default is 8 <br />
```bash
nohup sh $PIPELINE_DIR/runAltPipeline.sh -threads 25 -memory 16 -active_case_limit 8 > output_alt_pipeline_run.txt 2>&1&'
```

## Dependencies

There is a script file called download_dependencies.sh that will help you download all of these programs if running on ubuntu,
similar code for Red-hat is commented out which can be removed if necessary.  Please note that this file will NOT download GATK and Annovar as those programs have license agreements 
to launch the script in the terminal type, this simple script will download all dependencies in the directory it currently resides in
```bash
bash download_dependencies.sh
```

bash
-variables need to be set in ~./.bashrc file
-some of the code uses bash syntax so need to make sure bash installed on linux distribution

To install, proceed to install in the order below

git 
```bash
sudo yum install git #Red-hat
sudo apt-get install git #ubuntu
```
zip
```bash
sudo yum install unzip #Red-hat
sudo apt-get install unzip #ubuntu
```

java
```bash
sudo yum install java-1.8.0-openjdk-devel #Red-hat
sudo apt-get install openjdk-8-jdk #Ubuntu

```
wget
```bash
sudo yum install wget #Red-hat
sudo apt-get install wget #Ubuntu
```

gcc
```bash
sudo yum install gcc #red-hat
sudo apt-get install gcc #ubuntu
```

python-devel
```bash
sudo yum install python-devel #Red-hat
sudo yum install python-dev #Ubuntu
```
zlib
```bash
sudo yum install zlib-devel #Red-hat
sudo apt-get install zlib1g-dev #ubuntu
```
g++
```bash
sudo yum install gcc-c++ #red-hat
sudo apt-get install g++ #ubuntu
```

download the git repository
```bash
gitclone https://github.com/schneiderthomas/AltAmpDs
```

### Python Dependencies

#### pip
```bash
wget https://dl.fedoraproject.org/pub/epel/epel-release-latest-7.noarch.rpm #red-hat
sudo yum install epel-release-7.noarch.rpm #red-hat
sudo yum install python-pip #red-hat
sudo apt-get install python-pip #ubuntu
```
#### biopython (1.66)
```bash
sudo pip install biopython==1.66
```
#### pysam (0.8.4)
```bash
sudo pip install pysam==0.8.4
```
#### pyvcf (0.6.7)
```bash
sudo pip install pyvcf==0.6.7
```
#### Pandas (0.16.2)
```bash
sudo pip install pandas==0.16.2
```
#### regex (2015.3.18)
```bash 
sudo pip install regex==2015.3.18
```

### Linux/Shell Dependencies
#### Zenity
to display dialog boxes from shell script (to let tech know that processing is done)
```bash
sudo yum install zenity #red-hat
sudo apt-get install zenity #ubuntu
```

xterm
```bash
sudo yum install xterm #red-hat
sudo yum install xorg-x11-xauth.x86_64 xorg-x11-server-utils.x86_64 dbus-x11.x86_64 #red-hat
sudo apt-get install xterm xorg dbus #ubuntu

```

bcl2fastq (v2.17)
to convert files from bcl to fastq files
```bash
#optional, already provided as zip
#Red-hat
wget 'ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/software/bcl2fastq/bcl2fastq2-v2.17.1.14-Linux-x86_64.zip'
unzip bcl2fastq2-v2.17.1.14-Linux-x86_64.zip
yum localinstall bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm
#if unbuntu
sudo apt-get install alien dpkg-dev debhelper build-essential #needed for unbuntu
sudo alien bcl2fastq2-v2.17.1.14-Linux-x86_64.rpm
sudo dpkg -i bcl2fastq2-v2.17.1.14-Linux-x86_64.deb

```



## Major Program dependencies
these programs need to be downloaded and/or compiled and their resulting directories need to be placed in this 
directory, a more recent version may be used but there may be some compatibility issues with the pipeline as it is

GATK 3.5 (VERY IMPORTANT AT LEAST 3.5) <br />
annovar - 2014-11-12 <br />

freebayes v0.9.20 <br />
bcftools-1.2 <br />
FastQC v0.11.3 <br />
htslib-1.2.1 <br />
IGVTools 2.3.57 <br />
picard 2.10 <br />
samtools_1.2 <br />
snpeff 4.1g 2015-05-17 <br />
varscan v2.3.9 <br />
bwa 0.7.10 <br />
vcflib v.1.0.0 <br />
CoverageQC - for debugging <br />
bedtools2 -> Version 2.19.1 <br />
Trimmomatic 0.33 <br />



Please note annovar and GATK have license agreements must be accepted before you download them and therefore they cannot be downloaded using the above script.
Instructions to download these files are below:

#annovar
please download annovar, version  2014-11-12 was used originally (therefore is the preferred version to ensure compatibility), to download Annovar click [here](http://www.openbioinformatics.org/annovar/annovar_download_form.php).  After downloading annovar place the annovar folder entitled "annovar" in the current directory
Note: the original splicing threshold for annovar is to 2, this can be modified if one goes to file table_annovar.pl and modifies the line 
```python
$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $protocol -hgvs -outfile $tempfile.$protocol -exonsort $queryfile $dbloc";
```
to 
```python
$sc = "annotate_variation.pl -geneanno -buildver $buildver -dbtype $protocol -splicing_threshold 5 -hgvs -outfile $tempfile.$protocol -exonsort $queryfile $dbloc";
```

# GATK
version 3.5 is being used for this pipeline
get the latest software [here](https://software.broadinstitute.org/gatk/download/)
if download version higher than 3.5, need to change line 34 in amplicon_ds_pipeline.sh
as appropriate


### Extra

In this repository there is a folder called ART, in here you will find shell that can be used to create
artifical FASTQ files similar to an ampliconDS run.  ART version ChocolateCherryCake-03-19-2015 was used in these scripts.

Download the latest ART program [here](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/).



## Reference Files
### hg19
will be downloaded if use download_dependencies.sh script

## ANNOVAR reference files
will download_dependencies.sh install clinvar, cosmic, exac, snp and 1000g
in the annovar directory, see download_dependencies.sh if curious



#### NOTES ON MAJOR FILES


##### runAltPipeline.sh

-the shell script which runs through the current directory (unless given) and feeds files to the pipeline shell script (location can be specified with -s command but default parameters are at the 
top of the shell script which can be changed if one moves the directory <br />

# ASSUMPTIONS: <br />
- the directory has a folder structure <br />
 BaseDirectory -> Data -> Intensities -> BaseCalls <br />
 will exit if this is not seen <br />
  <br />
- there needs to be an even number of fastq files (not including the Undetermined FASTQ files) because there always be either two fastq files (or 8 when a NextSeq Folder with no lane splitting) in Amplicon DS pipeline, will exit if does not see this <br />
  <br />
- if no FASTQ files are present then there needs to be tiffs in BaseDirectory-> Images folder or bcl files in BaseFolder -> Data -> Intensities -> BaseCalls -> L001 & L002 & L003 & L004 so bcl2fastq can turn the images or bcl filtes to fastq files <br />
  <br />
- The BaseFolder name has to start to like 160518_N or 160518_M where the first part is a number and then there is an underscore and either an N letter or a M letter (this tells the script if it is dealing with a NextSeq or MiSeq folder), if the name does not start
  like this it should exit in an error <br />
 
 
