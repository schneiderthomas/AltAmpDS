#!/bin/sh

#####DEFAULT LOCATION PARAMETERS
#pipeline dir should be a local variable set on your computer
HOME_DIR=$PIPELINE_DIR



###EXAMPLE USAGE
#'sh /u01/tom_pipeline/final_CMP26pipeline/runAltPipeline -folder /u01/tom_pipeline/Test_Folder/ART_151_test/ -s /u01/tom_pipeline/Bioinformatics/amplicon_ds_pipeline.sh -sp TruSight_Tumor_26 -manifest_libA /u01/tom_pipeline/final_CMP26pipeline/Manifest_Folder/TruSightTumor-FPA-Manifest_RevB.txt -manifest_libB /u01/tom_pipeline/final_CMP26pipeline/Manifest_Folder/TruSightTumor-FPB-Manifest_RevB.txt > output_alt_pipeline_run.txt 2>&1&'

usage ()
{
  echo 'Usage: see shell script for default variables, all parameters are optional'
  echo 'Optional Parameters'
  echo '-folder <the NextSeq or MiSeq run folder to perform pipeline on; default: current directory in terminal> '
  echo 'Can specify parameters if desire -s <shell_file for pipeline> -f <name of Sample Project in Illumina SampleSheet for the cases you want to run, bcl2fastq will make a folder in directory with this name, will look for folder>'
  echo '-manifest_libA <the manifest file for library A> -manifest_libB <the manifest file for library B> '
  echo '-validation <boolean true/false - will perfrom more processing steps like seeing if there is a trusight tumor file to compare results to>'
  echo '-debugging <boolean true/false - if true no temporary files will be removed>'
  echo '-home_dir <directory location of where the pipeline folder is located>'
  echo '-active_case_limit <integer: number of cases to process at one time, default is 8>'
  echo '-threads <integer: number of threads to use when calling variables: default 25>'
  echo '-memory <integer: the amount of memory to specify for the java virtual manager to use: default 16>'
  
  exit
}


echo $PWD
#/usr/local/bin/bcl2fastq
NAME=$PWD

#default parameters, can be changed by specifying as a parameter
shell_file=$HOME_DIR/amplicon_ds_pipeline.sh
manifest_libA=$HOME_DIR/Manifest_Folder/TruSightTumor-FPA-Manifest_RevB.txt
manifest_libB=$HOME_DIR/Manifest_Folder/TruSightTumor-FPB-Manifest_RevB.txt
sample_project_folder=""
active_case_limit=8
THREADS=25
MEM=16

validation=false
debugging=false

while [ "$1" != "" ]; do
case $1 in
        -h )   usage
                    ;;
        -folder )   shift
                    NAME=$1
                    ;;
        -s )   shift
                    shell_file=$1
                    ;;
        -sp )   shift
                    sample_project_folder=$1
                    ;;            
        -manifest_libA )   shift
                    manifest_libA=$1
                    ;;
        -manifest_libB )   shift
                    manifest_libB=$1
                    ;;
        -validation )   shift
                    validation=$1
                    ;;
        -debugging )   shift
                    debugging=$1
                    ;;
        -home_dir )   shift
                    HOME_DIR=$1
                    ;;
        -active_case_limit )   shift
                    active_case_limit=$1
                    ;;
        -threads )   shift
                    THREADS=$1
                    ;;
        -memory )   shift
                    MEM=$1
                    ;;
                
    esac
    shift
done   


echo "$NAME"
NoSampleProjectSpecified=False
if [ -d "$NAME/Data/Intensities/BaseCalls/$sample_project_folder" ]
then
	echo "$NAME/Data/Intensities/BaseCalls/$sample_project_folder EXISTS!"
else
	echo "$NAME/Data/Intensities/BaseCalls/$sample_project_folder DOES NOT EXIST!"
	echo "Make sure you are in the directory of a run folder and the name has not been changed"
	echo "exiting"
	zenity --info --title="ALT PIPELINE ERROR" --text '<span foreground="black" font="10">/Data/Intensities/BaseCalls/ DOES NOT EXIST! Make sure you are in the directory of a run folder and the name has not been changed </span>'
	exit
fi
#currentFolderBase will output the base of the folder, ie if in /u01/tom_pipeline/Bioinformatics then currentFolderBase = Bioinformatics
currentFolderBase=$(basename "$NAME")
echo "The current folder base is "$currentFolderBase
sequencer=""
if [[ $currentFolderBase =~ [[:digit:]+]_M ]]
then
        echo "The Sequencer is assumed to be MiSeq"
        sequencer=MiSeq
elif  [[ $currentFolderBase =~ [[:digit:]+]_N ]]
then
        echo "The Sequencer is assumed to be NextSeq"
        sequencer=NextSeq
else
      echo "A sequencer cannot be determined, program will exit"
      echo "Make sure you are in the directory of a run folder and the name of the folder has not been changed"
      zenity --info --title="ALT PIPELINE ERROR" --text '<span foreground="black" font="10">Error! See output_alt_pipeline.txt file. A sequencer cannot be determined, program will exit. Make sure you are in the directory of a run folder and the name of the folder has not been changed. </span>'
      exit
  
fi

search_dir=$NAME/Data/Intensities/BaseCalls/$sample_project_folder

NoSampleProjectSpecified=False
find $search_dir -maxdepth 2 -name "*.fastq*" -not -path "*/Archer_Run*/*" -not -path "*/Sarcoma_Run*/*" -not -path "*/Alt_Alignment*/*" > tmpfilelist.txt
if [ -s "tmpfilelist.txt" ]
then
     echo "There are fastq files, therefore, bcl2fastq does not need to be run"
else
     echo "There are no fastq files, therefore, bcl2fastq needs to be run"
     #need to do lane splitting for basespace as they only accept files that are split by lane
     #/usr/local/bin/bcl2fastq --no-lane-splitting &
     input=$(zenity --title="PICK the CMP26 SampleSheet for this run" --file-selection)
     if [[ "$input" != "${input/ /}" ]]
    then
        zenity --info --title="ALT PIPELINE" --text '<span foreground="black" font="16">Sample Sheet Contains Spaces! Rename!</span>'
        exit
    fi
    
     /usr/local/bin/bcl2fastq --sample-sheet "${input}" &
     wait
     find $search_dir -maxdepth 2 -name "*.fastq*" -not -path "*/Archer_Run*/*" -not -path "*/Sarcoma_Run*/*" -not -path "*/Alt_Alignment*/*" > tmpfilelist.txt
fi

 

sort tmpfilelist.txt > tmpfilelistsorted.txt
rm tmpfilelist.txt



validfilesarray=()

lanesplit=False
foundlanesplit=False
foundmerged=False

while read current_File
do

   if [[ $current_File =~ .*Undetermined.* ]]
   then
       echo "made it in"
       echo "There is an Undetermined file"
   else 
	echo $current_File
	if [[ $current_File =~ [[:digit:]+]_L00[234].* ]]
	  then
	         echo "Made it in"
	         lanesplit=True
	         foundlanesplit=True
	elif [[ $current_File =~ [[:digit:]+]_L001.* ]]
	    then
	       echo "MiSeq file or NextSeq split file"
	         
	 else
	    foundmerged=True
	fi
	
	
        validfilesarray+=($current_File)
   fi


done < tmpfilelistsorted.txt

echo "The length of the array is below"
echo "${#validfilesarray[@]}"


if [ $foundmerged == True ] && [ $foundlanesplit == True ]
then
    echo "There is a mix of merged fastq and lane split fastq"
    echo "Please only have one type of fastq when processing"
    zenity --info --title="ALT PIPELINE ERROR" --text '<span foreground="black" font="10">There is a mix of merged fastq and lane split fastq. Please only have one type of fastq file when processing</span>'
    exit
    
fi



echo "Now doing array test"


#http://www.cyberciti.biz/faq/finding-bash-shell-array-length-elements/
validfilesarrayLength=${#validfilesarray[@]}
echo "$validfilesarrayLength"


if [ $lanesplit = False ]
then
remainder=`expr $validfilesarrayLength % 4`
echo "$remainder"

if [ $remainder == 0 ]
then
#now work only from here
echo "inside"
active_cases=0
for (( c=0; c<$validfilesarrayLength; c=c+4 ))
do
   echo "Welcome $c times"
   echo ${validfilesarray[c]}
   echo ${validfilesarray[c+1]}
   echo ${validfilesarray[c+2]}
   echo ${validfilesarray[c+3]}
   fqr1=$(basename ${validfilesarray[c]})
   if [[ $fqr1 =~(.*).gz ]]
    then
    fqr1=${BASH_REMATCH[1]}
   fi
   fqr2=$(basename ${validfilesarray[c+1]})
   if [[ $fqr2 =~(.*).gz ]]
    then
    fqr2=${BASH_REMATCH[1]}
   fi
   fqr3=$(basename ${validfilesarray[c+2]})
   if [[ $fqr3 =~(.*).gz ]]
    then
   fqr3=${BASH_REMATCH[1]}
   fi
   fqr4=$(basename ${validfilesarray[c+3]})
   if [[ $fqr4 =~(.*).gz ]]
    then
    fqr4=${BASH_REMATCH[1]}
   fi
   

   if [[ $fqr1 =~(.*)_S ]]
   then 
       currentFile=${BASH_REMATCH[1]}
       echo $currentFile
       echo "The base is : $base"
       mkdir $search_dir/Alt_Alignment
       
       #now need to check if there is a trusight vcf
       if [ -f $search_dir/Alignment/$currentFile.vcf ]
       then
          echo "There is a vcf file"
	  echo "$shell_file -fqrdir $search_dir -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentFile -vcforig $search_dir/Alignment/$currentFile.vcf -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging  -threads $THREADS -memory $MEM"
          sh $shell_file -fqrdir $search_dir -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentFile -vcforig $search_dir/Alignment/$currentFile.vcf -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging -threads $THREADS -memory $MEM &

       else
       
       bash $shell_file -fqrdir $search_dir -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentFile -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging -threads $THREADS -memory $MEM & 

       fi #end if there is a vcf file

       

   fi  # end if can match the name

 #don't want to be processing more than 8 cases at a time  
 #DEBUGGING CHANGE LATER
 active_cases=$(($active_cases+1))
 if [ $active_cases = $active_case_limit ]
 then
    echo $active_cases
    echo "Made it into the wait if statement"
    active_cases=0
    echo "The active cases have been set to zero: $active_cases"
    wait
 
 fi   
 
 
done #for loop looping through the files in the folder

else
   echo "The remainder was not zero. This means there are some unexpected files in $NAME/Data/Intensities/BaseCalls"
   echo "Please remove files in $NAME/Data/Intensities/BaseCalls and reprocess again"
   zenity --info --title="ALT PIPELINE ERROR" --text '<span foreground="black" font="10">The remainder was not zero. This means there are some unexpected files in current_folder/Data/Intensities/BaseCalls directory. Please remove files in $NAME/Data/Intensities/BaseCalls and reprocess again.</span>'
   exit
fi #$remainder !=0

fi #lanesplit = False
echo "The lanesplit is $lanesplit"
if [ $lanesplit = True ]
then
remainder=`expr $validfilesarrayLength % 16`
echo "$remainder"

if [ $remainder == 0 ]
then
#now work only from here
echo "inside"
active_cases=0
for (( c=0; c<$validfilesarrayLength; c=c+16 ))
do
   echo "Welcome $c times"
   echo ${validfilesarray[c]}
   echo ${validfilesarray[c+1]}
   echo ${validfilesarray[c+2]}
   echo ${validfilesarray[c+3]}
   echo ${validfilesarray[c+4]}
   echo ${validfilesarray[c+5]}
   echo ${validfilesarray[c+6]}
   echo ${validfilesarray[c+7]}
   if [[ ${validfilesarray[c]} =~(.*)/(.*)_S([[:digit:]]+).* ]]
   then
        currentCase=${BASH_REMATCH[2]}
        samplenumber1=${BASH_REMATCH[3]}
        fqr1=${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}_R1_001.fastq
        fqr2=${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}_R2_001.fastq
        
        
        
   else
      echo "Can't process files"
      exit
   fi
   if [[ ${validfilesarray[c+8]} =~(.*)/(.*)_S([[:digit:]]+).* ]]
   then
        currentCase=${BASH_REMATCH[2]}
        samplenumber2=${BASH_REMATCH[3]}
        fqr3=${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}_R1_001.fastq
        fqr4=${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}_R2_001.fastq
        
   else
      echo "Can't process files"
      exit
   fi
   #echo ${validfilesarray[c]}
   echo $currentCase
   
   mkdir $search_dir/Alt_Alignment
   #####AVOID DIRNAME , NO IDEA WHY THIS DOESN'T WORK
   #mkdir $search_dir/Alt_Alignment
   #directory= $(dirname ${validfilesarray[c]})
   #echo $directory
   #>/u01/cp26/160204_NS500796_0009_AH523YAFXX/Data/Intensities/BaseCalls/
   #ls ${directorykk}/Alt_Alignment
   #>ls: cannot access /Alt_Alignment: No such file or directory

   
   #http://www.cyberciti.biz/faq/linux-unix-bsd-xargs-construct-argument-lists-utility/
   #http://stackoverflow.com/questions/20269590/final-arguments-for-xargs
   #note: don't use the -0 option for xargs, doesn't work
   ls ${search_dir}/${currentCase}*.fastq.gz | xargs -I file cp file $search_dir/Alt_Alignment 

   gunzip ${search_dir}/Alt_Alignment/${currentCase}*.fastq*
   
   
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}*_R1*.fastq | xargs echo "cat"
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}*_R1*.fastq | xargs cat > $fqr1
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}*_R2*.fastq | xargs echo "cat"
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber1}*_R2*.fastq | xargs cat > $fqr2
   
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}*_R1*.fastq | xargs echo "cat"
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}*_R1*.fastq | xargs cat > $fqr3
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}*_R2*.fastq | xargs echo "cat"
   ls ${search_dir}/Alt_Alignment/${currentCase}_S${samplenumber2}*_R2*.fastq | xargs cat > $fqr4
   
   
   
   fqr1=$(basename ${fqr1})
   fqr2=$(basename ${fqr2})
   fqr3=$(basename ${fqr3})
   fqr4=$(basename ${fqr4})
   echo $fqr1
   echo $fqr2
   echo $fqr3
   echo $fqr4
   
   
    
      
  
       
       #now need to check if there is a trusight vcf
       if [ -f $search_dir/Alignment/$currentFile.vcf ]
       then
          echo "There is a vcf file"
	  echo "$shell_file -fqrdir $search_dir/Alt_Alignment -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentCase -vcforig $search_dir/Alignment/$currentFile.vcf -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging -threads $THREADS -memory $MEM"
          sh $shell_file -fqrdir $search_dir/Alt_Alignment -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentCase -vcforig $search_dir/Alignment/$currentFile.vcf -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging -threads $THREADS -memory $MEM &

       else
         echo "$shell_file -fqrdir $search_dir/Alt_Alignment -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentCase -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging & -threads $THREADS -memory $MEM"
         bash $shell_file -fqrdir $search_dir/Alt_Alignment -fqr1 $fqr1 -fqr2 $fqr2 -fqr3 $fqr3 -fqr4 $fqr4 -s $search_dir/Alt_Alignment -sname $currentCase -sequencer $sequencer -manifest_libA $manifest_libA -manifest_libB $manifest_libB -validation $validation -debugging $debugging -threads $THREADS -memory $MEM & 


       fi #end if there is a vcf file

       

   

 #don't want to be processing more than 8 cases at a time  
 #DEBUGGING CHANGE LATER
 active_cases=$(($active_cases+1))
 if [ $active_cases = $active_case_limit ]
 then
    echo $active_cases
    echo "Made it into the wait if statement"
    active_cases=0
    echo "The active cases have been set to zero: $active_cases"
    wait
 
 fi   
 
 
 
done #for loop looping through the files in the folder

else
   echo "The remainder was not zero. This means there are some unexpected files in $NAME/Data/Intensities/BaseCalls"
   echo "Please remove files in $NAME/Data/Intensities/BaseCalls and reprocess again"
   zenity --info --title="ALT PIPELINE ERROR" --text '<span foreground="black" font="10">The remainder was not zero. This means there are some unexpected files in current_folder/Data/Intensities/BaseCalls directory. I am assuming only FASTQ files cooresponding to AmpliconDS files are here.  Please remove files that do not coorespond to AmpliconDS Please remove unexpected files in $NAME/Data/Intensities/BaseCalls and reprocess again.</span>'
   exit
fi #$remainder !=0

fi #lanesplit = True


#waiting for computer to finish running through folder
wait
rm -f ${search_dir}/Alt_Alignment/*.fastq*

echo "ALL DONE WITH PROCESSING!"

#gui to let the tech know that processing is complete
zenity --info --title="ALT PIPELINE" --text '<span foreground="black" font="32">Processing Complete</span>'
