#!/usr/bin/python

import csv
import sys, getopt
#sample usage
#python create_bed_from_manifest.py -i /home/tom/smartgit/Bioinformatics/TruSightTumor-FPA-Manifest_RevB.txt -o /home/tom/workspaceLuna/testing/amplicons.bed  







def main(argv):
   inputmanifest = ''
   outputbed = '' 
   try:
      opts, args = getopt.getopt(argv, "hi:o:e:", ["ifile=", "ofile=", "enlarge="])
   except getopt.GetoptError:
      print 'create_bed_from_manifest.py -i <input_manifest> -o <output_bed_file>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
        print 'create_bed_from_manifest.py -i <input_manifest> -o <output_bed_file>'
        sys.exit()
      elif opt in ("-i", "--ifile"):
         inputmanifest = arg
      elif opt in ("-o", "--ofile"):
         outputbed = arg
      elif opt in ("-e", "--enlarge"):
         enlarge = arg

        

   #inputmanifest='/home/tom/smartgit/Bioinformatics/TruSightTumor-FPA-Manifest_RevB.txt'
   #bedname='/home/tom/workspaceLuna/testing/amplicons.bed'
   file_read = csv.reader(open(inputmanifest, 'r'), delimiter='\t');

   header=[]
   headernext=False
   inIntervalPortion=False

   new_file=open(outputbed,'w+')
   new_file_write=csv.writer(new_file,delimiter='\t',lineterminator='\n')
   subfactor=0
   addfactor=1
   if enlarge == '1':
       addfactor=300
       subfactor=300
   for row in file_read:
    
    if row[0]=='[Intervals]':
        headernext=True
        continue
    elif headernext==True:
        header=row
        headernext=False
        inIntervalPortion=True
        continue
    elif inIntervalPortion==True:
        #have to do the +/- 300 because there are regions that are covered by probes that fall outside the interval and samtools will throw those reads
        #I know Illumina uses the interval for variant calling and only calls variants in the interval, but I don't want that I also want to call 
        #what is covered by a probe hence also expanding because this is used by GATK
        targetStart=str(int(row[header.index('Target Start')])-subfactor)
        targetStop=str(int(row[header.index('Target Stop')])+addfactor)
        chr=row[header.index('Chromosome')]
        selectRow = [chr,targetStart,targetStop]
        #print selectRow
        new_file_write.writerow(selectRow)    
    else:
        continue
    

   new_file.close()




if __name__ == "__main__":
   main(sys.argv[1:])   
