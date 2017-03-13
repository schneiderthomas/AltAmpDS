#!/usr/bin/python

import regex
import csv
import sys
import getopt



def main(argv):

    inputfile=[]
    bedfile=[]
    outputfile=[]
    try:
        opts, args = getopt.getopt(argv,"hi:o:b:",["ifile=","ofile=","bfile="])
    except getopt.GetoptError:
        print 'append_bed_name.py -i <input_gatk_depthOfCoverage_file>  -o <appended_file> -b <bed_file_with_name_in_column_4>'
        sys.exit(2)
    for opt, arg in opts:
      if opt == '-h':
        print 'append_bed_name.py -i <input_gatk_depthOfCoverage_file>  -o <appended_file> -b <bed_file_with_name_in_column_4>'
        sys.exit()
      elif opt in ("-i", "--ifile"):
        inputfile = arg
      elif opt in ("-o", "--ofile"):
        outputfile = arg
        print outputfile
      elif opt in ("-b", "--bfile"):
        bedfile = arg
    print str(outputfile)
    print str(inputfile)
    f = open(outputfile,'wt')
    f2 = open(bedfile,'rb')
    bed_reader=csv.reader(f2,delimiter="\t")

    #with open(inputfile,'rb') as tsvfile:
     # reader=csv.DictReader(tsvfile,delimiter="\t")
      #fieldnames=['Name','Target','total_coverage','average_coverage','libA_total_cvg','libB_total_cvg',
	#	  'libA_mean_cvg','libB_mean_cvg','libA_%_above_250','libB_%_above_250'
	#	  ,'libA_%_above_500','libB_%_above_500'
	#	  ,'libA_%_above_1000','libB_%_above_1000']
      #writer= csv.DictWriter(f,delimiter="\t",fieldnames=fieldnames)
      #writer.writeheader()
      #iteration=0
      #for row in reader:
	  #print((row['Otherinfo']))
	  #if iteration!=0:
	  #    bed_row=bed_reader.next()
	#  bed_row=bed_reader.next()
	#  row['Name']=bed_row[3]
	#  writer.writerow({'Name': bed_row[3], 'Target': row['Target'], 'total_coverage': row['total_coverage'],'average_coverage': row['average_coverage'],
	#		  'libA_total_cvg':row['libA_total_cvg'],'libB_total_cvg':row['libB_total_cvg'],
	#		  'libA_mean_cvg': row['libA_mean_cvg'],'libB_mean_cvg':row['libB_mean_cvg']
	#		  ,'libA_%_above_250':row['libA_%_above_250'],'libB_%_above_250':row['libB_%_above_250']
	#	  ,'libA_%_above_500': row['libA_%_above_500'],'libB_%_above_500':row['libB_%_above_500']
	#	  ,'libA_%_above_1000': row['libA_%_above_1000'],'libB_%_above_1000':row['libB_%_above_1000']})
    #f.close()
    #f2.close()
    
    
    with open(inputfile,'rb') as tsvfile:
      reader=csv.DictReader(tsvfile,delimiter="\t")
      fieldnames=['Name','Target','total_coverage','average_coverage'
		  ,'total_%_above_500','total_%_above_1000']
      writer= csv.DictWriter(f,delimiter="\t",fieldnames=fieldnames)
      writer.writeheader()
      iteration=0
      for row in reader:
	  #print((row['Otherinfo']))
	  #if iteration==1:
	  bed_row=bed_reader.next()
	  #else: 
	  #    iteration=1
	  #    continue
	  writer.writerow({'Name': bed_row[3], 'Target': row['Target'], 'total_coverage': row['total_coverage'],'average_coverage': row['average_coverage']
		  ,'total_%_above_500': row['total_%_above_500'],'total_%_above_1000':row['total_%_above_1000']})
    f.close()
    f2.close()
    
    
    #with open(inputfile,'rb') as tsvfile:
    #    reader=csv.reader(tsvfile,delimiter="\t")
    #    writer= csv.writer(f,delimiter="\t")
    #    iteration=0
     #   for row in reader:
            #print((row['Otherinfo']))
    #        newrow=row[::-1]
            
    #        if iteration!=0:
     #          bed_row=bed_reader.next()
     #          newrow.append(bed_row[3])
     #       else:
     #           newrow.append('Name')
     #           iteration=1
    #        newrow=newrow[::-1]    
    #        writer.writerow(newrow)
    #f.close()
    #f2.close()



if __name__ == "__main__":
   main(sys.argv[1:])     
