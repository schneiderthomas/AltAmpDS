#!/usr/bin/python

import pysam
import csv
import Bio.Seq
import regex
import sys, getopt
# http://www.tutorialspoint.com/python/python_command_line_arguments.htm

# Returns "Start" if You need to look for the primer at the "start" of the read, returns "End" if you need to look for the primer at the end of the read
# parameters are passed by reference by default
#FOR REFERENCE NOT USING ANYMORE
def where_to_look_for_primer(selectedread):

       
     
           
        if selectedread.is_reverse:
            if selectedread.is_read1:
                    currentPrimerList = end_pos
            else:
                    currentPrimerList = end_neg
        else:
            if selectedread.is_read1:
                    currentPrimerList = start_neg
               
            else:
                    currentPrimerList = start_pos
        return currentPrimerList
            

           
#FOR REFERENCE NOT USING ANYMORE
def getmatch(selectedread, primerList, target_portion_read):
    count_end = len(primerList)-1
    count = 0
    endreadmatchlength = None
    returnlist = None
    foundPrimerLength = None
    while (count <= count_end):
        # CHECKING IF THE READ IS REVERSED, IF IT IS REVERSED THEN ENDING SEQUENCE CONTAINS THE PRIMER
        # IF THE READ IS NOT REVERSED THEN THE BEGINNED SEQUENCE CONTAINS THE PRIMER
        currentPrimer = primerList[count]
        if selectedread.is_reverse:
            endPrimer = currentPrimer[::-1]
            endtestprimer = "(" + endPrimer + ")" + "{e<=2}"
            end_reverse_read = target_portion_read[::-1]
            match = regex.match(endtestprimer, end_reverse_read)
                    
        else:
            # searching the beginning, with allowable error of 2 or less
            testprimer = "(" + currentPrimer + ")" + "{e<=2}"
            match = regex.match(testprimer, target_portion_read)
                    
        # index_of_primer_start = read.query_sequence.find(currentPrimer,0,40)
                
        # if index_of_primer != -1:
        if match != None:
            count = count_end + 1
            endreadmatchlength = match.span()[1]
            foundPrimerLength = len(currentPrimer)
            returnlist = [endreadmatchlength, foundPrimerLength]                                       
        else:
            if count != count_end:
                count = count + 1
            else:
                # THIS MEANS THERE IS NO MATCHES AND THEREFORE UNSURE WHAT PRIMER USED AND SO CUTTING OFF THE LENGTH OF THE LARGEST PRIMER WHICH IS THIRTY
                count = count + 1
    
    
    return returnlist


def getmatchall(readforward, readreverse, customPrimerList, usePosStrand, target_portion_read_forward, target_portion_read_reverse):
    # forward and reverse fasta files should be the same length
    currentChr=readforward.tid
    currentPos=readforward.pos
    count_end = len(customPrimerList)-1
    count = 0
    endreadmatchlengthforward = None
    endreadmatchlengthreverse = None
    returnlist = None
    foundPrimerLengthforward = None
    foundPrimerLengthreverse = None
    while (count <= count_end):
        
        if customPrimerList[count][0] == '-' and usePosStrand:
            count= count +1
            continue
        elif int(customPrimerList[count][1]) != currentChr:
            count= count +1
            continue
        elif abs(int(customPrimerList[count][3])-currentPos) > 20000:#note the tile is based on exon size, the ie the start and stop from the manifest cooresponds to exon size, the largest exon in the human genome is 18200 bp    
            count= count + 1
            continue
        # CHECKING IF THE READ IS REVERSED, IF IT IS REVERSED THEN ENDING SEQUENCE CONTAINS THE PRIMER
        # IF THE READ IS NOT REVERSED THEN THE BEGINNED SEQUENCE CONTAINS THE PRIMER
        currentPrimerforward = customPrimerList[count][4]
        currentPrimerreverse = customPrimerList[count][5]
       # if selectedread.is_reverse:
        endPrimerreverse = currentPrimerreverse[::-1]
        testprimerreverse = "(" + endPrimerreverse + ")" + "{e<=2}"
        end_reverse_read = target_portion_read_reverse[::-1]
        matchreverse = regex.match(testprimerreverse, end_reverse_read)
                    
        # else:
            # searching the beginning, with allowable error of 2 or less
        testprimerforward = "(" + currentPrimerforward + ")" + "{e<=2}"
        matchforward = regex.match(testprimerforward, target_portion_read_forward)
                 
        # index_of_primer_start = read.query_sequence.find(currentPrimer,0,40)
                
        # if index_of_primer != -1:
        if matchreverse != None and matchforward != None:
            count = count_end + 1
            endreadmatchlengthforward = matchforward.span()[1]
            endreadmatchlengthreverse = matchreverse.span()[1]
            foundPrimerLengthforward = len(currentPrimerforward)
            foundPrimerLengthreverse = len(currentPrimerreverse)
            returnlist = [endreadmatchlengthforward, foundPrimerLengthforward, endreadmatchlengthreverse, foundPrimerLengthreverse]                                       
        elif matchreverse != None:
	    #print "Reverse matched but forward did not"
	    #print "The read forward is below"
	    #print readforward
	    #print target_portion_read_forward
	    #print "The read reverse is below"
	    #print readreverse
	    #print target_portion_read_reverse
	    #print "The fasta file forward is below"
	    #print primerListforward.filename
	    #print "The fasta file reverse is below"
	    #print primerListreverse.filename
	    count = count + 1
	elif matchforward != None:    
	    #print "Forward matched but reverse did not"
	    #print "The read forward is below"
	    #print readforward
	    #print target_portion_read_forward
	    #print "The read reverse is below"
	    #print readreverse
	    #print target_portion_read_reverse
	    #print "The fasta file forward is below"
	    #print primerListforward.filename
	    #print "The fasta file reverse is below"
	    #print primerListreverse.filename
	    count = count + 1
        else:
            count = count + 1
    
    
    return returnlist


def change_cigar_either(booleanIsReverse,originalcigar, endreadmatchlength, numberlist, numberlistint, letterlist, expectedEndReadMatchLength, readLength, read):

# because changing read.position
    global forwardread
    global reverseread
    oldAlignmentStart = 0
    # print sum(numberlistint)  
    #if read.qname != "NS500796:7:H77VTAFXX:1:11107:3753:3856":	
     # return False
     
   # print 'Found NS500796:7:H77VTAFXX:1:11107:3753:3856'
    #  print 'The endreadmatchlength is ' + str(endreadmatchlength)
      
    if booleanIsReverse:
        numberlist = numberlist[::-1]
        numberlistint = numberlistint[::-1]
        letterlist = letterlist[::-1]
    else:
      oldAlignmentStart = forwardread.qstart
    # #THERE MAY BE A SITUATION IF THERE IS TANDEM DUPLICATION shortly after a probe and the duplicated entry gets a match even though THOUGH THERE ARE TOO MANY DELETIONS, NEED TO GET RID OF THIS SCENARIO BECAUSE ONLY EXPECTING 2 deletions at most
    
    # IF IT IS NOT LESS THEN OR EQUAL TO endreadmatchlength
    indexchanging = 0
    deletionPresent = False
    # numberToSubtract is dependent on whether or not there is a deletion present in the primer area, need to adjust position by the number of deletions present
    # (shouldn't be more than 2 since only 2 deletions are allowed, if more then read should be thrown out)
    # numberToSubtract and InsertionCount, only works if in the primer
    DeletionCount = 0 
    DeletionCountExceptEnd=0
    InsertionCount = 0
    insertionAtEnd = 0
    booleanDeletionAtEnd = False
    booleanInsertionAtEnd = False
    runningtotalBeforeDeletion=0
    if not sum(numberlistint) <= endreadmatchlength:
        
        runningTotal = 0
        for index in range(len(numberlist)):
           
            #
            if endreadmatchlength <= numberlistint[index] + runningTotal:
                
                indexchanging = index
                runningTotal = runningTotal + numberlistint[index]
                
                #if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		#  print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer'
		#  print 'The endreadmatchlength is ' + str(endreadmatchlength)
		#  print 'The current letter is ' + letterlist[index]
		#  print 'The runningTotal is ' + str(runningTotal)
		#  print 'The DeletionCount is ' + str(DeletionCount)
		  
                if letterlist[index] == 'D' or letterlist[index] == 'H' or letterlist[index] == 'N' or letterlist[index] == 'P':
                    #print 'I am trying to modify a read with a deletion at the end of the primer'
		    #print 'The read is' + str(read.qname)
                    deletionPresent = True
                    booleanDeletionAtEnd = True
                   # print 'The number to Subtract is ' + str(DeletionCount)
                    #print 'The runningTotal is ' + str(runningTotal)
                    #print 'The endereadmatchlength is ' + str(endreadmatchlength)
                    DeletionCountExceptEnd = DeletionCount
                    DeletionCount = DeletionCount + (runningTotal - endreadmatchlength)+1
                    runningtotalBeforeDeletion = runningTotal - numberlistint[index]
                    
                  
                if letterlist[index] == 'I':
		    # IF INSERTION IS AT END DON'T NEED TO MODIFY POSITION, SO DON'T ADD TO InsertionCount
		    # ie if 24M3I94M, and primer is 26 so want 26S, need 26S3I92M; it may not be perfectly aligned right, pathologist will need to adjust
		    insertionAtEnd = endreadmatchlength - (runningTotal - numberlistint[index])
		    #print 'There is an insertion at the end and the insertionAtEnd is' + str(insertionAtEnd)
		    #print 'The read name is ' + str(read.qname)
		    #print 'The original cigar is ' + str(originalcigar)
		    #print 'This is a reverse read ' + str(booleanIsReverse)
		    
                if DeletionCount + InsertionCount + insertionAtEnd >= 3:
			# there are too many deletions or insertions in the primer region, can't accept
			# ALSO CANT ACCEPT IF INSERTION AT END OF PRIMER LIKELY ERROR
			#if read.qname == "M01382:52:000000000-AAECF:1:2102:21256:8889":		
			 #   print 'throwing out M01382:52:000000000-AAECF:1:2102:21256:8889'
			  #  print 'InsertionCount is ' + str(InsertionCount)
			return True
		
		  
                break
		
            else:
                if letterlist[index] == 'D' or letterlist[index] == 'H' or letterlist[index] == 'N' or letterlist[index] == 'P':
                    deletionPresent = True
                    DeletionCount = DeletionCount + numberlistint[index]
                    #if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		     #    print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer where it should not be'
		     #    print 'The index is' + str(index)
		     #    print 'The endreadmatchlength is' + str(endreadmatchlength)
                elif letterlist[index] == 'I':
		    InsertionCount = InsertionCount + numberlistint[index]
		  
                runningTotal = runningTotal + numberlistint[index]
    else:
        return
    #if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
    #  print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 2'
    #  print 'The DeletionCount is ' + str(DeletionCount)            
    # check if deletion is present in the primer (if so need to adjust cigar so the cigar length is what is expected, need to substract from not the primer)
    # ie if 100M1D21M, need it to be 94M27S not 95M27S
    # if expected endreadmatchlength do not match that means there was either a 1 or 2 bp deletion or 1 or 2 bp insertion
    # need to subtract or add on to the last 
    # IF INSERTION NEED TO MOVE POSITION BACK INSTEAD OF FORWARD
    #print "The indexchanging is " + str(indexchanging)
    #print "The runningTotal is " + str(runningTotal)
    #print "The indexchanging letter is " + str(numberlistint[indexchanging])
    #print "The current letter is " + str(letterlist[indexchanging])
    #print "The endreadmatchlength is " + str(endreadmatchlength)
    #TO MAKE THINGS EASIER
    #numberToSubtractTotal = DeletionCount + InsertionCount + insertionAtEnd
    TotalInsertionsPresent = InsertionCount + insertionAtEnd
    
    # #PLEASENOTE: INSERTIONATEND DOES NOT CHANGE POSITION
    newcigar = ""
    newcigarnumberlist = []
    newcigarletterlist = []
    needtomodifynumber = False
    numberToSubtractTotal = 0
    for index in range(len(numberlist)):
        if index == 0:
            # newcigar=newcigar+str(endreadmatchlength)+'S'
            # newcigarnumberlist.append(str(endreadmatchlength-numberToSubtract))
            if booleanDeletionAtEnd == True:
                  newcigarnumberlist.append(str(runningtotalBeforeDeletion + InsertionCount -DeletionCountExceptEnd))   
            else:
                  newcigarnumberlist.append(str(endreadmatchlength + InsertionCount -DeletionCount))
            newcigarletterlist.append('S')
           
            if indexchanging == 0 and not (letterlist[indexchanging] == 'S' or letterlist[indexchanging] == 'H'):
                # newcigar=newcigar+str(runningTotal-endreadmatchlength)+letterlist[indexchanging]
                
		
		##TODO MAKE SUBFUNCTION AS CALLED TWICE
		newnumber= runningTotal - (endreadmatchlength + InsertionCount)  
		if newnumber > 0:
		   newcigarnumberlist.append(str(newnumber))
		   newcigarletterlist.append(letterlist[indexchanging])
		elif newnumber == 0:
		   continue
		else:
		 #  if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		  #    print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 3'
		   numberToSubtractTotal = abs(newnumber)
		   needtomodifynumber= True
                ###END TODO
                
                # else just keep the position the same
            else:
                if indexchanging == 0:
                    
                    # IF THE CURRENT S or H is larger then the primer then get rid of the primer length added to the newcigarnumberlist
                    del newcigarnumberlist[len(newcigarnumberlist) - 1]
                    newcigarnumberlist.append(str(numberlistint[0]))
                    
        elif index < indexchanging:
            continue
        elif index == indexchanging:
	       
            # ##DO I NEED ENDREADMATCHLENGTH?? WHY NOT USE EXPECTED??? -> NEED IN CASE DELETION
          	if letterlist[indexchanging] == 'D' or letterlist[indexchanging] == 'H':
		  #print 'I am trying to modify a read with a deletion at the end of the primer'
		  #print 'The read is ' + str(read.qname)
		  needtomodifynumber = True
		  booleanDeletionAtEnd = True
		  newnumber= runningTotal - (endreadmatchlength + InsertionCount)
		  #print 'The newnumber is ' + str(newnumber)
		  if newnumber > 0:
		    newcigarnumberlist.append(str(newnumber))
		    newcigarletterlist.append(letterlist[indexchanging])
		  else:
		    numberToSubtractTotal = abs(newnumber)
		    # DO NOTHING, THE DELETION LETTER WILL BE REMOVES (ONLY TWO SCENARIOS REMEMBER CAN'T HAVE MORE THAN 2 INS OR DEL
		elif letterlist[indexchanging] == 'I':  
		  #NOTE: WE ARE NOT MODIFYING IF THERE IS AN INSERTION AT THE END
		  needtomodifynumber = True
		  booleanInsertionAtEnd = True
		  numberToSubtractTotal = TotalInsertionsPresent
		  newcigarnumberlist.append(numberlist[indexchanging])
                  newcigarletterlist.append(letterlist[indexchanging])
		  
			
		else:
		 # if read.qname == "PTENex6.chr10.89711875_Del_right_after_End_of_Primer":		
		 #   print 'Found PTENex6.chr10.89711875_Del_right_after_End_of_Primer'
		 #   print 'The numberToSubtract is ' + str(numberToSubtractTotal)
		 #   print 'The number of the index being changed is ' + str(numberlist[indexchanging])
		 #   print 'The letter of the index being changed is ' + str(letterlist[indexchanging])
		 #   print 'The answer to the calculation is ' + str(runningTotal - endreadmatchlength - numberToSubtractTotal)
		  
		  ##TODO MAKE SUBFUNCTION AS CALLED TWICE
		  #newnumber= runningTotal - (endreadmatchlength + InsertionCount - (DeletionCount*2))
		  newnumber= runningTotal - (endreadmatchlength + InsertionCount)
		  if newnumber > 0:
		    newcigarnumberlist.append(str(newnumber))
		    newcigarletterlist.append(letterlist[indexchanging])
		  elif newnumber == 0:
		    continue
		  else:
		  #  if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		   #   print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 4'
		    #  print 'The number of the index being changed is ' + str(numberlist[indexchanging])
		     # print 'The letter of the index being changed is ' + str(letterlist[indexchanging])
		     # print 'runningTotal is ' + str(runningTotal)
		     # print 'numberToSubstract is' + str(numberToSubtractTotal)
		     # print 'endreadmatchlength is ' + str(endreadmatchlength)
		    numberToSubtractTotal = abs(newnumber)
		    needtomodifynumber= True
		  ##END TODO  
            # newcigarnumberlist.append(str(runningTotal-endreadmatchlength))
            # newcigarletterlist.append(letterlist[indexchanging])
       
            
        else:
            # newcigar=newcigar+numberlist[index]+letterlist[index]
           # if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		#   print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 1'
		 #  print 'The InsertionCount is' + str(InsertionCount)
		  # print 'The DeletionCount is' + str(DeletionCount)
		  # print 'The numberToSubtract is ' + str(numberToSubtractTotal)
		  # print 'Is there a deletion at the end' + str(booleanDeletionAtEnd)
		   
            if needtomodifynumber and letterlist[index] != 'D' and letterlist[index] !='I':
	      
	      #if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		#   print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 2FFFFF'
		 #  print 'The InsertionCount is' + str(InsertionCount)
		  # print 'The DeletionCount is' + str(DeletionCount)
		   #print 'The numberToSubtract is ' + str(numberToSubtractTotal)
		   #print 'Is there a deletion at the end' + str(booleanDeletionAtEnd)
		   #print 'The current number is ' + str(numberlistint[index])
		   #print 'The current letter is ' + str(letterlist[index])
	      newnumber = numberlistint[index] - numberToSubtractTotal
	      if newnumber > 0:
		newcigarnumberlist.append(str(newnumber))
		newcigarletterlist.append(letterlist[index])
		needtomodifynumber = False
	      else:
		numberToSubtractTotal = numberToSubtractTotal - numberlistint[index] + 1
		newcigarnumberlist.append(str(1))
		newcigarletterlist.append(letterlist[index])
	    elif needtomodifynumber and letterlist[index] == 'D' and InsertionCount >0:
	      #THIS IS A SPECIFIC SITUATION AND NEED TO ADD TO DELETIONCOUNT IN ORDER TO PROPERLY POSITION THE READ IF IT IS FORWARD
	      DeletionCount = DeletionCount + 1
	      #if read.qname == "PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer":		
		#   print 'Found PTENex6.chr10.89711875_Del_at_End_of_Primer_plus_Ins_in_primer 2kl;jl'
		#   print 'The InsertionCount is' + str(InsertionCount)
		#   print 'The DeletionCount is' + str(DeletionCount)
		#   print 'The numberToSubtract is ' + str(numberToSubtractTotal)
		#   print 'Is there a deletion at the end' + str(booleanDeletionAtEnd)
		#   print 'The current number is ' + str(numberlistint[index])
		#   print 'The current letter is ' + str(letterlist[index])
	      newnumber = numberlistint[index] - numberToSubtractTotal
	      if newnumber > 0:
		newcigarnumberlist.append(str(newnumber))
		newcigarletterlist.append(letterlist[index])
		needtomodifynumber = False
              elif newnumber == 0: 
		#KNOW THAT THIS IS A CERTAIN SITUATION IN WHICH WILL WANT TO REDUCE SIZE OF CLIP
		newcigarnumberlist[len(newcigarletterlist)-1]=str(int(newcigarnumberlist[len(newcigarletterlist)-1])-1)
		needtomodifynumber = False
	      else:
		numberToSubtractTotal = numberToSubtractTotal - numberlistint[index] + 1
		newcigarnumberlist.append(str(1))
		newcigarletterlist.append(letterlist[index])
	    else:
	      
	      
	      
	      newcigarnumberlist.append(numberlist[index])
	      newcigarletterlist.append(letterlist[index])
      
    
    if booleanIsReverse:
      newcigarnumberlist = newcigarnumberlist[::-1]
      newcigarletterlist = newcigarletterlist[::-1]
    
    
    for index in range(len(newcigarletterlist)):
	  ###TODO GET RID OF IF FOUND NOT TO HAVE ERRORS
	  ##if int(newcigarnumberlist[index]) <= 0:
	    ###  continue          
      newcigar = newcigar + newcigarnumberlist[index] + newcigarletterlist[index]     
      # Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
    
    if booleanIsReverse:
      reverseread.cigarstring = newcigar
    if not booleanIsReverse:
      forwardread.cigarstring = newcigar
      # need to adjust position if a deletion was present in the primer area so adding + numbertosubtract
      forwardread.pos = forwardread.pos + forwardread.qstart - oldAlignmentStart + DeletionCount - InsertionCount
    return False




def main(argv):
   global inputfile
   global outputfile
   global start_pos
   global end_pos
   global end_neg
   global start_neg
   global library_type
   global reverseread
   global forwardread
   global current
   
   inputfile = ''
   outputfile = ''
   library_type = '' 
   try:
      opts, args = getopt.getopt(argv, "hi:o:t:m:", ["ifile=", "ofile=", "librarytype=", "manifest="])
   except getopt.GetoptError:
      print 'clipPrimers.py -i <inputfile> -o <outputfile> -t <library_type> -m <Illumina manifest file>'
      sys.exit(2)
   for opt, arg in opts:
      if opt == '-h':
        print 'clipPrimers.py -i <inputfile> -o <outputfile> -t <library_type> -m <Illumina manifest file>'
        sys.exit()
      elif opt in ("-i", "--ifile"):
         inputfile = arg
      elif opt in ("-o", "--ofile"):
         outputfile = arg
      elif opt in ("-t", "--librarytype"):
         library_type = arg   
      elif opt in ("-m", "--manifest"):
         manifest = arg
        
   print 'The manifest file is "', manifest
   print 'Input file is "', inputfile
   print 'Output file "', outputfile
   #start pos is postive strand ULSO sequences
   start_pos=[]
   #start neg is negative strand DLSO sequences, note these are reverse complement of sequences in manifest
   start_neg=[]
   #end pos is positive strand DLSO sequences
   end_pos=[]
   #end neg is negative strand ULSO sequences, note these are reverse complement of sequences in manifest
   end_neg=[]
   header=[]
   #list where [strand,chr,startPosition,endPosition,currentStartPrimer,currentEndPrimer] 
   customPrimerList=[]
   rowcount=0
   file_read = csv.reader(open(manifest, 'r'), delimiter='\t');
   for row in file_read:
        if rowcount==1:
            #this is a header row
            header=row
            rowcount=rowcount+1
            continue
        elif row[0]=='[Targets]':
            break
        elif rowcount==0:
            rowcount=rowcount+1
            continue
        currentDLSO=row[header.index('DLSO Sequence')]
        currentULSO=row[header.index('ULSO Sequence')]
        startPosition=row[header.index('Start Position')]
        endPosition=row[header.index('End Position')]
        strand=row[header.index('Strand')]
        chr=row[header.index('Chromosome')][3:]
        currentStartPrimer=''
        currentEndPrimer=''
        if strand=='-':
            currentStartPrimer=Bio.Seq.reverse_complement(currentDLSO)
            currentEndPrimer=Bio.Seq.reverse_complement(currentULSO)
            start_neg.append(currentStartPrimer)
            end_neg.append(currentEndPrimer)
        else:
            currentStartPrimer=currentULSO
            currentEndPrimer=currentDLSO
            start_pos.append(currentStartPrimer)
            end_pos.append(currentEndPrimer)
        currentlist=[strand,chr,startPosition,endPosition,currentStartPrimer,currentEndPrimer]  
        customPrimerList.append((currentlist[:]))
        rowcount=rowcount+1
        
   samfile = pysam.AlignmentFile(inputfile)
   print "The length of start_pos is" + str(len(start_pos))
   print "The length of start_neg is" + str(len(start_neg))
   print "The length of start_pos is" + str(len(end_pos))
   print "The length of end_neg is" + str(len(end_neg))
   print "The first primer of start_pos is" + str(start_pos[0])
   print "The first primer of start_neg is" + str(start_neg[0])
   print "The first primer of end_pos is" + str(end_pos[0])
   print "The first primer of end_neg is" + str(end_neg[0])
   # template copies the header from the other alignment file
   newreads = pysam.AlignmentFile(outputfile, "wb", template=samfile)
   
    
   readnum = 0
   read1 = None
   read2 = None 
   for read in samfile:  # TMS KEY CHANGE (GET RID OF FETCH)
       
       
       # have to process reads1 and reads2
       
       
       # the read.rnext==read.rname is another double check they are on the same chromosome
       if read.is_read1 and read.rnext == read.rname:  # Store read1 of the pair and iterate the for loop.
            read1 = read
            
            continue
       elif read.is_read2 and read.rnext == read.rname:  # Store read2. Should now have read1 and read2 from a properly-paired read. Do a bit of QC-ing of the read-pair.
            read2 = read
            
            if read1 is None:
	        read1 = None
	        read2 = None
	        print "WARNING READ1 EMPTY WHILE READ2 WAS NOT, IF SEE ERROR MULTIPLE TIMES CHECK IF SORTED BY QUERYNAME"
	        #print "continue at elif read.is_read2 and read.rnext == read.rname:"
	        continue
            elif read2.qname != read1.qname:
                read1 = None
                read2 = None
                #print "continue at elif read2.qname != read1.qname"
                continue
            
       if read1 is None or read2 is None:
           #print "continue at if read1 is None or read2 is None"
           continue
       # if read1.tid != read2.tid then that means the mate is mapped to a different chromosome, so getting rid of
       #print read1.cigarstring
       
       if (read1.query_length <= 35) or (read2.query_length <= 35) or (not isinstance(read.cigarstring, str)) or (read1.tid != read2.tid):
            # newreads.write(read)
            #print "continue at if (read1.query_length <= 35) or (read2.query_length <= 35) or (not isinstance(read.cigarstring, str)) or (read1.tid != read2.tid):" 
            continue
       # ASSUMING READ1 AND READ2 are opposite in read.is_reverse
       
       #if read1.qname != "NS500796:7:H77VTAFXX:1:11107:3753:3856" or read1.qname != "NS500796:7:H77VTAFXX:1:11107:3753:3856":	
        #   continue
     
       p = regex.compile(ur'\d+')
       p2 = regex.compile(ur'[MIDNSHP=X]')
       # would rather know which read is reverse and which is forward, will find out later if one is read1 or read2
       #if read1reverse
       usePosStrand=True
       if read1.is_reverse:
           reverseread = read1
           forwardread = read2
           currentPrimerListF=start_pos
           currentPrimerListR=end_pos
           usePosStrand=True
       else:
           reverseread = read2
           forwardread = read1
           currentPrimerListF=start_neg
           currentPrimerListR=end_neg
           usePosStrand=False
       
       
       #if forwardread.cigarstring == '121M' and reverseread.cigarstring == '121M':
	#   print "The read is 121M on both sides, this read should be written."
	#   print "The read name is : " + read1.qname
       numberlistR = regex.findall(p, reverseread.cigarstring)
       letterlistR = regex.findall(p2, reverseread.cigarstring)
       numberlistF = regex.findall(p, forwardread.cigarstring)
       letterlistF = regex.findall(p2, forwardread.cigarstring)
       # http://stackoverflow.com/questions/7368789/convert-all-strings-in-a-list-to-int
       numberlistintR = map(int, numberlistR)
       numberlistintF = map(int, numberlistF)
      
       if letterlistR[len(letterlistR) - 1] == 'H':
           if numberlistintR[len(numberlistintR) - 1] > 2:
               # then it clipped the end and hence it won't match a primer, filtering out
               # newreads.write(read)
               #print "continue at if numberlistintR[len(numberlistintR) - 1] > 2:"
               #print "the reverseread is :" + str(reverseread.cigarstring)
               continue
       if letterlistF[0] == 'H':
           if numberlistintF[0] > 2:
               # then it clipped the start and hence it won't match a primer, filtering out
               #print "continue at if numberlistintF[0] > 2:"
               #print "the forwardread is :" + str(forwardread.cigarstring)
               continue

       # if the primer is already soft-clipped, likely probe mismatch, should get rid of, poor quality		
       if letterlistR[len(letterlistR) - 1] == 'S':
           if numberlistintR[len(numberlistintR) - 1] > 10:     
               #print "continue at if numberlistintR[len(numberlistintR) - 1] > 10:"
               #print "the reverseread is :" + str(reverseread.cigarstring)
               continue
       if letterlistF[0] == 'S':
           if numberlistintF[0] > 10:
	       #print "continue at if numberlistintF[0] > 10:"
	       #print "the forwardread is :" + str(forwardread.cigarstring)
               continue	
       #throw out if heavily soft clipped at end	     
       if letterlistR[0] == 'S':
           if numberlistintR[0] > 50:     
               #print "continue at if numberlistintR[0] > 50:"
               #print "the reverseread is :" + str(reverseread.cigarstring)
               continue
       if letterlistF[len(numberlistintF) - 1] == 'S':
           if numberlistintF[len(numberlistintF) - 1] > 50:
	       #print "continue at if numberlistintF[0] > 10:"
	       #print "the forwardread is :" + str(forwardread.cigarstring)
               continue	
       
       # get read to check on
       # get the last 35 reads
       target_portion_read_R = reverseread.query_sequence[reverseread.query_length - 35:]
       # need read start
       target_portion_read_F = forwardread.query_sequence[0:35]
      
       
       # currentSequence = read.query_sequence
       # returns a match object
       # print read.query_sequence
       
       # print read.query_alignment_end
       # print read.next_reference_start
       # print read.is_reverse
       #currentPrimerListR = ''
       #currentPrimerListR = where_to_look_for_primer(reverseread)
       #currentPrimerListF = ''
       #currentPrimerListF = where_to_look_for_primer(forwardread)
       #[currentPrimerListF,currentPrimerListR] = where_to_look_for_primer_both(forwardread)
      
       
       # First only matched forward then reverse, to increase sensitivity the forward and reverse match have to have same index (ie can't match primer of PTEN forward and TP53 reverse)
       matchAllList = getmatchall(forwardread, reverseread, customPrimerList, usePosStrand, target_portion_read_F, target_portion_read_R)
       
       # throwing out reads if the primer did not match, likely low quality reads
       if matchAllList is None:
           #print "continue at if matchAllList is None:"
           #print "The target portion read forward is below"
           #print target_portion_read_F
           continue
       else:
	   matchF = matchAllList[0]
           matchFexpectedLength = matchAllList[1] 
           matchR = matchAllList[2]
           matchRexpectedLength = matchAllList[3]
       
           
        
                   
       originalreversecigar = reverseread.cigarstring
       originalforwardcigar = forwardread.cigarstring
       	
       # cigar_end_error=change_cigar_end(originalreversecigar,matchR-1,numberlistR,numberlistintR,letterlistR,matchRexpectedLength-1,reverseread.query_length)
       #cigar_end_error = change_cigar_end(originalreversecigar, matchRexpectedLength - 1, numberlistR, numberlistintR, letterlistR, matchRexpectedLength - 1, reverseread.query_length, read)
       cigar_end_error = change_cigar_either(True,originalreversecigar, matchRexpectedLength - 1, numberlistR, numberlistintR, letterlistR, matchRexpectedLength - 1, reverseread.query_length, read)
       
       # cigar_start_error=change_cigar_start(originalforwardcigar,matchF-1,numberlistF,numberlistintF,letterlistF,matchFexpectedLength-1,forwardread.query_length)
       #cigar_start_error = change_cigar_start(originalforwardcigar, matchFexpectedLength - 1, numberlistF, numberlistintF, letterlistF, matchFexpectedLength - 1, forwardread.query_length, read)
       cigar_start_error = change_cigar_either(False,originalforwardcigar, matchFexpectedLength - 1, numberlistF, numberlistintF, letterlistF, matchFexpectedLength - 1, forwardread.query_length, read)
            

       if cigar_end_error == True and cigar_start_error == True:
		#print "There was an error in making both cigars"		
		read1 = None
                read2 = None
                #print "continue at if cigar_end_error == True and cigar_start_error == True:"
                continue
       if cigar_start_error == True:
		#print "There was an error in making the start cigar"		
		read1 = None
                read2 = None
                #print "continue at if cigar_start_error == True"
                continue
       if cigar_end_error == True:
		#print "There was an error in making the end cigar"		
		read1 = None
                read2 = None
                #print "continue at if cigar_end_error == True"
                continue

       if reverseread.cigarstring == "121M" or forwardread.cigarstring == "121M":
           #print "WARNING A READ DIDN'T MATCH A PRIMER EVEN THOUGH IT SHOULD HAVE"
           # should have match a primer because I already had a check for this
           #print "continue at if reverseread.cigarstring == 121M or forwardread.cigarstring == 121M"
           continue
           #sys.exit()
       # final double checking of CIGAR invalids

       #if read.qname == "PTENex6.chr10.89711875_Del_right_after_End_of_Primer":		
        #    print 'found PTENex6.chr10.89711875_Del_right_after_End_of_Primer' 
	 #   print 'The orignal reverse cigar is' + str(originalreversecigar)
 	  #  print 'The orignal forward cigar is' + str(originalforwardcigar)
	   # print 'The new reverse cigar is' + str(reverseread.cigarstring)
 	    #print 'The new forward cigar is' + str(forwardread.cigarstring)
            #print 'The target portion read R is' + str(target_portion_read_R)
            #print 'The target portion read F is' + str(target_portion_read_F)  
	    #print 'The match forward length is' + str(matchF)
	    #print 'The match forward expected length is' + str(matchFexpectedLength)
	    #print 'The match forward length is' + str(matchR)
	    #print 'The match forward expected length is' + str(matchRexpectedLength)

       RealOnly = regex.compile(ur'[MIDN]')
       RealLetterlistR = regex.findall(RealOnly, reverseread.cigarstring)
       RealLetterlistF = regex.findall(RealOnly, forwardread.cigarstring)
       p = regex.compile(ur'\d+')
       p2 = regex.compile(ur'[MIDNSHP=X]')
       # according to the same file specifications
       # Sum of lengths of the M/I/S/=/X operations shall equal the length of SEQ
       numberlistR = regex.findall(p, reverseread.cigarstring)
       letterlistR = regex.findall(p2, reverseread.cigarstring)
       numberlistF = regex.findall(p, forwardread.cigarstring)
       letterlistF = regex.findall(p2, forwardread.cigarstring)
       # http://stackoverflow.com/questions/7368789/convert-all-strings-in-a-list-to-int
       numberlistintR = map(int, numberlistR)
       numberlistintF = map(int, numberlistF)
       calCigarLengthR = 0
       calCigarLengthF = 0
       for eachindex in range(0, len(letterlistR)):
           if (letterlistR[eachindex] != 'D') and (letterlistR[eachindex] != 'H') and (letterlistR[eachindex] != 'N') and (letterlistR[eachindex] != 'P'):
                calCigarLengthR = numberlistintR[eachindex] + calCigarLengthR
                
       for eachindex in range(0, len(letterlistF)):
           if (letterlistF[eachindex] != 'D') and (letterlistF[eachindex] != 'H') and (letterlistF[eachindex] != 'N') and (letterlistF[eachindex] != 'P'):
                calCigarLengthF = numberlistintF[eachindex] + calCigarLengthF
                
       
       # if letterlist has D,H,N,P; then don't write it down
       if (calCigarLengthR != reverseread.query_length) or (calCigarLengthF != forwardread.query_length):
           print "Cigar length did not match"
           print reverseread.qname
           print forwardread.query_length
           print "the length of letterlistR is: " + str(len(letterlistR))
           print "the length of the reverse matched endreadmatchlength is:" + str(matchR)
           print "the length of the reverse expected matched primer is:" + str(matchRexpectedLength)
           print "the length of the forward matched endreadmatchlength is:" + str(matchF)
           print "the length of the forward expected primer is:" + str(matchFexpectedLength)
           print "the new calculated reverse Cigar Length is: " + str(calCigarLengthR)
           print "the reverse cigar length should be " + str(reverseread.query_length)
           print "the original reverse cigar is: " + str(originalreversecigar)
           print "the reverse cigarstring is: " + str(reverseread.cigarstring)
           print "the new forward Cigar Length is: " + str(calCigarLengthF)
           print "the forward cigar length should be " + str(forwardread.query_length)
           print "the original forward cigar is: " + str(originalforwardcigar)
           print "the forward cigarstring is: " + str(forwardread.cigarstring)
           print "continue at if (calCigarLengthR != reverseread.query_length) or (calCigarLengthF != forwardread.query_length)"
           continue
       if (not RealLetterlistF) or (not RealLetterlistR):
           #print "continue at if (not RealLetterlistF) or (not RealLetterlistR)"
           continue
	   
            
       # print "THE ENDING CIGAR STRING IS : " + read.cigarstring 
       reverseread.mpos = forwardread.pos
       #print "The read is " + str(reverseread.qname)
       if reverseread.is_read1:
           newreads.write(reverseread)
           newreads.write(forwardread)
       else:
	   newreads.write(forwardread)
           newreads.write(reverseread)
    # END FOR READ IN SAMFILE
   newreads.close()
   samfile.close()


   


if __name__ == "__main__":
   main(sys.argv[1:])   

