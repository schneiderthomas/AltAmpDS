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


def getmatchall(readforward, readreverse, PrimerDictList, usePosStrand, target_portion_read_forward, target_portion_read_reverse):
    # forward and reverse fasta files should be the same length
    currentChr=readforward.tid
    currentPos=readforward.pos
    count_end = len(PrimerDictList)-1
    count = 0
    endreadmatchlengthforward = None
    endreadmatchlengthreverse = None
    returnlist = None
    foundPrimerLengthforward = None
    foundPrimerLengthreverse = None
    while (count <= count_end):
        
        if PrimerDictList[count]["strand"] == '-' and usePosStrand:
            count= count +1
            continue
        elif int(PrimerDictList[count]["chr"]) != currentChr:
            count= count +1
            continue
        elif abs(int(PrimerDictList[count]["startPos"])-currentPos) > 400:    
            count= count + 1
            continue
        # CHECKING IF THE READ IS REVERSED, IF IT IS REVERSED THEN ENDING SEQUENCE CONTAINS THE PRIMER
        # IF THE READ IS NOT REVERSED THEN THE BEGINNED SEQUENCE CONTAINS THE PRIMER
        currentPrimerforward = PrimerDictList[count]["currentStartPrimer"]
        currentPrimerreverse = PrimerDictList[count]["currentEndPrimer"]
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
            targetLength= PrimerDictList[count]["targetLength"]
            startPos = PrimerDictList[count]["startPos"]
            endPos = PrimerDictList[count]["endPos"]
            count = count_end + 1
            endreadmatchlengthforward = matchforward.span()[1]
            endreadmatchlengthreverse = matchreverse.span()[1]
            foundPrimerLengthforward = len(currentPrimerforward)
            foundPrimerLengthreverse = len(currentPrimerreverse)
            
            returnlist = [endreadmatchlengthforward, foundPrimerLengthforward, endreadmatchlengthreverse, foundPrimerLengthreverse,targetLength,startPos,endPos]                                       
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


def change_cigar_either(booleanIsReverse,originalcigar, endreadmatchlength, numberlist, numberlistint, letterlist, expectedEndReadMatchLength, readLength, read,targetLength,startManifestPos,endManifestPos):

    #returns list of two indexes one says if there was an error making the cigar and other says if there is overlapping of reads
    overlappingReads=False
# because changing read.position
    global forwardread
    global reverseread
    oldAlignmentStart = 0
    # print sum(numberlistint)  
    ###DEBUGGING
    #if forwardread.qname != "NS500796:7:H77VTAFXX:1:11101:7184:4318":
	#return True
    #  print 'Found M01382:52:000000000-AAECF:1:2102:21256:8889'
    #  print 'The endreadmatchlength is ' + str(endreadmatchlength)
    
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
    runningTotal = 0
    DeletionPos = 0
    if booleanIsReverse:
        numberlist = numberlist[::-1]
        numberlistint = numberlistint[::-1]
        letterlist = letterlist[::-1]
        oldAlignmentStartR = reverseread.qstart
        oldPositionStart = reverseread.pos
        #need to do this for reverse
        
        
        
        #WILL DO POSITION CHECK LATER FOR REVERSE
        #if endManifestPos-oldPositionStart-reverseread.qlen > 0:
	#   print read.qname
	#   print "The position is " + str(oldPositionStart)
	#   print "The q length is " + str(reverseread.qlen)
	   #when there is a deletion (ie there is less of the probe because the position start is different)
	   #the oldPositionStart is actually smaller and hence it will be > 0 since normal is zero
	   
	   #print "The reverse read position is " + str(reverseread.pos)
	   #print "The reverse read position is " + str(reverseread.pos+reverseread.qlen)
	   #print "The reverse read qlen is " + str(reverseread.qlen)
	   #print "The reverse read qstart is " + str(reverseread.qstart)
	   #FOR REV
	#   DeletionCount= DeletionCount + abs(endManifestPos-oldPositionStart-reverseread.qlen)
	#   DeletionPos=DeletionCount
	#elif endManifestPos-oldPositionStart-reverseread.qlen < 0:
	 #  InsertionCount= InsertionCount + (endManifestPos-oldPositionStart-reverseread.qlen)
        #else:
	#  return True
	  #print str(endManifestPos-forwardread.pos)
       
    else:
      oldAlignmentStart = forwardread.qstart
      oldPositionStart = forwardread.pos
      softclipbeginning = 0
      if letterlist[0]=='S':
	   softclipbeginning = int(numberlist[0])
      forwardreadOffset = startManifestPos - forwardread.pos + softclipbeginning
      if forwardreadOffset != 1 and (forwardreadOffset > 4 and forwardreadOffset < -2):
       
	#ie maps position "early" because skips first read on probe (like a deletion) or maps late because it added on a read
	 #if too match then too much error in the probe
	 return [True,overlappingReads]
      
      if forwardreadOffset > 1:
	#print "Hello insertion!!!!!!!!!!!!!!!!!!!"
	#print str(startManifestPos - forwardread.pos)
	InsertionCount= InsertionCount + (startManifestPos - forwardread.pos - 1)
	
        #the expectedpostionforward - forwardread.pos should equal 1 expected, because forwardread.pos is base 0? if difference < 1 deletion, if > 1 insertion
    
      if forwardreadOffset < 1:
	#ie maps position "early" because skips first read on probe (like a deletion)
	#print "Hello"
	#print str(startManifestPos - forwardread.pos)
	#need to add on DeletionPos makes it work (if not get cigar error where 25S125M when should be 25S126M)
	#DeletionPos will not be used to modify the position
	DeletionCount= DeletionCount + abs(startManifestPos - forwardread.pos - 1)
	DeletionPos=DeletionCount
	
        #runningTotal=runningTotal+DeletionCount
    # #THERE MAY BE A SITUATION IF THERE IS TANDEM DUPLICATION shortly after a probe and the duplicated entry gets a match even though THOUGH THERE ARE TOO MANY DELETIONS, NEED TO GET RID OF THIS SCENARIO BECAUSE ONLY EXPECTING 2 deletions at most
    
    
    
    
    
     
   

    if not sum(numberlistint) <= expectedEndReadMatchLength:
        
        
        for index in range(len(numberlist)):
            
            #
            if expectedEndReadMatchLength <= numberlistint[index] + runningTotal:
                
                indexchanging = index
                runningTotal = runningTotal + numberlistint[index]
                
		  
                if letterlist[index] == 'D' or letterlist[index] == 'H' or letterlist[index] == 'N' or letterlist[index] == 'P':
                   
                    deletionPresent = True
                    booleanDeletionAtEnd = True
                 
                    DeletionCountExceptEnd = DeletionCount
                    DeletionCount = DeletionCount + (runningTotal - expectedEndReadMatchLength) + 1 
                    runningtotalBeforeDeletion = runningTotal - numberlistint[index]
                if letterlist[index] == 'I':
        		    insertionAtEnd = expectedEndReadMatchLength - (runningTotal - numberlistint[index])
        		   
		    
                if DeletionCount + InsertionCount + insertionAtEnd >= 3:
        			return [True,overlappingReads]
				  
                break
		
            else:
                if letterlist[index] == 'D' or letterlist[index] == 'H' or letterlist[index] == 'N' or letterlist[index] == 'P':
                    deletionPresent = True
                    DeletionCount = DeletionCount + numberlistint[index]

                elif letterlist[index] == 'I':
		    InsertionCount = InsertionCount + numberlistint[index]
		  
                runningTotal = runningTotal + numberlistint[index]
    else:
        return [True,overlappingReads]
    
    
    
    
    #The position modification is not valid if there is a deletion in the primer when a reverse read, hence need to check if a deletion is present in the reverse 
    #read
    if booleanIsReverse:
        #get total Deletions and total insertions
        totalDeletions=0
        totalInsertions=0
        totalSoftClip=0
        for index in range(len(numberlist)):
	     if letterlist[index] == 'D':
	       totalDeletions = totalDeletions + numberlistint[index]
	     elif letterlist[index] == 'I':
	       totalInsertions = totalInsertions + numberlistint[index]
	     elif letterlist[index] == 'S':
	       totalSoftClip == totalSoftClip + numberlistint[index]
	       
        reversereadOffset= endManifestPos-oldPositionStart-reverseread.qlen - totalDeletions + totalInsertions - totalSoftClip
        
        if abs(reversereadOffset) > 2:
	    #if too match then too much error in the probe
	     return [True,overlappingReads]
        if reversereadOffset > 0:
	   DeletionCount= DeletionCount + abs(reversereadOffset)
	   DeletionPos=abs(reversereadOffset)
	   #if reversereadOffset > 2:
	   #    print "Oh no!"
           #    print read.qname
           #    print read.cigarstring
           #    print "below is the results"
           #    print str(reversereadOffset)
           #    print "totalSoftClip is " + str(totalSoftClip)
           #    print "the chr location is " + str(read.tid)
           #    print "the location is " + str(endManifestPos)
	elif reversereadOffset < 0: 
           #print "We are in endManifestPos-oldPositionStart-reverseread.qlen - totalDeletions + totalInsertions"
           #print read.qname
           InsertionCount= InsertionCount + abs(endManifestPos-oldPositionStart-reverseread.qlen - totalDeletions + totalInsertions -totalSoftClip)
    
    
                
    # check if deletion is present in the primer (if so need to adjust cigar so the cigar length is what is expected, need to substract from not the primer)
    # ie if 100M1D21M, need it to be 94M27S not 95M27S
    # if expected endreadmatchlength do not match that means there was either a 1 or 2 bp deletion or 1 or 2 bp insertion
    # need to subtract or add on to the last 
    # IF INSERTION NEED TO MOVE POSITION BACK INSTEAD OF FORWARD
    
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
            # newcigar=newcigar+str(expectedEndReadMatchLength)+'S'
            # newcigarnumberlist.append(str(expectedEndReadMatchLength-numberToSubtract))
            if booleanDeletionAtEnd == True:
                  newcigarnumberlist.append(str(runningtotalBeforeDeletion + InsertionCount -DeletionCountExceptEnd)) 
                  
            else:
                  #newcigarnumberlist.append(str(expectedEndReadMatchLength + InsertionCount -DeletionCount))
                  newcigarnumberlist.append(str(expectedEndReadMatchLength + InsertionCount -DeletionCount))
            newcigarletterlist.append('S')
           
            if indexchanging == 0 and not (letterlist[indexchanging] == 'S' or letterlist[indexchanging] == 'H'):
                # newcigar=newcigar+str(runningTotal-expectedEndReadMatchLength)+letterlist[indexchanging]
                
		
		##TODO MAKE SUBFUNCTION AS CALLED TWICE
		newnumber= runningTotal - (expectedEndReadMatchLength + InsertionCount) + DeletionPos 
		if newnumber > 0:
		   newcigarnumberlist.append(str(newnumber))
		   newcigarletterlist.append(letterlist[indexchanging])
		elif newnumber == 0:
		   continue
		else:
		 
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
	       
            
          	if letterlist[indexchanging] == 'D' or letterlist[indexchanging] == 'H':
		  #print 'I am trying to modify a read with a deletion at the end of the primer'
		  #print 'The read is ' + str(read.qname)
		  needtomodifynumber = True
		  booleanDeletionAtEnd = True
		  newnumber= runningTotal - (expectedEndReadMatchLength + InsertionCount) + DeletionPos 
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
		
		  
		  ##TODO MAKE SUBFUNCTION AS CALLED TWICE
		  newnumber= runningTotal - (expectedEndReadMatchLength + InsertionCount) + DeletionPos 
		  
		  if newnumber > 0:
		    newcigarnumberlist.append(str(newnumber))
		    newcigarletterlist.append(letterlist[indexchanging])
		    #numberToSubtractTotal = abs(DeletionCount)
		    #needtomodifynumber= True
		  elif newnumber == 0:
		    continue
		  else:
		
		    numberToSubtractTotal = abs(newnumber)
		    needtomodifynumber= True
		  ##END TODO  
            
       
            
        else:
            
         
		   
            if needtomodifynumber and letterlist[index] != 'D' and letterlist[index] !='I':
	      
	      
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
      
    
    
    
    
    
  
   

    if newcigarletterlist[0] != 'S':
        print "there was an error and a cigar was not clipped"
        return [True,overlappingReads]
    dontClipEnd=False
    numbertoClip=0
    indexchangingAll=len(newcigarnumberlist)-1
    deletionEnd=0
    insertionEnd=0
    softEnd=0
    matchEnd=0
    foundIndexChanged=False
    #print "Started new end read clip"
    if sum(numberlistint) > expectedEndReadMatchLength + targetLength:
        
        numbertoClip = 0-targetLength
        #starting at one because not wanting to count the probe that was soft clipped
        
        for index in range(1,len(newcigarnumberlist)):
        
	    
            if numbertoClip < 0:
	      
	      if (newcigarletterlist[index] != 'D') and (newcigarletterlist[index] != 'I') and (newcigarletterlist[index] != 'S') and (newcigarletterlist[index] != 'H') and (newcigarletterlist[index] != 'N') and (newcigarletterlist[index] != 'P'):
                    numbertoClip=numbertoClip+int(newcigarnumberlist[index])
              #because less of target will be covered if there is an insertion, DO NOT DO ANYTHING IF INSERTION (THE NUMBER OF MATCHES IN CIGAR ALREADY LESS)
              #elif newcigarletterlist[index] == 'I':
	      #numbertoClip=numbertoClip
	      elif newcigarletterlist[index] == 'D':
		    #because more of target will be covered if there is an deletion
		    #DONT NEED TO MULTIPLY BY TWO BECAUSE THE NUMBER COVERED WILL BE LESS
		    #numbertoClip=numbertoClip+(int(2*newcigarnumberlist[index]))
		    numbertoClip=numbertoClip+(int(newcigarnumberlist[index]))
	      elif newcigarletterlist[index] == 'S':
		    #soft clips are only at the end and since for loop started at 1 this means at the end, if it hasn't
		    #made the target yet will set numbertoClip at zero so in the next if statement it will be set as softEnd
		    numbertoClip=0
		    
	      
		#has to be equal too since 0 represents it exactly reached the target point, using matchEnd, softEnd and insertionEnd after that  
            if numbertoClip >= 0:
	        
	        if foundIndexChanged==False:
		   foundIndexChanged=True
		   indexchangingAll=index
                if newcigarletterlist[index] == 'I':
                    #numbertoClip=numbertoClip+int(newcigarnumberlist[index])
                    insertionEnd=insertionEnd+int(newcigarnumberlist[index])
                    
                elif newcigarletterlist[index] == 'M' and index > indexchangingAll:
		    matchEnd= matchEnd + int(newcigarnumberlist[index])
                elif newcigarletterlist[index] == 'S':
		    softEnd= matchEnd + int(newcigarnumberlist[index])
                elif newcigarletterlist[index] == 'D':
		    if index==indexchangingAll:
		       #can't have indexchanging on a deletion, reseting
		       #indexchangingAll=len(newcigarnumberlist)-1
		       numbertoClip=0
                    deletionEnd = deletionEnd + int(newcigarnumberlist[index])
        #numbertoClip=runningTotal-targetLength
        
      
       
        
        if numbertoClip > 0 and dontClipEnd==False:
            newcigar2 = ""
            newcigarnumberlist2 = []
            newcigarletterlist2 = []
            
            
            for index in range(len(newcigarnumberlist)):
                
                if index ==indexchangingAll:
                    currentNumber=int(newcigarnumberlist[index])
                    currentLetter=newcigarletterlist[index]
                    if currentLetter =='M':
                        if currentNumber-numbertoClip > 0:
                            newcigarletterlist2.append('M')
                            newcigarnumberlist2.append(str(currentNumber-numbertoClip))
                            newcigarletterlist2.append('S')
                            newcigarnumberlist2.append(str(numbertoClip+insertionEnd+matchEnd+softEnd))
                        else:
                            newcigarletterlist2.append('S')
                            newcigarnumberlist2.append(str(numbertoClip+insertionEnd+matchEnd+softEnd))
                    else:
                        newcigarletterlist2.append('S')
                        newcigarnumberlist2.append(str(numbertoClip+insertionEnd+matchEnd+softEnd))
                    
                    
                elif index < indexchangingAll:
                    newcigarletterlist2.append(newcigarletterlist[index])
                    newcigarnumberlist2.append(str(newcigarnumberlist[index]))
                
                
        else:
            dontClipEnd=True
                
                
            
    else:
        dontClipEnd=True
      
      
    if dontClipEnd != True:
        overlappingReads=True
        newcigarletterlist=newcigarletterlist2
        newcigarnumberlist=newcigarnumberlist2
      
      
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
      if dontClipEnd != True:
          #####reverseread.pos = reverseread.pos + reverseread.qstart - oldAlignmentStart + DeletionCount - InsertionCount - insertionEnd + deletionEnd
          reverseread.pos = reverseread.pos + reverseread.qstart - oldAlignmentStartR - insertionEnd + deletionEnd
          ###
          #if abs(startManifestPos-oldAlignmentStartR) < abs(endManifestPos-oldAlignmentStartR):
          #        reverseread.pos = startManifestPos + expectedEndReadMatchLength
          #else:
          #        reverseread.qstart = EndManifestPos + expectedEndReadMatchLength
          
    if not booleanIsReverse:
      forwardread.cigarstring = newcigar
      # need to adjust position if a deletion was present in the primer area so adding + numbertosubtract
      if startManifestPos - oldPositionStart != 1:
	   
           forwardread.pos = forwardread.pos + forwardread.qstart - oldAlignmentStart + DeletionCount - InsertionCount + (startManifestPos - oldPositionStart - 1)
      else:
           forwardread.pos = forwardread.pos + forwardread.qstart - oldAlignmentStart + DeletionCount - InsertionCount
      ###forwardread.pos = readExpectedPostion + DeletionCount - InsertionCount
      #if abs(startManifestPos-oldAlignmentStart) < abs(endManifestPos-oldAlignmentStart):
      #            forwardread.qstart = startManifestPos + expectedEndReadMatchLength
      #else:
      #            forwardread.qstart = EndManifestPos + expectedEndReadMatchLength
      
      
    
    RealOnly = regex.compile(ur'[MIDN]')
    RealLetterlist = regex.findall(RealOnly, newcigar) 
    # http://stackoverflow.com/questions/7368789/convert-all-strings-in-a-list-to-int
    numberlistint = map(int, newcigarnumberlist)  
    calCigarLength = 0
    zeroLengthElement=False
    for eachindex in range(0, len(newcigarletterlist)):
	if (newcigarletterlist[eachindex] != 'D') and (newcigarletterlist[eachindex] != 'H') and (newcigarletterlist[eachindex] != 'N') and (newcigarletterlist[eachindex] != 'P'):
	    calCigarLength = abs(numberlistint[eachindex]) + calCigarLength
	  
    for eachindex in range(0, len(numberlistint)):
	if numberlistint[eachindex]== 0:
	    zeroLengthElement=True
	   
       
       
       
       
       # if letterlist has D,H,N,P; then don't write it down
   
   
    if (calCigarLength != readLength) or zeroLengthElement:
	print "Cigar length did not match or Zero Length Element"
	print read.qname
	print read.query_length
	if booleanIsReverse:
	   print "the deletion count is " + str(DeletionCount)
	   print "the insertion count is " + str(InsertionCount)
	   print "the endManifestPos is " + str(endManifestPos)
	   print "the oldPositionStart is " + str(oldPositionStart)
	   print "the totalDeletions is " + str(totalDeletions)
	   print "the totalInsertions is " + str(totalInsertions)
	   print "the reverseread.qlen is " + str(readLength)
	   print "the endManifestPos-oldPositionStart-reverseread.qlen - totalDeletions + totalInsertions is " + str(endManifestPos-oldPositionStart-readLength - totalDeletions + totalInsertions)
	
	print "is read revserse " + str(booleanIsReverse)
	print "the length of letterlistR is: " + str(len(letterlist))
	print "the length of the matched endreadmatchlength is:" + str(endreadmatchlength)
	print "the length of the expected matched primer is:" + str(expectedEndReadMatchLength)
	print "the new calculated Cigar Length is: " + str(calCigarLength)
	print "the cigar length should be " + str(read.query_length)
	print "the original cigar is: " + str(originalcigar)
	print "the new cigarstring is: " + str(newcigar)
	print "the target length is: " + str(targetLength)
	print "return True at if (calCigarLength != read.query_length)"
	return [True,overlappingReads]
    if (not RealLetterlist):
	print "continue at if (not RealLetterlistF) or (not RealLetterlistR)"
	return [True,overlappingReads]
    
    
 
    return [False,overlappingReads]




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
   #customPrimerList=[]
   rowcount=0
   file_read = csv.reader(open(manifest, 'r'), delimiter='\t');
#    for row in file_read:
#         if rowcount==1:
#             #this is a header row
#             header=row
#             rowcount=rowcount+1
#             continue
#         elif row[0]=='[Targets]':
#             break
#         elif rowcount==0:
#             rowcount=rowcount+1
#             continue
#         currentDLSO=row[header.index('DLSO Sequence')]
#         currentULSO=row[header.index('ULSO Sequence')]
#         startPosition=row[header.index('Start Position')]
#         endPosition=row[header.index('End Position')]
#         strand=row[header.index('Strand')]
#         chr=row[header.index('Chromosome')][3:]
#         currentStartPrimer=''
#         currentEndPrimer=''
#         if strand=='-':
#             currentStartPrimer=Bio.Seq.reverse_complement(currentDLSO)
#             currentEndPrimer=Bio.Seq.reverse_complement(currentULSO)
#             start_neg.append(currentStartPrimer)
#             end_neg.append(currentEndPrimer)
#         else:
#             currentStartPrimer=currentULSO
#             currentEndPrimer=currentDLSO
#             start_pos.append(currentStartPrimer)
#             end_pos.append(currentEndPrimer)
#         currentlist=[strand,chr,startPosition,endPosition,currentStartPrimer,currentEndPrimer]  
#         customPrimerList.append((currentlist[:]))
#         rowcount=rowcount+1
        
   header=[]
   headernext=False
   inTargetsPortion=False
   inProbesPortion=False
   PrimerDictList=[]
   for row in file_read:
 
    if row[0] == '[Probes]':
      headernext = True
      inProbesPortion = True
      inTargetsPortion = False
      continue
    elif row[0] == '[Targets]':
     headernext = True
     inTargetsPortion = True
     inProbesPortion = False
     inTargetInteger = 0
     continue
    elif row[0] == '[Intervals]':
     break
    elif headernext == True:
     header = row
     headernext = False
     continue
    elif inProbesPortion == True:
     currentULSO = str(row[header.index('ULSO Sequence')])
     currentDLSO = str(row[header.index('DLSO Sequence')])
     currentDict = {'ULSO Sequence':currentULSO, 'DLSO Sequence':currentDLSO}
     
     currentDLSO = row[header.index('DLSO Sequence')]
     currentULSO = row[header.index('ULSO Sequence')]
     # startPosition=row[header.index('Start Position')]
     # endPosition=row[header.index('End Position')]
     strand = row[header.index('Strand')]
     chr = row[header.index('Chromosome')][3:]
     currentStartPrimer = ''
     currentEndPrimer = ''
     if strand == '-':
        currentStartPrimer = Bio.Seq.reverse_complement(currentDLSO)
        currentEndPrimer = Bio.Seq.reverse_complement(currentULSO)

     else:
        currentStartPrimer = currentULSO
        currentEndPrimer = currentDLSO
     currentDict = {"strand":strand, "chr":chr, "startPos":0, "endPos":0, "currentStartPrimer":currentStartPrimer, "currentEndPrimer":currentEndPrimer,
                  "startPrimerLength":0, "endPrimerLength":0, "targetLength":0}  
     PrimerDictList.append(currentDict)
     continue
    elif inTargetsPortion == True:
     
     chr = str(row[header.index('Chromosome')])
     startPos = int(str(row[header.index('Start Position')]))
     endPos = int(str(row[header.index('End Position')]))
     PrimerDictList[inTargetInteger]["startPos"] = startPos
     PrimerDictList[inTargetInteger]["endPos"] = endPos
     startPrimerLength = len(PrimerDictList[inTargetInteger]["currentStartPrimer"])
     endPrimerLength = len(PrimerDictList[inTargetInteger]["currentEndPrimer"])
     targetLength = endPos - startPos - len(PrimerDictList[inTargetInteger]["currentStartPrimer"]) - len(PrimerDictList[inTargetInteger]["currentEndPrimer"])
     #adding plus 3 to match the clipping for the start probe (because actually clip before end of probe at the start, which is what illumina does)
     PrimerDictList[inTargetInteger]["targetLength"] = targetLength+3
     PrimerDictList[inTargetInteger]["startPrimerLength"] = startPrimerLength
     PrimerDictList[inTargetInteger]["endPrimerLength"] = endPrimerLength
     # print PrimerDictList[inTargetInteger]["targetLength"] 
     inTargetInteger = inTargetInteger + 1
    else:
     continue
   
   samfile = pysam.AlignmentFile(inputfile)
   
   # template copies the header from the other alignment file
   newreads = pysam.AlignmentFile(outputfile, "wb", template=samfile)
   
    
   readnum = 0
   read1 = None
   read2 = None
   for read in samfile:  # TMS KEY CHANGE (GET RID OF FETCH)
       
       
       # have to process reads1 and reads2
       
       
       # the read.rnext==read.rname is another double check they are on the same chromosome
       #hence I am filtering out secondary alignments
       if read.is_read1 and read.rnext == read.rname and not read.is_secondary:  # Store read1 of the pair and iterate the for loop.
            read1 = read
            
            continue
       elif read.is_read2 and read.rnext == read.rname and not read.is_secondary:  # Store read2. Should now have read1 and read2 from a properly-paired read. Do a bit of QC-ing of the read-pair.
            read2 = read
            
            if read1 is None:
	        read1 = None
	        read1 = None
	        print "WARNING READ1 EMPTY WHILE READ2 WAS NOT, IF SEE ERROR MULTIPLE TIMES CHECK IF SORTED BY QUERYNAME"
	        #print "continue at elif read.is_read2 and read.rnext == read.rname:"
	        continue
            elif read2.qname != read1.qname:
                read1 = None
                read1 = None
                #print "continue at elif read2.qname != read1.qname"
                continue
	    elif read1.is_reverse and read2.is_reverse:
	        #one of the reads has to be forward, filtering out
	        read1 = None
	        read2 = None
	        continue
	    elif not read1.is_reverse and not read2.is_reverse:
	        #one of the reads has to be reverse, filtering out
	        read1 = None
	        read2 = None
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
       numberlistR=[]
       letterlistR=[]
       numberlistF=[]
       letterlistF=[]
       try:
         numberlistR = regex.findall(p, reverseread.cigarstring)
         letterlistR = regex.findall(p2, reverseread.cigarstring)
         numberlistF = regex.findall(p, forwardread.cigarstring)
         letterlistF = regex.findall(p2, forwardread.cigarstring)
       except:
	 read1 = None
	 read2 = None
	 print "ERROR: in making numberlist and letterlist"
	 continue
	 
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
       matchAllList = getmatchall(forwardread, reverseread, PrimerDictList, usePosStrand, target_portion_read_F, target_portion_read_R)
       
       # throwing out reads if the primer did not match, likely low quality reads
       if matchAllList == None:
           #print "continue at if matchAllList == None:"
           #print "The target portion read forward is below"
           #print target_portion_read_F
           #print "The name is " + reverseread.qname
           continue
       else:
	   matchF = matchAllList[0]
           matchFexpectedLength = matchAllList[1] 
           matchR = matchAllList[2]
           matchRexpectedLength = matchAllList[3]
           targetLength=matchAllList[4]
           startPos=matchAllList[5]
           endPos=matchAllList[6]
       
           
        
                   
       originalreversecigar = reverseread.cigarstring
       originalforwardcigar = forwardread.cigarstring
       
       
       
   
       
       
       # cigar_end_error=change_cigar_end(originalreversecigar,matchR-1,numberlistR,numberlistintR,letterlistR,matchRexpectedLength-1,reverseread.query_length)
       #cigar_end_error = change_cigar_end(originalreversecigar, matchRexpectedLength - 1, numberlistR, numberlistintR, letterlistR, matchRexpectedLength - 1, reverseread.query_length, read)
       [cigar_end_error,overlappingReadsReverse] = change_cigar_either(True,originalreversecigar, matchRexpectedLength - 1, numberlistR, numberlistintR, letterlistR, matchRexpectedLength - 1, reverseread.query_length, reverseread, targetLength,startPos,endPos)
       
       # cigar_start_error=change_cigar_start(originalforwardcigar,matchF-1,numberlistF,numberlistintF,letterlistF,matchFexpectedLength-1,forwardread.query_length)
       #cigar_start_error = change_cigar_start(originalforwardcigar, matchFexpectedLength - 1, numberlistF, numberlistintF, letterlistF, matchFexpectedLength - 1, forwardread.query_length, read)
       [cigar_start_error,overlappingReadsForward] = change_cigar_either(False,originalforwardcigar, matchFexpectedLength - 1, numberlistF, numberlistintF, letterlistF, matchFexpectedLength - 1, forwardread.query_length, forwardread, targetLength,startPos,endPos)
            
            
     
       #from samtools specifications 
       #Bit 0x4 (read is unmapped) is the only reliable place to tell whether the read is unmapped. If 0x4 is set, no assumptions
       #can be made about RNAME, POS, CIGAR, MAPQ, and bits 0x2, 0x100, and 0x800
       #note need to comment out what I am doing to rname because even though it is in the specifications, it is giving GATK issues on validation
       if cigar_end_error == True and cigar_start_error == True:
	#TEMP	print "There was an error in making both cigars"		
		read1 = None
                read2 = None
                #print "in cigar_end_error == True and cigar_start_error == True:"
                #print "the original start cigar is" + str(originalforwardcigar)
                #print "the original end cigar is" + str(originalreversecigar)
                continue
       if cigar_start_error == True:
	        
	        #read1 = None
	        #read2 = None
	        #continue
	        #will unmap forward read if only cigar error in one of the reads
	        #print "in cigar_start_error == True"
	        #print "the original start cigar is" + str(originalforwardcigar)
	        if not (forwardread.is_unmapped):
		    tosubtract=0
		    if forwardread.is_secondary:
		      tosubtract=tosubtract+256
		    if forwardread.is_proper_pair:
		      tosubtract=tosubtract+2
		    if forwardread.is_supplementary:
		      tosubtract=tosubtract+2048
	            forwardread.flag = forwardread.flag + 4 - tosubtract
	            forwardread.mapping_quality = 0
	            forwardread.cigarstring=''
	            forwardread.pos=0
	            #forwardread.rname=0
	            reverseread.flag = reverseread.flag + 8
	        
       if cigar_end_error == True:
	        #read1 = None
	        #read2 = None
	        #continue
	        #will unmap reverse read if only cigar error in one of the reads
	        #print "in cigar_end_error == True"
	        #print "the original end cigar is" + str(originalreversecigar)
	        if not (reverseread.is_unmapped):
		    
		    tosubtract=0
		    if reverseread.is_secondary:
		      tosubtract=tosubtract+256
		    if reverseread.is_proper_pair:
		      tosubtract=tosubtract+2
		    if reverseread.is_supplementary:
		      tosubtract=tosubtract+2048
	            reverseread.flag = reverseread.flag + 4 - tosubtract
	            reverseread.mapping_quality = 0
	            reverseread.cigarstring=''
	            reverseread.pos=0
	            #reverseread.rname=0
	            forwardread.flag = forwardread.flag + 8

       
       # final double checking of CIGAR invalids

       #and requires that all of its parts in the expression evaluate to True for the whole expression to be True.
       #or is much less picky, as soon as any part of the expression evaluates to True the whole expression is True.

       
	        
       # print "THE ENDING CIGAR STRING IS : " + read.cigarstring 
       reverseread.mpos = forwardread.pos
       forwardread.mpos = reverseread.pos
       #print "The read is " + str(reverseread.qname)
       
       ###################NOW ADDING A CHECK WHERE I AM UNMAPPING A READ if there are > 10 mismatches or if reads overlap, drop read with more mismatches
       ###from sam  format specifications
       #MD ZString for mismatching positions. Regex : [0-9]+(([A-Z]|\^[A-Z]+)[0-9]+)*10
       #to only get mismatches regex should be [0-9]+([A-Z]+)
       #since don't want it where 6^ATTTTTTA is used because that represents a deletion
       
       
       
       if not (cigar_end_error or cigar_start_error):
	 
	           
	 try:
	     mismatchTagForword=forwardread.get_tag('MD');
             mismatchTagReverse=reverseread.get_tag('MD');
             
             
             p = regex.compile(ur'[0-9]+([A-Z]+)')
             mismatchLetterListForward=regex.findall(p,mismatchTagForword)
             mismatchLetterListReverse=regex.findall(p,mismatchTagReverse)
             numberofMismatchesF=len(mismatchLetterListForward)
             numberofMismatchesR=len(mismatchLetterListReverse)
             
             if numberofMismatchesF > 10 and not (forwardread.is_unmapped):
	       
	       tosubtract=0
	       if forwardread.is_secondary:
		  tosubtract=tosubtract+256
	       if forwardread.is_proper_pair:
		  tosubtract=tosubtract+2
	       if forwardread.is_supplementary:
		  tosubtract=tosubtract+2048
	       forwardread.flag = forwardread.flag + 4 - tosubtract
	       forwardread.mapping_quality = 0
	       forwardread.cigarstring=''
	       forwardread.pos=0
	       
	        #telling it that the mate is unmapped
               reverseread.flag = reverseread.flag + 8
               
               
	          
	          
	      
             elif overlappingReadsReverse and not overlappingReadsForward and (numberofMismatchesR < numberofMismatchesF):
	      
	       tosubtract=0
	       if forwardread.is_secondary:
		  tosubtract=tosubtract+256
	       if forwardread.is_proper_pair:
		  tosubtract=tosubtract+2
	       if forwardread.is_supplementary:
		  tosubtract=tosubtract+2048
	       forwardread.flag = forwardread.flag + 4 - tosubtract
	       forwardread.mapping_quality = 0
	       forwardread.cigarstring=''
	       forwardread.pos=0
	       
	        #telling it that the mate is unmapped
               reverseread.flag = reverseread.flag + 8 
               
             if numberofMismatchesR > 10 and not (reverseread.is_unmapped):
	    
               tosubtract=0
	       if reverseread.is_secondary:
		  tosubtract=tosubtract+256
	       if reverseread.is_proper_pair:
		  tosubtract=tosubtract+2
	       if reverseread.is_supplementary:
		  tosubtract=tosubtract+2048
	       reverseread.flag = reverseread.flag + 4 - tosubtract
	       reverseread.mapping_quality = 0
	       reverseread.cigarstring=''
	       reverseread.pos=0
	       
	       #telling it that the mate is unmapped
 	       forwardread.flag = forwardread.flag + 8
 	       
	    
             elif overlappingReadsForward and not overlappingReadsReverse and (numberofMismatchesF < numberofMismatchesR):
	        #suggests low quality reverse if it was clipped and still had more mismatches
	        tosubtract=0
	        if reverseread.is_secondary:
		  tosubtract=tosubtract+256
	        if reverseread.is_proper_pair:
		  tosubtract=tosubtract+2
	        if reverseread.is_supplementary:
		  tosubtract=tosubtract+2048
	        reverseread.flag = reverseread.flag + 4 - tosubtract
	        reverseread.mapping_quality = 0
	        reverseread.cigarstring=''
	        reverseread.pos=0
	        
	        #telling it that the mate is unmapped
	        forwardread.flag = forwardread.flag + 8 
	        
	 except:
	       pass
	       #Do nothing, ie there is no MD
	       
	      
	  
       
       ###########################
       if reverseread.mapping_quality == 0 and forwardread.mapping_quality == 0:
	   read1 = None
           read2 = None
           continue
       #if read.qname == 'M01382:152:000000000-AM7DD:1:1110:14900:4306':
	    #print "line 1230"
	    #print "In read M01382:152:000000000-AM7DD:1:1110:14900:4306"
	    #print "cigar_end_error is " + str(cigar_end_error)
	    #print "cigar_start_error is " + str(cigar_end_error)
	    #print "forward read qname is " + str(forwardread.qname)
	    #print "reverse read qname is " + str(reverseread.qname) 
	    #print str(reverseread.is_read1)
	    #print str(reverseread.flag)
	    #print str(forwardread.is_read1)
	    #print str(forwardread.flag)
       
       
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

