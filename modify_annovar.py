#!/usr/bin/python

import csv
import vcf
import regex
import math
import sys, getopt

#using pandas for ease of merging columns to tsv files
import pandas as pd
#from dns.tokenizer import DELIMITER



def main(argv):
    global inputfile
    global outputfile
    global vcffile
    global pipeline
    global consensus_refseq_csv
    try:
        opts, args = getopt.getopt(argv,"hi:o:v:p:c:",["ifile=","ofile=","vcffile=","pipelineused=","consensusfile="])
    except getopt.GetoptError:
        print 'modify_annovar.py -i <inputannovar> -o <outputannovar> -v <vcffile> -p <pipelineused> -c <consensus_refseq_file>'
        sys.exit(2)

    for opt, arg in opts:
      if opt == '-h':
        print 'modify_annovar.py -i <inputannovar> -o <outputannovar> -v <vcffile> -p <pipelineused> -c <consensus_refseq_file>'
        print 'The purpose of the simple python script is to add columns to annovar files from vcfs which give allele data, which aligner and variant caller was used, and allele data of each library, plus add filter data columns as well'
        print 'Please note annovar data must NOT have otherinfo because this program cannot handle empty headers' 
        print 'Note: if multisample - there must not be more than 2 samples; the multisample processing is specific to trusight tumor'
        print 'the consensus refseq file is a two column csv file with headers with first column gene name and the consensus transcript used, header names can be anything'
        sys.exit()
      elif opt in ("-i", "--ifile"):
        inputfile = arg
      elif opt in ("-o", "--ofile"):
        outputfile = arg
      elif opt in ("-v", "--vcffile"):
        vcffile = arg   
      elif opt in ("-p", "--pipelineused"):
        pipeline = arg
      elif opt in ("-c", "--consensusfile"):
        consensus_refseq_csv = arg
    
    #In [17]: s1 = Series(['X0', 'X1', 'X2', 'X3'], name='X')
    #In [18]: result = concat([df1, s1], axis=1)
    #http://pandas.pydata.org/pandas-docs/stable/merging.html
    vcf_reader = vcf.Reader(open(vcffile, 'r'))

    ##############################################################
    file_read = csv.reader(open(inputfile, 'r'), delimiter='\t')
    #assuming annovar has otherinfo
    #annovar doesn't label the otherinfo columns and thus need to add or pandas will crash
    row1=len(file_read.next())
    row2=len(file_read.next())
    difference=row2-row1
    newheader=[]
    for i in xrange(difference):
        if i==2:
            newheader.append('VCF.CHR')
        elif i==3:
            newheader.append('VCF.POS')
        elif i==9:
            newheader.append('VCF.INFO')
        else:
            newheader.append(i)  
    
    file_read = csv.reader(open(inputfile, 'r'), delimiter='\t')
    inputfile2=inputfile+"2.txt"
    new_file=open(inputfile2,'w+')
    new_file_write=csv.writer(new_file,delimiter='\t')      
    #print row_count_2
    #print row_count_1
    row_number=0
    for row in file_read:
        if row_number==0:
            new_file_write.writerow(row+newheader)
        else:
            new_file_write.writerow(row)
        row_number=row_number+1    
    new_file.close()
    ##########################################

    #file_read.close()
    annovar_original = pd.read_table(inputfile2,header=0)
    annovar_modified=annovar_original.copy()

    consensus_genes=pd.read_csv(consensus_refseq_csv,header=0, index_col=0)

    #outputs VCFv4.1
    #print vcf_reader.metadata['fileformat']  when using this as the file : EGFR_deletion_BRAF_snp_uniquify_lib_A.gatk.hg19_multianno_no_otherinfo


    numofrecords=0;
    numofsamples=len(vcf_reader.samples)
    newcolumns=pd.DataFrame(columns=('Read Depth', 'Alt Variant Freq', 'Alt Read Depth'))
    filtercolumn=pd.Series(name='Filters')
    libAcolumns=pd.DataFrame(columns=('Read Depth Library A', 'Alt Variant Freq Library A', 'Alt Read Depth Library A'))
    libBcolumns=pd.DataFrame(columns=('Read Depth Library B', 'Alt Variant Freq Library B', 'Alt Read Depth Library B'))
    parsecolumns=pd.DataFrame(columns=('cDNA_Variant', 'Protein_Variant', 'Exon','Transcript','Consensus.AAChange' ))
    rawvcfinfocolumns=pd.DataFrame(columns=('SNPEFF annotation','Info','Qual','VCF_Filter','Format','Sample1','Sample2'))




    if numofsamples>2:
        raise Exception('The vcf file has more than two samples, this program will not work more than 2 two samples')
    #end if numofsamples>2

    libraryAname='None'
    libraryBname='None'
    if numofsamples==2:
        libraryAname=vcf_reader.samples[0]
        libraryBname=vcf_reader.samples[1]

    #print annovar_modified.columns.values.tolist()
    #annovar_modified.rename(columns={'Start':'Coordinate'}, inplace=True);
    #annovar_modified.rename(columns={'End':'Coordinate_End'}, inplace=True);
    #annovar_modified.rename(columns={'ExonicFunc.refGene':'Type'}, inplace=True);
    #annovar_modified.rename(columns={'Gene.refGene':'Gene'}, inplace=True);
    #print annovar_modified.columns.values.tolist()
    for record in vcf_reader:
        
        for altrecord in xrange(len(record.ALT)):
            allele_depth_record=0
            currentrecord=numofrecords+altrecord 
            #setting up variables
            consensus_aachange=""
            transcript=""
            cDNA_variant=""
            protein_variant=""
            exon=""
            snpeff=""
            #print "The record Position is below"
            #print record.POS
            #print "The annovar vcf position is below"
            #print annovar_modified.loc[currentrecord].get('VCF.POS')
            
            if record.POS != annovar_modified.loc[currentrecord].get('VCF.POS'):
                raise Exception('ERROR: VCF AND ANNOVAR DO NOT MATCH!  CANNOT MODIFY ANNOVAR! ANNOVAR MIGHT HAVE THROWN OUT A RECORD')
            
            if record.INFO.has_key('EFF'):
                snpeff=str(','.join(record.INFO.get('EFF')))
            ##############now getting the current gene name
            current_gene_name = annovar_modified.loc[currentrecord].get('Gene.refGene');
            ###have to do this because ANNOVAR may put parenthesis around     
            p=regex.compile(r'[(]')
            gene_regex=regex.finditer(p,current_gene_name)
            try: 
                end_index = gene_regex.next().start()
            except StopIteration: 
                end_index = -1
            if end_index!=-1:
                #if it is a splicing variant won't put in AAchange, moving it so it does, technically not AAchange but it is where I am looking for annotation
                #if annovar_modified.loc[currentrecord]['Func.refGene'] =="splicing":
                #    new_aachange=current_gene_name[end_index+1:len(current_gene_name)-1] #+1  and -1 because it is in parenthesis in current annovar so getting rid of
                #    annovar_modified.loc[currentrecord,'AAChange.refGene']=new_aachange
                #    annovar_modified.loc[currentrecord,'Gene.refGene']=current_gene_name[:end_index]
                
                current_gene_name=current_gene_name[:end_index]
        ###have to do this because ANNOVAR may put semicolon around     
            p=regex.compile(r'[;]')
            gene_regex=regex.finditer(p,current_gene_name)
            try: 
                end_index = gene_regex.next().start()
            except StopIteration: 
                end_index = -1
            if end_index!=-1:
                #if it is a splicing variant won't put in AAchange, moving it so it does, technically not AAchange but it is where I am looking for annotation
                #if annovar_modified.loc[currentrecord]['Func.refGene'] =="splicing":
                #    new_aachange=current_gene_name[end_index+1:len(current_gene_name)-1] #+1  and -1 because it is in parenthesis in current annovar so getting rid of
                #    annovar_modified.loc[currentrecord,'AAChange.refGene']=new_aachange
                #    annovar_modified.loc[currentrecord,'Gene.refGene']=current_gene_name[:end_index]
                
                current_gene_name=current_gene_name[:end_index]
                annovar_modified.loc[currentrecord,'Gene.refGene']=current_gene_name
            #print current_gene_name
            annovar_modified.loc[currentrecord]['Func.refGene']    
            ################################
            ####finished getting current gene name, now getting the consensus refseq
            try:
                transcript=consensus_genes.loc[current_gene_name].get('RefSeq')
                if annovar_modified.loc[currentrecord]['Func.refGene'] =="splicing":
                    hgvsc=annovar_modified.loc[currentrecord].get('GeneDetail.refGene')
                else:
                    hgvsc=annovar_modified.loc[currentrecord].get('AAChange.refGene')
                #print hgvsc
                if isinstance(hgvsc,str): #if not a string then nan hence it is a float and nothing is there
                    current_aachange=hgvsc
                    #print current_aachange
                    current_aachangelist=current_aachange.split(',')
                    for x in range(0,len(current_aachangelist)):
                        if current_aachangelist[x].find(transcript) >=0:
                            consensus_aachange=current_aachangelist[x]
                            consensus_aachangelist=consensus_aachange.split(':')
                            for x in range(0,len(consensus_aachangelist)):
                                if consensus_aachangelist[x].find("c.") >=0:
                                    cDNA_variant=consensus_aachangelist[x]
                                    #print cDNA_variant
                                if consensus_aachangelist[x].find("p.") >=0:
                                    protein_variant=consensus_aachangelist[x]
                                if consensus_aachangelist[x].find("exon") >=0:
                                    exon=consensus_aachangelist[x] 
                        
            except:
                print "ERROR: Gene not in consensus, the current gene is below; okay if EGFR-AS1"
		print current_gene_name
            ############################
            parsecolumns.loc[currentrecord]=[cDNA_variant,protein_variant,exon,transcript,consensus_aachange]
            
            
            sample1str="."
            sample2str="."
            
            
        
            if record.num_called < numofsamples:
                filtercolumn.loc[currentrecord]='PB'
            #for some reason pyvcf doesn't seem to parse varscan and vardict well, even though the filter says pass, so if empty considering it passed    
            elif (record.FILTER==None or record.FILTER=="" or not record.FILTER):
                filtercolumn.loc[currentrecord]='PASS'
            else:
                filtercolumn.loc[currentrecord]=', '.join(record.FILTER)
                
            for sample in record.samples:
                #print hasattr(sample.data, 'DP')
                
                if numofsamples==2:
                    if libraryAname==sample.sample:
                        sample1str=str(sample.data)
                    else:
                        sample2str=str(sample.data)
                else:
                    sample1str=str(sample.data)
                
                #need this check because illumina vcf doesn't have DP here
                if hasattr(sample.data,'DP'):
		  if getattr(sample.data, 'DP')==None and numofsamples==2:
		      #means one of the 2 samples does not have a variant here, hence probe bias
		      if libraryAname==sample.sample:
			  #there is nothing in libA
			  sample1str=str(sample.data)
			  libAcolumns.loc[currentrecord]=[0,0,0]
		      elif libraryBname==sample.sample:
			  libBcolumns.loc[currentrecord]=[0,0,0]
			  sample2str=str(sample.data)
		      else:
			  raise Exception('The vcf file has more than two samples, this program will not work with more than 2 two samples')
		      continue
		     
                
                if hasattr(sample.data, 'AD'): 
                    #FOR ALL VARIANT CALLERS OTHER THAN VARSCAN, AD is a list [x,y] the first value in the list(ie x) is the number of reference reads, the second ie [1] or y is number of alt reads
                    if getattr(sample.data, 'AD')!=None:
                        currentAD=getattr(sample.data,'AD')
                        if isinstance(currentAD,list):
                            allele_depth_record=float(currentAD[1])
                            if not hasattr(sample.data,'DP'):
			      #need to do this for illumina vcf
			      read_depth_record=float(currentAD[1])+float(currentAD[0])
                        else:
                            #for varscan AD is not a list, AD is a single value representing the allele depth
                            allele_depth_record=float(currentAD)   
                elif hasattr(sample.data, 'AO'):           
                    if getattr(sample.data, 'AO')!=None:
                        currentAO=getattr(sample.data, 'AO')
                        if isinstance(currentAO, list):
                            #print 'Add_Allele_Data_From_VCF_To_Annovar:Warning Multiple Different Allele Calls at same spot, annovar separates into 2
                            #print currentAO
                            if currentAO[altrecord] is None:
                                allele_depth_record=0
                            else:
                                allele_depth_record=currentAO[altrecord]
                        else:    
                            allele_depth_record=float(getattr(sample.data,'AO'))
                #need this check because illumina vcf doesn't have DP here
                if hasattr(sample.data, 'DP'):
		  read_depth_record=float(getattr(sample.data, 'DP'))
                allele_frequency_record=((allele_depth_record/read_depth_record)*100)
                if numofsamples==1:
                    newcolumns.loc[currentrecord]=[read_depth_record,allele_frequency_record,allele_depth_record]
                elif libraryAname==sample.sample:
                    libAcolumns.loc[currentrecord]=[read_depth_record,allele_frequency_record,allele_depth_record]
                elif libraryBname==sample.sample:
                    libBcolumns.loc[currentrecord]=[read_depth_record,allele_frequency_record,allele_depth_record]
                else:
                    raise Exception('The vcf file has more than two samples, this program will not work with more than 2 two samples')
            #end for sample in record.samples
            
            rawvcfinfocolumns.loc[currentrecord]=[snpeff,str(record.INFO),str(record.QUAL),str(record.FILTER),str(record.FORMAT),sample1str,sample2str]
            
            
            
            #now making an average depth and allele data if more than one sample
            if numofsamples==2:
                read_depth_record=(libAcolumns.loc[currentrecord][0]+libBcolumns.loc[currentrecord][0])
                allele_depth_record=(libAcolumns.loc[currentrecord][2]+libBcolumns.loc[currentrecord][2])
                allele_frequency_record=(allele_depth_record/read_depth_record)*100
                newcolumns.loc[currentrecord]=[read_depth_record,allele_frequency_record,allele_depth_record]
                
                
            
            #end if numofsamples==2
        #end for altrecord in xrange(len(record.ALT))
        numofrecords=numofrecords+len(record.ALT)
    #end for record in vcf_reader

    #adding in pipeline info
    pipelineColumn=pd.Series(name='Pipeline')
    for i in range(numofrecords):
        pipelineColumn.loc[i]=pipeline

                
    result=pd.concat([annovar_modified,parsecolumns,newcolumns,pipelineColumn,filtercolumn],axis=1)
    if numofsamples==2:
        result=pd.concat([result,libAcolumns,libBcolumns],axis=1)
    result=pd.concat([result,rawvcfinfocolumns],axis=1)    
    #if don't put index false it will output row names that the pandas library gave, so to make it look like the original, set it to false
    #now rearranging the columns
    gene_col=result['Gene.refGene']
    func_col=result['Func.refGene']
    cDNA_col=result['cDNA_Variant']
    protein_col=result['Protein_Variant']
    exonicFunc_col=result['ExonicFunc.refGene']
    filters_col=result['Filters']
    pipeline_col=result['Pipeline']
    variant_freq_col=result['Alt Variant Freq']
    read_depth_col=result['Read Depth']
    alt_read_depth_col=result['Alt Read Depth']
    result.drop(labels=['Gene.refGene','Func.refGene','cDNA_Variant','Protein_Variant','ExonicFunc.refGene','Pipeline','Filters','Alt Variant Freq','Read Depth','Alt Read Depth'],axis=1,inplace=True)
    result.insert(0,'Gene.refGene',gene_col)
    result.insert(1,'Func.refGene',func_col)
    result.insert(2,'cDNA_Variant',cDNA_col)
    result.insert(3,'Protein_Variant',protein_col)
    result.insert(4,'ExonicFunc.refGene',exonicFunc_col)
    result.insert(5,'Pipeline',pipeline_col)
    result.insert(6,'Filters',filters_col)
    result.insert(7,'Alt Variant Freq',variant_freq_col)
    result.insert(8,'Read Depth',read_depth_col)
    result.insert(9,'Alt Read Depth',alt_read_depth_col)

    result.to_csv(outputfile, sep='\t',index=False)





#############end def main  


if __name__ == "__main__":
   main(sys.argv[1:])               
