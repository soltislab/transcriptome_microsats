#!/usr/bin/env python

########################################################################################################################
#
#	LocateSSRsandORFs.py -i output of get_orfs_or_cdss.py -r misa results file -o Output file -l length in bp 
#			of significant overlap -s summary information file
#
#	Written by: Matt Gitzendanner
#				University of Florida
#				Department of Biology
#				magitz@ufl.edu
#
# Script compare the location of an SSR and the longest ORF in a given scaffold.
########################################################################################################################

from Bio import SeqIO
import sys
import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file of amino acid ORFs with their coordinates, ie. the output of get_orfs_or_cdss.py")
parser.add_argument("-r", help="input the SSR results file from misa.pl or after filtering for repeated loci ith SSR_RepeatFilter.py")
parser.add_argument("-o", help="Output file SSR and ORF locations")
parser.add_argument("-l", help="Length in bp of overlap to be considered significant in counts, default=15", default=15)
parser.add_argument("-s", help="Summary info output file")

args = parser.parse_args()

infile= args.i
resultfile=args.r
outfile=args.o
MinOverlapLen=args.l
summaryfile=args.s


try:
	ResultFile=open(resultfile, 'r')
except IOError:
	print "Can't open file", resultfile

try:
	OutFile=open(outfile, 'w')
except IOError:
	print "Can't open file", outfile
	
try:
	SummaryFile=open(summaryfile, 'a')
except IOError:
	print "Can't open file", summaryfile


SSR_dict={}  # Dictionary for info about each locus
SSR_locus_count=0  #Counter for number of loci for summary info
SSR_type_count_dict={} #Dictionary to hold counters for different types of repeat loci (compound, mono, di, tri...)

#Get the first line so we can skip it
FirstLine=ResultFile.readline()

#Read through the rest of the result file and get the info about the SSR loci
for Line in ResultFile:
	Line=Line.strip('\n')
	Line_bits=re.split('\t', Line)
	scaffold=Line_bits[0]
	SSRstart=Line_bits[5]	#start	0
	SSRend=Line_bits[6]		#end	1 
	SSRtype=Line_bits[2]	#type	2
	SSR=Line_bits[3]		#sequence 3
	
	SSR_locus_count+=1
	if scaffold in SSR_dict: #Many scaffolds have multiple SSR loci, add them as extra sets of items in the list
		SSR_dict[scaffold].append(SSRstart) #.append can only do one thing at a time.
		SSR_dict[scaffold].append(SSRend)
		SSR_dict[scaffold].append(SSRtype)
		SSR_dict[scaffold].append(SSR) 
	else:
		SSR_dict[scaffold]=[SSRstart,SSRend,SSRtype,SSR]
	
	try:	#Increment count for the type if we've seen it before
		SSR_type_count_dict[SSRtype]+=1	
	except:	#else add a new type and set to 1
		SSR_type_count_dict[SSRtype]=1


SSR_no_overlap_count_dict={} #Dictionary to hold counts for non-overlapping ORFs and SSRs
SSR_overlap_count_dict={} #Dictionary to hold counts of overlapping ORFs and SSRs
SSR_sig_overlap_count_dict={} #Dictionary to hold counts of significant overlapping ORFs and SSRs
ORF_count=0
ORF_scaffold_dict={}  #dictionary to track start and end coordinates of ORFs on each scaffold--many have more than 1 (multiple ORFs of the same length)


def CompareSSRtoORF(scaffold):
	
	NumberOfLoci=int(len(SSR_dict[scaffold])/4) #How many loci are we looking at on this scaffold
	NumberOfORFs=int(len(ORF_scaffold_dict[scaffold])/2) #How many ORFs are we looking at on this scaffold
	
	for Locus in range(1, NumberOfLoci+1):
		
		ThisSSRStartindex=(int(Locus) - 1 ) * 4  #Since we are adding multiple loci in a single list, we need to offset for the 2nd, 3rd, etc locus.
		ThisSSREndIndex=ThisSSRStartindex+1 
		ThisSSRTypeIndex=ThisSSRStartindex+2
		ThisSSRSSRIndex=ThisSSRStartindex+3
	
		SSRrange=range(int(SSR_dict[scaffold][ThisSSRStartindex]),int(SSR_dict[scaffold][ThisSSREndIndex]))  
		
		MaxOverlap=0
		for ORF in range(1,NumberOfORFs+1):
			ThisORFstartIndex=(int(ORF)-1)*2 #As above but for multiple ORFs with coordinates listed in pairs
			ThisORFendIndex=ThisORFstartIndex+1
			
			ORFrange=range(int(ORF_scaffold_dict[scaffold][ThisORFstartIndex]),int(ORF_scaffold_dict[scaffold][ThisORFendIndex]))
		
			SSRrangeSet=set(SSRrange)
			ORFrangeSet=set(ORFrange)
		
			overlap=SSRrangeSet.intersection(ORFrangeSet) #Get the length of the overlap between SSR and ORF
			
			if len(overlap) > MaxOverlap :	#for this SSR, find the maximum overlap with all ORFs on the scaffold
				MaxOverlap=len(overlap)
			
			
		if int(MaxOverlap) > 0:
			try:
				SSR_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]+=1  #SSR_dict[scaffold][2] is the type of the SSR for this scaffold
			except:
				SSR_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]=1
		
			if MaxOverlap > MinOverlapLen:
				try:
					SSR_sig_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]+=1
				except:
					SSR_sig_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]=1
		else:
			try:
				SSR_no_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]+=1  
			except:
				SSR_no_overlap_count_dict[SSR_dict[scaffold][ThisSSRTypeIndex]]=1
			
		OutFile.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(record.id,SSR_dict[scaffold][ThisSSRTypeIndex],SSR_dict[scaffold][ThisSSRSSRIndex],SSR_dict[scaffold][ThisSSRStartindex],SSR_dict[scaffold][ThisSSREndIndex],ORFstart,ORFend,len(overlap)))


try:
	InFile=open(infile, 'r')
except IOError:
	print "Can't open file", infile



#Read through the ORF file and get ORF info and write stats
for record in SeqIO.parse(InFile, 'fasta'):
	ORF_coords=re.search('(\d+)\.\.(\d+)', record.description)
	ORFstart=ORF_coords.group(1)
	ORFend=ORF_coords.group(2)
	IDbits=re.split('\|',record.id)
	ORF_count+=1

	try:		#Some scaffolds have multiple ORFs, add their start and end coords as a list.
		ORF_scaffold_dict[IDbits[0]].append(ORFstart)
		ORF_scaffold_dict[IDbits[0]].append(ORFend)	# IDbits[0] is scaffold
	except:
		ORF_scaffold_dict[IDbits[0]]=[ORFstart,ORFend]
	
skip_key=0
for key in ORF_scaffold_dict:
	try:
		overlap=CompareSSRtoORF(key)
	except:
		#Some ORFs that were originally pulled have now been filtered as redundant, so we can skip them.
		skip_key+=1
		

#Write summary info:

SummaryFile.write('%s\t' %(infile)) #print name of sample (using infile)
SummaryFile.write('%s\t' %(SSR_locus_count)) #print number of loci
SummaryFile.write('%s\t' %(len(ORF_scaffold_dict))) #print number of scaffolds with ORFs

if float(SSR_locus_count) > 0: #Avoid div by 0 errors where no SSR loci
	SummaryFile.write('%s\t%s\t\t' %(ORF_count, (float(ORF_count)/float(SSR_locus_count)))) #print number of ORFs and fraction of SSR loci with ORFs--some loci don't have ORFs
else:
	SummaryFile.write('%s\tNA\t\t' %(ORF_count))
	
Sum=0
for type in ('c', 'c*', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6'):
	try:
		SummaryFile.write('%s\t' %(SSR_type_count_dict[type])) #print count for each type in list above Note: Cant use the SSR locus type counts as many scaffolds have multiple equal length ORFs
		Sum=Sum+int(SSR_type_count_dict[type])
	except:
		SummaryFile.write('0\t') #or if that type wasn't found print 0

SummaryFile.write('%s\t\t' %(Sum)) #Print Sum and leave a blank column

Sum=0
for type in ('c', 'c*', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6'):
	try:
		SummaryFile.write('%s\t' %(SSR_no_overlap_count_dict[type])) #print non overlapping count for each type in list above Note: Cant use the SSR locus type counts as many scaffolds have multiple equal length ORFs
		Sum=Sum+int(SSR_no_overlap_count_dict[type])
	except:
		SummaryFile.write('0\t') #or if that type wasn't found print 0
SummaryFile.write('%s\t\t' %(Sum)) #Print Sum and leave a blank column

Sum=0
for type in ('c', 'c*', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6'):
	try:
		SummaryFile.write('%s\t' %(SSR_overlap_count_dict[type])) #print total overlaps of any length for each type
		Sum=Sum+int(SSR_overlap_count_dict[type])
	except: 
		SummaryFile.write('0\t') #or if that type wasn't found print 0
SummaryFile.write('%s\t\t' %(Sum)) #Print Sum and leave a blank column

Sum=0
for type in ('c', 'c*', 'p1', 'p2', 'p3', 'p4', 'p5', 'p6'):
	try:
		SummaryFile.write('%s\t' %(SSR_sig_overlap_count_dict[type])) #print total overlaps greater than cutoff
		Sum=Sum+int(SSR_sig_overlap_count_dict[type])
	except: 
		SummaryFile.write('0\t') #or if that type wasn't found print 0
SummaryFile.write('%s\n' %(Sum)) #Print Sum and leave new line at the end
