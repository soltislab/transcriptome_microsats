#!/usr/bin/env python

########################################################################################################################
# compare_pcrs.py -a CODE -b CODE
# CODE is the 1KP 4-letter code for a sample. 
#
# Script to compare the ePCR results
#
# Note paths/names that need to be set manually below.
#
# Written by: Matt Gitzendanner
#				University of Florida
#				Department of Biology
#				magitz@ufl.edu
#
# Note: There may be some parts of this script that end up not being used. There were plans to do more here, but some 
#	parts were shifted to other scripts or left undeveloped. But the script does what it needs to...
#	Just don't be too confused if you see parts that don't end up getting used...
#
########################################################################################################################

from Bio import SeqIO
from Bio import pairwise2
from Bio import AlignIO
from Bio.Align import AlignInfo

import sys
import argparse
import re
import tempfile

parser = argparse.ArgumentParser()
parser.add_argument("-a", help="name of sample a")
parser.add_argument("-b", help="name of sample b")

args = parser.parse_args()

sampleA= args.a
sampleB= args.b

######################################################################
#
# Need to set some paths:
primer_file_path="primers/"	#Folder where primer files are
primer_file_suffix=".ePCR"	#Suffix added to name to get filename
ePCR_file_path="ePCR_out/" #Folder where ePCR results are located
out_path="ePCR_summary_160303/" #where to write result files

scaffold_path="../Assemblies/" #where the assemblies are located. This script needs to read all the assemblies into memory to quickly compare the amplicons. As such it will take a fair bit of RAM when it is run.
scaffold_suffix="-SOAPdenovo-Trans-assembly.fa" #Bit after the sampleA part

#Use the filtered versions here
misa_output_path="../FilteredSSRs/"
misa_output_suffix="-SOAPdenovo-Trans-assembly.fa.misa_SSR.results"

summary_locus_count= out_path + "Loci_amplified_summary.160303.txt"
#
######################################################################


#Assumption is that the ePCR result names are: sampleA.on.sampleB.ePCR.out
sampleA_primers_file=primer_file_path + str(sampleA) + primer_file_suffix
sampleB_primers_file=primer_file_path + str(sampleB) + primer_file_suffix

sampleA_scaffold_file= scaffold_path + sampleA + scaffold_suffix
sampleB_scaffold_file= scaffold_path + sampleB + scaffold_suffix

ePCR_result_file= ePCR_file_path + str(sampleA) + ".on." + str(sampleB) + ".ePCR.out"

out_file= out_path + str(sampleA) + ".on." + str(sampleB) + ".ePCR.summary"

amplification_summary= str(sampleA) + ".on." + str(sampleB)

sampleA_misa_file= misa_output_path + str(sampleA) + misa_output_suffix
sampleB_misa_file= misa_output_path + str(sampleB) + misa_output_suffix



try:
	ePCR_result=open(ePCR_result_file, 'r')
except IOError:
	print "Can't open file", ePCR_result_file


try:
	OutFile=open(out_file, 'w')
except IOError:
	print "Can't open file", out_file
	

try:
	SummaryFile=open(summary_locus_count, 'a')
except IOError:
	print "Can't open summary file: ", summary_locus_count
	

class SSRlocus:
	'A class to store information about an SSR locus--make it easier to add locus info to dictionary.'
	def __init__(self, SSR_type, SSR_seq, SSR_size, SSR_start, SSR_end,SSR_amplicon_start, SSR_amplicon_end):
		self.SSR_type = SSR_type
		self.SSR_seq = SSR_seq
		self.SSR_size = SSR_size
		self.SSR_start = SSR_start
		self.SSR_end = SSR_end
		self.SSR_amplicon_seq=""
		self.SSR_amplicon_start=SSR_amplicon_start
		self.SSR_amplicon_end=SSR_amplicon_end
	
def Read_primers_to_dict(file): 
	'Reads primer file into a dictionary with a list of expected amplicon sizes for each scaffold. There can be multiple loci on a scaffold, thus the need to use a list here'
	try:
		primers=open(file, 'r')
	except IOError:
		print "Can't open file", file
	
	#Get the first line so we can skip it
	FirstLine=primers.readline()
	num_loci=0 #Count the number of loci
	locus_count_dict={} #Keep track of which locus we are looking at for a scaffold. There are multiple, but we don't have them numbered anymore.

	primer_dict={}
	
	for Line in primers:
		Line=Line.strip('\n')
		Line_bits=re.split('\t', Line)
		num_loci+=1
		
		primer_scaffold=Line_bits[0]
		try:
			locus_count_dict[primer_scaffold]+=1 #Already seen this scaffold in sampleA, need to increment.
		except:
			locus_count_dict[primer_scaffold]=0 #This is the first time we've seen this scaffold in sampleA (index from 0)
	

		try:
			primer_dict[primer_scaffold].append(int(Line_bits[3])) # Append size of amplicon to list of amplicons for this scaffold  (Line_bits[3] = size of amplicon)
		except:
			primer_dict[primer_scaffold]=[int(Line_bits[3])] #If this is the first, just add it to a list.
			
	return(primer_dict, num_loci)
		
def Read_seqs_to_dict(file):
	try:
		sequences=open(file, 'r')
	except IOError:
		print "Can't open file", file
	
	record_dict = SeqIO.to_dict(SeqIO.parse(file, "fasta"))
	sequences.close()
	return(record_dict)
	
	
def Read_misa_seqs_to_dict(file):
	try:
		misa_info=open(file, 'r')
	except IOError:
		print "Can't open file", file
	#Get the first line so we can skip it
	FirstLine=misa_info.readline()
	
	SSR_dict={}
	for Line in misa_info:
		Line=Line.strip('\n')
		Line_bits=re.split('\t', Line)
	
		Locus_ID=Line_bits[0]
	
		#For each scaffold, create a dictionary that is a list of the SSRLocus objects found on that scaffold.
	
		try:		#Some scaffolds have multiple SSRs, add their start and end coords as a list.
			SSR_dict[Locus_ID].append(SSRlocus(Line_bits[2], Line_bits[3], Line_bits[4], Line_bits[5], Line_bits[6], Line_bits[14], Line_bits[15]))
		except:
									#SSR_type			SSR_seq			SSR_size	SSR_start		SSR_end		product start  product end
			SSR_dict[Locus_ID]=[SSRlocus(Line_bits[2], Line_bits[3], Line_bits[4], Line_bits[5], Line_bits[6], Line_bits[14], Line_bits[15])]

	return(SSR_dict)

 
#Read the primer files for each sample	
(sampleA_dict, sampleA_num_loci)=Read_primers_to_dict(sampleA_primers_file)
(sampleB_dict, sampleB_num_loci)=Read_primers_to_dict(sampleB_primers_file)


#Read the scaffolds for each file
sampleA_scaffold_dict=Read_seqs_to_dict(sampleA_scaffold_file)
sampleB_scaffold_dict=Read_seqs_to_dict(sampleB_scaffold_file)

#Read the SSR locus information from Misa results for sampleA
sampleA_SSR_locus_dict=Read_misa_seqs_to_dict(sampleA_misa_file)


ePCR_dict={}
locus_count_dict={} #Keep track of which locus we are looking at for a scaffold. There are multiple, but we don't have them numbered anymore.
sampleA_loci_amplifying=[]

for Line in ePCR_result:
	Line=Line.strip('\n')
	Line_bits=re.split('\t', Line)
	
	sampleA_scaffold=Line_bits[1]		#ePCR output format is name in B, then name in A
	sampleB_scaffold=Line_bits[0]
	
	try:
		locus_count_dict[sampleA_scaffold]+=1 #Already seen this scaffold in sampleA, need to increment.
	except:
		locus_count_dict[sampleA_scaffold]=0 #This is the first time we've seen this scaffold in sampleA (index from 0)
	
	
	sampleB_start=int(Line_bits[3])-1 #Need to do -1 this value for the correct sequence.
	sampleB_end=int(Line_bits[4])-1
	
	amplicon_size=(sampleB_end-sampleB_start)+1
	
	try:
		size_difference= amplicon_size - int(sampleA_dict[sampleA_scaffold][locus_count_dict[sampleA_scaffold]])
	
	except:
		print ("Can't get size for %s locus %d" %(sampleA_scaffold, locus_count_dict[sampleA_scaffold]))
		
		
# 	if sampleB_scaffold in sampleB_dict:		#Was the scaffold where the ePCR amplified in sample B also in the list of SSR loci for sample B?
# 		locus_in_sampleB=1
# 	else:
# 		locus_in_sampleB=0
	
	sampleB_amplicon=sampleB_scaffold_dict[sampleB_scaffold].seq[sampleB_start:sampleB_end]
	
	try:
		amplicon_start=int(sampleA_SSR_locus_dict[sampleA_scaffold][locus_count_dict[sampleA_scaffold]].SSR_amplicon_start) #Lookup the amplicon start, use locus_count_dict to see which SSR for the scaffold we are working with
	
		amplicon_end=int(sampleA_SSR_locus_dict[sampleA_scaffold][locus_count_dict[sampleA_scaffold]].SSR_amplicon_end)
		
		amplicon_motif=sampleA_SSR_locus_dict[sampleA_scaffold][locus_count_dict[sampleA_scaffold]].SSR_seq

		locus_name=sampleA_scaffold + "_" + str(locus_count_dict[sampleA_scaffold])
				
		if locus_name not in sampleA_loci_amplifying:
			sampleA_loci_amplifying.append(locus_name)

	except:
		print "Cant get sequence for %s locus %s" %(sampleA_scaffold, locus_count_dict[sampleA_scaffold])
		
	sampleA_amplicon=sampleA_scaffold_dict[sampleA_scaffold].seq[amplicon_start:amplicon_end]
	
	
	
	OutFile.write("%s\t%s\t%s\t%s\t%s\n" %(sampleA_scaffold, sampleB_scaffold, sampleA_amplicon, sampleB_amplicon, amplicon_motif) )
	

SummaryFile.write("%s\t%d\t%d\n" %(amplification_summary,sampleA_num_loci,len(sampleA_loci_amplifying)))	#Write the summary info for each amplification #of loci # amplifying.

	



