#!/usr/bin/env python

########################################################################################################################
# Add genetic distance to table of microsats
#
#  Written by: Matt Gitzendanner
#				University of Florida
#				Department of Biology
#				magitz@ufl.edu
#
# add_distances.py -i <summary table from DiffMicrosatsBatch.R> -g <a distance matrix> -a <summary file>
#		-s <file with 1KP code to species translation> -o <output filename> 
#
########################################################################################################################

import argparse
import re

parser = argparse.ArgumentParser()
parser.add_argument("-i", help="input file with the summary table of microsats")
parser.add_argument("-g", help="genetic distance file")
parser.add_argument("-a", help="Locus amplification summary data from compare_pcrs.py")
parser.add_argument("-s", help="species translation file to go from 4 letter code to species name, also has # of loci in that sample's SSR primer set")
parser.add_argument("-o", help="output file")

args = parser.parse_args()

infile = args.i
distance_file= args.g
ampSummary = args.a
species_file=args.s
outfile = args.o

try:
	IN=open(infile, 'r')
except IOError:
	print "Can't open file", infile

try:
	DISTANCES=open(distance_file, 'r')
except IOError:
	print "Can't open file", distance_file

try:
	AMPSUMM=open(ampSummary, 'r')
except IOError:
	print "Can't open file", ampSummary


try:
	SPECIES=open(species_file, 'r')
except IOError:
	print "Can't open file", species_file

try:
	OUT=open(outfile, 'w')
except IOError:
		print "Can't open file", outfile
		

#Dictionary for looking up species names with 4-letter 1KP code.
Species_dict={}
#Dictionary for storing the number of loci that sample had.
Species_SSRCount_dict={}

for Line in SPECIES:
	Line = Line.strip('\n')
	Line_bits=re.split("\t",Line)
	
	try:
		Species_dict[Line_bits[0]]=Line_bits[1]
		Species_SSRCount_dict[Line_bits[0]]=Line_bits[2]
	except:
		print ("Error, couldn't add %s to Species_dict with value %s" %(Line_bits[0],Line_bits[1]))



#Load the amplification summary data into dictionaries
TotalLociDict={}
AmplifyDict={}
for Line in AMPSUMM:
	Line = Line.strip('\n')			
	Line_bits=re.split("\t",Line)		#JKNQ.on.DSVQ	3663	726
										# PCR			Total	Amp
										
	PCR='"' + Line_bits[0] + '"' #Add quotes around the PCR name to match as it is in input summary file.
	TotalLociDict[PCR]=Line_bits[1]
	AmplifyDict[PCR]=Line_bits[2]


Distance_dict={} # Dictionary to store the rows of the distance matrix as a list of values.
Column_dict={} #Dictionary to store which column the species is in for the distance matrix.
Line_number=1


#Load the distance matrix into a dictionary.
for Line in DISTANCES:
	if Line_number == 1:
		#The first line has the order of the species columns in the matrix
		Line = Line.strip('\n')
		Line_bits=re.split(",",Line)

		for i in range(len(Line_bits)):
			Column_dict[Line_bits[i]]=i-1  #Subtract 1 since we remove the species name in the list, so columns are 1 off of what we have here
		Line_number+=1
	else:
		Line = Line.strip('\n')
		Line_bits=re.split(",",Line, maxsplit=1) #only want to split the species name in the 1st column from the rest of the data, which will remain as a list
		
		Distance_dict[Line_bits[0]]=Line_bits[1].split(",")


for Line in IN:
	Line = Line.strip('\n')
	Line_bits=re.split(",",Line, maxsplit=1) #only want to split the ePCR name from the rest of the list.
	
	PCR_string=Line_bits[0].replace('\"','')
	
	PCR_bits=re.split("\.", PCR_string)

	#Get species of PCR Primer source
	Source_species=Species_dict[PCR_bits[0]]
	
	#Get species of amplified target species
	Target_species=Species_dict[PCR_bits[2]]
	
	#Get the distance from the matrix for these species
	PairDist=Distance_dict[Source_species][Column_dict[Target_species]]
	
	try:
		OUT.write("%s,%s,%s,%s,%s,%s,%s\n" %(Line_bits[0], Source_species, Target_species, PairDist, AmplifyDict[Line_bits[0]], Species_SSRCount_dict[PCR_bits[0]], Line_bits[1]))
	except:
		print "No in for: ", Line_bits[0]
	
	
	

	
	