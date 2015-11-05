#!/usr/bin/env python

## This script identifies SSR loci in a query file that occur in 
## translated regions by comparing the query file to a subject file
## of regions annotated as "CDS"

## Script written by Richie Hodel and Xiaoxian Liu

## The two inputs are the Query file (in this case Query.txt) and the
## Subject file (in this case CDS.txt)

import sys

## Read in the input files

QueryFileName = sys.argv[1] #Query input file
SubjectFileName = sys.argv[2] #Subject input file
	
## The user can assign a meaningful prefix to the output file

file_index = sys.argv[3] #beginning of output file name

## Initializing counters

fail_counter=0.0
hit_counter=0.0
total_counter=0.0

## The output file will be 'user prefix' + '_in_coding.txt'
OutFileName1 = file_index + "_in_coding.txt"
	
WriteOutFile = True
InFile1 =open(QueryFileName, 'r')
InFile2 =open(SubjectFileName, 'r')


if WriteOutFile:
	OutFile1 = open(OutFileName1, 'w')
	
## Go through Query file one line at a time, and capture QueryID, ScaffoldID,
## and location of primer

	for Line in InFile1:
		Line=Line.strip('\n')
		QueryList = Line.split('\t')	
		QueryScaffoldID = float(QueryList[2])
		QueryLow = int(QueryList[9])
		QueryHi = int(QueryList[10])
		QueryID = str(QueryList[0])
		total_counter=total_counter+1
		

		InFile2 =open(SubjectFileName, 'r')		

## For each line of the Query file, this for loop is run for the Subject file

		for Lines in InFile2:
			Lines=Lines.strip('\n')
			ElementList2 = Lines.split('\t')
			SubjectScaffoldID = float(ElementList2[0])
			SubjectLow = int(ElementList2[1])
			SubjectHi= int(ElementList2[2])
	
## If statement determines if the scaffold IDs match, and if the primer falls
## in a translated region
			
			if (QueryScaffoldID == SubjectScaffoldID and 
			QueryLow >= SubjectLow and QueryHi <= SubjectHi):					
				OutFile1.write(Line+"\n")
				hit_counter=hit_counter+1
			else: 
				fail_counter=fail_counter+1

## At the end of the outfile, indicated the percent of translated loci

	In_coding_percent = 100*(hit_counter/total_counter)		
	OutFile1.write("The percentage of SSRs in coding regions is ")
	OutFile1.write(str(In_coding_percent))

				
			
			
