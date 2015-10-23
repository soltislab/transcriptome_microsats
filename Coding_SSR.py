#!/usr/bin/env python

## Match Query and Subject ##
## Script written by Richie Hodel and Xiaoxian Liu

import sys

QueryFileName = sys.argv[1] #Query input file
SubjectFileName = sys.argv[2] #Subject input file
	
file_index = sys.argv[3] #beginning of output file name

fail_counter=0.0
hit_counter=0.0
total_counter=0.0
	
OutFileName1 = file_index + "_in_coding.txt"
	
WriteOutFile = True
InFile1 =open(QueryFileName, 'r')
InFile2 =open(SubjectFileName, 'r')
	
if WriteOutFile:
	OutFile1 = open(OutFileName1, 'w')
	
	
	for Line in InFile1:
		Line=Line.strip('\n')
		QueryList = Line.split('\t')	
		QueryScaffoldID = float(QueryList[2])
		QueryLow = int(QueryList[9])
		QueryHi = int(QueryList[10])
		QueryID = str(QueryList[0])
		total_counter=total_counter+1
		

		InFile2 =open(SubjectFileName, 'r')		
				
		for Lines in InFile2:
			Lines=Lines.strip('\n')
			ElementList2 = Lines.split('\t')
			SubjectScaffoldID = float(ElementList2[0])
			SubjectLow = int(ElementList2[1])
			SubjectHi= int(ElementList2[2])
			
			#print "Subject ID", SubjectScaffoldID
			#print "S low", SubjectLow
			#print "S hi", SubjectHi				
			
			if (QueryScaffoldID == SubjectScaffoldID and 
			QueryLow >= SubjectLow and QueryHi <= SubjectHi):					
				OutFile1.write(Line+"\n")
				hit_counter=hit_counter+1
			#	print Lines
			else: 
				fail_counter=fail_counter+1
				
	In_coding_percent = 100*(hit_counter/total_counter)		
	OutFile1.write("The percentage of SSRs in coding regions is ")
	OutFile1.write(str(In_coding_percent))

				
			
			