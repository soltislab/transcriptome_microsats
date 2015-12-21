# transcriptome_microsats
Scripts used in evaluating microsatellite loci developed from transcriptomes

##Overview

These scripts were used to process the output of the microsatellite finding program MISA ([Dieringer and Schlötterer, 2003](http://onlinelibrary.wiley.com/doi/10.1046/j.1471-8286.2003.00351.x/abstract)). MISA was run on the transcriptome datasets from the [1,000 Plant Transcriptome Project](http://onekp.com) and the results summarized with these and previously released scripts. 


Please see the [Tutorial](Tutorial/README.md) for a walkthrough of developing microsatellite loci from a sample dataset and more complete usage and output information.


What follows is a brief description of the scripts available here.

##Scripts and their functions:

###LocateSSRsandORFs.py
A script to compare the location of an SSR and the longest ORF in a given scaffold.
Note that for the project, ORFs were identified using [get_orfs_or_cdss.py](https://github.com/peterjc/pico_galaxy/blob/master/tools/get_orfs_or_cdss/).

###getScaffolds.py
A script to pull sequences from a fasta file based on a list and write to a new file.

###SSR_RepeatFilter.py
Filters out microsatellite loci likely derived from inferred isoforms of the same locus.

###compare_pcrs.py
Script to compare the ePCR results of amplifying one sample with the primers developed from another sample.

###add_distances.py
Add genetic distance to table of microsats. 

###DiffMicrosatsBatch.R
This scrpt is summarizes the ePCR results. It generates a summary containing SSR and flanking region differences.

###Graphs10GeneticID.R
Generates figure ## of Hodel et al 2015.

###PAL_to_BLAST.R
This script takes a set of PALs (potentially amplifiable loci) from PAL_FINDER
and removes duplicates. It then outputs a fasta file of the loci that can be used in a BLAST search.

###BLAST_to_Coding_SSR.R
This script takes in a .gff annotated genome file from NCBI and BLAST results
and creates the two input files for the script Coding_SSR.py

###Coding_SSR.py
This script compares two input files--one ("query") containing the genomic
locations of SSR loci identified, and one ("subject") containing the locations of translated regions of the genome (CDS).

###Repeat_Numb_Search.R
This script uses the output from the script Coding_SSR.py and the program
PAL_FINDER to identify the distribution of repeat motifs among the PALs.

