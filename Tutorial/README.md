#Tutorial on running the various scripts included here.

This tutorial will walk you through the basic steps of designing microsatellites from transcriptome assemblies with the scripts developed for Hodel *et al.* (2016). The sample data are in the Tutorial folder and are a subset of the transcriptome of *Amborella trichopoda*, taking the 1st 1,000 scaffolds from the file available [here](http://mirrors.iplantcollaborative.org/browse/iplant/home/shared/onekp_pilot/taxa/AMBO-Amborella_trichopoda/assemblies/AMBO.faa) (Matasci *et al.* 2015).  

This tutorial assumes you have a linux like system and have Python 2.7 and Perl installed. Including the following Python packages:
 argparse

MISA can be downloaded from [here](http://pgrc.ipk-gatersleben.de/misa/). Modified versions of the primer3 scripts are in the MISA_changes directory of this repository.

For Step 7, where open reading frames (ORFs) are identified, you will need the script [get_orfs_or_cdss.py](https://github.com/peterjc/pico_galaxy/blob/master/tools/get_orfs_or_cdss/get_orfs_or_cdss.py)

## Step 1: Clone the repository

The first step is to get the scripts, data and other files associated with this repository.

>git clone https://github.com/soltislab/transcriptome_microsats.git

## Step 2: Change directories to the Tutorial directory

To run the tutorial, change directories in to the Tutorial directory your clone of the repository.

>cd transcriptome_microsats/Tutorial

## Step 3: Run MISA
Now we will run MISA. You will need to substitute PATH with the path to your install of MISA.
The misa.ini configuration file used for this study is in the Tutorial directory.

> perl PATH/misa.pl AMBO.1K.fna

This should generate two files. AMBO.1k.fna.misa with 691 lines (should be the sample as Tutorial.AMBO.1k.fna.misa) and AMBO.1K.fna.statistics which should start with the following information:

```
Specifications
==============

Sequence source file: "AMBO.1K.fna"

Definement of microsatellites (unit size / minimum number of repeats):
(1/10) (2/6) (3/5) (4/5) (5/5) (6/5) 

Maximal number of bases interrupting 2 SSRs in a compound microsatellite:  25



RESULTS OF MICROSATELLITE SEARCH
================================

Total number of sequences examined:              1000
Total size of examined sequences (bp):           1816420
Total number of identified SSRs:                 757
Number of SSR containing sequences:              421
Number of sequences containing more than 1 SSR:  185
Number of SSRs present in compound formation:    67

```

The whole file is 115 lines and should be identical to Tutorial.AMBO.1K.fna.statistics.


## Step 4: Run the modified script p3_in.MAG.pl
Run the modified version of the p3_in.pl script. This generates the input file needed for primer3.

>perl ../MISA_changes/p3_in.MAG.pl AMBO.1K.fna.misa 

The output should say:

    690 records created.

And the file AMBO.1K.fna.p3in will be created, which should have 4,140 lines and be identical to Tutorial.AMBO.1K.fna.p3in 

## Step 5: Run primer3

Now you will run primer3 using the input file generated above.

> primer3 < AMBO.1K.fna.p3in > AMBO.1K.fna.p3out

This will generate a primer3 output file with 70,208 lines, identical to Tuturial.AMBO.1K.fna.p3out.

## Step 6: Run the modified script p3_out.MAG.pl

Next you will take this primer3 output and use the modified version of p3_out.pl to summarize the primers for microsatellites.

> perl ../MISA_changes/p3_out.MAG.pl AMBO.1K.fna.p3out AMBO.1K.fna.misa 

The output should say:
    Primer modelling was successful for 576 sequences.
    Primer modelling failed for 114 sequences.

And two files should be created:
 1. AMBO.1K.fna.results, with 577 lines, identical to Tutorial.AMBO.1K.fna.results
 2. AMBO.1K.fna.goodSeqs, with 576 lines, identical to Tutorial.AMBO.1K.fna.goodSeqs
 
 The AMBO.1K.fna.results file has the primers designed for the 576 microsatellite loci. Some of these may have been designed from inferred isoforms in the transcriptome assembly, and are therefore really from the same locus.
 

## Step 7: Filter Loci with same priming site
While not a perfect filter, one easy way to identify loci that are developed from inferred isoforms is to filter out loci that use the same priming site as a previous locus. The script SSR_RepeatFilter.py does this. It takes as it's input the results file from step 6 and the a name, in this case, we keep the same AMBO.1kp.fna used throughout the tutorial.

>python ../SSR_RepeatFilter.py AMBO.1K.fna.results AMBO.1K.fna

This results in two files:
1. AMBO.1K.fna.repeat.txt, the information of loci that were identified as repeats and removed from the results. This file should have 62 lines and be identical to Tutorial.AMBO.1K.fna.repeat.txt.'
2. AMBO.1K.fna.misa_SSR.results, the retained loci. This file should have 516 lines and be identical to Tutorial.AMBO.1K.fna.misa_SSR.results. These are the 515 good microsatellite loci that could be used for primer testing and screening for variation. 

The next steps will explore where the microsatellite loci developed above lie relative to the longest open reading frame on the scaffold. This will allow you to filter out the loci that are likely to be in a coding region and are more likely to be under selection.


##Step 8: Locate the SSRs and determine their location relative to open reading frames (ORFs).

The script get_orfs_or_cdss.py takes the transcriptome and finds the longest ORF on each scaffold. First we need to trim the list of scaffolds to only include the ones with good microsatellites.

> python ../getScaffolds.py -i AMBO.1K.fna.goodSeqs -s AMBO.1K.fna -o AMBO.1K.fna.goodScaffolds

This should create a file with 13,530 lines identical to AMBO.1K.fna.goodScaffolds.

Next find the location of the longest ORF in these scaffolds using get_orfs_or_cdss.py (you need to install this):

> python get_orfs_or_cdss.py AMBO.1K.fna.goodScaffolds fasta 1 ORF open top 10 both AMBO.1K.ORFS.fna AMBO

Which should produce the output:
    Genetic code table 1
    Minimum length 10 aa
    Found 367 ORFs in 365 sequences

And the files AMBO.1K.ORFS.fna and AMBO.1K.ORFS.faa with 7,113 and 2,739 lines respectively and identical to Tutorial.AMBO.1K.ORFS.fna and Tutorial.AMBO.1K.ORFS.faa.

Now we can summarize the location of the microsatellites relative to the identified ORFs using LocateSSRsandORFs.py:

> python ../LocateSSRsandORFs.py -i AMBO.1K.ORFS.faa -r AMBO.1K.fna.misa_SSR.results -o AMBO.1K.Locations -s SSR_ORF_summary.txt

Note that this uses the default of 15 bp of overlap with an ORF as the cutoff for reporting. To change this value, add the -l flag with a number of bp to use. The script reports three sets of overlap results: no overlap, any overlap and overlap greater than the cutoff.

Which creates two files:
1. SSR_ORF_summary.txt (the script will append to this file, so you can run on multiple files and accumulate results). The file will look like:

    AMBO.1K.ORFS.faa	515	365	367	0.712621359223		42	2	95	327	40	1	2	6	515		35	1	83	274	13	1	2	2	411		7	1	12	53	27	0	0	4	104	23	13	0	0	3	49

The columns are:
```
1. AMBO.1K.ORFS.faa: Name of the input file
2. 515: Number of microsatellite loci
3. 365: Number of scaffolds with a detected ORF
4. 367: Number of ORFs (some have more than 1 ORF of same length)
5. 	0.712621359223: Fraction of SSR loci on scaffolds with ORFs
6. BLANK
7-15. 42	2	95	327	40	1	2	6	515: The numbers of loci of each type in order: c, c*, p1, p2, p3, p4, p5, p6
16. BLANK
17-25: 35	1	83	274	13	1	2	2	411: The numbers of loci of each type that do not overlap the ORF in order: c, c*, p1, p2, p3, p4, p5, p6
26: BLANK
7-34: 7	1	12	53	27	0	0	4: The numbers of loci of each type that have any overlap with an ORF in order: c, c*, p1, p2, p3, p4, p5, p6
35: BLANK
36-42: 104	23	13	0	0	3	49: The numbers of loci of each type that have more than the designated overlap with the ORF in order: c, c*, p1, p2, p3, p4, p5, p6
```

2. AMBO.1K.Locations with the columns:
  ```
1. Scaffold Name|ORF number

2. Type of microsatellite repeat (c, c*, p1, p2, p3, p4, p5, or p6)

3. Repeat sequence

4. BP start of the microsatellite

5. BP end of the microsatellite

6. BP stat of the ORF

7. BP end of the ORF

8. Length of overlap between microsatellite and ORF in BP
```


From these files you can now decide which microsatellite loci you wish to order primers for and attempt to amplify.



## Tutorial for evaluating if loci are in translated regions using Glycine

This tutorial will walk the user through generating a list of PALs (potentially 
amplifiable loci) for genomic (as opposed to transcriptomic) microsatellites, running 
scripts and a BLAST search on those PALs to determine if they are in translated regions 
of the genome, and obtaining the distribution of repeat number motifs for the loci.


This section uses the software PAL_FINDER, the annotated genome of Glycine max, 
and custom scripts described below.  

###1. Go to the main directory of your cloned repository

  >cd transcriptome_microsats/

###2. Download and install PAL_FINDER 
	(http://sourceforge.net/projects/palfinder/)
	
	The file "README.txt" will guide the user on how to use PAL_FINDER
	
	The required input for PAL_FINDER is a config file (default config.txt)
	
###3. Run PAL_FINDER
	
	We use the file "config_g1.txt" to run PAL_FINDER:
	
	>perl pal_finder_v0.02.04.pl config_g1.txt
	
	PAL_FINDER will output: a set of PALs in a .txt file
	The file "P_to_B_input.txt" is a sample file that would
	be the output of running PAL_FINDER, and is ready to use
	in the next step.
	
###4. Run the script PAL_to_BLAST.R
	You can use the provided file "P_to_B_input.txt" or use 
 	an output file from PAL_FINDER.  If you use a different file,
	you need to change the name of the file read into the variable
	'Input' in the script.
	If you run this script, it  will output two files: 'Forward.fasta' and 'Reverse.fasta'
	These two files are already included in the main directory of the repository.

###5. Perform a BLASTn search by querying the two files from the previous 
	step against the Glycine max annotated genome (NCBI Assembly Accession
	number GCF_000004515.3).  If you run BLASTn through the NCBI website, you
	may need to divide the two fasta files into several smaller files to
	ensure that the job will finish before it runs out of computing 
	resources.  Alternatively, you could build a local BLAST database to
	perform the search, or build a BLAST database on a computing resource
	such as galaxy, which would be able to perform a BLAST search on the
	entire file.  
		
	The files we obtained from the BLASTn search are "G_R_0764-0510-HitTable.csv"
	and "G_R_0803-0253-HitTable.csv", which are available at:
	https://gatorbox.rc.ufl.edu/index.php/s/361JoKPasE00yql

###6. Run BLAST_to_Coding_SSR.R (requires R module tidyr)
	This script requires three input files: the two BLAST hit tables from the 
	previous step, and the Glycine max annotated genome file:
	"GCF_000004515.3_V1.1_genomic.gff" available at 
	https://gatorbox.rc.ufl.edu/index.php/s/361JoKPasE00yql
	
	This script will output the files "CDS.txt" and "query.txt"

###7. Run Coding_SSR.py

	This script takes the files "CDS.txt" and "query.txt" as input.  The user will be
	asked to enter a meaningful prefix, and the output file will be called:
	"<prefix>_in_coding.txt".  For example, if your prefix is "test", then your 
	output file would be called "test_in_coding.txt"
	
	>python Coding_SSR.py Query.txt CDS.txt test

###8. Run Repeat_Numb_Search.R

  The input file for this script are "test_in_coding.txt", which is the output of the
  script "Coding_SSR.py", and "PAL_summary_input.txt", which is the output of running
  PAL_FINDER.  This script will produce an output file "Coding_motif_output.txt"





#Citations

Matasci, N., L.-H. Hung, Z. Yan, E. J. Carpenter, N. J. Wickett, S. Mirarab, N. Nguyen, T. Warnow, S. Ayyampalayam, M. Barker, J. G. Burleigh, M. A. Gitzendanner, E. Wafula, J. P. Der, C. W. dePamphilis, B. Roure, H. Philippe, B. R. Ruhfel, N. W. Miles, S. W. Graham, S. Mathews, B. Surek, M. Melkonian, D. E. Soltis, P. S. Soltis, C. Rothfels, L. Pokorny, J. A. Shaw, L. DeGironimo, D. W. Stevenson, J. C. Villarreal, T. Chen, T. M. Kutchan, M. Rolf, R. S. Baucom, M. K. Deyholos, R. Samudrala, Z. Tian, X. Wu, X. Sun, Y. Zhang, J. Wang, J. Leebens-Mack, and G. K.-S. Wong. 2014. Data access for the 1,000 Plants (1KP) project. GigaScience 3:17. [Link](http://www.gigasciencejournal.com/content/3/1/17/abstract)

