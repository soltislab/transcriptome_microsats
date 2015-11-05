
#############################################################################

## BLAST_to_Coding_SSR.R
## Script written by Richie Hodel

## This script uses a .gff file (annotated Glycine max genome from NCBI), and
## a BLAST report for SSR Loci blasted against the Glycine max genome to 
## prepare two files, which will be used in a subsequent script 
## (Coding_SSR.py) to determine which loci are in translated regions of
## the genome (i.e., regions that are annotated as "CDS")

## The output of this script is two files (Query.txt, CDS.txt).  Query.txt 
## contains the SSR loci identified from the BLAST search, with some
## unncessary columns and duplicates removed.  CDS.txt contains the regions 
## of the Glycine max genome that are annotated as "CDS".  These two files
## are the two input files used in the script Coding_SSR.py

#############################################################################


## Read in the Glycine max genome .gff file
## Another file name could be substituted here is the user has a different
## taxon of interest

Genome <- read.table("GCF_000004515.3_V1.1_genomic.gff", sep = "\t", skip=7,
                     quote="")

## Remove unnecessary columns

Trimmed <- Genome
Trimmed$V2 <- NULL
Trimmed$V6 <- NULL
Trimmed$V7 <- NULL
Trimmed$V8 <- NULL
Trimmed$V9 <- NULL

## Remove entries from the annotated genome that are not listed as "CDS"

CDS <- subset(Trimmed, Trimmed$V3 == "CDS")
CDS <- subset(CDS, !duplicated(CDS$V4))
CDS$V3 <- NULL


## Installing necessary packages

install.packages("tidyr")
require(tidyr)


## Trimming characters so that the scaffold names are purely numeric

CDS <- separate(CDS, V1, c('V2', 'V3'), sep="_")
CDS$V3 <- gsub("[^[:alnum:]///' ]", "", CDS$V3)
CDS$V2 <- NULL

## Write one of the output files

write.table(CDS, "CDS.txt", row.names=FALSE, quote=FALSE, 
            col.names=FALSE, sep="\t")

## Read in the BLAST results and trim duplicates

Search1 <- read.csv("G_R_0764-0510-HitTable.csv", header=FALSE)
List1 <- subset(Search1, !duplicated(Search1$V1))

Search2 <- read.csv("G_R_0803-0253-HitTable.csv", header=FALSE)
List2 <- subset(Search2, !duplicated(Search2$V1))

## Merge two BLAST hit tables.  It was sometimes necessary to subset the 
## data so the BLAST searches would run more efficiently.  The data are now
## being merged.

List <- rbind(List1, List2)

## Trimming characters so that the scaffold names are purely numeric

List <- separate(List, V2, c('V2_1', 'V2'), sep="_")
List$V2_1 <- NULL
List$V2 <- gsub("[^[:alnum:]///' ]", "", List$V2)

## Write the other output file

write.table(List, "Query.txt", quote=FALSE, 
            col.names=FALSE, sep="\t")


#############################################################################



