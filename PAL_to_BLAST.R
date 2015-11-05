
#############################################################################

##            PAL_to_BLAST.R     written by Richie Hodel

#############################################################################

##     This R script will take a .txt file of potential SSR loci,
##       remove duplicates, and output a fasta file for a BLAST search.
##     The input for this script is an output file from PAL_FINDER.
##     The output from PAL_FINDER is a tab-delimited .txt file
##       with 14 columns.


##     The default input file (P_to_B_input.txt) is in the line of code 
##       immediately following this comment. Replace this file name with
##       your input file.


Input <- read.table(file = "P_to_B_input.txt", header = TRUE, sep = "\t")


# This section removes duplicate loci from the data set

Input_1 <- na.omit(Input)
Input_2 <- subset(Input_1, !duplicated(Input_1$Forward.Primer))
Input_3 <- subset(Input_2, !duplicated(Input_2$Reverse.Primer))
Input_4 <- subset(Input_3, !duplicated(Input_3$SequenceID))
Filtered <- Input_4


## This section is sorting the different motif types

Loci2 <- which (Filtered$Repeat.Motif.Size == 2)
Loci3 <- which (Filtered$Repeat.Motif.Size == 3)
Loci4 <- which (Filtered$Repeat.Motif.Size == 4)
Loci5 <- which (Filtered$Repeat.Motif.Size == 5)
Loci6 <- which (Filtered$Repeat.Motif.Size == 6)

Total <- length(Loci2)+length(Loci3)+length(Loci4)+length(Loci5)+length(Loci6)

## This section calculates the distribution of motifs and writes a file.
## The default file is called "Motif_output.txt"

Prop2 <- 100*(length(Loci2)/Total)
Prop3 <- 100*(length(Loci3)/Total)
Prop4 <- 100*(length(Loci4)/Total)
Prop5 <- 100*(length(Loci5)/Total)
Prop6 <- 100*(length(Loci6)/Total)

Motifs <- data.frame(row.names=NULL, Motif=c("2", "3", "4", "5", "6"),
                     Count=c(length(Loci2), length(Loci3), length(Loci4), 
                             length(Loci5), length(Loci6)), Percent=c(Prop2,
                             Prop3, Prop4, Prop5, Prop6))


write.table(Motifs, file="Motif_output.txt", row.names=FALSE, quote=FALSE,
            sep="\t")

# This section removes columns that are not going to be needed in the FASTA output

Filtered$Occurances.of.Forward.Primer.in.Reads <- NULL
Filtered$Repeat.Motif <- NULL
Filtered$Repeat.Motif.Size <- NULL
Filtered$Total.Repeats.In.Amplicon   <- NULL
Filtered$Primer.Designed..1.y.0.n. <- NULL
Filtered$Number.Tandem.Repeats <- NULL
Filtered$Occurances.of.Reverse.Primer.in.Reads <- NULL
Filtered$Occurances.of.Amplifiable.Primer.Pair.in.Reads <- NULL
Filtered$Occurances.of.Amplifiable.Primer.Pair.in.PALs <- NULL
Filtered$row.names   <- NULL
Filtered$SequenceID <- NULL

# This section produces one fasta file for the forward primers, and 
# one for the reverse primers

DataF <- Filtered
DataR <- Filtered

DataF$R.Primer.Name <- NULL
DataF$Reverse.Primer <- NULL

DataR$F.Primer.Name <- NULL
DataR$Forward.Primer <- NULL

DataF$F.Primer.Name <- paste(">", DataF$F.Primer.Name)
DataR$R.Primer.Name <- paste(">", DataR$R.Primer.Name)


## The following steps are when the fasta files are being written.  
## You can change the file name(s) as desired.  

write.table(DataF, file="Forward.fasta", row.names=FALSE, sep="\n",
            col.names=FALSE, quote=FALSE, eol="\n")

write.table(DataR, file="Reverse.fasta", row.names=FALSE, sep="\n",
            col.names=FALSE, quote=FALSE, eol="\n")

#############################################################################


