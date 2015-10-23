
#############################################################################

##     This R script takes a list of loci known to be in coding regions
##       (based on a BLAST search), and identifies the distribution of
##       repeat motifs using the output from PAL_FINDER.
##     The input for this script are output files from CodinSSR.py and 
##       PAL_FINDER.
##     The output from CodinSSR.py is a tab-delimited .txt file
##       with 13 columns and one final line that is trimmed by this script.
##     The output from PAL_FINDER is a tab-delimited .txt file
##       with 14 columns.

##############  This section determines the motif type 

Loci <- read.table("coding_loci_input.txt", fill=TRUE)
Loci <- head(Loci, -1)

All_PAL <- read.csv("PAL_summary_input.csv")

All_PAL <- na.omit(All_PAL)

Loci$V1 <- NULL
Loci$V3 <- NULL
Loci$V4 <- NULL
Loci$V5 <- NULL
Loci$V6 <- NULL
Loci$V7 <- NULL
Loci$V8 <- NULL
Loci$V9 <- NULL
Loci$V10 <- NULL
Loci$V11 <- NULL
Loci$V12 <- NULL
Loci$V13 <- NULL

Loci$R.Primer.Name <- Loci$V2

Joined <- merge(Loci, All_PAL)
nrow(Joined)

Dinuc <- subset(Joined, Joined$Repeat.Motif.Size == 2)
Trinuc <- subset(Joined, Joined$Repeat.Motif.Size == 3)
Tetranuc <- subset(Joined, Joined$Repeat.Motif.Size == 4)
Pentanuc <- subset(Joined, Joined$Repeat.Motif.Size == 5)
Hexanuc <- subset(Joined, Joined$Repeat.Motif.Size == 6)

Sum <- nrow(Dinuc)+nrow(Trinuc)+nrow(Tetranuc)+nrow(Pentanuc)+nrow(Hexanuc)

Prop2 <- 100*(nrow(Dinuc)/Sum)
Prop3 <- 100*(nrow(Trinuc)/Sum)
Prop4 <- 100*(nrow(Tetranuc)/Sum)
Prop5 <- 100*(nrow(Pentanuc)/Sum)
Prop6 <- 100*(nrow(Hexanuc)/Sum)

Coding_Motifs <- data.frame(row.names=NULL, Motif=c("2", "3", "4", "5", "6"),
                            Count=c(nrow(Dinuc), nrow(Trinuc), nrow(Tetranuc), 
                                    nrow(Pentanuc), nrow(Hexanuc)), Percent=c(Prop2,
                                    Prop3, Prop4, Prop5, Prop6))

write.table(Coding_Motifs, file="Coding_motif_output.txt", row.names=FALSE, quote=FALSE,
            sep="\t")


#############################################################################


