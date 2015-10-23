
#############################################################################


Genome <- read.table("GCF_000004515.3_V1.1_genomic.gff", sep = "\t", skip=7,
                     quote="")


Trimmed <- Genome

Trimmed$V2 <- NULL
Trimmed$V6 <- NULL
Trimmed$V7 <- NULL
Trimmed$V8 <- NULL
Trimmed$V9 <- NULL


CDS <- subset(Trimmed, Trimmed$V3 == "CDS")
CDS <- subset(CDS, !duplicated(CDS$V4))
CDS$V3 <- NULL

install.packages("tidyr")
require(tidyr)

CDS <- separate(CDS, V1, c('V2', 'V3'), sep="_")
CDS$V3 <- gsub("[^[:alnum:]///' ]", "", CDS$V3)
CDS$V2 <- NULL
write.table(CDS, "CDS.txt", row.names=FALSE, quote=FALSE, 
            col.names=FALSE, sep="\t")


Search1 <- read.csv("G_R_0764-0510-HitTable.csv", header=FALSE)
List1 <- subset(Search1, !duplicated(Search1$V1))

Search2 <- read.csv("G_R_0803-0253-HitTable.csv", header=FALSE)
List2 <- subset(Search2, !duplicated(Search2$V1))

List <- rbind(List1, List2)

List <- separate(List, V2, c('V2_1', 'V2'), sep="_")
List$V2_1 <- NULL
List$V2 <- gsub("[^[:alnum:]///' ]", "", List$V2)

write.table(List, "Query.txt", quote=FALSE, 
            col.names=FALSE, sep="\t")


#############################################################################



