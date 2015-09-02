## This scrpt is summarizing the ePCR summaries containing SSRs and flanking regions
## the args has the script run through all files listed in filenames, obtained by ls ePCR_summary/ > filenames in the folder /bio/soltis/1kp/SSRs
## authors: Charlotte Germain-Aubrey, FLMNH, with Francois Michonneau, FLMNH.
## September 2015

library(Hmisc)
library(Biostrings)

args=(commandArgs(TRUE))

FinalTable <- matrix(data=NA, nrow = 0, ncol = 15)
colnames(FinalTable) <- c("PCR_name", "microsat_length_same", "microsat_length_diff", "flanking1_length_same",  "flanking1_length_diff", "flanking1_same", "flanking1_diff", "flanking2_length_same",  "flanking2_length_diff", "flanking2_same", "flanking2_diff","all_sequence_length_same", "all_sequence_length_diff", "all_sequence_same", "all_sequence_diff")

#table <- read.table("AJBK.on.AQXA.ePCR.summary", header=FALSE, sep="\t")
file_name <- args()
table <- read.table(file_name, header=FALSE, sep="\t")
colnames(table) <- c("name1", "name2", "seq1", "seq2", "micro")
table$micro <- as.character(table$micro)
table$seq1 <- as.character(table$seq1)
table$seq2 <- as.character(table$seq2)

#extract name of ePCR
parts <- strsplit(as.character(file_name[1]), split=".")
one <- unlist(lapply(parts,"[",1))
two <- unlist(lapply(parts, "[",2))
three <- unlist(lapply(parts, "[",3))
PCR_name <- paste(one, two, three, sep=".")
###########################
## my code

seq1 <- as.character(table$seq1)

## extract motif by replacing everthing that is not ATCG
motif <- gsub("[^actg]", "", table$micro, ignore.case = TRUE)
table <- cbind(table, motif)


n_repeat <- mapply(function(m, p){
                       ## identify where the motif occurs at least 2 times
                       pos <- gregexpr(paste0("(", m, "){2,}"), p)
                       ## where does it produce the longest repeat?
                       correct_pos <- which.max(attr(pos[[1]], "match.length"))
                       to_repl <- regmatches(p, pos)[[1]][correct_pos]
                       ## flank the longest repeat with |
                       gsub(paste0("(", to_repl, ")"), "|\\1|", p)
                   }, motif, seq1)

## cut each sequence where the | occur
cut_seq <- strsplit(n_repeat, "|", fixed = TRUE)

## a little better:
table(sapply(cut_seq, length))

## add parts to table:
seq1a <- unlist(lapply(cut_seq, "[", 1))
micro1 <- unlist(lapply(cut_seq, "[", 2))
seq1b <- unlist(lapply(cut_seq, "[", 3))

## do the same for sequence of 2nd individual
seq2 <- as.character(table$seq2)
n_repeat2 <- mapply(function(m, p){
  pos <- gregexpr(paste0("(", m, "){2,}"), p)
  correct_pos <- which.max(attr(pos[[1]], "match.length"))
  to_repl <- regmatches(p, pos)[[1]][correct_pos]
  gsub(paste0("(", to_repl, ")"), "|\\1|", p)
}, motif, seq2)
cut_seq2 <- strsplit(n_repeat2, "|", fixed = TRUE)
table(sapply(cut_seq2, length))
seq2a <- unlist(lapply(cut_seq2, "[", 1))
micro2 <- unlist(lapply(cut_seq2, "[", 2))
seq2b <- unlist(lapply(cut_seq2, "[", 3))


## put everything together in table2
table2 <- cbind(table, seq1a, micro1, seq1b, seq2a, micro2, seq2b)
## get rid of NA rows
table2 <- na.omit(table2)
table3 <- table2[!table2$micro2 == "NA" , ]

micro1 <- table3$micro1
micro2 <- table3$micro2
motif <- table3$motif
micro1 <- as.character(micro1)
micro2 <- as.character(micro2)
motif<- as.character(motif)

  
z <- nrow(table3)
r <- table(1,0)
d <- table(1,0)
## now compare the microsat repeats
for (i in 1:z){
## how many times is the motif in micro1 ?
lr <- nchar(motif[i])
r1<- nchar(micro1[i])/lr
r2 <- nchar(micro2[i])/lr

## are the microsats the same length ? 
n <- r1==r2
q <- cbind(q, n)
}
##how many are the same length ? 

mt <- length(which(q == "TRUE"))
mf <- length(which(q == "FALSE"))

## Now do the same with sequences

seq1a <- table3$seq1a
seq1a <- as.character(seq1a)
seq1b <- table3$seq1b
seq1b <- as.character(seq1b)
seq2a <- table3$seq2a
seq2a <- as.character(seq2a)
seq2b <- table3$seq2b
seq2b <- as.character(seq2b)


l <- list()
s <- list()
for (i in 1:z){
  ## for the first franking region
  s1 <- seq1a[i]
  l1 <- nchar(seq1a[i])
  s2 <- seq2a[i]
  l2 <- nchar(seq2a[i])
 # are the sequences the same length ? 
  m <- l1 == l2
  ## are the sequences the same  ?
  n <- s1 == s2
  s <- cbind(s, n)
  l <- cbind(l, m)
}
##how many are the same  ? 
st <- length(which(s == "TRUE"))
sf <- length(which(s == "FALSE"))
## how many are the same length ? 
lt <- length(which(l == "TRUE"))
lf <- length(which(l == "FALSE"))

## Same thing for flanking sequence 2 (5'end)
v <- list()
w <- list()
for (i in 1:z){
  ## for the 5' franking region
  s3 <- seq1b[i]
  l3 <- nchar(seq1b[i])
  s4 <- seq2b[i]
  l4 <- nchar(seq2b[i])
  ## are the sequences the same  ? 
  x <- s3 == s4
  y <- l3 == l4
  v <- cbind(v, x)
  w <- cbind(w, y)
}
##how many are the same  ? 
kst <- length(which(v == "TRUE"))
ksf <- length(which(v == "FALSE"))
klt <- length(which(w == "TRUE"))
klf <- length(which(w == "FALSE"))

## check if whole sequences are the same for any individual)

seq1 <- table3$seq1
seq1 <- as.character(seq1)
seq2 <- table3$seq2
seq2 <- as.character(seq2)

a <- list()
b <- list()
for (i in 1:z){
  s5 <- seq1[i]
  l5 <- nchar(seq1[i])
  s6 <- seq2[i]
  l6 <- nchar(seq2[i])
  ## are the sequences the same  ? 
  c <- s5 == s6
  d <- l5 == l6
  a <- cbind(a, c)
  b <- cbind(b, d)
}
##how many are the same  sequence? 
at <- length(which(a == "TRUE"))
af <- length(which(a == "FALSE"))
## how many are the same length ? 
bt <- length(which(d == "TRUE"))
bf <- length(which(d == "FALSE"))

##now fill in table
FinalTable[1,1] <- PCR_name
FinalTable[1,2] <- mt
FinalTable[1,3] <- mf
FinalTable[1,4] <-lt
FinalTable[1,5] <-lf
FinalTable[1,6] <-st
FinalTable[1,7] <-sf
FinalTable[1,8] <-klt
FinalTable[1,9] <-klf
FinalTable[1,10] <-kst
FinalTable[1,11] <-ksf
FinalTable[1,12] <- bt
FinalTable[1,13] <- bf
FinalTable[1,14] <- at
FinalTable[1,15] <- af

table_name <- paste("output/", PCR_name, "_out", sep="" )
write.csv(FinalTable, table_name, sep = "," )