### This code will analyze the output from the ePCR from Oenothera and plot regression lines from the different columns
### Author: Charlotte Germain-Aubrey

library(scatterplot3d)
library(Rcmdr)
library(plotrix)

#load file and set working directory
setwd("~/Desktop/Microsat_review_paper/MyTrials/summary/")
dat <- read.csv("Oenothera_10geneDistances.csv", header=TRUE)

#here we determine if the ePCR was on the same species of a different one
same <- subset(dat, dat$Primer_Source == dat$PCR_Target)
diff <- subset(dat, dat$Primer_Source != dat$PCR_Target)

#proportion of  regions being the same
z <- ((dat$flanking1_same/(dat$flanking1_same + dat$flanking1_diff)) + (dat$flanking2_same/(dat$flanking2_same + dat$flanking2_diff))/2) 
fl <- ((dat$flanking1_length_same/(dat$flanking1_same + dat$flanking1_diff)) + (dat$flanking2_length_same/(dat$flanking2_same + dat$flanking2_diff))/2) 
m <- (dat$microsat_length_same/(dat$microsat_length_same + dat$microsat_length_diff)) 
all.length <- (dat$all_sequence_length_same/(dat$all_sequence_length_same + dat$all_sequence_length_diff)) 
all <- (dat$all_sequence_same/(dat$all_sequence_same + dat$all_sequence_diff))


ID <- dat$Genetic_Identity
total <- dat$Number.Loci.from.Source
amp <- (dat$microsat_length_same+dat$microsat_length_diff)
prop <- amp/total

#plot lines of similarity by the genetic ID (distance)
pdf("TrendLines10.pdf")
par(xpd=FALSE)
plot(ID,prop,pch='.', xlab="Genetic Identity", ylab="Proportion", col="grey", cex=0.5)
abline(lm(z ~ ID), col = "Purple")
abline(lm(fl ~ ID), col = "Blue")
abline(lm(m ~ ID), col = "Green")
abline(lm(all ~ ID), col = "Yellow")
abline(lm(all.length ~ ID), col = "Red")
legend('topleft', bty = 'n', c("Proportion of loci amplifying", "Flanking regions same", "Flanking regions length same", "Microsatellite same", "PCR sequence same", "PCR sequence length same"), pch = c(15,  NA, NA, NA, NA, NA), col = c("grey", "Purple", "Blue", "Green", "Yellow", "Red"), border=NA, lty = c(NA, 1, 1, 1, 1, 1, 1), x.intersp = 0.8, y.intersp = 1, cex=1)
dev.off()

