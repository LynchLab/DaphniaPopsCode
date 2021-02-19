setwd("/scratch/rraborn/")
#load("/data/LynchLabCME/Daphnia/DaphPopsTSRs/binaries/PdSTRIPE_complete.rdata")
library(dplyr)
library(stringr)
library(reshape2)

experimentName <- PdSTRIPE
tagCountThreshold <- 2

tss.start <- PdSTRIPE@tssCountData[[19]]
tssRepsCounts$test2 <- PdSTRIPE@tssCountData[[20]]
tssRepsCounts$test3 <- PdSTRIPE@tssCountData[[21]]
tssRepsCounts$test4 <- PdSTRIPE@tssCountData[[22]]
tssRepsCounts$test5 <- PdSTRIPE@tssCountData[[23]]
tssRepsCounts$test6 <- PdSTRIPE@tssCountData[[24]]
tssRepsCounts$test7 <- PdSTRIPE@tssCountData[[25]]

#tsr.df <- experimentName@tsrDataMerged[[tsrSet]]
tss.df <- tss.start
 
tss.df <-
  tss.df[tss.df$nTAGs >= tagCountThreshold, ] #filtering the tss.start dataset

this.tssSet <- tssRepsCounts[[1]] #priming the loop to come
this.vec <- paste(this.tssSet$seq, this.tssSet$TSS, this.tssSet$strand, sep="_")

for (j in 2:length(tssRepsCounts)) {
  #need to add iterative merging
  this.tssSet <- tssRepsCounts[[j]]
  this.vec2 <- paste(this.tssSet$seq, this.tssSet$TSS, this.tssSet$strand, sep=" ")
  this.vec <- c(this.vec, this.vec2)
  }

tss.df <- colsplit(this.vec, " ", c("seq", "TSS", "strand"))

mergeCols <- c("seq", "TSS", "strand")

for (j in 1:length(tssRepsCounts)) { 
  print(j)
  this.tssSet <- tssRepsCounts[[j]]
  #...  we are discarding counts below the tag count threshold tagCountThreshold:
  this.tssSet <-
    this.tssSet[this.tssSet$nTAGs >= tagCountThreshold, ]
  
tss.df <-  right_join(tss.df, this.tssSet, by=mergeCols)

}

colnames(tss.df) <- c("seq", "TSS", "strand", "nTAGs_1", "isTRUE_1", "nTAGs_2", "isTRUE_2", "nTAGs_3", "isTRUE_3","nTAGs_4", "isTRUE_4","nTAGs_5", "isTRUE_5","nTAGs_6", "isTRUE_6", "nTAGs_7", "isTRUE_7")
  
