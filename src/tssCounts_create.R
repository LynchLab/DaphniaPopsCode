setwd("/scratch/rraborn/")
#load("/data/LynchLabCME/Daphnia/DaphPopsTSRs/binaries/PdSTRIPE_complete.rdata")
library(dplyr)
library(stringr)
library(reshape2)
library(GenomicRanges)

experimentName <- PdSTRIPE
tagCountThreshold <- 2

#starting with W17
tss.start <- PdSTRIPE@tssCountData[[19]]
tssRepsCounts$test2 <- PdSTRIPE@tssCountData[[19]]
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

#length(this.vec)
this.vec.new <- unique(this.vec)
#length(this.vec.new)

tss.df <- colsplit(this.vec.new, " ", c("seq", "TSS", "strand"))

mergeCols <- c("seq", "TSS", "strand")

for (j in 1:length(tssRepsCounts)) { 
  print(j)
  this.tssSet <- tssRepsCounts[[j]]
  #...  we are discarding counts below the tag count threshold tagCountThreshold:
#  this.tssSet <-
#    this.tssSet[this.tssSet$nTAGs >= tagCountThreshold, ]
  
tss.df <-  right_join(tss.df, this.tssSet, by=mergeCols)
}

colnames(tss.df) <- c("seq", "TSS", "strand", "nTAGs_1", "isTRUE_1", "nTAGs_2", "isTRUE_2", "nTAGs_3", "isTRUE_3","nTAGs_4", "isTRUE_4","nTAGs_5", "isTRUE_5","nTAGs_6", "isTRUE_6", "nTAGs_7", "isTRUE_7")
tss.df.W17 <- cbind(tss.df[,1:3], tss.df$nTAGs_1, tss.df$nTAGs_2, tss.df$nTAGs_3, tss.df$nTAGs_4, tss.df$nTAGs_5, tss.df$nTAGs_6, tss.df$nTAGs_7)
colnames(tss.df.W17) <- c("seq", "TSS", "strand", "W17_r1", "W17_r2", "W17_r3", "W17_r4", "W17_r5", "W17_r6", "W17_r7")
tss.df.W17.new <- tss.df.W17 %>% 
  arrange(seq,TSS, strand)

####### doing the same for TEX36

tssRepsCounts <- NULL
tss.start <- PdSTRIPE@tssCountData[[16]]
tssRepsCounts$test1 <- PdSTRIPE@tssCountData[[16]]
tssRepsCounts$test2 <- PdSTRIPE@tssCountData[[17]]
tssRepsCounts$test3 <- PdSTRIPE@tssCountData[[18]]

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

length(this.vec)
this.vec.new <- unique(this.vec)
length(this.vec.new)

tss.df <- colsplit(this.vec.new, " ", c("seq", "TSS", "strand"))

mergeCols <- c("seq", "TSS", "strand")

for (j in 1:length(tssRepsCounts)) { 
  print(j)
  this.tssSet <- tssRepsCounts[[j]]
  #...  we are discarding counts below the tag count threshold tagCountThreshold:
  #  this.tssSet <-
  #    this.tssSet[this.tssSet$nTAGs >= tagCountThreshold, ]
  
  tss.df <-  right_join(tss.df, this.tssSet, by=mergeCols)
}

colnames(tss.df) <- c("seq", "TSS", "strand", "nTAGs_1", "isTRUE_1", "nTAGs_2", "isTRUE_2", "nTAGs_3", "isTRUE_3")
tss.df.TEX36 <- cbind(tss.df[,1:3], tss.df$nTAGs_1, tss.df$nTAGs_2, tss.df$nTAGs_3) 
colnames(tss.df.TEX36) <- c("seq","TSS","strand", "TEX36_r1", "TEX36_r2", "TEX36_r3")
tss.df.TEX36.new <- tss.df.TEX36 %>% 
  arrange(seq,TSS, strand)

#replacing NAs with 0 for the counts file
tss.df.W17.new2 <- tss.df.W17.new %>% replace(is.na(.), 0) 
tss.df.TEX36.new2 <- tss.df.TEX36.new %>% replace(is.na(.), 0)

#####
tss.df.W17.new2.p <- tss.df.W17.new2[tss.df.W17.new2$strand=="+",]
tss.df.W17.new2.m <- tss.df.W17.new2[tss.df.W17.new2$strand=="-",]
tss.df.TEX36.new2.p <- tss.df.TEX36.new2[tss.df.TEX36.new2$strand=="+",]
tss.df.TEX36.new2.m <- tss.df.TEX36.new2[tss.df.TEX36.new2$strand=="-",]

tss.df.TEX36.new2.p.start <- tss.df.TEX36.new2.p$TSS -1
tss.df.TEX36.new2.p.end <- tss.df.TEX36.new2.p$TSS

tss.df.TEX36.new2.m.start <- tss.df.TEX36.new2.m$TSS+1
tss.df.TEX36.new2.m.end <- tss.df.TEX36.new2.m$TSS

tss.df.TEX36.p.se <- cbind(tss.df.TEX36.new2.p$seq, tss.df.TEX36.new2.p.start, tss.df.TEX36.new2.p.end, tss.df.TEX36.new2.p$strand, tss.df.TEX36.new2.p[,4:6])
tss.df.TEX36.m.se <- cbind(tss.df.TEX36.new2.m$seq, tss.df.TEX36.new2.m.start, tss.df.TEX36.new2.m.end, tss.df.TEX36.new2.m$strand, tss.df.TEX36.new2.m[,4:6])
colnames(tss.df.TEX36.p.se) <- c("seq","start","end","strand", "TEX36_r1", "TEX36_r2", "TEX36_r3")
colnames(tss.df.TEX36.m.se) <- c("seq","start","end", "strand", "TEX36_r1", "TEX36_r2", "TEX36_r3")

tss.df.TEX36.se <- rbind(tss.df.TEX36.p.se, tss.df.TEX36.m.se)
tss.df.TEX36.se.sorted <- tss.df.TEX36.se %>% 
  arrange(seq, start, strand)

tss.TEX36.gr <- makeGRangesFromDataFrame(tss.df.TEX36.se.sorted,
                         keep.extra.columns=TRUE,
                         ignore.strand=FALSE,
                         seqinfo=NULL,
                         seqnames.field=c("seq"),
                         start.field="start",
                         end.field=c("end", "stop"),
                         strand.field="strand",
                         starts.in.df.are.0based=FALSE)


tsr.TEX36 <- PdSTRIPE@tsrData[[9]] #in progress 

#writing the data to a text file
write.table(tss.df.W17.new2, file="W17_all_reps_TSS_counts.txt", sep="\t")
write.table(tss.df.TEX36.new2, file="TEX36_all_reps_TSS_counts.txt", sep = "\t")

