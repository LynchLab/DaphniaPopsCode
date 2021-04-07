setwd("/scratch/rraborn/")
load("/data/LynchLabCME/Daphnia/DaphPopsTSRs/binaries/PdSTRIPE_complete.rdata")
library(TSRchitect)
library(TSRexploreR)
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

###### doing the same for OA85

tssRepsCounts <- NULL
tss.start <- PdSTRIPE@tssCountData[[10]]
tssRepsCounts$test1 <- PdSTRIPE@tssCountData[[10]]
tssRepsCounts$test2 <- PdSTRIPE@tssCountData[[11]]

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

colnames(tss.df) <- c("seq", "TSS", "strand", "nTAGs_1", "isTRUE_1", "nTAGs_2", "isTRUE_2")
tss.df.OA85 <- cbind(tss.df[,1:3], tss.df$nTAGs_1, tss.df$nTAGs_2) 
colnames(tss.df.OA85) <- c("seq","TSS","strand", "OA85_r1", "OA85_r2")
tss.df.OA85.new <- tss.df.OA85 %>% 
  arrange(seq,TSS, strand)



#### Now with LPB

tssRepsCounts <- NULL
tss.start <- PdSTRIPE@tssCountData[[4]]
tssRepsCounts$test1 <- PdSTRIPE@tssCountData[[4]]
tssRepsCounts$test2 <- PdSTRIPE@tssCountData[[4]]

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

colnames(tss.df) <- c("seq", "TSS", "strand", "nTAGs_1", "isTRUE_1", "nTAGs_2", "isTRUE_2")
tss.df.LPB <- cbind(tss.df[,1:3], tss.df$nTAGs_1, tss.df$nTAGs_2) 
colnames(tss.df.LPB) <- c("seq","TSS","strand", "LPB_r1", "LPB_r2")
tss.df.LPB.new <- tss.df.LPB %>% 
  arrange(seq,TSS, strand)

######

#replacing NAs with 0 for the counts file
tss.df.W17.new2 <- tss.df.W17.new %>% replace(is.na(.), 0) 
tss.df.TEX36.new2 <- tss.df.TEX36.new %>% replace(is.na(.), 0)
tss.df.OA85.new2 <- tss.df.OA85.new %>% replace(is.na(.), 0) 
tss.df.LPB.new2 <- tss.df.LPB.new %>% replace(is.na(.), 0) 

#####
tss.df.W17.new2.p <- tss.df.W17.new2[tss.df.W17.new2$strand=="+",]
tss.df.W17.new2.m <- tss.df.W17.new2[tss.df.W17.new2$strand=="-",]
tss.df.TEX36.new2.p <- tss.df.TEX36.new2[tss.df.TEX36.new2$strand=="+",]
tss.df.TEX36.new2.m <- tss.df.TEX36.new2[tss.df.TEX36.new2$strand=="-",]
tss.df.OA85.new2.p <- tss.df.OA85.new2[tss.df.TEX36.new2$strand=="+",]
tss.df.OA85.new2.m <- tss.df.OA85.new2[tss.df.TEX36.new2$strand=="-",]
tss.df.LPB.new2.p <- tss.df.LPB.new2[tss.df.TEX36.new2$strand=="+",]
tss.df.LPB.new2.m <- tss.df.LPB.new2[tss.df.TEX36.new2$strand=="-",]

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

#tsr.index <- which(PdSTRIPE@replicateIDs=="9") #grabbing the TEX36 TSR indices
#tsrs.TEX36 <- PdSTRIPE@tsrData[[tsr.index]] 

tsr.combined <- makeGRangesFromTSR(PdSTRIPE, "merged", 1)

TEX36_tss_tsr_ol <- findOverlaps(tss.TEX36.gr, tsr.combined)
TEX36_tss_tsr.ind <- as.data.frame(TEX36_tss_tsr_ol)
head(TEX36_tss_tsr.ind)

TEX36.OL.tss <- tss.df.TEX36.se.sorted[TEX36_tss_tsr.ind$queryHits,]
TEX36.OL.tsr <- as.data.frame(tsr.combined)
TEX36.OL.tsr <- TEX36.OL.tsr[TEX36_tss_tsr.ind$subjectHits,]
#rownames(TEX36.OL.tss) <- rownames(TEX36.OL.tsr) #adding tsr names to TSRs
#tex36_tsr_names <- paste("tsr",TEX36_tss_tsr.ind$subjectHits, sep="_") #adding new tsr names to TSRs
tex36_tsr_names <- as.numeric(TEX36_tss_tsr.ind$subjectHits)
TEX36.OL.tss <- TEX36.OL.tss[,-3:-4]
TEX36.OL.tss <- TEX36.OL.tss[,-1]
TEX36.tss.df <- cbind(tex36_tsr_names, TEX36.OL.tss)
colnames(TEX36.tss.df) <- c("ID", "start", "TEX36_r1", "TEX36_r2", "TEX36_r3")

#writing the TSSs to a file
write.table(TEX36.tss.df, file="TEX36sampleTSSs.txt", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

#writing the data to a text file
write.table(tss.df.W17.new2, file="W17_all_reps_TSS_counts.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(tss.df.TEX36.new2, file="TEX36_all_reps_TSS_counts.txt", sep = "\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

W17.sum.vec <- rowSums(tss.df.W17.new2[, c(4:10)])
TEX36.sum.vec <- rowSums(tss.df.TEX36.new2[,c(4:6)])
OA85.sum.vec <- rowSums(tss.df.OA85.new2[, c(4:5)])
LPB.sum.vec <- rowSums(tss.df.LPB.new2[, c(4:5)])
W17.tss.combined <- cbind(tss.df.W17.new2[,1:3], W17.sum.vec)
TEX36.tss.combined <- cbind(tss.df.TEX36.new2[,1:3], TEX36.sum.vec)
OA85.tss.combined <- cbind(tss.df.OA85.new2[,1:3], OA85.sum.vec)
LPB.tss.combined <- cbind(tss.df.LPB.new2[,1:3], LPB.sum.vec)


tss.out <-  inner_join(W17.tss.combined, TEX36.tss.combined, by=mergeCols)
tss.out2 <- inner_join(tss.out, OA85.tss.combined, by=mergeCols)
tss.out3 <- inner_join(tss.out2, LPB.tss.combined, by=mergeCols)
colnames(tss.out3) <- c("seq", "TSS", "strand", "W17_counts", "TEX36_counts", "OA85_counts", "LPB_counts")
write.table(tss.out3, file="DaphPopstssOutFile.csv", row.names=TRUE, col.names=TRUE, quote=FALSE)
