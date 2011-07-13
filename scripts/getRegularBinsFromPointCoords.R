#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################


#PRE-REQUISITE
#source("getTssFromCanonicalGenesTable.R")
### could use overlap select but might be faster now I have the data to re-create the bin co-ordinates from the TSS locations. 
# e.g. chr1  14409  14409 is in ch1	14000	15000
# need to then check for bins that appear more than once. 
head(simpleCanonTesBed.allProt)


#eliminate bins with duplicated TSS or TES
#TSS
tssInBins <- simpleCanonTssBed.allProt
tssInBins$bin.start <- floor(tssInBins$start/1000)*1000
tssInBins$bin.end <- ceiling(tssInBins$end/1000)*1000
tssInBins$binRef <- paste(tssInBins$chr,tssInBins$bin.start, tssInBins$bin.end)
head(tssInBins)
#table(table(tssInBins$binRef))
#singleTssBins <- names(table(tssInBins$binRef)==1)
binCounts <- table(tssInBins$binRef)
tssInBins$count <- as.numeric(binCounts[tssInBins$binRef])
nrow(tssInBins)
tssInBins <- subset(tssInBins, tssInBins$count==1)
nrow(tssInBins)
protCounts <- table(tssInBins$name)
tssInBins$protCount <- as.numeric(protCounts[tssInBins$name])
nrow(tssInBins)
tssInBins <- subset(tssInBins, tssInBins$protCount==1)
nrow(tssInBins)
rm(binCounts, protCounts)

#TES
tesInBins <- simpleCanonTesBed.allProt
tesInBins$bin.start <- floor(tesInBins$start/1000)*1000
tesInBins$bin.end <- ceiling(tesInBins$end/1000)*1000
tesInBins$binRef <- paste(tesInBins$chr,tesInBins$bin.start, tesInBins$bin.end)
head(tesInBins)
#table(table(tesInBins$binRef))
#singletesBins <- names(table(tesInBins$binRef)==1)
binCounts <- table(tesInBins$binRef)
tesInBins$count <- as.numeric(binCounts[tesInBins$binRef])
nrow(tesInBins)
tesInBins <- subset(tesInBins, tesInBins$count==1)
nrow(tesInBins)
protCounts <- table(tesInBins$name)
tesInBins$protCount <- as.numeric(protCounts[tesInBins$name])
nrow(tesInBins)
tesInBins <- subset(tesInBins, tesInBins$protCount==1)
nrow(tesInBins)
rm(binCounts, protCounts)

## these lists now only contain bins with a single TSS or TES. 
#remove non-standard chroms
tssInBins <- tssInBins[-grep("_", tssInBins$chr),]
tesInBins <- tesInBins[-grep("_", tesInBins$chr),]
if(output) {write.table(tssInBins, file="tssInBins.tab", row.names=F, quote=F, sep="\t")}
if(output) {write.table(tesInBins, file="tesInBins.tab", row.names=F, quote=F, sep="\t")}

#might want to remove all the 'canonical' objects
# rm(list=ls()[grep(ls(), pattern="Canon", ignore.case=T)])




