#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################


library(adephylo)


setwd("C:/Users/dave/HalfStarted/PhenotypeTree/")

old.o <- options(scipen=999)           # need to disable scientific notation for large genome coordinates. Restore options at end of script.



## load data as matrix
allCounts <- read.delim("binCountsTest/genomeCountsAllMarks.tab", header=T)
sampleInfo <- read.delim("binCountsTest/sampleInfo.tab", header=T)
sampleInfo$validName <-make.names(sampleInfo$fileName)  # required to deal with hyphens etc.
sampleInfo$niceName <- paste(sampleInfo$hMark, sampleInfo$Tissue, sampleInfo$individual, sep=".")
names(allCounts)[-c(1:3)] <- sampleInfo$niceName[match(names(allCounts)[-c(1:3)], sampleInfo$validName)]


poorReads <- c("H3K4me3.Skeletal_Muscle.62", "H3K4me3.Adult_Liver.5", "H3K27me3.Adult_Liver.3","H3K27me3.Adult_Liver.5")
poorReadIndex <- match(poorReads, names(allCounts))
counts.Filtered <- allCounts[,-(poorReadIndex)]
ucscCoordNames <- paste(counts.Filtered[,1],counts.Filtered[,2],counts.Filtered[,3])
rm(allCounts)
row.names(counts.Filtered) <- make.names(ucscCoordNames, unique=T)
rm(ucscCoordNames)

# TEMP filter on tss containing bins for speed
tssInBins <- read.delim("output/tssInBins.tab", header=T)
tss.counts <- counts.Filtered[match(make.names(tssInBins$binRef), row.names(counts.Filtered)),]

tss.counts <- tss.counts[,-c(1:3)]


#countAcross <- rowSums(counts.Filtered[,-c(1:3)] > 0 )
#minTotalCounts <- 1000
#counts.Filtered.H3K4me3 <- counts.Filtered[countAcross > minTotalCounts ,grep('H3K4me3',names(counts.Filtered))]
#counts.Filtered.H3K4me3.sample <- counts.Filtered.H3K4me3[sample(1:nrow(counts.Filtered.H3K4me3),250),]
#countMatrix <- t(counts.Filtered.H3K4me3.sample)

sumAcross <- rowSums(tss.counts   )
minSum <- 1000
tss.counts.H3K4me3 <- tss.counts[ sumAcross > minSum,grep('H3K4me3',names(tss.counts))]
countMatrix <- t(tss.counts.H3K4me3)


## load tree
tissueTreeText <- "((((H3K4me3.Skeletal_Muscle.19,H3K4me3.Skeletal_Muscle.63)SKEL_MUSCLE,(H3K4me3.Stomach_Smooth_Muscle.28,H3K4me3.Rectal_Smooth_Muscle.30)SMOO_MUSCLE)MUSCLE,H3K4me3.Fetal_Kidney.UW_H-22676)MESODERM,(H3K4me3.Pancreatic_Islets.NA, H3K4me3.Adult_Liver.3)ENDODERM, H3K4me3.Fetal_Brain.HuFNSC-T)Root;"
tissueTree <- read.tree(text =tissueTreeText)
plot(tissueTree)
W <- proxTips(tissueTree, met = "Abouheif")


# create phylo4 object

tissue.H3K4me3.phylo4d <- phylo4d(tissueTree, countMatrix)
#table.phylo4d(tissue.H3K4me3.phylo4d ) # too many traits


# run pPCA

tissue.ppca <- ppca(tissue.H3K4me3.phylo4d,scale=FALSE,scannf=FALSE,nfposi=3,nfnega=3, method="Abouheif")

barplot(tissue.ppca$eig,main='pPCA eigenvalues',cex.main=1.8)

plot(tissue.ppca,ratio.tree=.7)
dotchart(tissue.ppca$c1[,1],lab=rownames(tissue.ppca$c1),main="Global principal component 1")


#screeplot(tissue.ppca)

scatter(tissue.ppca, useLag=TRUE)
#plot(tissue.ppca, useLag=TRUE)
#table.phylo4d(tissue.ppca, ratio=.7, var.lab=c("1st global PC", "1st local PC"), tip.label=myLab, box=FALSE,cex.lab=1.4, cex.sym=1.2, show.node.label=TRUE)


### select traits by PC
pc.index <- 3
a <- tissue.ppca$c1[,pc.index] # loadings 
names(a) <- row.names(tissue.ppca$c1)
highContrib <- a[a< quantile(a,0.01) | a>quantile(a,0.99)]
datSel <- cbind.data.frame(countMatrix[, names(highContrib)], tissue.ppca$li)
temp <- phylo4d(tissueTree, datSel)
table.phylo4d(temp) # plot of most structured traits
head(highContrib[(order(highContrib))])	# show the high contributing traits in order.
tail(highContrib[(order(highContrib))])


## PHYLOGENETIC AUTOCORRELATION TESTS FOR THESE TRAITS
prox <- proxTips(LargeTissueTree, method="Abouheif")
abouheif.moran(luscombe.tree.data[, names(highContrib)], prox)



