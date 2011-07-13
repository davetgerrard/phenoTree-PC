
old.o <- options(scipen=999)           # need to disable scientific notation for large genome coordinates. Restore options at end of script.

### test with chr11 data only 
### H3K4me3, H3K27me3. 8 tissue samples (some duplicates)
### pdx1 (pancreatic fate determination) is on chr13 (NOT chr11!)


setwd("C:/Users/dave/HalfStarted/PhenotypeTree/binCountsTest/")

countFiles <- dir("chr11", full.names=T)
#chrom specific files have names like:-
# GSM669625_UCSF-UBC.Fetal_Brain.H3K27me3.HuFNSC-T.chr11.1000.binCounts
# but these are not perfectly structured. Cannot use strsplit()
#sampleInfo <- data.frame(rawName=basename(countFiles))
#sampleInfo$h_mark <- NA
#write.table(basename(countFiles),file="sampleNames.txt", quote=F, row.names=F, col.names=F)

# I have made a file with sampleInfo for this test
sampleInfo <- read.delim("sampleInfo.tab", header=T)
sampleInfo$niceName <- paste(sampleInfo$hMark, sampleInfo$Tissue, sampleInfo$individual, sep=".")




#thisCount <- read.delim(countFiles[1], header=F, as.is=T)
#thisCount.2  <- read.delim(countFiles[2], header=F)


allCounts <- data.frame()
for (i in 1:length(countFiles))  {
	thisCount <- read.delim(countFiles[i], header=F)
	if(i==1) {
		prototype <- thisCount[,1:3]
		allCounts <- prototype
		names(allCounts) <- c("chr", "start", "end")
	} else {
		# chr start and end of each file should be identical 
		stopifnot(identical(thisCount[,1:3], prototype[,1:3]))
	}
	allCounts[,basename(countFiles[i])] <- thisCount[,4]
}
names(allCounts)[4:19] <- sampleInfo$niceName



## What happens if take out rows with all zeros?
# Not a lot. Removes about 5000 of 135,000 rows.
#allCounts$rowTotal <- rowSums(allCounts[, -(1:3)])
#hist(allCounts)
#hist(allCounts$rowTotal)
#allCounts.noZero <- subset(allCounts, rowTotal > 0)
#nrow(allCounts.noZero)
#nrow(allCounts)
#cor.test(allCounts[,6], allCounts[,5])
#cor.test(allCounts.noZero[,6], allCounts.noZero[,5])
#cor.test(allCounts.noZero[,4], allCounts.noZero[,5])

# The pancreatic islet samples have less reads in total. 300k (H3K4me3) and 600k (H3K27me3) vs 1-1.5 million for others.
colSums(allCounts[,4:19])



##### heatmap with euclidean distance
numbRandBins <- 500
randBinIndex <- sample(1:nrow(allCounts),numbRandBins )
titleText <- paste(numbRandBins, "random bins. Euclidean distance")
heatmap(t(as.matrix(allCounts[randBinIndex , 4:19])), main=titleText )
## splits H3K4me3 and H3K27me3 with exception of Pancreatic islets H3K4me3, 
## which sit as close outgroup to H3K27me3 group.

##### heatmap with Binary distance
numbRandBins <- 500
randBinIndex <- sample(1:nrow(allCounts),numbRandBins )
titleText <- paste(numbRandBins, "random bins. Binary distance")
heatmap(t(as.matrix(allCounts[randBinIndex , 4:19])), main=titleText, distfun=function(c) dist(c,method="binary") )
## No obvious grouping at all. 


##Try a pca
randomSet.pca <- princomp(allCounts[randBinIndex , 4:19])  # v.quick

fullChromSet.pca <- princomp(allCounts[, 4:19])  #still pretty quick

fullChromSet.pca
summary(fullChromSet.pca)

loadings(fullChromSet.pca)
loadings(fullChromSet.pca)[]

pdf("fullChromSet.pca.chr11.1000bins.pdf")
par(mfrow=c(2,3))
plot(fullChromSet.pca) 	# see relative sizes of PC variances.
#biplot(fullChromSet.pca, choices=c(1,2))

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,2], xlim=c(0,1), xlab="PC1", ylab="PC2")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,2], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,3], xlim=c(0,1), xlab="PC1", ylab="PC3")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,3], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,4], xlim=c(0,1), xlab="PC1", ylab="PC4")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,4], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)


plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,5], xlim=c(0,1), xlab="PC1", ylab="PC5")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,5], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,6], xlim=c(0,1), xlab="PC1", ylab="PC6")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,6], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
dev.off()

# 2 vs 3 
plot(loadings(fullChromSet.pca)[,2], loadings(fullChromSet.pca)[,3], xlim=c(-0.3,1), xlab="PC2", ylab="PC3")
text(loadings(fullChromSet.pca)[,2], loadings(fullChromSet.pca)[,3], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

# 2 vs 6 (both vary along H3K27me3)
plot(loadings(fullChromSet.pca)[,2], loadings(fullChromSet.pca)[,6], xlim=c(-0.3,1), xlab="PC2", ylab="PC6")
text(loadings(fullChromSet.pca)[,2], loadings(fullChromSet.pca)[,6], 
	sampleInfo$niceName, cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)


fullChromSet.pca.6.scores <- data.frame(fullChromSet.pca$scores[ , 1:6])
ucscCoordNames <- paste(allCounts[,1],allCounts[,2],allCounts[,3])
row.names(fullChromSet.pca.6.scores) <- ucscCoordNames

write.table(fullChromSet.pca.6.scores, file="fullChromSet.pca.6.scores.tab", sep="\t", quote=F)

#not sure if this is useful, but shows where the zero count bins fall on the pc scores. 
hist(fullChromSet.pca.6.scores$Comp.1, breaks=5000,xlim=c(-100,100))
hist(fullChromSet.pca.6.scores$Comp.2, breaks=5000,xlim=c(-100,100))
hist(fullChromSet.pca.6.scores$Comp.3, breaks=5000,xlim=c(-100,100))
hist(fullChromSet.pca.6.scores$Comp.4, breaks=5000,xlim=c(-100,100))
hist(fullChromSet.pca.6.scores$Comp.5, breaks=5000,xlim=c(-100,100))
hist(fullChromSet.pca.6.scores$Comp.6, breaks=5000,xlim=c(-100,100))

allCounts[head(order(fullChromSet.pca.6.scores$Comp.1)),]
allCounts[head(order(fullChromSet.pca.6.scores$Comp.2)),]
allCounts[head(order(fullChromSet.pca.6.scores$Comp.3)),]
allCounts[head(order(fullChromSet.pca.6.scores$Comp.4)),]
allCounts[head(order(fullChromSet.pca.6.scores$Comp.5)),]
allCounts[head(order(fullChromSet.pca.6.scores$Comp.6)),]

allCounts[tail(order(fullChromSet.pca.6.scores$Comp.1)),]
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.2)),]
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.3)),]
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.4)),]
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.5)),]
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.6)),]

## this came up high on PC2
allCounts[tail(order(fullChromSet.pca.6.scores$Comp.2)),]
fullChromSet.pca.6.scores[31830,]
# compare with score distributions
summary(fullChromSet.pca.6.scores)


## look at the top and bottom 1000 bins in each PC to get profiles of scores.
pdf("boxPlotsOfPrincCompsWholeGenome.pdf",width=12, height=8)
par(mar=c(3,20,3,2),mfrow=c(2,1))
numbToUse <- 1000
colours.hMark <- ifelse(sampleInfo$hMark == "H3K4me3", "green", "red")
for(i in 1:6)  {
	thisComp <- paste("Comp",i,sep=".")
	titleText <- paste("Prin.",thisComp, numbToUse,"highest scoring bins")
	boxplot(allCounts[tail(order(fullChromSet.pca.6.scores[,thisComp]),n=1000),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)
	titleText <- paste("Prin.",thisComp, numbToUse,"lowest scoring bins")
	boxplot(allCounts[head(order(fullChromSet.pca.6.scores[,thisComp]),n=1000),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)


}
dev.off()







##### TIDY UP


options(old.o) #restore options.
