#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

old.o <- options(scipen=999)  

orderedChroms <- c(paste("chr", 1:22, sep=""),"chrX", "chrY", "chrM")

sampleInfo <- read.delim("sampleInfo.tab", header=T)
sampleInfo$validName <-make.names(sampleInfo$fileName)  # required to deal with hyphens etc.
sampleInfo$niceName <- paste(sampleInfo$hMark, sampleInfo$Tissue, sampleInfo$individual, sep=".")



allCounts <- read.delim("genomeCountsAllMarks.tab", header=T)
names(allCounts)[-c(1:3)] <- sampleInfo$niceName[match(names(allCounts)[-c(1:3)], sampleInfo$validName)]


# was thinking of some crazy reordering index to get chroms in 'nice' order.
#temp <- allCounts$chr
#temp <- sub( "chr", "", temp)
#temp <- sub("X", "23",temp)
#temp <- sub("Y", "24",temp)
#temp <- as.numeric(sub("M", "25",temp))
## MUCH FASTER to order the summary tables using e.g. chromProps[match( orderedChroms,row.names(chromProps) ),]


##plot some summaries per sample, per chrom. Sex chrom ratios? 

chromCounts <- rowsum(allCounts[,-c(1:3)],allCounts$chr,reorder=FALSE)
chromCounts <- chromCounts[match( orderedChroms,row.names(chromCounts) ),]
heatmap(as.matrix(rowsum(allCounts[,-c(1:3)],allCounts$chr,reorder=FALSE)), mar=c(10,10))

# is their a gender effect in the reads from chrX and chrY? Yes. 
barplot(as.numeric(chromCounts["chrY",]/chromCounts["chrX",]), names.arg=sampleInfo$sex[match(names(chromCounts),sampleInfo$niceName)], horiz=T, las=1)
barplot(as.numeric(chromCounts))

sample.readTotals <- colSums(chromCounts)
chromProps <- t(t(chromCounts)/sample.readTotals)
#chromProps <- chromProps[match( orderedChroms,row.names(chromProps) ),]	#not needed if chromCounts already ordered.

###### NICE PLOT OF CHROM DATA
par(mar=c(3,20,2,2), mfrow=c(2,1))
barplot(sample.readTotals, horiz=T, las=1, main="Total Reads")
barplot(as.matrix(chromProps), horiz=T, las=1, main="Reads per chromosome as proportion")
text(cumsum(chromProps[,1])[1:5],(ncol(chromProps) + ((ncol(chromProps)-1)*0.2)),labels=names(chromProps[1:5,1]),
	pos=2, col="dark grey")


# complex index to reorder levels but not results of groupTable.
#	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))



## Run pca on full table
fullChromSet.pca <- princomp(allCounts[, -c(1:3)]) 

fullChromSet.pca
summary(fullChromSet.pca)

loadings(fullChromSet.pca)
loadings(fullChromSet.pca)[]

pdf("fullChromSet.pca.1000bins.pdf")
par(mfrow=c(2,3))
plot(fullChromSet.pca) 	# see relative sizes of PC variances.
#biplot(fullChromSet.pca, choices=c(1,2))

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,2], xlim=c(-1,1),  xlab="PC1", ylab="PC2")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,2], 
	row.names(loadings(fullChromSet.pca)[]), cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,3], xlim=c(-1,1), xlab="PC1", ylab="PC3")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,3], 
	row.names(loadings(fullChromSet.pca)[]), cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,4], xlim=c(-1,1), xlab="PC1", ylab="PC4")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,4], 
	row.names(loadings(fullChromSet.pca)[]), cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)


plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,5], xlim=c(-1,1), xlab="PC1", ylab="PC5")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,5], 
	row.names(loadings(fullChromSet.pca)[]), cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,6], xlim=c(-1,1), xlab="PC1", ylab="PC6")
text(loadings(fullChromSet.pca)[,1], loadings(fullChromSet.pca)[,6], 
	row.names(loadings(fullChromSet.pca)[]), cex=1, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)
dev.off()


fullChromSet.pca.6.scores <- data.frame(fullChromSet.pca$scores[ , 1:6])
ucscCoordNames <- paste(allCounts[,1],allCounts[,2],allCounts[,3])
row.names(fullChromSet.pca.6.scores) <- ucscCoordNames

#write.table(fullChromSet.pca.6.scores, file="fullChromSet.pca.6.scores.whole.genome.tab", sep="\t", quote=F)




## look at the top and bottom 1000 bins in each PC to get profiles of scores.
pdf("boxPlotsOfPrincCompsWholeGenome.pdf",width=12, height=8)
par(mar=c(3,20,3,2),mfrow=c(2,1))
numbToUse <- 1000
colours.hMark <- ifelse(sampleInfo$hMark[match(names(allCounts)[-c(1:3)], sampleInfo$niceName)] == "H3K4me3", "green", "red")
for(i in 1:6)  {
	thisComp <- paste("Comp",i,sep=".")
	titleText <- paste("Prin.",thisComp, numbToUse,"highest scoring bins")
	boxplot(allCounts[tail(order(fullChromSet.pca.6.scores[,thisComp]),n=numbToUse ),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)
	titleText <- paste("Prin.",thisComp, numbToUse,"lowest scoring bins")
	boxplot(allCounts[head(order(fullChromSet.pca.6.scores[,thisComp]),n=numbToUse ),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)


}
dev.off()









