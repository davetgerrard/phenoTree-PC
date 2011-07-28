




poorReads <- c("H3K4me3.Skeletal_Muscle.62", "H3K4me3.Adult_Liver.5", "H3K27me3.Adult_Liver.3","H3K27me3.Adult_Liver.5")
poorReadIndex <- match(poorReads, names(allCounts))

counts.Filtered <- allCounts[,-(poorReadIndex)]
rm(allCounts,fullChromSet.pca) 
filteredChromSet.pca <- princomp(counts.Filtered[, -c(1:3)]) 

filteredChromSet.pca
summary(filteredChromSet.pca)

loadings(filteredChromSet.pca)
loadings(filteredChromSet.pca)[]

pdf("filteredChromSet.pca.1000bins.pdf")
#par(mfrow=c(3,2), mar=c(4,4,3,1))
plot(filteredChromSet.pca) 	# see relative sizes of PC variances.
#biplot(filteredChromSet.pca, choices=c(1,2))

plot(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,2], xlim=c(-1,1),  xlab="PC1", ylab="PC2")
text(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,2], 
	row.names(loadings(filteredChromSet.pca)[]), cex=0.7, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,3], xlim=c(-1,1), xlab="PC1", ylab="PC3")
text(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,3], 
	row.names(loadings(filteredChromSet.pca)[]), cex=0.7, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,4], xlim=c(-1,1), xlab="PC1", ylab="PC4")
text(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,4], 
	row.names(loadings(filteredChromSet.pca)[]), cex=0.7, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)


plot(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,5], xlim=c(-1,1), xlab="PC1", ylab="PC5")
text(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,5], 
	row.names(loadings(filteredChromSet.pca)[]), cex=0.7, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

plot(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,6], xlim=c(-1,1), xlab="PC1", ylab="PC6")
text(loadings(filteredChromSet.pca)[,1], loadings(filteredChromSet.pca)[,6], 
	row.names(loadings(filteredChromSet.pca)[]), cex=0.7, pos=4) ;abline(v=0,lty=2) ; abline(h=0,lty=2)

shortNames <- row.names(loadings(filteredChromSet.pca)[])
pchValue <- rep(15, length(shortNames))
pchValue[grep("H3K27",row.names(loadings(filteredChromSet.pca)[]))] <- 1
shortNames <- sub("H3K4me3.", "", shortNames)
shortNames <- sub("H3K27me3.", "", shortNames)

plot(loadings(filteredChromSet.pca)[,4], loadings(filteredChromSet.pca)[,5], xlim=c(-1,1), xlab="PC4", ylab="PC5", pch= pchValue)
mtext("A", 3, 1, adj = 0)
legend("bottomleft", c("H3K4me3", "H3K27me3"), pch=c(15,1))
text(loadings(filteredChromSet.pca)[,4], loadings(filteredChromSet.pca)[,5], 
	shortNames, cex=0.8, pos=4) 
abline(v=0,lty=2); abline(h=0,lty=2)

dev.off()


filteredChromSet.pca.6.scores <- data.frame(filteredChromSet.pca$scores[ , 1:6])
ucscCoordNames <- paste(counts.Filtered[,1],counts.Filtered[,2],counts.Filtered[,3])
row.names(filteredChromSet.pca.6.scores) <- ucscCoordNames


#write.table(filteredChromSet.pca.6.scores, file="filteredChromSet.pca.6.scores.whole.genome.tab", sep="\t", quote=F)

## look at the top and bottom 1000 bins in each PC to get profiles of scores.
pdf("boxPlotsOfFilteredPrincCompsWholeGenome.pdf",width=12, height=8)
par(mar=c(3,20,3,2),mfrow=c(2,1))
numbToUse <- 1000
colours.hMark <- ifelse(sampleInfo$hMark[match(names(counts.Filtered)[-c(1:3)], sampleInfo$niceName)] == "H3K4me3", "green", "red")
for(i in 1:6)  {
	thisComp <- paste("Comp",i,sep=".")
	titleText <- paste("Prin.",thisComp, numbToUse,"highest scoring bins")
	boxplot(counts.Filtered[tail(order(filteredChromSet.pca.6.scores[,thisComp]),n=numbToUse ),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)
	titleText <- paste("Prin.",thisComp, numbToUse,"lowest scoring bins")
	boxplot(counts.Filtered[head(order(filteredChromSet.pca.6.scores[,thisComp]),n=numbToUse ),-c(1:3)], horizontal=T, las=1, col=colours.hMark, main=titleText)


}
dev.off()


