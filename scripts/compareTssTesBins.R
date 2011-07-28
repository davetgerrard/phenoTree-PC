#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################
old.o <- options(scipen=999)  

setwd("C:/Users/dave/HalfStarted/PhenotypeTree/")

#previous sequence


# getTssFromCanonicalGenesTable.R
# getRegularBinsFromPointCoords.R	# quick fix to re-create 1000bp bins around point coords.
# workWithGenomePcaScore.R	# load ALL the pca scores. 
# compareTssTesBins.R	# THIS FILE


#source(getRegul...
tssInBins <- read.delim("output/tssInBins.tab", header=T)
tesInBins <- read.delim("output/tesInBins.tab", header=T)

# compare pca scores for tss and tes bins
# might be better to use raw counts scores for this?

#tssInBins,tesInBins

match(tssInBins$binRef, row.names(filteredChromSet.pca.6.scores))

tss.pca.scores <- filteredChromSet.pca.6.scores[match(tssInBins$binRef, row.names(filteredChromSet.pca.6.scores)),]
tes.pca.scores <- filteredChromSet.pca.6.scores[match(tesInBins$binRef, row.names(filteredChromSet.pca.6.scores)),]
random.pca.scores <- filteredChromSet.pca.6.scores[sample(1:nrow(filteredChromSet.pca.6.scores),20000), ]
rm(filteredChromSet.pca.6.scores)  # shouldn't need this anymore.

write.table(tss.pca.scores, file="output/tss.pca.score.tab", row.names=T, quote=F, sep="\t")
write.table(tes.pca.scores, file="output/tes.pca.score.tab", row.names=T, quote=F, sep="\t")
write.table(random.pca.scores, file="output/random.pca.score.tab", row.names=T, quote=F, sep="\t")

par(mfrow=c(2,3))
for ( i in 1:6)  {
	thisComp <- paste("Comp" , i, sep=".")
	#thisComp <- "Comp.6"
	y_limits <- c(	min(min(tss.pca.scores[,thisComp]),min(tes.pca.scores[,thisComp]),min(random.pca.scores[,thisComp])), 
				max(max(tss.pca.scores[,thisComp]),max(tes.pca.scores[,thisComp]),max(random.pca.scores[,thisComp])))
	label.names=c("tss","tes","random")
	boxplot(tss.pca.scores[,thisComp],xlim=c(0,4), ylim=y_limits, main=thisComp)
	boxplot(tes.pca.scores[,thisComp], add=T, at =2)
	boxplot(random.pca.scores[,thisComp], add=T, at =3, names="random")
	mtext(label.names,1,1, at=1:3)
}




