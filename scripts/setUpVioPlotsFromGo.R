addGroupScoresToGroupTable <- function(groupTable,groupIds,groupName,scoreTable,idColumn,scoreColum) {
	index <- na.omit(match(groupIds , scoreTable[,idColumn]))
	scores <- scoreTable[index ,scoreColum]
	groupTable <- rbind(groupTable,data.frame(score=scores,group=groupName))
	
}



#pc.i <- 1

library(lattice)
library(gplots)

bwPlotsAsPdf <- function(pca.scores,orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15)  {
	pdfName <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,"pdf",sep=".")

	pdf(pdfName, paper="a4")
	frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results",sep=" ")
	textplot(frontPageText)
	for(pc.i in 1:length(summmaryPcResultList)) {
		plotPCsFromSummary(pca.scores=pca.scores, orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
	}
	dev.off()
}

bwPlotsAsFigs <- function(pca.scores, orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, figType="tiff")  {
	#frontPageText <- paste(pdfName,"\nUsing", orderBy.PC, "below", pcTableSigThreshold,"\nLimited to",numbGraphResults,"results",sep=" ")
	#textplot(frontPageText)
	#op <- par(cex=1.5)
	for(pc.i in 1:length(summmaryPcResultList)) {
		plotDir <- paste("vioPlotsByPC",orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep="_")
		if(!file.exists(plotDir))  {dir.create(plotDir)}
		fileName <- paste("vioPlotsByPC",pc.i,orderBy.PC,"Below",pcTableSigThreshold,"Max",numbGraphResults,figType,sep=".")
		fileName <- paste(plotDir,fileName,sep="/")
		tiff(fileName, compression="lzw",width=180, height=480,units="mm",res=300)	
		plotPCsFromSummary(pca.scores=pca.scores, orderBy.PC=orderBy.PC, pcTableSigThreshold=pcTableSigThreshold, numbGraphResults=numbGraphResults, pc.i=pc.i)
		dev.off()	
	}
	#par(op)
}


plotPCsFromSummary <- function(pca.scores, orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =15, pc.i=1)  {

	sigTable <- subset(summmaryPcResultList[[pc.i]],summmaryPcResultList[[pc.i]][,orderBy.PC] < pcTableSigThreshold)
	sigTable <- sigTable[order(sigTable[,orderBy.PC]),]
	sigTable$outputRank <- rank(sigTable[,orderBy.PC])

	sigClusters <- levels(as.factor(as.character(sigTable$bestCluster)))
	for(thisCluster in sigClusters)  {
		if(is.na(thisCluster)) {break}
		rankToSet <- min(na.omit(sigTable[sigTable$bestCluster == thisCluster,"outputRank"]))
		sigTable[which(sigTable$bestCluster == thisCluster),"outputRank"] <- rankToSet
	}	
	sigTable <- sigTable[order(sigTable$outputRank),]

	## need to count each cluster only once and each independent go term once. 
	
	plotTable <- sigTable[match(unique(sigTable$outputRank),sigTable$outputRank),]
	numbGraphResultsLocal <- min(nrow(plotTable),numbGraphResults)
	if(numbGraphResultsLocal < 1)  {
		warningMess <- paste("No sig results for PC", pc.i,"\n using", orderBy.PC, "below", pcTableSigThreshold)
		textplot(warningMess)
		return(NULL)
	}
	plotTable <- plotTable[1:numbGraphResultsLocal,]

	# could try strwrp or substr
	compHead <- paste("Comp.",pc.i,sep="")
	baseScore <- data.frame(score=pca.scores[,compHead],group="                            All Proteins")
	groupTable <- baseScore
	plotList <- list()
	for(i in 1:nrow(plotTable))  {
		if(is.na(plotTable$bestCluster[i]))  {
			plotList[[i]] <- intersect(listProtsInGo(goTerm=plotTable[i,"goTerm"],ontology=plotTable[i,"ontology"]),validProts)
			goID <- paste(paste(strwrap(plotTable[i,"description"],width=40),collapse="\n"),plotTable[i,"goTerm"],sep="\n")
		}	else {
			plotList[[i]] <- seedList[[plotTable$bestCluster[i]]]
			goID <- paste("CLUSTER", plotTable$bestCluster[i])
		}
		groupTable <- addGroupScoresToGroupTable(groupTable=groupTable,groupIds=plotList[[i]],
					groupName=goID,scoreTable=pca.scores,
					idColumn="spAccession",scoreColum=compHead)
	}


	bwTitle <- compHead
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	bwTitle <- compHead
	# complex index to reorder levels but not results of groupTable.
	bymedian <- with(groupTable, reorder(group, rev(as.numeric(row.names(groupTable))), min))	
	plot(bwplot(bymedian ~ score, groupTable, main=bwTitle,  aspect="xy", cex.main=2, cex.lab=2,
		par.settings = list(layout.widths = list(axis.left = 0, ylab.axis.padding = 40)),
	       panel = function(..., box.ratio) {
		   panel.violin(..., col = "darkgrey",
				varwidth = FALSE, box.ratio = box.ratio)
		   panel.bwplot(..., fill = "lightgrey", box.ratio = .1)
	       } )
	)
}




##################################


bwPlotsAsPdf(pca.scores=tss.pca.score.annot, orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10)
bwPlotsAsPdf(pca.scores=tss.pca.score.annot, orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10)
bwPlotsAsPdf(pca.scores=tss.pca.score.annot, orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10)

bwPlotsAsPdf(pca.scores=tss.pca.score.annot, orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-05, numbGraphResults =11)


bwPlotsAsFigs(pca.scores=tss.pca.score.annot, orderBy.PC="elimWilcox", pcTableSigThreshold=1.0e-04, numbGraphResults =10)
bwPlotsAsFigs(pca.scores=tss.pca.score.annot, orderBy.PC="elimAbsWilcox", pcTableSigThreshold=1.0e-03, numbGraphResults =10)
bwPlotsAsFigs(pca.scores=tss.pca.score.annot, orderBy.PC="elimKS", pcTableSigThreshold=1.0e-05, numbGraphResults =10)
