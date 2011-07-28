#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################
old.o <- options(scipen=999)  

#previous sequence


# getTssFromCanonicalGenesTable.R
# getRegularBinsFromPointCoords.R	# quick fix to re-create 1000bp bins around point coords.
# workWithGenomePcaScore.R	# load ALL the pca scores. 
# compareTssTesBins.R	
# 

# run Gene set enrichment on tss.pca.scores

setwd("C:/Users/dave/HalfStarted/PhenotypeTree/")
tss.pca.scores <- read.delim("output/tss.pca.score.tab", header=T)




tss.pca.score.annot  <- tss.pca.scores
tss.pca.score.annot$binRef  <- row.names(tss.pca.score.annot)
tss.pca.score.annot <- merge(tss.pca.score.annot,tssInBins, by="binRef")
head(tss.pca.score.annot)
#have a quick look at muscle related genes. 
tss.pca.score.annot[grep('MYO', tss.pca.score.annot$name),]

# split genename, uniprot ID and uniprot accession from 'name'
tss.pca.score.annot$spAccession <- matrix(unlist(strsplit(as.character(tss.pca.score.annot$name),";",fixed=T)),ncol=3,byrow=T)[,2]
tss.pca.score.annot$spId <- matrix(unlist(strsplit(as.character(tss.pca.score.annot$name),";",fixed=T)),ncol=3,byrow=T)[,3]
tss.pca.score.annot$spAccession <- sub('-[0-9]+','',tss.pca.score.annot$spAccession)	# remove protein isoform numbers

## hypothesise that PC4 should be enriched for muscle development genes
## pc5 should be enriched for 

plot(tss.pca.score.annot$Comp.4, tss.pca.score.annot$Comp.5)


# TASKS:-
# load topGo and protein to go map
# create 'geneList' of PC scores. Conc. on PCs 4 and 5. (2 & 3 should be brain though).
# run GO analysis on scores. Find enriched terms. 

#############


source("scripts/setUpProteinTopGo.R")


#### list proteins for a GO term, accounting for ontology
listProtsInGo <- function(goTerm,ontology)  {
	if(ontology == "BP") {
		genesInTerm(GOdata.BP,goTerm)[[1]]
	} else if(ontology == "MF") {
		genesInTerm(GOdata.MF,goTerm)[[1]]
	} else if(ontology == "CC") {
		genesInTerm(GOdata.CC,goTerm)[[1]]
	}
}

geneList <- tss.pca.score.annot$Comp.1
names(geneList) <- tss.pca.score.annot$spAccession
#nodeSize

GOdata.BP <- new("topGOdata",
  description =  "Histone mods data set",
              ontology = "BP" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenesReturnAll ,
              nodeSize = topGo.nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.MF <- new("topGOdata",
  description =  "Histone mods data set",
              ontology = "MF" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenesReturnAll ,
              nodeSize = topGo.nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )

GOdata.CC <- new("topGOdata",
  description =  "Histone mods data set",
              ontology = "CC" ,
              allGenes = geneList,
  geneSelectionFun = topDiffGenesReturnAll ,
              nodeSize = topGo.nodeSizeValue ,
              annot = annFUN.GO2genes,
GO2genes=go2prot 
               )



#multiple runs. cycle through PCs, ontology terms.
pcs <- 1:4


#ubi.pca.5.scores
#####testing params
#i<-4
#thisGOgraph <- "BP"
#thisGOgraph <- "MF"
thisGOgraph <- "CC"



runGoForAll <- function()  {

summmaryPcResultList <- list()	# gonna have a separate data frame for each PC.
for (i in 1:6)  {
	summmaryPcResultList[[i]] <- data.frame()
}


for(thisGOgraph in goGraphs )  {

#initialise a topGO object, takes a little while for 17,000 genes!
geneList <- tss.pca.score.annot$Comp.1
names(geneList) <- tss.pca.score.annot$spAccession
GOdata <- new("topGOdata",
		  description =  "Liver proteins data set",
              ontology = thisGOgraph ,
              allGenes = geneList,
		  geneSelectionFun = topDiffGenesReturnAll,
              nodeSize = nodeSizeValue ,
              annot = annFUN.GO2genes,
		GO2genes=go2prot 
		               )





for (i in 1:6)  {


compHead <- paste("Comp.",i,sep="")
geneList <- tss.pca.score.annot[,compHead]
names(geneList) <- tss.pca.score.annot$spAccession
GOdata <- updateGenes(GOdata,geneList,topDiffGenesReturnAll)



test.stat <- new("classicScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox tests") 
resultWilcox <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTest2Sided, name = "Wilcox test", cutOff = elimCutOff)    
resultElimWilcox <- getSigGroups(GOdata, test.stat)


test.stat <- new("classicScore", testStatistic = GOWilcoxTest1Sided , scoreOrder="increasing", name = "Wilcox tests") 
resultWilcoxIncreasing <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTest1Sided , scoreOrder="increasing", name = "Wilcox test", cutOff = elimCutOff)    #,alternative="greater"
resultElimWilcoxIncreasing <- getSigGroups(GOdata, test.stat)


test.stat <- new("classicScore", testStatistic = GOWilcoxTest1Sided, scoreOrder="decreasing", name = "Wilcox tests") 
resultWilcoxDecreasing <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTest1Sided, scoreOrder="decreasing", name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcoxDecreasing <- getSigGroups(GOdata, test.stat)


#resultElimWilcoxAdj <- resultElimWilcox 
#score(resultElimWilcoxAdj) <- qvalue(score(resultElimWilcox))$qvalue

test.stat <- new("classicScore", testStatistic = GOKSTest, name = "KS tests")     # ,alternative="less"  
resultKS <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOKSTest, name = "KS test", cutOff = elimCutOff)    #,alternative="less"
resultElimKS <- getSigGroups(GOdata, test.stat)

# convert scores to absolute values for absWilcox tests. 
geneList <- abs(tss.pca.score.annot[,compHead])
names(geneList) <- tss.pca.score.annot$spAccession
GOdata <- updateGenes(GOdata,geneList,topDiffGenesReturnAll)
test.stat <- new("classicScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox tests") 
resultWilcoxAbs <- getSigGroups(GOdata, test.stat)

test.stat <- new("elimScore", testStatistic = GOWilcoxTestGreater, name = "Wilcox test", cutOff = elimCutOff)    #,alternative="less"
resultElimWilcoxAbs <- getSigGroups(GOdata, test.stat)


manual.Wilcox <- data.frame(Wilcox=score(resultWilcox),goTerm=names(score(resultWilcox)))
manual.elimWilcox <- data.frame(elimWilcox=score(resultElimWilcox),goTerm=names(score(resultElimWilcox)))
manual.absWilcox <- data.frame(absWilcox =score(resultWilcoxAbs),goTerm=names(score(resultWilcoxAbs)))
manual.elimAbsWilcox <- data.frame(elimAbsWilcox=score(resultElimWilcoxAbs),goTerm=names(score(resultElimWilcoxAbs)))
manual.WilcoxIncreasing <- data.frame(WilcoxIncreasing=score(resultWilcoxIncreasing),goTerm=names(score(resultWilcoxIncreasing)))
manual.elimWilcoxIncreasing <- data.frame(elimWilcoxIncreasing =score(resultElimWilcoxIncreasing),goTerm=names(score(resultElimWilcoxIncreasing)))
manual.WilcoxDecreasing <- data.frame(WilcoxDecreasing=score(resultWilcoxDecreasing),goTerm=names(score(resultWilcoxDecreasing)))
manual.elimWilcoxDecreasing <- data.frame(elimWilcoxDecreasing=score(resultElimWilcoxDecreasing),goTerm=names(score(resultElimWilcoxDecreasing)))
#resultWilcoxIncreasing
#resultElimWilcoxIncreasing
#resultWilcoxDecreasing
#resultElimWilcoxDecreasing
manual.KS <- data.frame(KS=score(resultKS),goTerm=names(score(resultKS)))
manual.elimKS <- data.frame(elimKS=score(resultElimKS),goTerm=names(score(resultElimKS)))
manual.termStats <- termStat(GOdata,names(score(resultWilcox)))
manual.termStats <- subset(manual.termStats, select=-c(Significant,Expected))
manual.termStats$goTerm <- row.names(manual.termStats)

thisGoResults <- NULL
thisGoResults <- merge(manual.termStats, manual.Wilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.absWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimAbsWilcox, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.WilcoxIncreasing, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimWilcoxIncreasing, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.WilcoxDecreasing, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimWilcoxDecreasing, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.KS, by="goTerm")
thisGoResults <- merge(thisGoResults ,manual.elimKS, by="goTerm")

thisGoResults$ontology <- thisGOgraph
thisGoResults$description <-  as.character(topGO:::.getTermsDefinition(as.character(thisGoResults$goTerm), ontology(GOdata),numChar=200))
goGroup <- as.character(thisGoResults$goTerm)
#thisGoResults$number <- unlist(lapply(goGroup, FUN=function(x) length(genesInTerm(GOdata,x)[[1]])))

summmaryPcResultList[[i]] <- rbind(summmaryPcResultList[[i]],thisGoResults)

}	#end of PC
}	#end of GO tree

return(summmaryPcResultList)
}	# end of runGoForAll

############

#summmaryPcResultList <- list()	# gonna have a separate data frame for each PC.
#for (i in 1:6)  {#
#	summmaryPcResultList[[i]] <- data.frame()
#}






### save the results.

for(i in 1:length(summmaryPcResultList)) {
	fileName <- paste("GoSummaryByPc",i,"tab",sep=".")
	write.table(summmaryPcResultList[[i]], file=fileName, sep="\t",row.names=F,quote=F)
}



#"GO:0007420"	# brain development. Highly sig for PC3 in elimWilcox, elimAbsWilcox but not elimKS





######################################### FUNC CLUSTERING

#### need to re-do functional clustering. 

protTerms <- na.omit(prot2go[c(tss.pca.score.annot$spAccession)])







