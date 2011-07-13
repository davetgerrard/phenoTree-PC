#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################




setwd("C:/Users/dave/HalfStarted/PhenotypeTree/binCountsTest/")



filteredChromSet.pca.6.scores <- read.delim("filteredChromSet.p2ca.6.scores.whole.genome.tab", header=T)




# smooth muscle over skeletal muscle
tail(filteredChromSet.pca.6.scores[order(filteredChromSet.pca.6.scores$Comp.5),], n=300)

# skeletal muscle over smooth muscle
head(filteredChromSet.pca.6.scores[order(filteredChromSet.pca.6.scorefils$Comp.5),], n=300)


# pull out a specific region
 filteredChromSet.pca.6.scores[grep('chr11 17740000', row.names(filteredChromSet.pca.6.scores)),]
