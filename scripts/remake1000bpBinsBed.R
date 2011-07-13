#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################




## needed to quickly reconstruct the 1000bp bins bed file.
bins1000 <- row.names(filteredChromSet.pca.6.scores)
write(bins1000, file="bins1000.bed")
bins1000 <- read.delim("bins1000.bed", sep=' ', header=F)
write.table(bins1000, file="bins1000.bed", sep="\t",  row.names=F,quote=F, col.names=F)

## could easily have used 'cut' on a binCount file on Templar. 
# do I need names for each bin?