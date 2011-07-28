#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################


pc5.distinct.genes <- tss.pca.score.annot[(tss.pca.score.annot$Comp.4 > -200) & (tss.pca.score.annot$Comp.5 < -150),]$spAccession
write(pc5.distinct.genes, file="output/pc5.dinstinct.genes.txt")
pc5.distinct.genes.binRef <- tss.pca.score.annot[(tss.pca.score.annot$Comp.4 > -200) & (tss.pca.score.annot$Comp.5 < -150),]$binRef
write(pc5.distinct.genes.binRef, file="output/pc5.dinstinct.genes.binRef.txt")

pc4.5.topLeft.distinct.genes <- tss.pca.score.annot[(tss.pca.score.annot$Comp.4 < -150) & (tss.pca.score.annot$Comp.5 > 100),]$spAccession
write(pc4.5.topLeft.distinct.genes, file="output/pc4.5.topLeft.distinct.genes.txt")
pc4.5.topLeft.distinct.genes.binRef <- tss.pca.score.annot[(tss.pca.score.annot$Comp.4 < -150) & (tss.pca.score.annot$Comp.5 > 100),]$binRef
write(pc4.5.topLeft.distinct.genes.binRef, file="output/pc4.5.topLeft.distinct.genes.binRef.txt")


plot(tss.pca.score.annot$Comp.4, tss.pca.score.annot$Comp.5)
abline(h=-100)
abline(v=-150)
abline(h=100)
abline(h=0, lty=2) ;abline(v=0, lty=2)

listProtsInGo("GO:0051145", "BP")	# "smooth muscle cell differentiation"
focal_index <- match(listProtsInGo("GO:0051145", "BP"),tss.pca.score.annot$spAccession)
points(tss.pca.score.annot$Comp.4[focal_index],tss.pca.score.annot$Comp.5[focal_index], col="red", pch=22, bg="red")

listProtsInGo("GO:0030017", "CC")	# "sarcomere"
focal_index <- match(listProtsInGo("GO:0030017", "CC"),tss.pca.score.annot$spAccession)
points(tss.pca.score.annot$Comp.4[focal_index],tss.pca.score.annot$Comp.5[focal_index], col="blue", pch=23, bg="blue")

legend("topright", c("sarcomere", "smooth muscle diff."), col=c("blue", "red"), pch=c(23,22), pt.bg=c("blue","red"))


listProtsInGo("GO:0007608", "BP")	# olfactory
focal_index <- match(listProtsInGo("GO:0007608", "BP"),tss.pca.score.annot$spAccession)
tss.pca.score.annot[focal_index,]

# not useful.
#listProtsInGo("GO:0016568", "BP")	# chromatin remodelling
#focal_index <- match(listProtsInGo("GO:0016568", "BP"),tss.pca.score.annot$spAccession)
#tss.pca.score.annot[focal_index,]
#points(tss.pca.score.annot$Comp.4[focal_index],tss.pca.score.annot$Comp.5[focal_index], col="yellow", pch=24, bg="yellow")








