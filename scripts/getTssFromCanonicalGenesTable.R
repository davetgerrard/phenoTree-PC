#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################
old.o <- options(scipen=999)  

output <- FALSE

canon.genes <- read.delim("C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_geneStrandStartEndUniprotSymbolCanonicalLinked.Edited.tab", header=T)
nrow(canon.genes)


canon.genes <- na.omit(canon.genes)

nrow(canon.genes)
head(canon.genes)

# store single points as 1-based so adjust bed start by +1.  
canon.genes$hg19.canonical.tss <- ifelse(canon.genes$hg19.knownGene.strand == '+', canon.genes$hg19.knownGene.txStart + 1, canon.genes$hg19.knownGene.txEnd)
canon.genes$hg19.canonical.tes <- ifelse(canon.genes$hg19.knownGene.strand == '-', canon.genes$hg19.knownGene.txStart + 1, canon.genes$hg19.knownGene.txEnd)

if(output) {write.table(canon.genes, file="C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_canonicalGenesTssUniprot.tab", row.names=F, quote=F)}

#simpleCanonTssBed <- cbind(canon.genes$hg19.knownGene.chrom, canon.genes$hg19.canonical.tss ,canon.genes$hg19.canonical.tss )
simpleCanonTssBed <- data.frame(	chr=canon.genes$hg19.knownGene.chrom, 
						start=(canon.genes$hg19.canonical.tss -1),
						end=canon.genes$hg19.canonical.tss,
						name=paste(canon.genes$hg19.kgXref.geneSymbol,canon.genes$hg19.kgXref.spID,canon.genes$hg19.kgXref.spDisplayID,sep=";") )
head(simpleCanonTssBed, n=100)
if(output) write.table(simpleCanonTssBed, file="C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_canonicalTssOnly.tab", row.names=F, quote=F)


# now wnat only genes with a swissprot id
canon.genes.allProt <- subset(canon.genes, hg19.kgXref.spDisplayID != "")
if(output) {write.table(canon.genes.allProt, file="C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_canonicalGenesTssUniprotProtOnly.tab", row.names=F, quote=F)}
simpleCanonTssBed.allProt <- data.frame(	chr=canon.genes.allProt$hg19.knownGene.chrom, 
						start=(canon.genes.allProt$hg19.canonical.tss - 1) ,
						end=canon.genes.allProt$hg19.canonical.tss,
						name=paste(canon.genes.allProt$hg19.kgXref.geneSymbol,canon.genes.allProt$hg19.kgXref.spID,canon.genes.allProt$hg19.kgXref.spDisplayID,sep=";") )
if(output) {write.table(simpleCanonTssBed.allProt, file="C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_canonicalTssOnlyWithProtOnly.bed", row.names=F, quote=F)}

# for comparison of enrichment, I also want an equivelent table of TESs (transcription End sites)
#canon.genes.allProt$hg19.canonical.tes <- ifelse(canon.genes.allProt$hg19.knownGene.strand == '-', canon.genes.allProt$hg19.knownGene.txStart, canon.genes.allProt$hg19.knownGene.txEnd)
simpleCanonTesBed.allProt <- data.frame(	chr=canon.genes.allProt$hg19.knownGene.chrom, 
						start=(canon.genes.allProt$hg19.canonical.tes - 1) ,
						end=canon.genes.allProt$hg19.canonical.tes,
						name=paste(canon.genes.allProt$hg19.kgXref.geneSymbol,canon.genes.allProt$hg19.kgXref.spID,canon.genes.allProt$hg19.kgXref.spDisplayID,sep=";") )
if(output) {write.table(simpleCanonTesBed.allProt, file="C:/Users/dave/HalfStarted/PhenotypeTree/data/hg19_canonicalTesOnlyWithProtOnly.bed", row.names=F, quote=F)}








