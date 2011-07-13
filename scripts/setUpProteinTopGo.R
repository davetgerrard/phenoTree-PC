#!/usr/bin/Rscript

#######################################################
#
#       Dave Gerrard, University of Manchester
#       2011
#
######################################################

#######################LIBRARIES
library(topGO)
#require(gplots)
#library(qvalue)



######################## FUNCTIONS

GOWilcoxTestGreater <- function (object,alternativeType="greater") 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}



GOWilcoxTest2Sided <- function (object,alternativeType="two.sided") 
{
    N <- numAllMembers(object)
    na <- numMembers(object)
    if (na == 0 || na == N) 
        return(1)
    x.a <- rankMembers(object)
    return(wilcox.test(x.a, seq_len(N)[-x.a], alternative = alternativeType)$p.value)
}

topDiffGenesReturnAll <- function(allScore) {
  return(allScore )
}

######################### MAIN

prot2go <-  readMappings("C:/Users/dave/LiverProteins/data/go2prot.map")
go2prot <- inverseList(prot2go)

goGraphs <- c("BP","MF","CC")
topGo.nodeSizeValue <- 10
topGo.topNodesValue <- 20
topGo.elimCutOff <- 0.01   # 0.001 is too harsh
topGo.topTerms <- 50







