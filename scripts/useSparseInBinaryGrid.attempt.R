


fullProtList <- names(prot2go)
fullTermList <- unique(as.character(unlist(prot2go)))


list.1 <- fullProtList[1:1000]
list.2 <- fullTermList[1:1000]


binaryGrid.sparse <- Matrix(0,length(list.1),length(list.2), dimnames=list(c(list.1), c(list.2)), sparse=T)

str(binaryGrid.sparse)


setSparseValuesFromList <- function()  {
	
}

pairedList <- prot2go 

testVec <- prot2go["A0AV02"]
binaryGrid.sparse[names(testVec), as.character(unlist(testVec))] <- 1 


binaryGrid.sparse[head(list.1),]
sapply(head(fullProtList) ,  FUN= function(x) binaryGrid.sparse[x, as.character(unlist(prot2go[x]))] <- 1)

for(i in 1:length(list.1))  {
	binaryGrid.sparse[list.1[i], as.character(unlist(prot2go[list.1[i]]))] <- 1 

}
