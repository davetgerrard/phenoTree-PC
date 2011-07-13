this_chrom <- 'chr11'
binSize <- 1000
chr_length <- 135006516

options("scipen"= 15)
starts <- seq(from=0, to=chr_length, by=binSize )
ends <-  starts + binSize
ends[length(ends)] <- chr_length # shorten the last bin to end of chromosome
table <- cbind(this_chrom, starts, ends)

#table[,2:3] <- format(table[,2:3], scientific = FALSE)


write.table(table, file=paste(this_chrom,binSize,'bed',sep="."), quote=F,row.names=F, col.names=F,sep="\t") 
	
