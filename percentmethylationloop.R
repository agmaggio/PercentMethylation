#load your libraries 
library(dplyr)

#assign variables that aren't specific to your sample name ($samples may need to be changed)
colnames <- c("chromosome", "position", "strand", "countmethylated", "countunmethylated","ccontext","trinuccontext")
filenames <- list.files(path=".", pattern="_R_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt.gz", full.names=TRUE)
samples <- c("AM381","AM105","AM183","AM249","AM331","AM340","AM382","AM421","AM453","AM524","AM554","AM569","AM631","AM632","AM732","AM745","AM135","AM171","AM229","AM285","AM288","AM407","AM613","AM666","AM850","AM858","AM16","AM84","AM140","AM273","AM279","AM293","AM338","AM397","AM403","AM523","AM558","AM693","AM705","AM708","AM790","AM875","AM133","AM159","AM269","AM320","AM341","AM385","AM488","AM553","AM578","AM676")
results <- data.frame(matrix(nrow = 3, ncol = 2))
colnames(results) <- c("ID","AVGPM") 
DMR <- "ESRRG" 
chrom <- "chr1" 
start <- 216721310
stop <- 216721816 

#assign functions you will use
percent.methylation <- function(x,y) {
  x/(x+y)
}

#start your loop (check that you dont need to change variable names)
for (i in samples) { 
  infile<- paste(i,"_R_1_val_1_bismark_bt2_pe.deduplicated.bismark.cov.gz.CpG_report.txt.gz",sep="")
  outfile<- paste(i,DMR,".csv",sep="")
  a <- paste(i)
  b <-paste(i,chrom,sep="")
  c <-paste(i,DMR,sep="")
  
  #read cov file
  message("Reading file...",i)
  a <- read.table(infile,header=FALSE, col.names=colnames)
  
  #subset by chrom 
  message("Subsetting by chromosome...")
  b <- subset(a, chromosome==chrom)
  
  #subset by position on chrom
  message("Subsetting by position on chromosome...")
  c <- subset(b, position >= start & position <= stop)
  dim(c)
  
  #add percent methylation and cov columns for each row
  message("Add %methylation and cov columns...")
  d <- mutate(c, PM=percent.methylation(c[,4],c[,5]), cov = (c[,4]+c[,5]))
	
  #filter by coverage 
  e <- subset(d, cov >= 3) 

  #mean methylation plus and minus strand 
  f <- mean(e$PM,na.rm=TRUE)
  results[i,1] <- i
  results[i,2] <- f
  
  #write output csv 
  message("Writing output csv for sample\ ",i)
  write.csv(e, file= outfile)

}
  #saveRData 
  save.image(file="avgpm.RData")

  #write out averages for each file 
  avgfilename <- paste(DMR,"_avgpm.csv",sep="")
  write.csv(results, file=avgfilename)
  message("DONE")
