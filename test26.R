
setwd("/home/abayega/R/tests/test26")


homol <- read.table("Bo.all_proteins_UP000000803_blastp_sorted.txt.new",header=T,sep = '\t',quote = "",row.names = 1 )

homol_cov <- function(homol) {
  perc = c()
  num = c()
  k = homol$qcov2
  for (i in 1:100) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  plot(rev(perc),rev(num), type = "l",xlab = "Query (blue) or target (green) percentage coverage", col="blue",
       ylab = "number of genes", main = "percentage coverage across B. oleae Dmel homologues")
  #lines(perc,num,type="l",col="blue")
  perc = c()
  num = c()
  k = homol$scov2
  for (i in 1:100) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  #par(new=F)
  #plot(perc,num,pch=23,col="green",add=TRUE)
  lines(rev(perc),rev(num),type="l",col="green")
}
homol_cov(homol)


homol_cov <- function(homol) {
  perc = c()
  num = c()
  k = homol$qcov2
  for (i in 1:100) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  plot(perc,num, axes=F, type = "l",xlab = "Query (blue) or target (green) percentage coverage", col="blue",
       ylab = "number of genes", main = "percentage coverage across B. oleae Dmel homologues")
  axis(1, at = c(0,20,40,60,80,100), labels = c(100,80,60,40,20,0), tick = TRUE)
  perc = c()
  num = c()
  k = homol$scov2
  for (i in 1:100) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  #par(new=F)
  #plot(perc,num,pch=23,col="green",add=TRUE)
  lines(perc,num,type="l",col="green")
}
homol_cov(homol)



homol_cov <- function(homol) {
  n=100
  perc = c()
  num = c()
  k = homol$qcov
  for (i in c(100,90,80,70,60,50)) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  plot(perc,num)
}

homol_cov(homol)
dim(homol)
head(homol)
homol$qcov2 <- 100*homol$aln_length/homol$qlen
homol$scov2 <- 100*homol$aln_length/homol$slen

homol_cov <- function(homol) {
  n=100
  perc = c()
  num = c()
  k = homol$qcov2
  for (i in c(100,90,80,70,60,50)) {
    perc <- c(perc,i)
    num <- c(num,length(which(k>=i)))
  }
  plot(perc,num)
}

homol_cov(homol)
colnames(homol)<-c('sseqid',	'pident',	'aln_length',	'qstart',	'qend',	'sstart',	'send',	'qlen2',	'qlen','slen',	'qcov',	'scov',	'evalue','qcov2','scov2')
