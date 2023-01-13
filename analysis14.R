rm(list = ls())
setwd("~/R/tests/test14")

library(ggplot2)
###Trial
# create a data frame
variety=rep(LETTERS[1:7], each=40)
treatment=rep(c("high","low"),each=20)
note=seq(1:280)+sample(1:150, 280, replace=T)
data=data.frame(variety, treatment ,  note)

# grouped boxplot
ggplot(data, aes(x=variety, y=note, fill=treatment)) + 
  geom_boxplot() 

#Now load the data
Bo.1H <- read.table("positive_strand_genes_common_positive_Bo_E_1H_C010_10_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
Bo.2H <- read.table("positive_strand_genes_common_positive_Bo_E_2H_C010_09_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
Bo.3H <- read.table("positive_strand_genes_common_positive_Bo_E_3H_C010_08_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
Bo.4H <- read.table("positive_strand_genes_common_positive_Bo_E_4H_C010_06_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
Bo.5H <- read.table("positive_strand_genes_common_positive_Bo.E.5H_all_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
Bo.6H <- read.table("positive_strand_genes_common_positive_Bo_E_6H_C010_07_bedtools_coverage.out.bed_gene_coverage", header = T, row.names = 1, sep = '\t')
head(Bo.1H)

#ERCC rows 
nBo.1H.ERCC <- 69
nBo.2H.ERCC <- 73
nBo.3H.ERCC <- 62
nBo.4H.ERCC <- 65
nBo.5H.ERCC <- 64
nBo.6H.ERCC <- 71 

Bo.1H <-subset.data.frame(Bo.1H, rat_cov_5 < 2.1)
Bo.1H <-subset.data.frame(Bo.1H, rat_cov_3 < 2.1)
Bo.2H <-subset.data.frame(Bo.2H, rat_cov_5 < 2.1)
Bo.2H <-subset.data.frame(Bo.2H, rat_cov_3 < 2.1)
Bo.3H <-subset.data.frame(Bo.3H, rat_cov_5 < 2.1)
Bo.3H <-subset.data.frame(Bo.3H, rat_cov_3 < 2.1)
Bo.4H <-subset.data.frame(Bo.4H, rat_cov_5 < 2.1)
Bo.4H <-subset.data.frame(Bo.4H, rat_cov_3 < 2.1)
Bo.5H <-subset.data.frame(Bo.5H, rat_cov_5 < 2.1)
Bo.5H <-subset.data.frame(Bo.5H, rat_cov_3 < 2.1)
Bo.6H <-subset.data.frame(Bo.6H, rat_cov_5 < 2.1)
Bo.6H <-subset.data.frame(Bo.6H, rat_cov_3 < 2.1)


#Grouped boxPlot coverage box
Bo.1H.genes_5 <- Bo.1H[(nBo.1H.ERCC+1):length(rownames(Bo.1H)),5]
Bo.1H.ERCC_5 <- Bo.1H[(1:nBo.1H.ERCC),5]
Bo.1H.genes_3 <- Bo.1H[(nBo.1H.ERCC+1):length(rownames(Bo.1H)),6]
Bo.1H.ERCC_3 <- Bo.1H[(1:nBo.1H.ERCC),6]

Bo.2H.genes_5 <- Bo.2H[(nBo.2H.ERCC+1):length(rownames(Bo.2H)),5]
Bo.2H.ERCC_5 <- Bo.2H[(1:nBo.2H.ERCC),5]
Bo.2H.genes_3 <- Bo.2H[(nBo.2H.ERCC+1):length(rownames(Bo.2H)),6]
Bo.2H.ERCC_3 <- Bo.2H[(1:nBo.2H.ERCC),6]

Bo.3H.genes_5 <- Bo.3H[(nBo.3H.ERCC+1):length(rownames(Bo.3H)),5]
Bo.3H.ERCC_5 <- Bo.3H[(1:nBo.3H.ERCC),5]
Bo.3H.genes_3 <- Bo.3H[(nBo.3H.ERCC+1):length(rownames(Bo.3H)),6]
Bo.3H.ERCC_3 <- Bo.3H[(1:nBo.3H.ERCC),6]

Bo.4H.genes_5 <- Bo.4H[(nBo.4H.ERCC+1):length(rownames(Bo.4H)),5]
Bo.4H.ERCC_5 <- Bo.4H[(1:nBo.4H.ERCC),5]
Bo.4H.genes_3 <- Bo.4H[(nBo.4H.ERCC+1):length(rownames(Bo.4H)),6]
Bo.4H.ERCC_3 <- Bo.4H[(1:nBo.4H.ERCC),6]

Bo.5H.genes_5 <- Bo.5H[(nBo.5H.ERCC+1):length(rownames(Bo.5H)),5]
Bo.5H.ERCC_5 <- Bo.5H[(1:nBo.5H.ERCC),5]
Bo.5H.genes_3 <- Bo.5H[(nBo.5H.ERCC+1):length(rownames(Bo.5H)),6]
Bo.5H.ERCC_3 <- Bo.5H[(1:nBo.5H.ERCC),6]

Bo.6H.genes_5 <- Bo.6H[(nBo.6H.ERCC+1):length(rownames(Bo.6H)),5]
Bo.6H.ERCC_5 <- Bo.6H[(1:nBo.6H.ERCC),5]
Bo.6H.genes_3 <- Bo.6H[(nBo.6H.ERCC+1):length(rownames(Bo.6H)),6]
Bo.6H.ERCC_3 <- Bo.6H[(1:nBo.6H.ERCC),6]

coverage.df <- cbind(Bo.1H.genes_5,Bo.1H.ERCC_5,Bo.1H.genes_3,Bo.1H.ERCC_3, Bo.2H.genes_5,Bo.2H.ERCC_5,Bo.2H.genes_3,Bo.2H.ERCC_3,
                                   Bo.3H.genes_5,Bo.3H.ERCC_5,Bo.3H.genes_3,Bo.3H.ERCC_3, Bo.4H.genes_5,Bo.4H.ERCC_5,Bo.4H.genes_3,Bo.4H.ERCC_3,
                                   Bo.5H.genes_5,Bo.5H.ERCC_5,Bo.5H.genes_3,Bo.5H.ERCC_3, Bo.6H.genes_5,Bo.6H.ERCC_5,Bo.6H.genes_3,Bo.6H.ERCC_3)
boxplot.matrix(coverage.df, use.cols = T, ylim=c(0,2))


time.points <- c(rep(c("Bo.1H","Bo.2H","Bo.3H","Bo.4H","Bo.5H","Bo.6H"), c((length(rownames(Bo.1H))*2),(length(rownames(Bo.2H))*2),(length(rownames(Bo.3H))*2),
                                                                           (length(rownames(Bo.4H))*2),(length(rownames(Bo.5H))*2),(length(rownames(Bo.6H))*2))))

#genes.ercc <- c(rep(c("Bo.1H.genes_5","Bo.1H.ERCC_5","Bo.1H.genes_3","Bo.1H.ERCC_3", 
#                      "Bo.2H.genes_5","Bo.2H.ERCC_5","Bo.2H.genes_3","Bo.2H.ERCC_3",
#                      "Bo.3H.genes_5","Bo.3H.ERCC_5","Bo.3H.genes_3","Bo.3H.ERCC_3",
#                      "Bo.4H.genes_5","Bo.4H.ERCC_5","Bo.4H.genes_3","Bo.4H.ERCC_3",
#                      "Bo.5H.genes_5","Bo.5H.ERCC_5","Bo.5H.genes_3","Bo.5H.ERCC_3",
#                      "Bo.6H.genes_5","Bo.6H.ERCC_5","Bo.6H.genes_3","Bo.6H.ERCC_3"), c(length(Bo.1H.genes_5),length(Bo.1H.ERCC_5),length(Bo.1H.genes_3),length(Bo.1H.ERCC_3),
#                                                                                             length(Bo.2H.genes_5),length(Bo.2H.ERCC_5),length(Bo.2H.genes_3),length(Bo.2H.ERCC_3),
#                                                                                             length(Bo.3H.genes_5),length(Bo.3H.ERCC_5),length(Bo.3H.genes_3),length(Bo.3H.ERCC_3),
#                                                                                             length(Bo.4H.genes_5),length(Bo.4H.ERCC_5),length(Bo.4H.genes_3),length(Bo.4H.ERCC_3),
#                                                                                             length(Bo.5H.genes_5),length(Bo.5H.ERCC_5),length(Bo.5H.genes_3),length(Bo.5H.ERCC_3),
#                                                                                             length(Bo.6H.genes_5),length(Bo.6H.ERCC_5),length(Bo.6H.genes_3),length(Bo.6H.ERCC_3))))

genes.ercc <- c(rep(c("genes_5","ERCC_5","genes_3","ERCC_3", 
                      "genes_5","ERCC_5","genes_3","ERCC_3",
                      "genes_5","ERCC_5","genes_3","ERCC_3",
                      "genes_5","ERCC_5","genes_3","ERCC_3",
                      "genes_5","ERCC_5","genes_3","ERCC_3",
                      "genes_5","ERCC_5","genes_3","ERCC_3"), c(length(Bo.1H.genes_5),length(Bo.1H.ERCC_5),length(Bo.1H.genes_3),length(Bo.1H.ERCC_3),
                                                                length(Bo.2H.genes_5),length(Bo.2H.ERCC_5),length(Bo.2H.genes_3),length(Bo.2H.ERCC_3),
                                                                length(Bo.3H.genes_5),length(Bo.3H.ERCC_5),length(Bo.3H.genes_3),length(Bo.3H.ERCC_3),
                                                                length(Bo.4H.genes_5),length(Bo.4H.ERCC_5),length(Bo.4H.genes_3),length(Bo.4H.ERCC_3),
                                                                length(Bo.5H.genes_5),length(Bo.5H.ERCC_5),length(Bo.5H.genes_3),length(Bo.5H.ERCC_3),
                                                                length(Bo.6H.genes_5),length(Bo.6H.ERCC_5),length(Bo.6H.genes_3),length(Bo.6H.ERCC_3))))


coverage <- c(Bo.1H.genes_5,Bo.1H.ERCC_5,Bo.1H.genes_3,Bo.1H.ERCC_3,
              Bo.2H.genes_5,Bo.2H.ERCC_5,Bo.2H.genes_3,Bo.2H.ERCC_3,
              Bo.3H.genes_5,Bo.3H.ERCC_5,Bo.3H.genes_3,Bo.3H.ERCC_3,
              Bo.4H.genes_5,Bo.4H.ERCC_5,Bo.4H.genes_3,Bo.4H.ERCC_3,
              Bo.5H.genes_5,Bo.5H.ERCC_5,Bo.5H.genes_3,Bo.5H.ERCC_3,
              Bo.6H.genes_5,Bo.6H.ERCC_5,Bo.6H.genes_3,Bo.6H.ERCC_3)

cov.data=data.frame(time.points,genes.ercc,coverage)

ggplot(data.frame(cov.data), aes(x=time.points, fill=genes.ercc)) + 
  geom_bar() 

ggplot(cov.data, aes(x=time.points, y=coverage, fill=genes.ercc)) + 
  geom_boxplot() #+ scale_y_continuous(limits=c(0,10)) #scale_y_discrete(limits=c("0", "2"))


##Doing some analysis on Bo. Gene expression
Bo.all.gex <- read.table("Bo.gene.expression", header = T, row.names = 1, sep = '\t')
head(Bo.all.gex)

#Sum expression of genes accross timepoints
k <- apply(Bo.all.gex,1,sum)
l <- apply(Bo.all.gex[,1:6],2,sum)
plot(l)

head(k)
typeof(k)
dim(k)
o <- c(1,2,3)
dim(o)

#Find genes with expression
(100*length(which(k>0)))/length(k)

Bo.all.gex$sum <- k
head(Bo.all.gex)

