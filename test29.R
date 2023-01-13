k <- read.table("/home/abayega/R/tests/test29/all_tags2017_readlengths.txt",header=F, sep = '\t')
boxplot(k$V2, ylim=c(0,3000), main="read length distribution of tagging sequences")
colnames(k) <- c("seq_name", "length")

library(ggplot2)

time.points <- c(rep(c("read_length"), each=(nrow(k))))

ftr_mat <- c(k$length)

cov.data=data.frame(time.points,ftr_mat)

p <- ggplot(cov.data, aes(x=time.points, y=ftr_mat)) + #, color=fac)
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  ggtitle("read length distribution of tagging sequences") +
  labs(x = "tags") +
  labs(fill = '') +
  #ylim(c(0,3000)) +
  labs(y = "read length (bases)") + 
  theme(
    plot.title = element_text(color="black", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, face="bold")
  )
p


k <- read.table("/home/abayega/R/tests/test29/contig_assignmentbyprobes_plus_genes",header=T, sep = '\t')
k <- read.table("/home/abayega/R/tests/test29/contig_assignmentbyblast_plus_genes",header=T, sep = '\t')
head(k)

chr1 = k[k$chromosome=="I",]
chr1 = sum(chr1$No._genes)


chr2 = k[k$chromosome=="II",]
chr2 = sum(chr2$No._genes)

chr3 = k[k$chromosome=="III",]
chr3 = sum(chr3$No._genes)

chr4 = k[k$chromosome=="IV",]
chr4 = sum(chr4$No._genes)

chr5 = k[k$chromosome=="V",]
chr5 = sum(chr5$No._genes)

chrUn = k[k$chromosome=="Un",]
chrUn = sum(chrUn$No._genes)

chrom = c(chr1,chr2,chr3,chr4,chr5,chrUn)
chrom2 = c(chr1,chr2,chr3,chr4,chr5,chrUn)
chromb = cbind.data.frame(chrom,chrom2)

barplot(chromb, main = "Number of genes per chromosome for probe assignment", names.arg = c("I","II","III","IV","V","Un"), xlab = 'chromosomes', ylab='Total No. of genes')

        