rm(list = ls())

setwd('/home/abayega/R/tests/test48')

classn <- read.table("C1-15H_pass_FL_canu.flair.sqanti_classification.txt",header=T,sep = '\t')


head(classn)
dim(classn)

#Creating dmel zygotic dataframe
#make sure the 2 vectors u are creating look exactly the same. One of them has spaces around which u can remove using gsub
min_cov_less3 <- subset(classn, classn$min_cov <3)
dim(min_cov_less3)

boxplot(min_cov_less3$gene_exp)
summary(min_cov_less3$gene_exp)
length(which(min_cov_less3$gene_exp < 1))

min_cov_less3_nogeneExp <- subset(classn, classn$min_cov <3 & classn$gene_exp < 1)
write.table(min_cov_less3_nogeneExp, file="min_cov_less3_nogeneExp", row.names = F, quote = F, sep = "\t")

#Try ggplot violin plots
library(ggplot2)
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# Gene length
#I will only focus on D Renzis zygotic early category and leave dmel_zyg$gene_length, 
geneLen <- c(novel_gene_abund$count_fl)
Species <- c(rep('Novel', length(novel_gene_abund$count_fl)))
finalTable <- data.frame(geneLen, Species)

colnames(finalTable) <- c('reads', 'Novel')
finalTable$species <- as.factor(finalTable$Novel)

length(which(finalTable$reads>=5))
conf_FL <- subset(gene_abundances, as.character(gene_abundances$pbid) %in% as.character(novelgeneIDs$V1) & gene_abundances$count_fl>=5)
dim(conf_FL)

p <- ggplot(finalTable, aes(x=Novel, y=reads, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'species') +
  labs(y = "Number of supporting reads") + 
  ylim(0,30) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

summary(finalTable$reads)

# Let's try find novel gene abundances from Illumina data

novel_gene_abund_ill <- subset(illumina_gene_abund, as.character(illumina_gene_abund$gene_id) %in% as.character(novelgeneIDs$V1))
novel_gene_abund_ill <- subset(illumina_gene_abund, as.character(illumina_gene_abund$gene_id) %in% as.character(conf_FL$pbid) & illumina_gene_abund$Sum >= 0.2)

head(finalTable_ill)
dim(novel_gene_abund_ill)

geneLen_ill <- c(novel_gene_abund_ill$Sum)
Species_ill <- c(rep('Novel', length(novel_gene_abund_ill$Sum)))
finalTable_ill <- data.frame(geneLen_ill, Species_ill)

colnames(finalTable_ill) <- c('TPM', 'Novel')
finalTable_ill$species <- as.factor(finalTable_ill$Novel)
head(finalTable_ill)

p <- ggplot(finalTable_ill, aes(x=Novel, y=TPM, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'TPM') +
  labs(y = "TPM") + 
  ylim(0,30) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

summary(finalTable_ill$TPM)
length(which(finalTable_ill$TPM >=0.2))

write.table(novel_gene_abund_ill, file="High_quality_novel_genes", row.names = F, quote = F, sep = "\t")
