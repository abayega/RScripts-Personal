
rm(list = ls())
setwd("/home/abayega/R/tests/test30")
library('gplots') #has the heatmap.2 function

exp_cor_ex <- read.table("illumina_exon_TPM_edited",header=T, row.names = 1, sep = '\t')
exp_cor_in <- read.table("illumina_intron_TPM_edited",header=T, row.names = 1, sep = '\t')
exp_cor <- read.table("illumina_exon_intron_TPM_embryos_only_edited",header=T, row.names = 1, sep = '\t')
exp_cor_zyg <- read.table("illumina_zygotic_genes_edited",header=T, row.names = 1, sep = '\t')

  
exp_cor <- cbind(exp_cor_ex$std_1H,exp_cor_in$std_1H,exp_cor_ex$std_2H,exp_cor_in$std_2H,exp_cor_ex$std_3H,exp_cor_in$std_3H,
                 exp_cor_ex$std_4H,exp_cor_in$std_4H,exp_cor_ex$std_5H,exp_cor_in$std_5H,exp_cor_ex$std_6H,exp_cor_in$std_6H,
                 exp_cor_ex$std_fem,exp_cor_in$std_fem,exp_cor_ex$std_mal,exp_cor_in$std_mal)


dim(exp_cor)

heatmap.2(as.matrix(exp_cor[,3:14]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE, key.xlab = 'adjusted_p-value',lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 ))

#main = "embryo intron retention by illumina"
exp_cor2 <- head(exp_cor, 1000)
heatmap.2(as.matrix(exp_cor2[,3:14]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE, key.xlab = 'adjusted_p-value',lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 ))

#main = "zygotic embryo intron retention by illumina"
heatmap.2(as.matrix(exp_cor_zyg[,3:14]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE, key.xlab = 'adjusted_p-value',lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 ))

#Now trying Nanopore
exp_cor_ont <- read.table("ont_introns_exons_mandalorion_edited",header=T, row.names = 1, sep = '\t') 
exp_cor_ont_zyg <- read.table("ont_zygotic_genes_edited2",header=T, row.names = 1, sep = '\t')

head(exp_cor_ont)
dim(exp_cor_ont)
#main = "embryo intron retention by nanopore"
heatmap.2(as.matrix(head(exp_cor_ont,3668)[,5:20]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE,lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 ))

#main = "zygotic embryo intron retention by nanopore"
heatmap.2(as.matrix(exp_cor_ont_zyg[,5:20]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE,lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 ))

##070919 many months later
exp_cor_ont <- read.table("ont_introns_exons_mandalorion_edited.cluster11",header=T, row.names = 1, sep = '\t') 
heatmap.2(as.matrix(exp_cor_ont[,5:20]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE,lwid = c(.01,.99),lhei = c(.01,.99),margins = c(10,20 )) #Title is ont_introns_exons_mandalorion_edited.cluster11-21-28-41
