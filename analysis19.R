rm(list = ls())

setwd("/home/abayega/R/tests/test19")

raw_reads <- read.table("ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_raw_reads_for_isoforms_minimap_with_alignments.txt", header = T, sep = '\t')
canu_cor <- read.table("ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_canu_corrected_reads_for_isoforms_minimap_with_alignments.txt", header = T, sep = '\t')
lordec1 <- read.table("x1_minimap2_with_alignments.txt", header = T, sep = '\t')
lordec <- read.table("x1_minimap2_with_alignments_e2.txt", header = T, sep = '\t')
canu_strand <- read.table("canu3_x1_minimap2_with_alignments.txt", header = T, sep = '\t')
tapis_cor <- read.table("ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_tapis_corrected_reads_for_isoforms_minimap_with_alignments.txt", header = T, sep = '\t')
sqanti_cor <- read.table("ERCC92_104688_ref_gapfilled_joined_lt9474.gt500.covgt10_chrMT_and_UN_refseq_IDs_sqanti_corrected_reads_for_isoforms_minimap_with_alignments.txt", header = T, sep = '\t')



algn_id = cbind(raw_reads$alignment_identity, canu_cor$alignment_identity, canu_strand$alignment_identity,lordec1$alignment_identity, lordec$alignment_identity, tapis_cor$alignment_identity, sqanti_cor$alignment_identity)

boxplot.matrix(algn_id, use.cols = T, ylab="align Identity (%)",col = c('red', 'green', 'yellow','grey','purple', 'blue', 'yellow'), names=c("raw_reads","canu_cor","canu_cor_strand","lordec_R1_only","lordec_all_rds","tapis_cor","sqanti_cor"), main = 'Alignment Identity and corection')
#lordec_all_rds=lordec refers to performing error correction using 4 read sets, that is read 1 and read 2 of 5 hour timepoint and sex_organs_testes_after
#lordec_R1_only=lordec1 refers to using just the read1 of the 5 hour timepoint. There seems to be no added benefit in using so many data sets
#There is no benefit of stranding reads before canu correction (although I used just the first 100,000 reads for the stranded)

summary(lordec$alignment_identity)


##09th/Feb/2019
#I want to redraw the graphs with violin plots

library(ggplot2)
#using ggplot2
expn <- c(raw_reads$alignment_identity, canu_cor$alignment_identity, lordec1$alignment_identity)
fac <- c(rep('raw_reads', length(raw_reads$alignment_identity)), rep('canu_cor', length(canu_cor$alignment_identity)), rep('canu_lordec_cor', length(lordec1$alignment_identity)))
expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('align identity', 'reads')
expn_mat$fac <- as.factor(expn_mat$fac)
expn_mat$fac <- factor(expn_mat$fac, levels = c('raw_reads','canu_cor', 'canu_lordec_cor'),ordered = TRUE)
head(expn_mat)
p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  labs(x = "reads") +
  labs(fill = '') +
  ylim(c(75,100)) +
  labs(y = "alignment identity (%)") + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, face="bold")
  )
p
##And it works, yeeyi...