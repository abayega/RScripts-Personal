rm(list = ls())
setwd("/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_5H/combined/minimap2_all_ncbi2/process")

file1="104688_ref_gapfilled_uncorrected_minimap2_with_alignments.sorted.txt"
file2="104688_ref_gapfilled_canu_corrected_minimap2_with_alignments.sorted.txt"
file3="104688_ref_gapfilled_canu_corr_uncorr_minimap2_with_alignments.sorted.txt"
file4="104688_ref_gapfilled_b4_canu_correct_minimap2_with_alignments.sorted.txt"

outdir <- "/mnt/parallel_scratch_mp2_wipe_on_december_2018/ioannisr/banthony/analysis/nanopore/olivefly/rna_seq2/Bo_E_5H/combined/minimap2_all_ncbi2/process"
outprefix <- "out"
pngpdfdir <- sprintf("%s/pngpdf", outdir)

if (! file.exists(pngpdfdir)){
    dir.create(pngpdfdir)
}

library(ggplot2)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

#Load the data
uncor <- read.table(file1, header = T, sep = '\t')
canu_cor <- read.table(file2, header = T, sep = '\t')
uncor_canu_cor <- read.table(file3, header = T, sep = '\t')
rds_b4_cor <- read.table(file4, header = T, sep = '\t')

head(rds_b4_cor)

#ERCC rows 
uncor_ERCC <- a <- 90471
canu_cor_ERCC <- b <- 69205
uncor_canu_cor_ERCC <- c <- 75359
rds_b4_cor_ERCC <- d <- 68730

#ploting ERCC correction parameters
time.points <- c(rep(c("Aln-identity","Insertions","Deletions","Perc_aln"), each=(uncor_ERCC+canu_cor_ERCC+uncor_canu_cor_ERCC+rds_b4_cor_ERCC)))

correction <- c(rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)))

percs <- c(uncor$alignment_identity[1:a],canu_cor$alignment_identity[1:b],uncor_canu_cor$alignment_identity[1:c],rds_b4_cor$alignment_identity[1:d],
           uncor$insertion[1:a],canu_cor$insertion[1:b],uncor_canu_cor$insertion[1:c],rds_b4_cor$insertion[1:d],
           uncor$deletion[1:a],canu_cor$deletion[1:b],uncor_canu_cor$deletion[1:c],rds_b4_cor$deletion[1:d],
           uncor$perc_read_aligned[1:a],canu_cor$perc_read_aligned[1:b],uncor_canu_cor$perc_read_aligned[1:c],rds_b4_cor$perc_read_aligned[1:d])

cov.data=data.frame(time.points,correction,percs)

pngpath <- sprintf('%s/%s_Assessing correction of ERCC reads using Canu.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
p <- ggplot(cov.data, aes(x=time.points, y=percs, fill=correction)) + 
  geom_boxplot() +
  ggtitle("Assessing correction of ERCC reads using Canu") +
  xlab("Alignment parameters") +
  ylab("Percentage") +
  my.theme
ggsave(pngpath)

#ploting gene correction parameters
#Gene rows 
uncor_gene <- nrow(uncor) - 90471
canu_cor_gene <- nrow(canu_cor) - 69205
uncor_canu_cor_gene <- nrow(uncor_canu_cor) - 75359
rds_b4_cor_gene <- nrow(rds_b4_cor) - 68730

time.points <- c(rep(c("Aln-identity","Insertions","Deletions","Perc_aln"), each=(uncor_gene+canu_cor_gene+uncor_canu_cor_gene+rds_b4_cor_gene)))

correction <- c(rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)))

percs <- c(uncor$alignment_identity[(a+1):nrow(uncor)],canu_cor$alignment_identity[(b+1):nrow(canu_cor)],uncor_canu_cor$alignment_identity[(c+1):nrow(uncor_canu_cor)],rds_b4_cor$alignment_identity[(d+1):nrow(rds_b4_cor)],
           uncor$insertion[(a+1):nrow(uncor)],canu_cor$insertion[(b+1):nrow(canu_cor)],uncor_canu_cor$insertion[(c+1):nrow(uncor_canu_cor)],rds_b4_cor$insertion[(d+1):nrow(rds_b4_cor)],
           uncor$deletion[(a+1):nrow(uncor)],canu_cor$deletion[(b+1):nrow(canu_cor)],uncor_canu_cor$deletion[(c+1):nrow(uncor_canu_cor)],rds_b4_cor$deletion[(d+1):nrow(rds_b4_cor)],
           uncor$perc_read_aligned[(a+1):nrow(uncor)],canu_cor$perc_read_aligned[(b+1):nrow(canu_cor)],uncor_canu_cor$perc_read_aligned[(c+1):nrow(uncor_canu_cor)],rds_b4_cor$perc_read_aligned[(d+1):nrow(rds_b4_cor)])

cov.data=data.frame(time.points,correction,percs)

pngpath <- sprintf('%s/%s_Assessing correction of gene reads using Canu.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
p <- ggplot(cov.data, aes(x=time.points, y=percs, fill=correction)) + 
  geom_boxplot() +
  ggtitle("Assessing correction of gene reads using Canu") +
  xlab("Alignment parameters") +
  ylab("Percentage") +
  my.theme
ggsave(pngpath)



#library(ROCR)
#library(caret)
#library(e1071)
#library(randomForest)
#library(partykit)
#library(ipred)
#library(rpart.plot)
#library(ROSE)
#library(pROC)
#library(MLmetrics)
