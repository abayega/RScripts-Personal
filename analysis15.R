rm(list = ls())
setwd("~/R/tests/test15")

library(ggplot2)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

#Load the data
uncor <- read.table("104688_ref_gapfilled_uncorrected_minimap2_with_alignments.sorted.txt", header = T, sep = '\t')
canu_cor <- read.table("104688_ref_gapfilled_canu_corrected_minimap2_with_alignments.sorted.txt", header = T, sep = '\t')
uncor_canu_cor <- read.table("104688_ref_gapfilled_canu_correctednUncorrected_minimap2_with_alignments.sorted.txt", header = T, sep = '\t')
rds_b4_cor <- read.table("104688_ref_gapfilled_reads_names_b4_correction_2_minimap2_with_alignments.sorted.txt", header = T, sep = '\t')

head(rds_b4_cor)

#ERCC rows 
uncor_ERCC <- 95250
canu_cor_ERCC <- 57654
uncor_canu_cor_ERCC <- 74837
rds_b4_cor_ERCC <- 59206

#ploting ERCC correction parameters
time.points <- c(rep(c("Aln-identity","Insertions","Deletions","Perc_aln"), each=(uncor_ERCC+canu_cor_ERCC+uncor_canu_cor_ERCC+rds_b4_cor_ERCC)))

correction <- c(rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_ERCC,canu_cor_ERCC,uncor_canu_cor_ERCC,rds_b4_cor_ERCC)))

percs <- c(uncor$alignment_identity[1:95250],canu_cor$alignment_identity[1:57654],uncor_canu_cor$alignment_identity[1:74837],rds_b4_cor$alignment_identity[1:59206],
           uncor$insertion[1:95250],canu_cor$insertion[1:57654],uncor_canu_cor$insertion[1:74837],rds_b4_cor$insertion[1:59206],
           uncor$deletion[1:95250],canu_cor$deletion[1:57654],uncor_canu_cor$deletion[1:74837],rds_b4_cor$deletion[1:59206],
           uncor$perc_read_aligned[1:95250],canu_cor$perc_read_aligned[1:57654],uncor_canu_cor$perc_read_aligned[1:74837],rds_b4_cor$perc_read_aligned[1:59206])

cov.data=data.frame(time.points,correction,percs)

ggplot(cov.data, aes(x=time.points, y=percs, fill=correction)) + 
  geom_boxplot() +
  ggtitle("Assessing correction of ERCC reads using Canu") +
  xlab("Alignment parameters") +
  ylab("Percentage") +
  my.theme


#ploting gene correction parameters
#Gene rows 
uncor_gene <- nrow(uncor) - 95250
canu_cor_gene <- nrow(canu_cor) - 57654
uncor_canu_cor_gene<- nrow(uncor_canu_cor) - 74837
rds_b4_cor_gene <- nrow(rds_b4_cor) - 59206

time.points <- c(rep(c("Aln-identity","Insertions","Deletions","Perc_aln"), each=(uncor_gene+canu_cor_gene+uncor_canu_cor_gene+rds_b4_cor_gene)))

correction <- c(rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)),
                rep(c("uncor","canu_cor","cor_uncor","b4_cor"), c(uncor_gene,canu_cor_gene,uncor_canu_cor_gene,rds_b4_cor_gene)))

percs <- c(uncor$alignment_identity[95251:nrow(uncor)],canu_cor$alignment_identity[57655:nrow(canu_cor)],uncor_canu_cor$alignment_identity[74838:nrow(uncor_canu_cor)],rds_b4_cor$alignment_identity[59207:nrow(rds_b4_cor)],
           uncor$insertion[95251:nrow(uncor)],canu_cor$insertion[57655:nrow(canu_cor)],uncor_canu_cor$insertion[74838:nrow(uncor_canu_cor)],rds_b4_cor$insertion[59207:nrow(rds_b4_cor)],
           uncor$deletion[95251:nrow(uncor)],canu_cor$deletion[57655:nrow(canu_cor)],uncor_canu_cor$deletion[74838:nrow(uncor_canu_cor)],rds_b4_cor$deletion[59207:nrow(rds_b4_cor)],
           uncor$perc_read_aligned[95251:nrow(uncor)],canu_cor$perc_read_aligned[57655:nrow(canu_cor)],uncor_canu_cor$perc_read_aligned[74838:nrow(uncor_canu_cor)],rds_b4_cor$perc_read_aligned[59207:nrow(rds_b4_cor)])

cov.data=data.frame(time.points,correction,percs)

ggplot(cov.data, aes(x=time.points, y=percs, fill=correction)) + 
  geom_boxplot() +
  ggtitle("Assessing correction of gene reads using Canu") +
  xlab("Alignment parameters") +
  ylab("Percentage") +
  my.theme




library(ROCR)
library(caret)
library(e1071)
library(randomForest)
library(partykit)
library(ipred)
library(rpart.plot)
library(ROSE)
library(pROC)
library(MLmetrics)
