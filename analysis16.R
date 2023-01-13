rm(list = ls())

setwd("/home/abayega/R/tests/test16")

gene_exp <- read.table("gene_RPG10", header = T, row.names = 1, sep = '\t')
ercc_exp <- read.table("ERCC_expression", header = T, row.names = 1, sep = '\t')
head(ercc_exp)

exp_1 <- length(which(gene_exp$X1H_TPM > 0))
exp_2 <- length(which(gene_exp$X2H_TPM > 0))
exp_3 <- length(which(gene_exp$X3H_TPM > 0))
exp_4 <- length(which(gene_exp$X4H_TPM > 0))
exp_5 <- length(which(gene_exp$X5H_TPM > 0))
exp_6 <- length(which(gene_exp$X6H_TPM > 0))

barplot(height = c(exp_1,exp_2,exp_3,exp_4,exp_5,exp_6), ylim = c(0,10000), names.arg = c("Bo.E.1H","Bo.E.2H","Bo.E.3H","Bo.E.4H","Bo.E.5H","Bo.E.6H"),
        main = "Number of genes expressed at different time points", col = "lightblue", ylab = "Number of genes expressed", xlab = "Time points")

exp_1 <- sum(gene_exp$X1H_TPM)
exp_2 <- sum(gene_exp$X2H_TPM)
exp_3 <- sum(gene_exp$X3H_TPM)
exp_4 <- sum(gene_exp$X4H_TPM)
exp_5 <- sum(gene_exp$X5H_TPM)
exp_6 <- sum(gene_exp$X6H_TPM)

barplot(height = c(exp_1,exp_2,exp_3,exp_4,exp_5,exp_6), ylim = c(0,10000), names.arg = c("1H","2H","3H","4H","5H","6H"),
        main = "Gene expression at different time points", col = "lightblue", ylab = "Total RPG10K", xlab = "Time points", cex.lab=1.5, cex.names = 1.5, cex.axis = 1.25)
barplot(height = c(exp_1,exp_2,exp_3,exp_4,exp_5,exp_6), ylim = c(0,10000), names.arg = c("1H","2H","3H","4H","5H","6H"),
        col = "lightblue", ylab = "Total RPG10K", xlab = "Time points", cex.lab=1.5, cex.names = 1.5, cex.axis = 1.25)

barplot(height = c(1258982557.9952,	636393756.498989,	1411143919.55606,	1324591841.61091,	1261153748.98927,	876890181.100687), names.arg = c("Bo.E.1H","Bo.E.2H","Bo.E.3H","Bo.E.4H","Bo.E.5H","Bo.E.6H"),
        main = "Gene expression at different time points", col = "lightblue", ylab = "Transcripts per embryo ", xlab = "Time points")

par(mgp=c(2.2,1,0),mar=c(4,4,4,2)+0.1)
barplot(height = c(12.6,	6.4,	14.1,	13.2,	12.6,	8.8), names.arg = c("1H","2H","3H","4H","5H","6H"),
        main = "Gene expression at different time points", col = "lightblue", ylab = expression(paste("Transcripts per embryo (x 10"^8*")",sep = '')), xlab = "Time points", cex.lab=1.5, cex.names = 1.5, cex.axis = 1.5)
barplot(height = c(12.6,	6.4,	14.1,	13.2,	12.6,	8.8), names.arg = c("1H","2H","3H","4H","5H","6H"),
        col = "lightblue", ylab = expression(paste("Transcripts per embryo (x 10"^8*")",sep = '')), xlab = "Time points", cex.lab=1.5, cex.names = 1.5, cex.axis = 1.5)

#Looking at ERCC expression
ercc_exp_1 <- sum(ercc_exp$X1H_TPM)
ercc_exp_2 <- sum(ercc_exp$X2H_TPM)
ercc_exp_3 <- sum(ercc_exp$X3H_TPM)
ercc_exp_4 <- sum(ercc_exp$X4H_TPM)
ercc_exp_5 <- sum(ercc_exp$X5H_TPM)
ercc_exp_6 <- sum(ercc_exp$X6H_TPM)

barplot(height = c(ercc_exp_1,	ercc_exp_2,	ercc_exp_3,	ercc_exp_4,	ercc_exp_5,	ercc_exp_6), names.arg = c("Bo.E.1H","Bo.E.2H","Bo.E.3H","Bo.E.4H","Bo.E.5H","Bo.E.6H"),
        main = "ERCC expression at different time points", col = "lightblue", ylab = "RPG10K", xlab = "Time points")

ercc_exp_1 <- sum(ercc_exp$X1H_abs_emb)
ercc_exp_2 <- sum(ercc_exp$X2H_abs_emb)
ercc_exp_3 <- sum(ercc_exp$X3H_abs_emb)
ercc_exp_4 <- sum(ercc_exp$X4H_abs_emb)
ercc_exp_5 <- sum(ercc_exp$X5H_abs_emb)
ercc_exp_6 <- sum(ercc_exp$X6H_abs_emb)
barplot(height = c(ercc_exp_1,	ercc_exp_2,	ercc_exp_3,	ercc_exp_4,	ercc_exp_5,	ercc_exp_6), names.arg = c("Bo.E.1H","Bo.E.2H","Bo.E.3H","Bo.E.4H","Bo.E.5H","Bo.E.6H"),
        main = "ERCC expression per embryo at different time points", col = "lightblue", ylab = "Transcripts per embryo", xlab = "Time points")
