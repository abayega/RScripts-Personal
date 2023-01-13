
rm(list = ls())
setwd("/home/abayega/R/tests/test46")
library('gplots') #has the heatmap.2 function

exp_cor_ex <- read.table("AllvAll_combined",header=T, row.names = 1, sep = '\t')

head(exp_cor_ex_TPE)

#Take only the transcripts per embryo and remove ERCC
exp_cor_ex_TPE <- exp_cor_ex[93:nrow(exp_cor_ex),24:29]

#Let us make low expressed genes, defined as those with TPE below detectable limits, equal to 0
exp_cor_ex_TPE$TPE_1H[exp_cor_ex_TPE$TPE_1H < 1110] = 0
exp_cor_ex_TPE$TPE_2H[exp_cor_ex_TPE$TPE_2H < 441] = 0
exp_cor_ex_TPE$TPE_3H[exp_cor_ex_TPE$TPE_3H < 922] = 0
exp_cor_ex_TPE$TPE_4H[exp_cor_ex_TPE$TPE_4H < 1371] = 0
exp_cor_ex_TPE$TPE_5H[exp_cor_ex_TPE$TPE_5H < 1546] = 0
exp_cor_ex_TPE$TPE_6H[exp_cor_ex_TPE$TPE_6H < 841] = 0
#exp_cor_ex_TPE[exp_cor_ex_TPE<500] = 0

#get the sums
exp_cor_ex_TPE$RowSum <- apply(exp_cor_ex_TPE,1,sum)

#remove rows genes that are not expressed
exp_cor_ex_TPE <- subset(exp_cor_ex_TPE, RowSum>0)

#let us try to log the values
exp_cor_ex <- log(exp_cor_ex_TPE+1)

#let us try to find most variable genes via standard deviation and coefficient of variation
exp_cor_ex$stdev <- apply(exp_cor_ex, 1, sd )
exp_cor_ex$mean <- apply(exp_cor_ex, 1, mean )
exp_cor_ex$covar <- exp_cor_ex$stdev/exp_cor_ex$mean

#try to get zscores
logExp_cor_ex <- exp_cor_ex
logExp_cor_ex$TPE_1H_z <- (exp_cor_ex$TPE_1H - exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$TPE_2H_z <- (exp_cor_ex$TPE_2H - exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$TPE_3H_z <- (exp_cor_ex$TPE_3H - exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$TPE_4H_z <- (exp_cor_ex$TPE_4H - exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$TPE_5H_z <- (exp_cor_ex$TPE_5H - exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$TPE_6H_z <- (exp_cor_ex$TPE_6H - exp_cor_ex$mean)/exp_cor_ex$stdev

#now, based on the z-score try to get the most differentially expressed genes via standard deviation and coefficient of variation
logExp_cor_ex <- logExp_cor_ex[,11:ncol(logExp_cor_ex)]
logExp_cor_ex$stdev <- apply(logExp_cor_ex, 1, sd )
logExp_cor_ex$mean <- apply(logExp_cor_ex, 1, mean )
logExp_cor_ex$covar <- logExp_cor_ex$stdev/logExp_cor_ex$mean

#Let's order rows with coefficient of variation from top to lowest
logExp_cor_ex2 <- logExp_cor_ex[order(logExp_cor_ex$covar, decreasing = F),]

#Let us write the results to file
write.table(logExp_cor_ex2, file = "Gene_expression_ordered_according_to_variability_in_TPE", sep = "\t", quote = F)

##010620 
#We try to find number of protein conding genes that are expressed at 1 hour AEL
exp_cor_ex <- read.table("boleae_1_hour_AEL_TPE",header=T, sep = '\t')
exp_cor_ex2 <- read.table("Boleae_NCBI_predicted_proteinCodingGenesNames",header=F)
exp_cor_ex3 <- read.table("boleae_1-6_hours_AEL_TPE",header=T, sep = '\t')
head(exp_cor_ex3)

length(which(exp_cor_ex$Transcripts_per_embryo>=1110 & exp_cor_ex$Genes %in% exp_cor_ex2$V1))

#Number of genes expressed at detection limit
length(which(exp_cor_ex3$reference != 'ERCC' & exp_cor_ex3$TPE_1H >= 1110))
length(which(exp_cor_ex3$reference != 'ERCC' & exp_cor_ex3$TPE_2H >= 441))
length(which(exp_cor_ex3$reference != 'ERCC' & exp_cor_ex3$TPE_2H >= 922))
