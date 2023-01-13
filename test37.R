#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test37")
library('gplots') #has the heatmap.2 function

exp_cor_ex <- read.table("boleae_final_transcriptsperembryo",header=T, row.names = 1, sep = '\t')

#remove adult heads
exp_cor_ex_1 <- exp_cor_ex[,1:6]
head(exp_cor_ex_1)

#set lower cut off to 1000 transcripts per embryo, this corresponds to ~0.01 RPG10K
exp_cor_ex_1[exp_cor_ex_1 <= 1000] = 0
head(exp_cor_ex_1)

#Remove all genes with no expression across the 6 timepoints
exp_cor_ex_2 <- exp_cor_ex_1[rowSums(exp_cor_ex_1) != 0,]
head(exp_cor_ex_2)

#Now let us try to get z-scores for the genes (x-mean)/SD
exp_cor_ex_2$mean <-apply(exp_cor_ex_2[,1:6],1,mean)
#exp_cor_ex_2$mean <- rowMeans(exp_cor_ex_2)

#standard deviation
exp_cor_ex_2$StdDev <-apply(exp_cor_ex_2[,1:6],1,sd)

#Now get zscores
exp_cor_ex_2$abs_emb_1H_z <- ((exp_cor_ex_2$abs_emb_1H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)
exp_cor_ex_2$abs_emb_2H_z <- ((exp_cor_ex_2$abs_emb_2H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)
exp_cor_ex_2$abs_emb_3H_z <- ((exp_cor_ex_2$abs_emb_3H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)
exp_cor_ex_2$abs_emb_4H_z <- ((exp_cor_ex_2$abs_emb_4H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)
exp_cor_ex_2$abs_emb_5H_z <- ((exp_cor_ex_2$abs_emb_5H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)
exp_cor_ex_2$abs_emb_6H_z <- ((exp_cor_ex_2$abs_emb_6H - exp_cor_ex_2$mean)/exp_cor_ex_2$StdDev)

heatmap.2(as.matrix(exp_cor_ex_2[,9:14]), Colv=FALSE, dendrogram='none') #,key.xlab = 'z-score', trace='none')

#Now let us try to get genes 
mat <- read.table("mat")
head(mat) 
mat <- mat$V2

mat_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% mat)

filetowrite <- mat_zscores[,9:14]
write.table(filetowrite, file = "mat_zscores.csv", sep = ',', quote = F)

#Now get genes spiking at 3 hrs AEL
three <- read.table("three")
head(three) 
three <- three$V2

three_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% three)

filetowrite <- three_zscores[,9:14]
write.table(filetowrite, file = "three_zscores.csv", sep = ',', quote = F)

#Now get genes spiking at 4 hrs AEL
four <- read.table("four")
head(four) 
four <- four$V2

four_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% four)

filetowrite <- four_zscores[,9:14]
write.table(filetowrite, file = "four_zscores.csv", sep = ',', quote = F)

#Now get genes spiking at 5 hrs AEL
five <- read.table("five")
head(five) 
five <- five$V2

five_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% five)

filetowrite <- five_zscores[,9:14]
write.table(filetowrite, file = "five_zscores.csv", sep = ',', quote = F)

#Now get genes spiking at 6 hrs AEL
six <- read.table("six")
head(six) 
six <- six$V2

six_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% six)

filetowrite <- six_zscores[,9:14]
write.table(filetowrite, file = "six_zscores.csv", sep = ',', quote = F)

#Now get zygotic genes
zyg <- read.table("zyg")
head(zyg) 
zyg <- zyg$V2

zyg_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% zyg)

filetowrite <- zyg_zscores[,9:14]
write.table(filetowrite, file = "zyg_zscores.csv", sep = ',', quote = F)
