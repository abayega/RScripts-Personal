rm(list = ls())

setwd("/home/abayega/R/tests/test17")

library(data.table)
canu_corr <- fread("Bo_E_all_pass_edited.correctedReads.names4")

canu_corr <-as.data.frame(canu_corr)

#Set the column names
colnames(canu_corr) <- c("read.names","after.canu.corr","b4.canu.corr","diff")

head(canu_corr)
dim(canu_corr)


canu_corr$V4 <- canu_corr$V3 - canu_corr$V2 
boxplot(canu_corr$V4)

sum(canu_corr$V4 == 0)
sum(canu_corr$V4 > 0)
sum(canu_corr$V4 > 10)
sum(canu_corr$V4 >= -10 & canu_corr$V4 <= 10)
sum(canu_corr$V4 >= -20 & canu_corr$V4 <= 20)
sum(canu_corr$V4 >= -30 & canu_corr$V4 <= 30)
sum(canu_corr$V4 >= -40 & canu_corr$V4 <= 40)
sum(canu_corr$V4 >= -50 & canu_corr$V4 <= 50)
sum(canu_corr$V4 >= -60 & canu_corr$V4 <= 60)
sum(canu_corr$V4 >= -70 & canu_corr$V4 <= 70)
sum(canu_corr$V4 >= -80 & canu_corr$V4 <= 80)
sum(canu_corr$V4 >= -90 & canu_corr$V4 <= 90)
sum(canu_corr$V4 >= -100 & canu_corr$V4 <= 100)

sum(canu_corr$V4 >= -100 & canu_corr$V4 <= 200)
k <- canu_corr[which(canu_corr$diff == 200)]

plot(canu_corr$V3, canu_corr$V2)
