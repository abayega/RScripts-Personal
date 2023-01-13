#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test49")
library('gplots') #has the heatmap.2 function
library('ggplot2')

exp_cor_ex <- read.table("C1-15H_counts_matrix.tsv.tpm_edited.tsv",header=T, row.names = 1, sep = '\t')
delete <- read.table("delete",header=T, row.names = 1, sep = '\t')

dim(exp_cor_ex)
#Remove all rows that have 0 expression
#delete1 <- delete[ delete != 0,]

#exp_cor_ex = subset(exp_cor_ex, apply(exp_cor_ex, 1, function(row) all(row !=0 )))
exp_cor_ex = subset(exp_cor_ex, rowSums(exp_cor_ex[2:ncol(exp_cor_ex)]) != 0)

#I want to do PCA to see how data looks

#let us try to log the values
logExp_cor_ex <- log(exp_cor_ex[1:ncol(exp_cor_ex)]+1)

#let us try to find most variable genes via standard deviation and coefficient of variation
logExp_cor_ex1 <- logExp_cor_ex
logExp_cor_ex1$stdev <- apply(logExp_cor_ex, 1, sd )
logExp_cor_ex1$mean <- apply(logExp_cor_ex, 1, mean )
logExp_cor_ex1$covar <- logExp_cor_ex1$stdev/logExp_cor_ex1$mean

#Let's order rows with coefficient of variation from top to lowest
logExp_cor_ex_noERCC <- logExp_cor_ex1[order(logExp_cor_ex1$covar, decreasing = T),]

#Take the most varrying transcripts
logExp_cor_ex_noERCC_4k <- head(logExp_cor_ex_noERCC, n=4000) 

#try to get zscores
logExp_cor_ex_noERCC$stdev <- apply(logExp_cor_ex_noERCC, 1, sd )
logExp_cor_ex_noERCC$mean <- apply(logExp_cor_ex_noERCC, 1, mean )
logExp_cor_ex_noERCC_zsc <- apply(logExp_cor_ex_noERCC, 2, function(a) ((a- rowMeans(logExp_cor_ex_noERCC[,1:150]))/logExp_cor_ex_noERCC$stdev))

#logExp_cor_ex_noERCC_4k_zsc <- logExp_cor_ex_noERCC_4k[,1:(ncol(logExp_cor_ex_noERCC_4k)-3)]
#logExp_cor_ex_noERCC_4k_zsc <- apply(logExp_cor_ex_noERCC_4k_zsc, 2, function(a) ((a- rowMeans(logExp_cor_ex_noERCC_4k_zsc))/logExp_cor_ex_noERCC_4k$stdev))
#logExp_cor_ex_noERCC_4k_zsc <- as.data.frame(logExp_cor_ex_noERCC_4k_zsc)

#now, based on the z-score try to get the most differentially expressed genes via standard deviation and coefficient of variation
logExp_cor_ex_noERCC_4k_zsc1 <- as.data.frame(logExp_cor_ex_noERCC_zsc)
#logExp_cor_ex_noERCC_4k_zsc1$stdev <- apply(logExp_cor_ex_noERCC_4k_zsc1, 1, sd )
#logExp_cor_ex_noERCC_4k_zsc1$mean <- apply(logExp_cor_ex_noERCC_4k_zsc1, 1, mean )
#logExp_cor_ex_noERCC_4k_zsc1$covar <- logExp_cor_ex_noERCC_4k_zsc1$stdev/logExp_cor_ex_noERCC_4k_zsc1$mean

#Lets order rows with coefficient of variation from top to lowest
logExp_cor_ex_noERCC_4k_zsc1 <- logExp_cor_ex_noERCC_4k_zsc1[order(logExp_cor_ex_noERCC_4k_zsc1$covar, decreasing = T),]
logExp_cor_ex_noERCC_4k_zsc1 = subset(logExp_cor_ex_noERCC_4k_zsc1, rowSums(logExp_cor_ex_noERCC_4k_zsc1[1:ncol(logExp_cor_ex_noERCC_4k_zsc1)]) != 0 & stdev != 1)

##try PCA
#set.seed(42)
#pcout <- princomp(matrix(rnorm(1000), 100, 10))
#biplot(pcout)
#biplot(pcout, xlabs=rep(".", dim(pcout$scores)[1])) # small symbol
#biplot(pcout, xlabs=rep("", dim(pcout$scores)[1])) # no symbol 

library('stringr')
prin_comp <- prcomp(t(logExp_cor_ex_noERCC_zsc[1:1000,1:150]), scale. = T)
edit_colnames <- function(setx){
  new_sample_ID = c()
  new_sex = c()
  for (i in rownames(setx)){
    samplex = str_extract(i, regex("C.+H"))
    sex = strsplit(i,"")[[1]][length(strsplit(i,"")[[1]])]
    new_sample_ID = c(new_sample_ID,samplex)
    new_sex = c(new_sex,sex)
    #sex = str_extract(setx, regex("H."))
    #splitx = strsplit(i,"C")[[1]][2]
    #for (x in splitx){
    #  if (x = "F" || x = "M")
    #    sex = x
    #  }
  }
  return(c(new_sample_ID,new_sex))
}

new_prin_comp = as.data.frame(prin_comp$x)
new_prin_comp$sample_id = edit_colnames(new_prin_comp)[1:nrow(new_prin_comp)]
new_prin_comp$sex = edit_colnames(new_prin_comp)[(nrow(new_prin_comp)+1):(nrow(new_prin_comp)*2)]

##get some details about the results
summary(prin_comp)
names(prin_comp)
prin_comp$x[1:5,1:5]
# Eigenvalues
eig <- (prin_comp$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.cumvar <- data.frame(eig = eig, variance = variance,
                         cumvariance = cumvar)

#ggplot(as.data.frame(prin_comp$x), aes(x=PC1,y=PC2)) + geom_point(size=4) +
#ggplot(as.data.frame(new_prin_comp), aes(x=PC1,y=PC2)) + geom_point(size=4) +
ggplot(new_prin_comp, aes(x=PC1,y=PC2, label=sample_id, color=sex)) +
  geom_label(aes(fill = sample_id), colour = "white", fontface = "bold")+
  theme_bw(base_size=32) +
  labs(x=paste0("PC1: ",round(variance[1],1),"%"),
       y=paste0("PC2: ",round(variance[2],1),"%")) +
  theme(legend.position="top")

biplot(prin_comp, scale = 0, xlabs=rep(".", dim(prin_comp$x)[1]), main="PCA of 1100 top differentially expressed; log then z-score")



##Scree plot http://www.sthda.com/english/wiki/print.php?id=207
#The importance of princpal components (PCs) can be visualized with a scree plot.

#Scree plot using base graphics :
  
barplot(eig.cumvar[, 2], names.arg=1:nrow(eig.cumvar),
        ylim = c(0,40),
        main = "Variances",
        xlab = "Principal Components",
        ylab = "Percentage of variances",
        col ="steelblue")
# Add connected line segments to the plot
lines(x = 1:nrow(eig.cumvar), 
      eig.cumvar[, 2], 
      type="b", pch=19, col = "red")

#Try heatmap
heatmap.2(as.matrix(logExp_cor_ex3),key.xlab = 'z-score', trace='none')
heatmap.2(as.matrix(logExp_cor_ex_noERCC_4k_zsc1[1:1000,1:150]), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=FALSE, key.xlab = 'adjusted_p-value',lwid = c(.01,.99),lhei = c(.01,.99),margins = c(7,15 ))
dev.off()

####Unloged analysis
#Let's order rows with standard deviation from top to lowest
exp_cor_ex2 <- exp_cor_ex[order(exp_cor_ex$covar, decreasing = T),]
exp_cor_ex3 <- exp_cor_ex2[1:100,1:4]

##try PCA
prin_comp <- prcomp(exp_cor_ex3, scale. = T)
biplot(prin_comp, scale = 0, xlabs=rep(".", dim(prin_comp$x)[1]), main="PCA of 100 top genes with highest standard deviation")

#Try heatmap
heatmap.2(as.matrix(exp_cor_ex3))


###Trying to print MoY expression
#female
moy_female_exp = read.table("delete3", header = F, row.names = 1)
moy_female_exp = as.data.frame(t(moy_female_exp))

time.points <- moy_female_exp$


ggplot(moy_female_exp, aes(x=ids, y=coverage, fill=genes.ercc)) + 
  geom_boxplot()
