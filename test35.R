#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test35")
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
zygotic <- read.table("zygotic")
head(zygotic) 
zygotic <- zygotic$V2

zygotic_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% zygotic)

filetowrite <- zygotic_zscores[,9:14]
write.table(filetowrite, file = "zygotic_zscores.csv", sep = ',', quote = F)

#Now get maternal genes
maternal <- read.table("maternal")
head(maternal) 
maternal <- maternal$V2

maternal_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% maternal)

filetowrite <- maternal_zscores[,9:14]
write.table(filetowrite, file = "maternal_zscores.csv", sep = ',', quote = F)

#Now get maternal-maternal genes
maternal <- read.table("maternal-maternal")
head(maternal) 
maternal <- maternal$V2

maternal_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% maternal)

filetowrite <- maternal_zscores[,9:14]
write.table(filetowrite, file = "maternal-maternal_zscores.csv", sep = ',', quote = F)

#Now get persistent genes
persistent <- read.table("persistent")
head(persistent) 
persistent <- persistent$V2

persistent_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% persistent)

filetowrite <- persistent_zscores[,9:14]
write.table(filetowrite, file = "persistent_zscores.csv", sep = ',', quote = F)

#Now get spiking genes
spiking <- read.table("spiking")
head(spiking) 
spiking <- spiking$V2

spiking_zscores <- subset(exp_cor_ex_2, rownames(exp_cor_ex_2) %in% spiking)

filetowrite <- spiking_zscores[,9:14]
write.table(filetowrite, file = "spiking_zscores.csv", sep = ',', quote = F)



#Now combine zygotic, maternal, persistent
maternaZygoticPersistentSpiking_zs <- rbind(zygotic_zscores,maternal_zscores,persistent_zscores,spiking_zscores)

heatmap.2(as.matrix(maternaZygotic_zs[,9:14]), trace='none', Colv=FALSE, dendrogram='none') #,key.xlab = 'z-score', trace='none')
filetowrite <- maternaZygoticPersistent_zs[,9:14]
write.table(filetowrite, file = "maternaZygoticPersistentSpiking_zs.csv", sep = ',', quote = F)
















#let us try to log the values
exp_cor_ex <- log(exp_cor_ex+1)

#let us try to find most variable genes via standard deviation and coefficient of variation
exp_cor_ex$stdev <- apply(exp_cor_ex, 1, sd )
exp_cor_ex$mean <- apply(exp_cor_ex, 1, mean )
exp_cor_ex$covar <- exp_cor_ex$stdev/exp_cor_ex$mean

#try to get zscores
logExp_cor_ex <- exp_cor_ex
logExp_cor_ex$egg_z <- (exp_cor_ex$egg-exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$larvae_z <- (exp_cor_ex$larvae -exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$pupae_z <- (exp_cor_ex$pupae-exp_cor_ex$mean)/exp_cor_ex$stdev
logExp_cor_ex$adult_z <- (exp_cor_ex$adult-exp_cor_ex$mean)/exp_cor_ex$stdev

#now, based on the z-score try to get the most differentially expressed genes via standard deviation and coefficient of variation
logExp_cor_ex <- logExp_cor_ex[,8:ncol(logExp_cor_ex)]
logExp_cor_ex$stdev <- apply(logExp_cor_ex, 1, sd )
logExp_cor_ex$mean <- apply(logExp_cor_ex, 1, mean )
logExp_cor_ex$covar <- logExp_cor_ex$stdev/logExp_cor_ex$mean

#Let's order rows with coefficient of variation from top to lowest
logExp_cor_ex2 <- logExp_cor_ex[order(logExp_cor_ex$covar, decreasing = T),]
logExp_cor_ex3 <- logExp_cor_ex2[1:1100,1:4]
colnames(logExp_cor_ex3) <- c("egg","larvae","pupae","adult")

dim(logExp_cor_ex2)
head(logExp_cor_ex2)

##try PCA
#set.seed(42)
#pcout <- princomp(matrix(rnorm(1000), 100, 10))
#biplot(pcout)
#biplot(pcout, xlabs=rep(".", dim(pcout$scores)[1])) # small symbol
#biplot(pcout, xlabs=rep("", dim(pcout$scores)[1])) # no symbol 

prin_comp <- prcomp(logExp_cor_ex3, scale. = T)
biplot(prin_comp, scale = 0, xlabs=rep(".", dim(prin_comp$x)[1]), main="PCA of 1100 top differentially expressed; log then z-score")

##get some details about the results
summary(prin_comp)
names(prin_comp)
# Eigenvalues
eig <- (prin_comp$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.cumvar <- data.frame(eig = eig, variance = variance,
              cumvariance = cumvar)

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
heatmap.2(as.matrix(logExp_cor_ex3), Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
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