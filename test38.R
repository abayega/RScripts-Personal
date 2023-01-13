rm(list = ls())
setwd("~/R/tests/test38")

library("ggplot2")
library('gplots') #has the heatmap.2 function

boleae <- read.table("boleae_sex_genes", header = T, row.names = 1, sep = "\t")
bdorsalis <- read.table("bdorsalis_sex_genes", header = T, row.names = 1, sep = "\t")
dmel <- read.table("dmel_sex_genes", header = T, row.names = 1, sep = "\t")

head(boleae)
colnames(boleae) <- c("1H",	"2H",	"3H", "4H",	"5H",	"6H")
colnames(bdorsalis) <- c("0-1h",	"2-4h",	"5-8h")
colnames(dmel) <- c("0-2h",	"2-4h",	"4-6h", "6-8h", "8-10h", "10-12h", "12-14h", "14-16h", "16-18h", "18-20h")

boleae <- as.data.frame(t(boleae))
bdorsalis <- as.data.frame(t(bdorsalis))
dmel <- as.data.frame(t(dmel))

boleae.colnames <- colnames(boleae)
bdorsalis.colnames <- colnames(bdorsalis)
dmel.colnames <- colnames(dmel)

plot(boleae[,1],type = "o", col = "red", xlab = "timepoint", ylab = "Transcripts per embryo", main = c(expression(italic("B oleae")),colnames(boleae)[1]), xaxt = "n")
axis(1, at=1:6, labels=rownames(boleae))
plot(bdorsalis[,1],type = "o", col = "red", xlab = "", ylab = "FPKM", main = colnames(bdorsalis)[1], xaxt = "n")
axis(1, at=1:3, labels=rownames(bdorsalis), las=2)
plot(dmel[,1],type = "o", col = "red", xlab = "timepoint", ylab = "FPKM", main = colnames(dmel)[1], xaxt = "n")
axis(1, at=1:10, labels=rownames(dmel),srt = 45)


par(mfrow=c(1,3), omi = rep(0,4), mar = c(4,2,2,0.3), mgp=c(1,0.5,0), cex.lab=1.2)

pdf("BoleaeDmel.pdf")
#par(mfrow=c(1,3)) #Not necessary for the pdf
for (i in (1:length(colnames(boleae)))){
  #jpeg(paste(colnames(boleae)[i],"jpg", sep = '.'))
  par(mfrow=c(1,2))
  plot(boleae[,i],type = "o", col = "red", xlab = "timepoint", ylab = "Transcripts per embryo", main = paste("B oleae",colnames(boleae)[i], sep = ', '), xaxt = "n")
  axis(1, at=1:6, labels=rownames(boleae))
  #plot(bdorsalis[,i],type = "o", col = "blue", xlab = "timepoint", ylab = "FPKM", main = paste("B dorsalis",colnames(bdorsalis)[i], sep = ', '), xaxt = "n")
  #axis(1, at=1:3, labels=rownames(bdorsalis))
  plot(dmel[,i],type = "o", col = "green", xlab = "", ylab = "FPKM", main = paste("D mel",colnames(dmel)[i], sep = ', '), xaxt = "n")
  axis(1, at=1:10, labels=rownames(dmel),las = 2)
  #dev.off()
}
dev.off()

###
#Trying to plot profile of Dmel genes in interesting Boleae clusters
cluster_orthos <- read.table("b_oleae_expression_over_1500_out_optimal_clustering.txt.dmelOrthos", sep = '\t', header = F, quote = "\"") #
clusters <- cluster_orthos$V1
orthos <- cluster_orthos$V3
dmel_embryo_fpkm <- read.table("dmel_embryo_FPKM", sep = '\t', header = T, quote = "\"", row.names = 1)
colnames(dmel_embryo_fpkm) <- c("Chromosome","Start","Stop","0-2hr","2-4hr","4-6hr","6-8hr","8-10hr",
                                "10-12hr","12-14hr","14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")

ptn_name_C18 <- subset(cluster_orthos, cluster_orthos$V1 == "18")
ptn_name_C18 <- ptn_name_C18$V3
head(ptn_name_C18)

dmel_embryo_fpkm_C18 <- subset(dmel_embryo_fpkm, rownames(dmel_embryo_fpkm) %in% ptn_name_C18)
dmel_embryo_fpkm_C18 <- dmel_embryo_fpkm_C18[,4:ncol(dmel_embryo_fpkm_C18)]

##remove genes with no expression
dmel_embryo_fpkm_C18 <- dmel_embryo_fpkm_C18[rowSums(dmel_embryo_fpkm_C18) != 0,]

#Now let us try to get z-scores for the genes (x-mean)/SD
dmel_embryo_fpkm_C18$mean <-apply(dmel_embryo_fpkm_C18,1,mean)

#standard deviation
dmel_embryo_fpkm_C18$StdDev <-apply(dmel_embryo_fpkm_C18,1,sd)

#Now get zscores
for (i in colnames(dmel_embryo_fpkm_C18[,1:12])){
  print(i)
}
for (i in colnames(dmel_embryo_fpkm_C18[,1:12])){
  colNamei <- paste(i,"zscore",sep = '.')
  dmel_embryo_fpkm_C18$colNamei <- ((dmel_embryo_fpkm_C18$i - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
}
head(dmel_embryo_fpkm_C18)

dmel_embryo_fpkm_C18$"0-2hr.zscore" <- ((dmel_embryo_fpkm_C18$`0-2hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"2-4hr.zscore" <- ((dmel_embryo_fpkm_C18$`2-4hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"4-6hr.zscore" <- ((dmel_embryo_fpkm_C18$`4-6hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"6-8hr.zscore" <- ((dmel_embryo_fpkm_C18$`6-8hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"8-10hr.zscore" <- ((dmel_embryo_fpkm_C18$`8-10hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"10-12hr.zscore" <- ((dmel_embryo_fpkm_C18$`10-12hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"12-14hr.zscore" <- ((dmel_embryo_fpkm_C18$`12-14hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"14-16hr.zscore" <- ((dmel_embryo_fpkm_C18$`14-16hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"16-18hr.zscore" <- ((dmel_embryo_fpkm_C18$`16-18hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"18-20hr.zscore" <- ((dmel_embryo_fpkm_C18$`18-20hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"20-22hr.zscore" <- ((dmel_embryo_fpkm_C18$`20-22hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)
dmel_embryo_fpkm_C18$"22-24hr.zscore" <- ((dmel_embryo_fpkm_C18$`22-24hr` - dmel_embryo_fpkm_C18$mean)/dmel_embryo_fpkm_C18$StdDev)

dmel_embryo_fpkm_C18.new <- dmel_embryo_fpkm_C18[,15:ncol(dmel_embryo_fpkm_C18)]
colnames(dmel_embryo_fpkm_C18.new) <- c("0-2hr","2-4hr","4-6hr","6-8hr","8-10hr","10-12hr","12-14hr",
                                        "14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")
write.table(dmel_embryo_fpkm_C18.new, file = "dmel_embryo_fpkm_C18.zscores", quote = F, sep = ',')
heatmap.2(as.matrix(dmel_embryo_fpkm_C18.new), trace='none')
##Python3
#plt.show(sns.clustermap(pd.read_csv(dr+'/'+'dmel_embryo_fpkm_C18.zscores'), cmap="coolwarm"))

##Cluster 24
ptn_name_C24 <- subset(cluster_orthos, cluster_orthos$V1 == "24")
ptn_name_C24 <- ptn_name_C24$V3
head(ptn_name_C24)

dmel_embryo_fpkm_C24 <- subset(dmel_embryo_fpkm, rownames(dmel_embryo_fpkm) %in% ptn_name_C24)
dmel_embryo_fpkm_C24 <- dmel_embryo_fpkm_C24[,4:ncol(dmel_embryo_fpkm_C24)]

##remove genes with no expression
dmel_embryo_fpkm_C24 <- dmel_embryo_fpkm_C24[rowSums(dmel_embryo_fpkm_C24) != 0,]

#Now let us try to get z-scores for the genes (x-mean)/SD
dmel_embryo_fpkm_C24$mean <-apply(dmel_embryo_fpkm_C24,1,mean)

#standard deviation
dmel_embryo_fpkm_C24$StdDev <-apply(dmel_embryo_fpkm_C24,1,sd)

#Now get zscores
dmel_embryo_fpkm_C24$"0-2hr.zscore" <- ((dmel_embryo_fpkm_C24$`0-2hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"2-4hr.zscore" <- ((dmel_embryo_fpkm_C24$`2-4hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"4-6hr.zscore" <- ((dmel_embryo_fpkm_C24$`4-6hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"6-8hr.zscore" <- ((dmel_embryo_fpkm_C24$`6-8hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"8-10hr.zscore" <- ((dmel_embryo_fpkm_C24$`8-10hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"10-12hr.zscore" <- ((dmel_embryo_fpkm_C24$`10-12hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"12-14hr.zscore" <- ((dmel_embryo_fpkm_C24$`12-14hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"14-16hr.zscore" <- ((dmel_embryo_fpkm_C24$`14-16hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"16-18hr.zscore" <- ((dmel_embryo_fpkm_C24$`16-18hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"18-20hr.zscore" <- ((dmel_embryo_fpkm_C24$`18-20hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"20-22hr.zscore" <- ((dmel_embryo_fpkm_C24$`20-22hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)
dmel_embryo_fpkm_C24$"22-24hr.zscore" <- ((dmel_embryo_fpkm_C24$`22-24hr` - dmel_embryo_fpkm_C24$mean)/dmel_embryo_fpkm_C24$StdDev)

dmel_embryo_fpkm_C24.new <- dmel_embryo_fpkm_C24[,15:ncol(dmel_embryo_fpkm_C24)]
colnames(dmel_embryo_fpkm_C24.new) <- c("0-2hr","2-4hr","4-6hr","6-8hr","8-10hr","10-12hr","12-14hr",
                                        "14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")
write.table(dmel_embryo_fpkm_C24.new, file = "dmel_embryo_fpkm_C24.zscores", quote = F, sep = ',')
heatmap.2(as.matrix(dmel_embryo_fpkm_C24.new), trace='none')


##Cluster 31
ptn_name_C31 <- subset(cluster_orthos, cluster_orthos$V1 == "31")
ptn_name_C31 <- ptn_name_C31$V3

dmel_embryo_fpkm_C31 <- subset(dmel_embryo_fpkm, rownames(dmel_embryo_fpkm) %in% ptn_name_C31)
dmel_embryo_fpkm_C31 <- dmel_embryo_fpkm_C31[,4:ncol(dmel_embryo_fpkm_C31)]

##remove genes with no expression
dmel_embryo_fpkm_C31 <- dmel_embryo_fpkm_C31[rowSums(dmel_embryo_fpkm_C31) != 0,]

#Now let us try to get z-scores for the genes (x-mean)/SD
dmel_embryo_fpkm_C31$mean <-apply(dmel_embryo_fpkm_C31,1,mean)

#standard deviation
dmel_embryo_fpkm_C31$StdDev <-apply(dmel_embryo_fpkm_C31,1,sd)

#Now get zscores
dmel_embryo_fpkm_C31$"0-2hr.zscore" <- ((dmel_embryo_fpkm_C31$`0-2hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"2-4hr.zscore" <- ((dmel_embryo_fpkm_C31$`2-4hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"4-6hr.zscore" <- ((dmel_embryo_fpkm_C31$`4-6hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"6-8hr.zscore" <- ((dmel_embryo_fpkm_C31$`6-8hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"8-10hr.zscore" <- ((dmel_embryo_fpkm_C31$`8-10hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"10-12hr.zscore" <- ((dmel_embryo_fpkm_C31$`10-12hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"12-14hr.zscore" <- ((dmel_embryo_fpkm_C31$`12-14hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"14-16hr.zscore" <- ((dmel_embryo_fpkm_C31$`14-16hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"16-18hr.zscore" <- ((dmel_embryo_fpkm_C31$`16-18hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"18-20hr.zscore" <- ((dmel_embryo_fpkm_C31$`18-20hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"20-22hr.zscore" <- ((dmel_embryo_fpkm_C31$`20-22hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)
dmel_embryo_fpkm_C31$"22-24hr.zscore" <- ((dmel_embryo_fpkm_C31$`22-24hr` - dmel_embryo_fpkm_C31$mean)/dmel_embryo_fpkm_C31$StdDev)

dmel_embryo_fpkm_C31.new <- dmel_embryo_fpkm_C31[,15:ncol(dmel_embryo_fpkm_C31)]
colnames(dmel_embryo_fpkm_C31.new) <- c("0-2hr","2-4hr","4-6hr","6-8hr","8-10hr","10-12hr","12-14hr",
                                        "14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")
write.table(dmel_embryo_fpkm_C31.new, file = "dmel_embryo_fpkm_C31.zscores", quote = F, sep = ',')
heatmap.2(as.matrix(dmel_embryo_fpkm_C31.new), trace='none')


##
##Cluster 1
ptn_name_C1 <- subset(cluster_orthos, cluster_orthos$V1 == "1")
ptn_name_C1 <- ptn_name_C1$V3

dmel_embryo_fpkm_C1 <- subset(dmel_embryo_fpkm, rownames(dmel_embryo_fpkm) %in% ptn_name_C1)
dmel_embryo_fpkm_C1 <- dmel_embryo_fpkm_C1[,4:ncol(dmel_embryo_fpkm_C1)]

##remove genes with no expression
dmel_embryo_fpkm_C1 <- dmel_embryo_fpkm_C1[rowSums(dmel_embryo_fpkm_C1) != 0,]

#Now let us try to get z-scores for the genes (x-mean)/SD
dmel_embryo_fpkm_C1$mean <-apply(dmel_embryo_fpkm_C1,1,mean)

#standard deviation
dmel_embryo_fpkm_C1$StdDev <-apply(dmel_embryo_fpkm_C1,1,sd)

#Now get zscores
dmel_embryo_fpkm_C1$"0-2hr.zscore" <- ((dmel_embryo_fpkm_C1$`0-2hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"2-4hr.zscore" <- ((dmel_embryo_fpkm_C1$`2-4hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"4-6hr.zscore" <- ((dmel_embryo_fpkm_C1$`4-6hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"6-8hr.zscore" <- ((dmel_embryo_fpkm_C1$`6-8hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"8-10hr.zscore" <- ((dmel_embryo_fpkm_C1$`8-10hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"10-12hr.zscore" <- ((dmel_embryo_fpkm_C1$`10-12hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"12-14hr.zscore" <- ((dmel_embryo_fpkm_C1$`12-14hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"14-16hr.zscore" <- ((dmel_embryo_fpkm_C1$`14-16hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"16-18hr.zscore" <- ((dmel_embryo_fpkm_C1$`16-18hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"18-20hr.zscore" <- ((dmel_embryo_fpkm_C1$`18-20hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"20-22hr.zscore" <- ((dmel_embryo_fpkm_C1$`20-22hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"22-24hr.zscore" <- ((dmel_embryo_fpkm_C1$`22-24hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)

dmel_embryo_fpkm_C1.new <- dmel_embryo_fpkm_C1[,15:ncol(dmel_embryo_fpkm_C1)]
colnames(dmel_embryo_fpkm_C1.new) <- c("0-2hr","2-4hr","4-6hr","6-8hr","8-10hr","10-12hr","12-14hr",
                                        "14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")
write.table(dmel_embryo_fpkm_C1.new, file = "dmel_embryo_fpkm_C1.zscores", quote = F, sep = ',')
heatmap.2(as.matrix(dmel_embryo_fpkm_C1.new), trace='none')

##
##Cluster 49, 84 = maternal
ptn_name_Matn1 <- subset(cluster_orthos, cluster_orthos$V1 == "49") # || cluster_orthos$V1 == "84")
ptn_name_Matn2 <- subset(cluster_orthos, cluster_orthos$V1 == "84") # || cluster_orthos$V1 == "84")
ptn_name_Matn <- ptn_name_Matn1$V3
ptn_name_Matn <- c(ptn_name_Matn1$V3, ptn_name_Matn2$V3)

dmel_embryo_fpkm_C1 <- subset(dmel_embryo_fpkm, rownames(dmel_embryo_fpkm) %in% ptn_name_Matn)
dmel_embryo_fpkm_C1 <- dmel_embryo_fpkm_C1[,4:ncol(dmel_embryo_fpkm_C1)]

##remove genes with no expression
dmel_embryo_fpkm_C1 <- dmel_embryo_fpkm_C1[rowSums(dmel_embryo_fpkm_C1) != 0,]

#Now let us try to get z-scores for the genes (x-mean)/SD
dmel_embryo_fpkm_C1$mean <-apply(dmel_embryo_fpkm_C1,1,mean)

#standard deviation
dmel_embryo_fpkm_C1$StdDev <-apply(dmel_embryo_fpkm_C1,1,sd)

#Now get zscores
dmel_embryo_fpkm_C1$"0-2hr.zscore" <- ((dmel_embryo_fpkm_C1$`0-2hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"2-4hr.zscore" <- ((dmel_embryo_fpkm_C1$`2-4hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"4-6hr.zscore" <- ((dmel_embryo_fpkm_C1$`4-6hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"6-8hr.zscore" <- ((dmel_embryo_fpkm_C1$`6-8hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"8-10hr.zscore" <- ((dmel_embryo_fpkm_C1$`8-10hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"10-12hr.zscore" <- ((dmel_embryo_fpkm_C1$`10-12hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"12-14hr.zscore" <- ((dmel_embryo_fpkm_C1$`12-14hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"14-16hr.zscore" <- ((dmel_embryo_fpkm_C1$`14-16hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"16-18hr.zscore" <- ((dmel_embryo_fpkm_C1$`16-18hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"18-20hr.zscore" <- ((dmel_embryo_fpkm_C1$`18-20hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"20-22hr.zscore" <- ((dmel_embryo_fpkm_C1$`20-22hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)
dmel_embryo_fpkm_C1$"22-24hr.zscore" <- ((dmel_embryo_fpkm_C1$`22-24hr` - dmel_embryo_fpkm_C1$mean)/dmel_embryo_fpkm_C1$StdDev)

dmel_embryo_fpkm_C1.new <- dmel_embryo_fpkm_C1[,15:ncol(dmel_embryo_fpkm_C1)]
colnames(dmel_embryo_fpkm_C1.new) <- c("0-2hr","2-4hr","4-6hr","6-8hr","8-10hr","10-12hr","12-14hr",
                                       "14-16hr","16-18hr","18-20hr","20-22hr","22-24hr")
write.table(dmel_embryo_fpkm_C1.new, file = "dmel_embryo_fpkm_Matn.zscores", quote = F, sep = ',')
heatmap.2(as.matrix(dmel_embryo_fpkm_C1.new), trace='none')
