rm(list = ls())
setwd("/home/abayega/R/tests/test35/gsea")

##Get GO enrichment
library(gProfileR)

#enriched GO terms in maternal0
a<-read.table("maternal_zscores.csv.dmel_homologue.0")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="maternal0.GO", quote = F, sep = '\t')
dframe1 <- read.table("maternal0.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal0 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal0 genes")
rm(op)

#enriched GO terms in maternal1
a<-read.table("maternal_zscores.csv.dmel_homologue.1")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="maternal1.GO", quote = F, sep = '\t')
dframe1 <- read.table("maternal1.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal1 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal1 genes")
rm(op)

#enriched GO terms in maternal2
a<-read.table("maternal_zscores.csv.dmel_homologue.2")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="maternal2.GO", quote = F, sep = '\t')
dframe1 <- read.table("maternal2.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal2 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal2 genes")
rm(op)

#enriched GO terms in maternal3
a<-read.table("maternal_zscores.csv.dmel_homologue.3")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="maternal3.GO", quote = F, sep = '\t')
dframe1 <- read.table("maternal3.GO",header = T,sep='\t')
dframe1 = head(dframe1, 50)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal3 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal3 genes")
rm(op)

#enriched GO terms in maternal-maternal
a<-read.table("maternal-maternal.dmelOrthos")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="maternal-maternal.GO", quote = F, sep = '\t')
dframe1 <- read.table("maternal-maternal.GO",header = T,sep='\t')
dframe1 = head(dframe1, 50)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal-maternal genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal3 genes")
rm(op)



##Now zygotes
#enriched GO terms in zygote0
a<-read.table("zygotic_zscores.csv.dmel_homologue.0")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="zygote0.GO", quote = F, sep = '\t')
dframe1 <- read.table("zygote0.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in zygotic genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in zygotic genes")
rm(op)



##Now persistent
#enriched GO terms in persistent0
a<-read.table("persistent_zscores.csv.dmel_homologue.0")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="persistent0.GO", quote = F, sep = '\t')
dframe1 <- read.table("persistent0.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in persistent0 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in persistent0 genes")
rm(op)

#enriched GO terms in persistent1
a<-read.table("persistent_zscores.csv.dmel_homologue.1")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="persistent1.GO", quote = F, sep = '\t')
dframe1 <- read.table("persistent1.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in persistent1 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in persistent1 genes")
rm(op)




##Now spiking/transient
#enriched GO terms in spiking0
a<-read.table("spiking_zscores.csv.dmel_homologue.0")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="spiking0.GO", quote = F, sep = '\t')
dframe1 <- read.table("spiking0.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in spiking0 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in spiking0 genes")
rm(op)

#enriched GO terms in spiking1
a<-read.table("spiking_zscores.csv.dmel_homologue.1")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="spiking1.GO", quote = F, sep = '\t')
dframe1 <- read.table("spiking1.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in spiking1 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in spiking1 genes")
rm(op)

#enriched GO terms in spiking2
a<-read.table("spiking_zscores.csv.dmel_homologue.2")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
write.table(dframe, file="spiking2.GO", quote = F, sep = '\t')
dframe1 <- read.table("spiking2.GO",header = T,sep='\t')
dframe1 = head(dframe1, 100)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in spiking2 genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in spiking2 genes")
rm(op)





