rm(list = ls())
setwd("/home/abayega/R/tests/test31/get_orthos")

##Get gene set enrichment
library(gProfileR)

#This is the final one I used. I am finding the gene set enrichment among the most highly expressed genes at 3 hpo
a<-read.table("egg.dmelOrthos")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 100) #dframe #
write.table(dframe1, file="egg.dmelOrthos.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("egg.dmelOrthos.dframe1",header = T,sep='\t')
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in egg genes")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in egg genes")
rm(op)

#finding the gene set enrichment among the maternal degraded expressed genes, 2 hpo
a<-read.table("larvae.dmelOrthos")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe2 = dframe #head(dframe, 100)
write.table(dframe2, file="larvae.dmelOrthos.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("larvae.dmelOrthos.dframe1",header = T,sep='\t')
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in larvae genes")
barplot(dframe1$overlap.perc[dframe1$domain == "MF"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "MF"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched molecular function in larvae genes")
op <- par(mar=c(17,4,4,2))
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal degraded genes")
rm(op)


#finding the gene set enrichment among the zygotic genes
a<-read.table("pupae.dmelOrthos")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe3 = head(dframe, 100)
write.table(dframe3, file="pupae.dmelOrthos.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("pupae.dmelOrthos.dframe1",header = T,sep='\t')
op <- par(mar=c(22,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", cex.names=0.9,ylim = c(0,40),ylab = "percentage of genes",main = "significantly enriched biological processes in pupae genes")
barplot(dframe1$overlap.perc[dframe1$domain == "MF"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "MF"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", cex.names=0.9,ylim = c(0,40),ylab = "percentage of genes",main = "significantly enriched molecular function in pupae genes")
barplot(dframe1$overlap.perc[dframe1$domain == "tf"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", cex.names=0.9,ylim = c(0,40),ylab = "percentage of genes",main = "significantly enriched transcription factors in pupae genes")
op <- par(mar=c(17,4,4,2))
rm(op)

#finding the gene set enrichment among the maternal stable/upregulated genes
a<-read.table("adult.dmelOrthos")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe4 = head(dframe, 100)
write.table(dframe4, file="adult.dmelOrthos.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("adult.dmelOrthos.dframe1",header = T,sep='\t')
op <- par(mar=c(22,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, cex.lab=1.5,las=2, col="skyblue", cex.names=0.9,ylim = c(0,40),ylab = "percentage of genes",main = "significantly enriched biological processes in adult genes")
op <- par(mar=c(15,4,4,2))
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",ylim = c(0,30),cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal stable/upregulated genes")
rm(op)
###Useful code ends here. The stuff below was just trial

