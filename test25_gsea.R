rm(list = ls())
setwd("/home/abayega/R/tests/test25/gsea")

##Get gene set enrichment
library(gProfileR)

a<-read.table("top150")
#b<-gprofiler(a,organism = "dmelanogaster",significant = F,src_filter=("GO:BP"),correction_method="fdr")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
#dframe = cbind(b$term.id,b$term.name,b$p.value,adjust_p_value)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
#dframe <- data.frame(dframe)
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
#write.table(dframe, file=paste("drosophila_enriched_GO_BP_terms_g_profiler_fdr0.05_cluster_",i,sep = ""), quote = F, sep = "\t")
write.table(dframe, file="top150_drosophila_enriched_GO_BP_terms_g_profiler", quote = F, sep = '\t')

#new_df <- subset(dframe, domain == "BP", select = c(overlap.perc,term.name))
#k=c(paste(new_df$term.name, sep = " ", collapse = NULL))
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe$overlap.perc[dframe$domain == "BP"], names.arg=c(paste(dframe$term.name[dframe$domain == "BP"],sep = " ", collapse = NULL)),horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 150 genes at 3 hpo")
rm(op)

#This is the final one I used. I am finding the gene set enrichment among the most highly expressed genes at 3 hpo
a<-read.table("top1_400")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 100)
write.table(dframe1, file="top1_400.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("top1_400.dframe1",header = T,sep='\t')
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 1_400 genes at 3 hpo")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in top 1_400 genes at 3 hpo")
rm(op)

#finding the gene set enrichment among the maternal degraded expressed genes, 2 hpo
a<-read.table("maternal_degraded")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 100)
write.table(dframe1, file="maternal_degraded.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("maternal_degraded.dframe1",header = T,sep='\t')
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in maternal degraded genes")
op <- par(mar=c(17,4,4,2))
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal degraded genes")
rm(op)


#finding the gene set enrichment among the zygotic genes
a<-read.table("zygotic_early_genes.Bo_homologue.proteins")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 100)
write.table(dframe1, file="zygotic_early_genes.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("zygotic_early_genes.dframe1",header = T,sep='\t')
op <- par(mar=c(22,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", cex.names=0.9,ylim = c(0,40),ylab = "percentage of genes",main = "significantly enriched biological processes in zygotic genes")
op <- par(mar=c(17,4,4,2))
rm(op)

#finding the gene set enrichment among the maternal stable/upregulated genes
a<-read.table("z.maternal_stable_upregulated.Bo_homologue.proteins")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 41)
write.table(dframe1, file="maternal_stable_upregulated.dframe1", quote = F, sep = '\t')
dframe1 <- read.table("maternal_stable_upregulated.dframe1",header = T,sep='\t')
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, cex.lab=1.5,las=2, col="skyblue", cex.names=1.5,ylim = c(0,60),ylab = "percentage of genes",main = "significantly enriched biological processes in maternal stable/upregulated genes")
op <- par(mar=c(15,4,4,2))
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",ylim = c(0,30),cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in maternal stable/upregulated genes")
rm(op)
###Useful code ends here. The stuff below was just trial



a<-read.table("top50")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 47)
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 50 genes at 3 hpo")
rm(op)

a<-read.table("top50_100")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 47)
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 50_100 genes at 3 hpo")
rm(op)

a<-read.table("top100_200")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 47)
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 100_200 genes at 3 hpo")
rm(op)

a<-read.table("top200_400")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 47)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 200_400 genes at 3 hpo")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched tf in top 200_400 genes at 3 hpo")
rm(op)


a<-read.table("top400_600")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
#BP  CC  hp keg  MF rea  tf 
#44  22   0   0   2   0   0 
dframe1 = head(dframe, 47)
op <- par(mar=c(18,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in next 300 genes at 3 hpo")
rm(op)
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched tf in next 300 genes at 3 hpo")

#Code used to get enrichment in DP_GP clusters
i = 44
get_gsea <- function(clus_dir) {
  pb <- txtProgressBar(min = 0, max = 87, style = 3)
  for (i in 1:3){
    print(paste("working on cluster",i))
    a<-read.table(paste(clus_dir,"/cluster_",i, sep = ""))
    b<-gprofiler(a,organism = "dmelanogaster",significant = F,src_filter=("GO:BP"),correction_method="fdr")
    adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
    dframe = cbind.data.frame(b$term.id,b$term.name,b$p.value,adjust_p_value)
    colnames(dframe) <- c('term.id','term.name','p.value','adjust_p_value')
    dframe = subset(dframe, p.value < 0.05)
    dframe = dframe[order(dframe$p.value),]
    write.table(dframe, file=paste("go_terms/drosophila_enriched_GO_BP_terms_g_profiler_fdr0.05_cluster_",i,sep = ""), row.names = F, quote = F, sep = "\t")
    Sys.sleep(0.1)
    setTxtProgressBar(pb, i)
  }
  close(pb)
}

get_gsea("clusters")

#Code used to get enrichment in combined DP_GP clusters
setwd('/home/abayega/R/tests/test23/combined_clusters')
get_gsea <- function(clus_dir) {
  pb <- txtProgressBar(min = 0, max = 3, style = 3)
  files <- list.files()
  files <- c(files[2],files[4],files[5])
  k = 1
  for (i in files){
    print(paste("working on cluster",i))
    a<-read.table(i)
    b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
    adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
    overlap.perc <- round((100*(b$overlap.size/b$query.size)),2)
    dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
    colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
    dframe = subset(dframe, p.value < 0.05)
    dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
    print(table(dframe$domain))
    dframe1 = head(dframe, 40)
    op <- par(mfrow=c(1,1),mar=c(20,4,4,2))
    barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
            horiz=F, las=2,col="skyblue", ylab = "percentage of genes",main = paste("significantly enriched biological processes in ",i,sep=""))
    #barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
            #horiz=F, las=2,col="skyblue", ylab = "percentage of genes", main = paste("significantly enriched tf in ",i,sep=""),cex.lab=1.2,cex.names=0.9)
    write.table(dframe, file=paste("go_terms/drosophila_enriched_terms_g_profiler_fdr0.05_cluster_",i,sep = ""), row.names = F, quote = F, sep = "\t")
    Sys.sleep(0.1)
    k = k + 1
    setTxtProgressBar(pb, k)
  }
  close(pb)
}
get_gsea("clusters")

a<-read.table("top1_400")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
overlap.perc <- round((100*(b$overlap.size/b$query.size)),2) 
dframe = cbind.data.frame(b$term.id,b$domain,b$term.name,b$p.value,adjust_p_value,overlap.perc)
colnames(dframe) <- c('term.id','domain','term.name','p.value','adjust_p_value','overlap.perc')
dframe = subset(dframe, adjust_p_value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$domain,dframe$adjust_p_value),]
table(dframe$domain)
dframe1 = head(dframe, 40)
op <- par(mar=c(20,4,4,2)) # 'mar' A numerical vector of the form 'c(bottom, left, top, right)' the 10 allows the names.arg below the barplot
barplot(dframe1$overlap.perc[dframe1$domain == "BP"], names.arg=c(paste(dframe1$term.name[dframe1$domain == "BP"],sep = " ", collapse = NULL)),
        horiz=F, cex = 1.2,las=2,col="skyblue", ylab = "percentage of genes",main = "significantly enriched biological processes in top 1_400 genes at 3 hpo")
barplot(dframe$overlap.perc[dframe$domain == "tf"], names.arg=c(paste(dframe$term.name[dframe$domain == "tf"],sep = " ", collapse = NULL)),
        horiz=F, las=2,col="skyblue", ylab = "percentage of genes",cex.lab=1.2,cex.names=0.7,main = "significantly enriched tf in top 1_400 genes at 3 hpo")
rm(op)

##tryint DESeq
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
install.packages("ggplot2")