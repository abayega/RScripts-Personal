rm(list = ls())
setwd("/home/abayega/R/tests/test23")

##Get gene set enrichment
library(gProfileR)
i = 44
a<-read.table("clusters/cluster_44")
b<-gprofiler(a,organism = "dmelanogaster",significant = F,src_filter=("GO:BP"),correction_method="fdr")
adjust_p_value<-p.adjust(b$p.value, method = c("fdr"), n = length(b$p.value))
dframe = cbind.data.frame(b$term.id,b$term.name,b$p.value,adjust_p_value)
#dframe = cbind(b$term.id,b$term.name,b$p.value,adjust_p_value)
colnames(dframe) <- c('term.id','term.name','p.value','adjust_p_value')
#dframe <- data.frame(dframe)
dframe = subset(dframe, p.value < 0.05) #   dframe[dframe$p.value < 0.05,]
dframe = dframe[order(dframe$p.value),]
write.table(dframe, file=paste("drosophila_enriched_GO_BP_terms_g_profiler_fdr0.05_cluster_",i,sep = ""), quote = F, sep = "\t")
#write.table(file="drosophila_enriched_GO_BP_terms_g_profiler_cluster-44",paste(b$term.id,b$term.name,b$p.value,adjust_p_value,sep ="\t"))

get_gsea <- function(clus_dir) {
  pb <- txtProgressBar(min = 0, max = 87, style = 3)
  for (i in 1:87){
    #print(paste("working on cluster",i))
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


##tryint DESeq
source("https://bioconductor.org/biocLite.R")
biocLite("DESeq")
install.packages("ggplot2")