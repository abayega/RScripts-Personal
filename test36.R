#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test36")
library('gplots') #has the heatmap.2 function
library(ggplot2)

exp_cor_ex <- read.table("Bo_EnHeads_tofu_sqanti_classification1.txt",header=T, row.names = 1, sep = '\t')

#print number of exons in novel genes compared to others
exons_novelty <- subset(exp_cor_ex,select = c(exons, native))
head(exons_novelty)
exons_novelty$fac <- as.factor(exons_novelty$native)

p <- ggplot(exons_novelty, aes(x=native, y=exons, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Transcript type') +
  labs(y = "Number of exons") + 
  #labs(y = "Number of exons") + 
  ylim(0,10) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#find percentage of novel genes with only 1 exon
100*(length(which(exp_cor_ex$native == "novel" & exp_cor_ex$exons==1))/length(which(exp_cor_ex$native == "novel")))

#find percentage of annotated genes with only 1 exon
100*(length(which(exp_cor_ex$native != "novel" & exp_cor_ex$exons==1))/length(which(exp_cor_ex$native != "novel")))
