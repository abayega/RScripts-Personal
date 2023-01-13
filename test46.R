#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test46")
library('gplots') #has the heatmap.2 function
library(ggplot2)

b4corr <- read.table("b4corr_alignment_stats",header=T, sep = '\t')
corr <- read.table("corr_alignment_stats",header=T, sep = '\t')

head(b4corr)

expn <- c(b4corr$alignment_identity_., corr$read_aligned_.)
fac <- c(rep('corrected',length(b4corr$alignment_identity_.)), rep('b4_corr', length(corr$read_aligned_.)))

expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('expn', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)

colnames(expn_mat) <- c('identity(%)', 'correction')
expn_mat$correction <- as.factor(expn_mat$correction)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac) #ggplot(expn_mat, aes(x='correction', y='identity(%)', fill='correction')) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "correction status") +
  labs(y = "alignment identity (%)") +
  labs(fill = 'key') +
  #scale_y_continuous(expression("alignment identity (%)")) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

# violin plot with median points
p + stat_summary(fun.y=median, geom="point", size=2, color="red")
p + geom_boxplot(width=0.1)
p + stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )