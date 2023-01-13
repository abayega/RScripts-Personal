setwd('/home/abayega/R/tests/test34')
raw_countNON <- read.table("zygotic_early_genes_coordinates_edited.NONzygotic",header=T, sep = '\t')
raw_countZygo <- read.table("zygotic_early_genes_coordinates_edited.zygotic",header=T, sep = '\t')
head(raw_countNON)


genelen <- matrix(raw_countNON$geneLength, raw_countZygo$geneLength)
No.exons <- matrix(raw_countNON$No.of.exons, raw_countZygo$No.of.exons)
intronLen <- matrix(raw_countNON$totalIntronLength, raw_countZygo$totalIntronLength)
exonLen <- matrix(raw_countNON$totalExonLength, raw_countZygo$totalExonLength)

#Try ggplot violin plots
library(ggplot2)
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# Gene length

expn <- c(raw_countNON$geneLength, raw_countZygo$geneLength)
fac <- c(rep('Non-zygotic',length(raw_countNON$geneLength)), rep('zygotic', length(raw_countZygo$geneLength)))
expn_mat <- data.frame(expn, fac)

colnames(expn_mat) <- c('geneLength', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Gene type') +
  labs(y = "Gene length (bases)") + 
  ylim(0,40000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p
wilcox.test(raw_countNON$geneLength, raw_countZygo$geneLength,paired = FALSE,conf.level = 0.99)

# No. of exons

expn <- c(raw_countNON$No.of.exons, raw_countZygo$No.of.exons)
fac <- c(rep('Non-zygotic',length(raw_countNON$No.of.exons)), rep('zygotic', length(raw_countZygo$No.of.exons)))
expn_mat <- data.frame(expn, fac)

colnames(expn_mat) <- c('No.of.exons', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Gene type') +
  labs(y = "Number of exons") + 
  ylim(0,20) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

# No. of introns
raw_countNON$No.of.introns <- raw_countNON$No.of.exons - 1
raw_countZygo$No.of.introns <- raw_countZygo$No.of.exons - 1
head(raw_countNON)

expn <- c(raw_countNON$No.of.introns, raw_countZygo$No.of.introns)
fac <- c(rep('Non-zygotic',length(raw_countNON$No.of.introns)), rep('zygotic', length(raw_countZygo$No.of.introns)))
expn_mat <- data.frame(expn, fac)

colnames(expn_mat) <- c('No.of.introns', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Gene type') +
  labs(y = "Number of introns") + 
  ylim(0,20) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

# totalIntronLength

expn <- c(raw_countNON$totalIntronLength, raw_countZygo$totalIntronLength)
fac <- c(rep('Non-zygotic',length(raw_countNON$totalIntronLength)), rep('zygotic', length(raw_countZygo$totalIntronLength)))
expn_mat <- data.frame(expn, fac)

colnames(expn_mat) <- c('totalIntronLength', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Gene type') +
  labs(y = "total Intron Length (bases)") + 
  ylim(0,2000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

# totalExonLength

expn <- c(raw_countNON$totalExonLength, raw_countZygo$totalExonLength)
fac <- c(rep('Non-zygotic',length(raw_countNON$totalExonLength)), rep('zygotic', length(raw_countZygo$totalExonLength)))
expn_mat <- data.frame(expn, fac)

colnames(expn_mat) <- c('totalExonLength', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Gene type') +
  labs(y = "total Exon Length (bases)") + 
  ylim(0,2000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p
wilcox.test(raw_countNON$totalExonLength[which(raw_countNON$totalExonLength <= 2000)], raw_countZygo$totalExonLength[which(raw_countZygo$totalExonLength<= 2000)],paired = FALSE,conf.level = 0.99)
