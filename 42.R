rm(list = ls())
setwd('/home/abayega/R/tests/test42')
#These files only contain a single genes and the transcript reported is that with the longest isoform
dmel_geneInfo <- read.table("dmel-all-r6.32.gtf.gene_mrna_info_singleGene",header=F,sep = '\t')
boleae_geneInfo <- read.table("boleae_singleGene_mrna_info",header=F, sep = '\t')
dmel_zygotic_pure <- read.table("zygotic_pure_new_ids.de_renzis", header = T, sep = '\t')
boleae_zygotic_pure <- read.table("Gfold_zygotic_genes_boleae", header = F, sep = '\t')
dmel_early_zygotic <- read.table("zygotic_early_new_ids",header = T, sep = '\t')

#dmel_geneInfo <- data.frame(dmel_geneInfo, stringsAsFactors=FALSE)
#boleae_geneInfo <- data.frame(boleae_geneInfo,stringsAsFactors = F)
#dmel_zygotic_pure <- data.frame(dmel_zygotic_pure, stringsAsFactors=FALSE)
#boleae_zygotic_pure <- data.frame(boleae_zygotic_pure, stringsAsFactors=FALSE)
head(dmel_early_zygotic)

colnames(dmel_zygotic_pure) <- c("submitted_ID","current_ID","converted_ID","related_record")
colnames(dmel_geneInfo) <- c("transcript_name","gene_name","protein","gene_length","transcript_name","longest_transcript_length", "strand")
colnames(boleae_geneInfo) <- c("transcript_name","gene_name","gene_ID","gene_length","longest_transcript_length", "strand")
colnames(boleae_zygotic_pure) <- c("gene_name")


#Creating dmel zygotic dataframe
#make sure the 2 vectors u are creating look exactly the same. One of them has spaces around which u can remove using gsub
dmel_zyg <- subset(dmel_geneInfo, as.character(dmel_geneInfo$gene_name) %in% gsub(" ","",as.character(dmel_zygotic_pure$Current.ID)))

#Now work on boleae dataframe
boleae_zyg <- subset(boleae_geneInfo, as.character(boleae_geneInfo$gene_name) %in% gsub(" ","",as.character(boleae_zygotic_pure$gene_name)))

#Now work on dmel early zygotic dataframe
dmel_early_zyg <- subset(dmel_geneInfo, as.character(dmel_geneInfo$gene_name) %in% gsub(" ","",as.character(dmel_early_zygotic$Current.ID)))

#Try ggplot violin plots
library(ggplot2)
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# Gene length
#I will only focus on D Renzis zygotic early category and leave dmel_zyg$gene_length, 
geneLen <- c(boleae_zyg$gene_length, dmel_early_zyg$gene_length)
Species <- c(rep('B.oleae', length(boleae_zyg$gene_length)), rep('D.melanogaster',length(dmel_early_zyg$gene_length)))
finalTable <- data.frame(geneLen, Species)

colnames(finalTable) <- c('Gene_length', 'species')
finalTable$species <- as.factor(finalTable$species)
head(finalTable)

p <- ggplot(finalTable, aes(x=species, y=Gene_length, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'species') +
  labs(y = "Gene length (bases)") + 
  ylim(0,5000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p #This results shows that Dmelanogaster pure zygotic genes are longer than early boleae zygotic genes


# Let's try to compare longest transcript length

geneLen <- c(dmel_zyg$longest_transcript_length, boleae_zyg$longest_transcript_length)
Species <- c(rep('D.melanogaster',length(dmel_zyg$gene_length)), rep('B.oleae', length(boleae_zyg$gene_length)))
finalTable <- data.frame(geneLen, Species)

colnames(finalTable) <- c('Transcript_length', 'species')
finalTable$species <- as.factor(finalTable$species)
head(finalTable)

p <- ggplot(finalTable, aes(x=species, y=Transcript_length, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'species') +
  labs(y = "Transcript length (bases)") + 
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


# Let's try to compare gene length for all genes

geneLen <- c(dmel_geneInfo$gene_length, boleae_geneInfo$gene_length)
Species <- c(rep('D.melanogaster',length(dmel_geneInfo$gene_length)),rep('Boleae',length(boleae_geneInfo$gene_length)))
finalTable <- data.frame(geneLen, Species)

colnames(finalTable) <- c('Gene_length', 'species')
finalTable$species <- as.factor(finalTable$species)
head(finalTable)

p <- ggplot(finalTable, aes(x=species, y=Gene_length, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'species') +
  labs(y = "Gene length (bases)") + 
  ylim(0,10000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#Let us directly compare the summaries of the early zygotic genes in Dmel and Boleae and compare that to the whole gene list
summary(dmel_early_zyg$gene_length)
summary(boleae_zyg$gene_length)

summary(dmel_geneInfo$gene_length)
summary(boleae_geneInfo$gene_length)

#So, Dmel and Boleae, surprisingly, have similar gene lengths with Dmel having higher
#Early zygotic genes are also longer in Boleae than Dmel.

#I just want to combine early zygotic and all genes to one graph
geneLen <- c(boleae_zyg$gene_length, dmel_early_zyg$gene_length,
             dmel_geneInfo$gene_length, boleae_geneInfo$gene_length)
Species <- c(rep('B.oleae, zygotic', length(boleae_zyg$gene_length)), rep('D.mel, zygotic',length(dmel_early_zyg$gene_length)),
             rep('D.mel, All genes',length(dmel_geneInfo$gene_length)),rep('B.oleae, All genes',length(boleae_geneInfo$gene_length)))
finalTable <- data.frame(geneLen, Species)

colnames(finalTable) <- c('Gene_length', 'species')
finalTable$species <- as.factor(finalTable$species)
head(finalTable)

p <- ggplot(finalTable, aes(x=species, y=Gene_length, fill=species)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'species') +
  labs(y = "Gene length (bases)") + 
  ylim(0,5000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14, face="bold"),
    axis.text.y = element_text(size=14, face="bold"),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

wilcox.test(boleae_zyg$gene_length,boleae_geneInfo$gene_length,paired = FALSE,conf.level = 0.99)
wilcox.test(dmel_early_zyg$gene_length,dmel_geneInfo$gene_length,paired = FALSE,conf.level = 0.99)
