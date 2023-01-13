rm(list = ls())

library(RColorBrewer)
library(ggplots)
library(DESeq2)
hmcol = colorRampPalette(brewer.pal(9, "GnBu"))(100)

source("https://bioconductor.org/biocLite.R")
biocLite("DESeq2")
install.packages("ggplot2")
install.packages("Rcpp")
install.packages("plyr")
install.packages("drat",repos="http://cran.rstudio.com" )
install.packages('devtools')
install.packages('xml2')
devtools::install_github('rstudio/shiny')
install.packages("tidyverse")
install.packages("factoextra")
devtools::install_github('rstudio/shiny')
devtools::install_github("tidyverse/readr")
install.packages("readr")

library(Rcpp)
library(tidyverse)  # data manipulation
library(cluster)    # clustering algorithms
library(factoextra) # clustering algorithms & visualization
library(gridExtra)


# Set working directory
setwd("/home/abayega/R/tests/test22")


# Counts-Matrix
raw_count <- read.table("b_oleae_expression",header=T,row.names=1)
raw_count <- read.table("combined_mandalorion_absolute",header=T,row.names=1)
raw_count2 <- read.table("combined_mandalorion_absolute",header=T,row.names=1)
#raw_count <- read.table("expression_over_1500",header=T,row.names=1)
#raw_count <- raw_count[,1:6]

#countdata1 <- raw_count[93:nrow(raw_count),17:24]
countdata1 <- round(raw_count[93:nrow(raw_count),17:24])

#countdata1 <- raw_count
#write.csv(countdata1, "final_filtered_data - Copy - Copy.csv")
head(countdata1)

# remove uninformative genes keep only genes that are expressed with at least 1200 absolute counts in at least one of the samples
countdata1[countdata1 < 1200 ] <- 0
countdata1 <- countdata1[rowSums(countdata1 > 1200) >= 1,]

# Remove all gene which has 0 value in all sample
#all <- apply(countdata1, 1, function(x) all(x==0) )
#newdata <- countdata1[!all,]
#write.csv(newdata, file = "count_0_filter.csv")
#head (newdata)

# remove uninformative genes keep only genes that are expressed with at least 1500 absolute counts in at least one of the samples
#dat <- newdata[rowSums(newdata > 1200) >= 1,]
#dat[dat < 1200 ] <- 0
#write.table(dat, file = "expression_over_1500_2", quote = F, sep = "\t")


# Convert to matrix just for the 
countdata2 <- as.matrix(countdata1)
colnames(countdata2) <- c('emb_1H','emb_2H','emb_3H','emb_4H','emb_5H','emb_6H','female','male')

#heatmap
heatmap(countdata2) #, scale = "column")

#let us try kmeans clustering
clusters <- kmeans(countdata,centers = 10, nstart = 25)
str(clusters)
fviz_cluster()

#heatmap
dat <- dat[!all,]
heatmap(as.matrix(dat), Rowv=NA, Colv=NA, scale="column") #, scale = "column")

###250718
setwd("/home/abayega/R/tests/test22")
raw_count <- data.frame(read.table("diff_expression",header=T,row.names=1))
dat2 = as.data.frame(dat)
# Create new sorted data frame. This is analogous to the second data frame
# by which you want to sort the first one.
raw_count.sorted = raw_count[order(rownames(raw_count)),]

raw_count.sorted.genes <- raw_count.sorted[71:length(rownames(raw_count.sorted)),]
all <- apply(raw_count.sorted.genes, 1, function(x) all(x==0) )
raw_count.sorted.genes <- raw_count.sorted.genes[!all,]

#Find genes for which 1H > 1200 and 0<=2H and 3H >= 1200
raw_count.sorted.genes.no1h <- raw_count.sorted.genes[which(raw_count.sorted.genes$X1H_abs_emb == 0),]

#Find genes for which 1H > 1200 and 0<=2H and 3H >= 1200
raw_count.sorted.genes.no3h <- raw_count.sorted.genes[which(raw_count.sorted.genes$X3H_abs_emb == 0),]

raw_count.sorted.genes$log2_1v2 <- log2(raw_count.sorted.genes$X1H_abs_emb/raw_count.sorted.genes$X2H_abs_emb)
raw_count.sorted.genes$log2_3v2 <- log2(raw_count.sorted.genes$X3H_abs_emb/raw_count.sorted.genes$X2H_abs_emb)


#This is how I created the table I used to find De Renzis' B oleae homologous genes
raw_count <- read.table("combined_mandalorion_absolute",header=T,row.names=1)
raw_count <- raw_count[,1:8]
raw_count[raw_count < 1200] <- 0
raw_count$gene_length <- raw_count2$gene_length
raw_count$gene_gc <- raw_count2$gene_gc
raw_count$transcript_length <- raw_count2$transcript_length
raw_count$transcript_gc <- raw_count2$transcript_gc
raw_count$gene_ID <- raw_count2$gene_ID
raw_count$transcript_name <- raw_count2$transcript_name
write.table(raw_count, file="combined_mandalorion_absolute_ed", quote = F, sep = "\t")