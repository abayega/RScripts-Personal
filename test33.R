
rm(list = ls())
setwd("/home/abayega/R/tests/test33")
library('gplots') #has the heatmap.2 function
library('ggplot2')
library(plyr)
library(plotly)
Rcount1 <- read.csv("10X_sampleing_stats_more_extra.csv",header=T)
head(Rcount1)

Rcount1$partitions <- as.factor(Rcount1$NBARCODES)
levels(Rcount1$partitions)
Rcount1$NG50 <- Rcount1$GN50.length/1000000

Rcount14k <- subset(Rcount1,NBARCODES==1400000)

my.theme <- theme(axis.text = element_text(colour="black", size=16),
                  text = element_text(size=12),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

ggplot(Rcount1, aes(x=Coverage, y=NG50, colour=partitions)) + 
  geom_line() + 
  geom_point() +
  ylab("NG50 (Mega bases)") +
  ylim(c(0.5,4)) +
  my.theme

ggplot(Rcount1, aes(x=Reads.per.partition, y=NG50, colour=partitions)) + 
  geom_line() + 
  geom_point() +
  xlab("Reads per partition") +
  ylab("NG50 (Mega bases)") +
  ylim(c(0.5,4)) +
  my.theme

ggplot(Rcount1, aes(x=Coverage, y=GN50.number, colour=partitions)) + 
  geom_line() + 
  geom_point() +
  xlab("Coverage") +
  ylab("LG50") +
  my.theme

###Now blast2GO
Rcount1 <- read.table("blast2GO.txt",header=F, sep = '\t')
head(Rcount1)
Rcount1$GO.Type <-as.factor(Rcount1$V4)
ggplot(data=Rcount1, aes(colour=GO.Type)) +
  xName="GO Type" +
  faceting=TRUE, facetingVarNames="GO.Type", 
                facetingDirection="horizontal") 


ggplot(Rcount1, aes(x=Coverage, y=GN50.number, colour=partitions)) + 
  geom_line() + 
  geom_point() +
  xlab("Coverage") +
  ylab("LG50") +
  my.theme

install.packages("devtools")
