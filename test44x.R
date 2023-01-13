# Set working directory
setwd("/home/banthony/R/tests/test44/repeat_seq")

read_stats <- read.table("repeat_reads_stats_2",header=T, sep = '\t')
head(read_stats)
dim(read_stats)

new_read_stats <- subset(read_stats, quality=="pass")

library("ggplot2")
require(data.table)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

#let's try to plot pass reads
species = factor(new_read_stats$sample, levels = c( 'C13H', 'C15H', 'CB1', 'CB2'),ordered = TRUE)
condition = new_read_stats$barcode
value = new_read_stats$No_reads
new_read_stats.df = data.frame(species,condition,value)

ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Assessing pass reads distribution among barcodes for capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Proportion of reads") +
  #ylim(c(0,1)) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00")) + #999999=blue,#56B4E9=grey,#E69F00=yellow not sure though
  #scale_fill_manual(values=c("#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999")) +
  my.theme

#let's try to plot fail reads
fail_read_stats <- subset(read_stats, quality=="fail")
species = factor(fail_read_stats$sample, levels = c( 'C13H', 'C15H', 'CB1', 'CB2'),ordered = TRUE)
condition = fail_read_stats$barcode
value = fail_read_stats$No_reads
fail_read_stats.df = data.frame(species,condition,value)

ggplot(fail_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Assessing fail reads distribution among barcodes for capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Proportion of reads") +
  my.theme

#Let's try to plot the total yields regardless of the barcodes
CB1p <- sum(read_stats$No_reads[which(read_stats$sample == 'CB1' & read_stats$quality=='pass')])/1000000
CB1f <- sum(read_stats$No_reads[which(read_stats$sample == 'CB1' & read_stats$quality=='fail')])/1000000

CB2p <- sum(read_stats$No_reads[which(read_stats$sample == 'CB2' & read_stats$quality=='pass')])/1000000
CB2f <- sum(read_stats$No_reads[which(read_stats$sample == 'CB2' & read_stats$quality=='fail')])/1000000

C13Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C13H' & read_stats$quality=='pass')])/1000000
C13Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C13H' & read_stats$quality=='fail')])/1000000

C15Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C15H' & read_stats$quality=='pass')])/1000000
C15Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C15H' & read_stats$quality=='fail')])/1000000

species = c( 'C13H', 'C13H', 'C15H', 'C15H', 'CB1', 'CB1', 'CB2', 'CB2')
species = factor(species, levels = c( 'C13H', 'C15H', 'CB1', 'CB2'),ordered = TRUE)
condition = c(rep(c('pass','fail'), 4))
value = c(C13Hp,C13Hf,C15Hp,C15Hf,CB1p,CB1f,CB2p,CB2f)
single_read_stats.df = data.frame(species,condition,value)

ggplot(single_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Assessing proportion of pass and fail reads among capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Proportion of reads") +
  my.theme

ggplot(single_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

#Now we try to see how barcodes are distributed across reads for each sample
C13H <- subset(read_stats, sample=='C13H')
species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
             rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(rep(c('fail','pass'), 13))
value = c((C13H$No_reads)/1000000)
C13H_read_stats.df = data.frame(species,condition,value)
ggplot(C13H_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among barcodes for capitata 13H samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

C15H <- subset(read_stats, sample=='C15H')
species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
             rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(rep(c('fail','pass'), 13))
value = c((C15H$No_reads)/1000000)
C15H_read_stats.df = data.frame(species,condition,value)
ggplot(C15H_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among barcodes for capitata 15H samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

CB1 <- subset(read_stats, sample=='CB1')
species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
             rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(rep(c('fail','pass'), 13))
value = c((CB1$No_reads)/1000000)
CB1_read_stats.df = data.frame(species,condition,value)
ggplot(CB1_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among barcodes for capitata CB1 samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

CB2 <- subset(read_stats, sample=='CB2')
species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
             rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(rep(c('fail','pass'), 13))
value = c((CB2$No_reads)/1000000)
CB2_read_stats.df = data.frame(species,condition,value)
ggplot(CB2_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among barcodes for capitata CB2 samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme


#Now let us try to get the stats of the aligned reads
align.stats <- fread("barcode_references_sorted_combined",header = F, sep = ' ')
head(align.stats)
colnames(align.stats) <- c('No_reads','sample','barcode','aligned')
species = factor(align.stats$barcode, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(align.stats$aligned)
value = c((align.stats$No_reads)/1000000)
new_read_stats.df = data.frame(species,condition,value)
new_read_stats.df$condition <- factor(new_read_stats.df$condition, levels = c("ERCC","MIT","NOG","GEN"))
ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing aligning partern of reads among barcodes for capitata 1H samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme
#In reverse order
#new_read_stats.df$condition <- factor(new_read_stats.df$condition, levels = rev(c("ERCC","MIT","NOG","GEN")))
#Alternatively we can use the function fct_rev library(forcats) like so: geom_bar(aes(fill = fct_rev(condition)),
#If we'd like to reverse the stacked order but but keeping the order of the legend, we use the argument position_stack(reverse = TRUE)
#geom_bar(aes(fill = fct_rev(quality)), stat = "identity", position = position_stack(reverse = TRUE)) +

loop.samples <- function(table){
  for(i in c('CB1','CB2','C13','C15')){
    samplex = i #paste('C',i,'H', sep = '')
    C1H <- subset(table, sample==samplex)
    species = factor(C1H$barcode, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
    condition = c(C1H$aligned)
    value = c((C1H$No_reads)/1000000)
    #print(head(value))
    new_read_stats.df = data.frame(species,condition,value)
    new_read_stats.df$condition <- factor(new_read_stats.df$condition, levels = c("ERCC","MIT","NOG","GEN"))
    p <- ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
      geom_bar(position="stack", stat="identity") +
      ggtitle(paste("Assessing aligning partern of reads among barcodes for capitata ",i,"H samples", sep = '')) +
      xlab("Time points (hours AEL)") +
      ylab("Number of reads (millions)") +
      my.theme
    print(p)
  }
}
loop.samples(align.stats)

