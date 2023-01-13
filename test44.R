# Set working directory
setwd("/home/banthony/R/tests/test44")

read_stats <- read.table("capitata_reads_stats",header=T, sep = '\t')
head(new_read_stats)
dim(new_read_stats)

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
species = factor(new_read_stats$sample, levels = c( 'C1H', 'C2H', 'C3H', 'C4H', 'C5H', 'C6H', 'C7H', 'C8H', 'C9H',
                                                    'C10H', 'C11H', 'C12H', 'C13H', 'C14H', 'C15H'),ordered = TRUE)
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
new_read_stats <- subset(read_stats, quality=="fail")
species = factor(new_read_stats$sample, levels = c( 'C1H', 'C2H', 'C3H', 'C4H', 'C5H', 'C6H', 'C7H', 'C8H', 'C9H',
                                                    'C10H', 'C11H', 'C12H', 'C13H', 'C14H', 'C15H'),ordered = TRUE)
condition = new_read_stats$barcode
value = new_read_stats$No_reads
new_read_stats.df = data.frame(species,condition,value)

ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Assessing fail reads distribution among barcodes for capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Proportion of reads") +
  my.theme

#Let's try to plot the total yields regardless of the barcodes
C1Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C1H' & read_stats$quality=='pass')])/1000000
C1Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C1H' & read_stats$quality=='fail')])/1000000

C2Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C2H' & read_stats$quality=='pass')])/1000000
C2Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C2H' & read_stats$quality=='fail')])/1000000

C3Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C3H' & read_stats$quality=='pass')])/1000000
C3Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C3H' & read_stats$quality=='fail')])/1000000

C4Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C4H' & read_stats$quality=='pass')])/1000000
C4Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C4H' & read_stats$quality=='fail')])/1000000

C5Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C5H' & read_stats$quality=='pass')])/1000000
C5Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C5H' & read_stats$quality=='fail')])/1000000

C6Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C6H' & read_stats$quality=='pass')])/1000000
C6Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C6H' & read_stats$quality=='fail')])/1000000

C7Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C7H' & read_stats$quality=='pass')])/1000000
C7Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C7H' & read_stats$quality=='fail')])/1000000

C8Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C8H' & read_stats$quality=='pass')])/1000000
C8Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C8H' & read_stats$quality=='fail')])/1000000

C9Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C9H' & read_stats$quality=='pass')])/1000000
C9Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C9H' & read_stats$quality=='fail')])/1000000

C10Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C10H' & read_stats$quality=='pass')])/1000000
C10Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C10H' & read_stats$quality=='fail')])/1000000

C11Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C11H' & read_stats$quality=='pass')])/1000000
C11Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C11H' & read_stats$quality=='fail')])/1000000

C12Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C12H' & read_stats$quality=='pass')])/1000000
C12Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C12H' & read_stats$quality=='fail')])/1000000

C13Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C13H' & read_stats$quality=='pass')])/1000000
C13Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C13H' & read_stats$quality=='fail')])/1000000

C14Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C14H' & read_stats$quality=='pass')])/1000000
C14Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C14H' & read_stats$quality=='fail')])/1000000

C15Hp <- sum(read_stats$No_reads[which(read_stats$sample == 'C15H' & read_stats$quality=='pass')])/1000000
C15Hf <- sum(read_stats$No_reads[which(read_stats$sample == 'C15H' & read_stats$quality=='fail')])/1000000

species = c( 'C1H', 'C1H', 'C2H', 'C2H', 'C3H', 'C3H', 'C4H', 'C4H', 'C5H','C5H', 'C6H', 'C6H', 'C7H', 'C7H', 'C8H','C8H',
             'C9H', 'C9H', 'C10H', 'C10H', 'C11H', 'C11H', 'C12H', 'C12H', 'C13H','C13H', 'C14H', 'C14H', 'C15H', 'C15H')
species = factor(species, levels = c( 'C1H', 'C2H', 'C3H', 'C4H', 'C5H', 'C6H', 'C7H', 'C8H', 'C9H','C10H', 'C11H', 'C12H', 'C13H', 'C14H', 'C15H'),ordered = TRUE)
condition = c(rep(c('pass','fail'), 15))
value = c(C1Hp,C1Hf,C2Hp,C2Hf,C3Hp,C3Hf,C4Hp,C4Hf,C5Hp,C5Hf,C6Hp,C6Hf,C7Hp,C7Hf,C8Hp,C8Hf,
          C9Hp,C9Hf,C10Hp,C10Hf,C11Hp,C11Hf,C12Hp,C12Hf,C13Hp,C13Hf,C14Hp,C14Hf,C15Hp,C15Hf)
new_read_stats.df = data.frame(species,condition,value)

ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Assessing proportion of pass and fail reads among capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Proportion of reads") +
  my.theme

ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among capitata samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

#Now we try to see how barcodes are distributed across reads for each sample
C1H <- subset(read_stats, sample=='C1H')
species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
             rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
condition = c(rep(c('fail','pass'), 13))
value = c((C1H$No_reads)/1000000)
new_read_stats.df = data.frame(species,condition,value)
ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
  geom_bar(position="stack", stat="identity") +
  ggtitle("Assessing number of pass and fail reads among barcodes for capitata 1H samples") +
  xlab("Time points (hours AEL)") +
  ylab("Number of reads (millions)") +
  my.theme

loop.samples <- function(table){
  for(i in 1:15){
    samplex = paste('C',i,'H', sep = '')
    C1H <- subset(read_stats, sample==samplex)
    species = c( rep('B01',2),rep('B02',2),rep('B03',2),rep('B04',2),rep('B05',2),rep('B06',2),rep('B07',2),rep('B08',2),
                 rep('B09',2),rep('B10',2),rep('B11',2),rep('B12',2),rep('NON',2) )
    species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','NON'), ordered = TRUE)
    condition = c(rep(c('fail','pass'), 13))
    value = c((C1H$No_reads)/1000000)
    #print(head(value))
    new_read_stats.df = data.frame(species,condition,value)
    p <- ggplot(new_read_stats.df, aes(x=species, y=value, fill=condition)) + 
      geom_bar(position="stack", stat="identity") +
      ggtitle(paste("Assessing number of pass and fail reads among barcodes for capitata ",i,"H samples", sep = '')) +
      xlab("Time points (hours AEL)") +
      ylab("Number of reads (millions)") +
      my.theme
    print(p)
  }
}
loop.samples(read_stats)

#Now let us try to get the stats of the aligned reads
align.stats <- fread("barcode_references_sorted_combined",header = F, sep = ' ')
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
  for(i in 1:15){
    samplex = paste('C',i,'H', sep = '')
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

#Let us try to get correlation of read lengths and alignment stats
readLength <- fread("read_lengths_identity",header = T, sep = '\t')
colnames(readLength) <- c('read_length','alignment_identity','read_aligned')

plot(readLength$read_length,readLength$`alignment_identity_%`)

lm_eqn1 = function(df1){
  m = lm(alignment_identity ~ read_length, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
lm_eqn2 = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(df1$alignment_identity,df1$read_length, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}
ggplot(data = readLength, aes(x = read_length, y = alignment_identity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Correlation of read length with alignment identity") +
  scale_x_continuous(expression("read length (bp)")) +
  scale_y_continuous(expression("alignment identity (%)")) +
  annotate(geom = "text", x = 2000, y = 10, label = lm_eqn1(readLength), parse = T, color = "red") +
  annotate(geom = "text", x = 2000, y = 20, label = lm_eqn2(readLength), parse = T, color = "red") +
  my.theme

##Plotting correlation of read length with percentage of read aligned
lm_eqn1 = function(df1){
  m = lm(read_aligned ~ read_length, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
lm_eqn2 = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(df1$read_aligned,df1$read_length, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}
ggplot(data = readLength, aes(x = read_length, y = read_aligned)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Correlation of read length with percentage of read aligned") +
  scale_x_continuous(expression("read length (bp)")) +
  scale_y_continuous(expression("Percentage of read aligned (%)")) +
  ylim(0,100) +
  annotate(geom = "text", x = 5000, y = 10, label = lm_eqn1(readLength), parse = T, color = "red") +
  annotate(geom = "text", x = 5000, y = 20, label = lm_eqn2(readLength), parse = T, color = "red") +
  my.theme
##Plotting correlation of alignment identity with percentage of read aligned
lm_eqn1 = function(df1){
  m = lm(alignment_identity ~ read_aligned, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                         b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
lm_eqn2 = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(df1$alignment_identity,df1$read_aligned, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}
ggplot(data = readLength, aes(x = read_aligned, y = alignment_identity)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("Correlation of alignment identity with percentage of read aligned") +
  scale_x_continuous(expression("Percentage of read aligned (%)")) +
  scale_y_continuous(expression("alignment identity (%)")) +
  xlim(0,100) +
  ylim(0,100) +
  annotate(geom = "text", x = 20, y = 10, label = lm_eqn1(readLength), parse = T, color = "red") +
  annotate(geom = "text", x = 20, y = 20, label = lm_eqn2(readLength), parse = T, color = "red") +
  my.theme

#Now, we shall try to see if read length is different for:
#1: fail and pass reads
#2: classified and unclassified reads
#3: aligned and unaligned reads
main_dir = '/media/banthony/f7f2c359-0d93-4e31-ac8e-888a669f866f/reads/nanopore/capitata/rna_seq/1/'
runame = 'fastq_runid_42a2db9e01c2eeb64918dc3fe576154d8be18349_renamed_readlengths_GC.txt'

combine_tables <- function(){
  C1H.fail <- fread(paste(main_dir,"fail/unclassified/",runame, sep = ''), header = F)
  colnames(C1H.fail) <- c('UN_length','UN_GC')
  for(i in 1:12){
    table <- fread(paste(main_dir,"fail/",i,"/",runame, sep = ''), header = F)
    if(i < 10){
      colnames(table) <- c(paste('B0',i,'_length',sep = ''),paste0('B0',i,'_GC'))
    }
    else{
      colnames(table) <- c(paste('B',i,'_length',sep = ''),paste0('B',i,'_GC'))
    }
    C1H.fail <- cbind(C1H.fail,table)
  }
  return(C1H.fail)
}
C1H.fail <- combine_tables()

combine_tables <- function(){
  C1H.pass <- fread(paste(main_dir,"pass/unclassified/",runame, sep = ''), header = F)
  colnames(C1H.pass) <- c('UN_length','UN_GC')
  for(i in 1:12){
    table <- fread(paste(main_dir,"pass/",i,"/",runame, sep = ''), header = F)
    if(i < 10){
      colnames(table) <- c(paste('B0',i,'_length',sep = ''),paste0('B0',i,'_GC'))
    }
    else{
      colnames(table) <- c(paste('B',i,'_length',sep = ''),paste0('B',i,'_GC'))
    }
    C1H.pass <- cbind(C1H.pass,table)
  }
  return(C1H.pass)
}
C1H.pass <- combine_tables()

species <- c(rep('B01',(length(C1H.fail$B01_length)+length(C1H.pass$B01_length))),rep('B02',(length(C1H.fail$B02_length)+length(C1H.pass$B02_length))),
             rep('B03',(length(C1H.fail$B03_length)+length(C1H.pass$B03_length))),rep('B04',(length(C1H.fail$B04_length)+length(C1H.pass$B04_length))),
             rep('B05',(length(C1H.fail$B05_length)+length(C1H.pass$B05_length))),rep('B06',(length(C1H.fail$B06_length)+length(C1H.pass$B06_length))),
             rep('B07',(length(C1H.fail$B07_length)+length(C1H.pass$B07_length))),rep('B08',(length(C1H.fail$B08_length)+length(C1H.pass$B08_length))),
             rep('B09',(length(C1H.fail$B09_length)+length(C1H.pass$B09_length))),rep('B10',(length(C1H.fail$B10_length)+length(C1H.pass$B10_length))),
             rep('B11',(length(C1H.fail$B11_length)+length(C1H.pass$B11_length))),rep('B12',(length(C1H.fail$B12_length)+length(C1H.pass$B12_length))),
             rep('UN',(length(C1H.fail$UN_length)+length(C1H.pass$UN_length)))
             )
species = factor(species, levels = c('B01','B02','B03','B04','B05','B06','B07','B08','B09','B10','B11','B12','UN'), ordered = TRUE)
condition <- c(rep('fail',length(C1H.fail$B01_length)),rep('pass',length(C1H.pass$B01_length)),rep('fail',length(C1H.fail$B02_length)),rep('pass',length(C1H.pass$B02_length)),
               rep('fail',length(C1H.fail$B03_length)),rep('pass',length(C1H.pass$B03_length)),rep('fail',length(C1H.fail$B04_length)),rep('pass',length(C1H.pass$B04_length)),
               rep('fail',length(C1H.fail$B05_length)),rep('pass',length(C1H.pass$B05_length)),rep('fail',length(C1H.fail$B06_length)),rep('pass',length(C1H.pass$B06_length)),
               rep('fail',length(C1H.fail$B07_length)),rep('pass',length(C1H.pass$B07_length)),rep('fail',length(C1H.fail$B08_length)),rep('pass',length(C1H.pass$B08_length)),
               rep('fail',length(C1H.fail$B09_length)),rep('pass',length(C1H.pass$B09_length)),rep('fail',length(C1H.fail$B10_length)),rep('pass',length(C1H.pass$B10_length)),
               rep('fail',length(C1H.fail$B11_length)),rep('pass',length(C1H.pass$B11_length)),rep('fail',length(C1H.fail$B12_length)),rep('pass',length(C1H.pass$B12_length)),
               rep('fail',length(C1H.fail$UN_length)),rep('pass',length(C1H.pass$UN_length))
                 )
value <- c(C1H.fail$B01_length,C1H.pass$B01_length,C1H.fail$B02_length,C1H.pass$B02_length,C1H.fail$B03_length,C1H.pass$B03_length,
           C1H.fail$B04_length,C1H.pass$B04_length,C1H.fail$B05_length,C1H.pass$B05_length,C1H.fail$B06_length,C1H.pass$B06_length,
           C1H.fail$B07_length,C1H.pass$B07_length,C1H.fail$B08_length,C1H.pass$B08_length,C1H.fail$B09_length,C1H.pass$B09_length,
           C1H.fail$B10_length,C1H.pass$B10_length,C1H.fail$B11_length,C1H.pass$B11_length,C1H.fail$B12_length,C1H.pass$B12_length,
           C1H.fail$UN_length,C1H.pass$UN_length)
new_read_stats.df = data.frame(species,condition,value)
new_read_stats.df2 <- subset(new_read_stats.df, new_read_stats.df[ , 3] < 1000)

new_read_stats.df2 <- subset(new_read_stats.df, new_read_stats.df[ , 3] > 10 & new_read_stats.df[ , 1] != "B03" & new_read_stats.df[ , 1] != "B08")

p <- ggplot(new_read_stats.df2, aes(x=species, y=value, fill=condition)) + #, color=fac)
  geom_violin(position=position_dodge(1)) +
  #geom_violin(fill='#A4A4A4', color="darkred", position=position_dodge(1)) +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "barcode") +
  labs(fill = 'quality') +
  labs(y = "read length (bp)") + 
  ylim(0,1000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

unmappedReads <- fread("unmapped_readlengths", header = F)
mappedReads <- fread("read_lengths_identity", header = T)
unmappedReads <- unmappedReads[1:1000000,]
species <- c(rep('unmapped',length(unmappedReads$V1)),rep('mapped',length(mappedReads$read_length)))
species = factor(species, levels = c('unmapped','mapped'), ordered = TRUE)
value <- c(unmappedReads$V1,mappedReads$read_length)

new_read_stats.df = data.frame(species,value)
p <- ggplot(new_read_stats.df, aes(x=species, y=value)) +
  geom_violin(fill='#A4A4A4', color="darkred") +
  #geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "barcode") +
  labs(fill = 'quality') +
  labs(y = "read length (bp)") + 
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

##Let's try alignment identity profile
readLength <- fread("mapped_read_lengths_identity_all",header = T, sep = '\t')
colnames(readLength) <- c('read_length','alignment_identity','read_aligned')

ggplot(readLength, aes(x=alignment_identity)) +
  geom_histogram(binwidth=2, fill='white', color="black") +
  #geom_histogram(aes(y=..density..), binwidth=2, fill='white', color="black") +
  #geom_density(alpha=.2, fill = "#FF6666") +
  geom_vline(aes(xintercept=mean(alignment_identity)),color="blue", linetype="dashed", size=1) +
  geom_vline(aes(xintercept=median(alignment_identity)),color="red", linetype="dashed", size=1) +
  labs(x = "read identity (%)") +
  labs(y = "count") + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
