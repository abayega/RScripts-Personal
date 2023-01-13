rm(list = ls())
setwd("~/R/tests/test39")

library("ggplot2")
library('gplots') #has the heatmap.2 function

boleae <- read.table("final_result.P_VSC_UK", header = F) #, row.names = 1, sep = "\t")
boleae2 <- read.table("final_result.P_VSC_UK2", header = F) #, row.names = 1, sep = "\t")
boleae3 <- read.table("final_result.P_VSC_UK_unvalidated", header = F)
boleae4 <- read.table("final_result.P_VSC_UK3", header = F)
boleae5 <- read.table("final_result.P_VSC_UK4", header = F) #VSC_K > 20
ygsResult <- read.table("final_result.edit_male_contigs_readlengths_GC.txt", row.names = 1, header = F)
ygsResult <- read.table("final_result.edit_male_contigs_readlengths_GC.txt", header = F)
head(ygsResult)

plot(boleae$V2, boleae$V1, ylim = c(0,100))

hist(boleae2$V1, ylim = c(0,3000))
hist(boleae3$V1)

colnames(boleae4) <- c('VSC_UK','P_VSC_UK')
plot(boleae4$V2,boleae4$V1, ylim = c(0,1000))

colnames(boleae5) <- c('VSC_UK','P_VSC_UK')
plot(boleae5$P_VSC_UK,boleae5$VSC_UK, ylim = c(0,1000))

colnames(ygsResult) <- c('contig1','length','GC','contig','NUM','MAX_K','K','UK','SC_K','SC_UK','P_SC_UK','VSC_K','VSC_UK','P_VSC_UK')

#histogram of P_VSC_UK
ygsResult_ed <- data.frame(lapply(ygsResult[which(ygsResult$P_VSC_UK != "."),], as.character), stringsAsFactors=FALSE)

hist(as.numeric(ygsResult_ed$P_VSC_UK), main = "Percentage of validated single copy unmatched", 
     xlab = "% of contig unmatched by female reads", ylab = "Number of scaffolds")

yValidatedAbove80 <- ygsResult_ed$contig[which(as.numeric(ygsResult_ed$P_VSC_UK) >= 80 & as.numeric(ygsResult_ed$VSC_K) > 20)]
write(yValidatedAbove80, file = "YGS_validated_above80")

yValidatedAbove90 <- ygsResult_ed$contig[which(as.numeric(ygsResult_ed$P_VSC_UK) >= 90 & as.numeric(ygsResult_ed$VSC_K) > 20)]
write(yValidatedAbove90, file = "YGS_validated_above90")

yValidatedAbove100 <- ygsResult_ed$contig[which(as.numeric(ygsResult_ed$P_VSC_UK) >= 100 & as.numeric(ygsResult_ed$VSC_K) > 20)]
write(yValidatedAbove100, file = "YGS_validated_above100")

ygsResult_ed2 <- subset(ygsResult_ed, as.numeric(P_VSC_UK) >= 90 & as.numeric(VSC_K) > 20)
head(ygsResult_ed2)
sum(as.numeric(ygsResult_ed2$length))

###Very clean down here to identify the total length of scaffolds that are common between YGS and CQ
cq_ygs <- c('Contig10795','Contig2397','Contig26517','Contig2933','Contig4856','Contig4896','Contig528',
            'Contig640','Contig7188','Contig719','Contig7205','Contig7589','Contig9401','Contig99') #These are the common scaffolds

ygsResult2 <- read.table("final_result.edit_male_contigs_readlengths_GC.txt", header = F)
colnames(ygsResult2) <- c('scaffolds','length','GC','contig','NUM','MAX_K','K','UK','SC_K','SC_UK','P_SC_UK','VSC_K','VSC_UK','P_VSC_UK')

sum(ygsResult2$length[which(ygsResult2$scaffolds %in% cq_ygs)])
#end

#090220
#Trying to find scaffolds that are common between CQ and YGS
cq_readLengths <- read.table("pbjelly_10X_v1_1.corrvalid_group_gt10_500000bc_74X_seed0.plus.all_LouisL_MP.sspace.final.scaffolds.sort-batch.fin_ONTFQ_PACBIO.jelly.out.repeatmasked.fillrepeat.Y_readlengths_GC.txt", header = F, sep = '\t')
colnames(cq_readLengths) <- c("contig_name", 'length', 'GC_percentage')

#number of common scaffolds
length(which(cq_readLengths$contig_name %in% ygsResult_ed2$contig1))

#Length of common scaffolds
sum(ygsResult_ed2$length(which(cq_readLengths$contig_name %in% ygsResult_ed2$contig1)))

CQ_YGS_intesection <- subset(ygsResult_ed2, ygsResult_ed2$contig1 %in% cq_readLengths$contig_name)
100*sum(as.numeric(CQ_YGS_intesection$length))/sum(as.numeric(ygsResult_ed2$length))

#So, 68% of sequence between YGS and CQ agrees that it is Y chromosome. They disagree on 32%

#So, both methods came up with the same amount of sequence that is Y chromosome, 3.9 Mb
sum(cq_readLengths$length)


#Trying to get scaffold lengths
lengths_table <- read.table('LGAM02_contigs_conversion_table.tsv', header = T, sep = '\t')
names_34 <- read.table('34_PCR_confirmed', header = F)
names_CQ <- read.table('CQ_mthod', header = F)
names_YGS <- read.table('YGS_method', header = F)
names_55 <- read.table('55_PCR_unconfirmed', header = F)

names_34_sub <- subset(lengths_table, lengths_table$accession %in% names_34$V1)
names_34_sub <- as.matrix(cbind(as.character(names_34_sub$accession),names_34_sub$length))
colnames(names_34_sub) <- c('accession','length')
write(names_34_sub, file = "names_34_sub", sep = '\t')

names_CQ_sub <- subset(lengths_table, lengths_table$accession %in% names_CQ$V1)
names_CQ_sub <- as.matrix(cbind(as.character(names_CQ_sub$accession),names_CQ_sub$length))
colnames(names_CQ_sub) <- c('accession','length')
write(names_CQ_sub, file = "names_CQ_sub")

names_YGS_sub <- subset(lengths_table, lengths_table$accession %in% names_YGS$V1)
names_YGS_sub <- as.matrix(cbind(as.character(names_YGS_sub$accession),names_YGS_sub$length))
colnames(names_YGS_sub) <- c('accession','length')
write(names_YGS_sub, file = "names_YGS_sub")

names_55_sub <- subset(lengths_table, lengths_table$accession %in% names_55$V1)
names_55_sub <- as.matrix(cbind(as.character(names_55_sub$accession),names_55_sub$length))
colnames(names_55_sub) <- c('accession','length')
write(names_55_sub, file = "names_55_sub", sep = '\t')
