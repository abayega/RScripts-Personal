rm(list = ls())
setwd("~/R/tests/test11")

file1 <- read.table('ERCC92_sl_ONT_raw_reads_ed_with_alignments.txt', header = T, sep = '\t')
file2 <- read.table('ERCC92_sl_ONT_canu_corrected_reads_with_alignments.txt', header = T, sep = '\t')
file3 <- read.table('ERCC92_sl_ONT_LSC_corrected_reads_ed_with_alignments.txt', header = T, sep = '\t')
file4 <- read.table('ERCC92_sl_Illumina_4525.004.R1_ed_with_alignments.txt', header = T, sep = '\t')
file5 <- read.table('saturation', header = T, sep = '\t')
file6 <- read.table('ERCC92_sl_ercc_mapped_corrected2.correctedReads_canu2_with_alignments.txt', header = T, sep = '\t')
file7 <- read.table('Bo_Y_contigs_olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_sl2_male_contig_read_names3_sorted_2_with_alignments.txt', header = T, sep = '\t')
file8 <- read.table('Bo_Y_contigs_olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_sl2_male_y_mapped_corrected.correctedReads_with_alignments.txt', header = T, sep = '\t')
file9 <- read.table('Bo_Y_contigs_olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_sl2_male_y_contigs_canu_lsc_corrected_ed_with_alignments.txt', header = T, sep = '\t')
file10 <- read.table('Bo_Y_contigs_olivefly_supernova_v1p1_corrvalid_group_gt10_500000bc_74X_seed0.pseudohap.1_sl2_male_y_contigs_canu_lsc_corrected_round2_ed_with_alignments.txt', header = T, sep = '\t')
head(file1)

a <- as.numeric(file1$align_identity)
b <- as.numeric(file2$align_identity)
c <- as.numeric(file3$align_identity)
d <- as.numeric(file4$align_identity)
e <- as.numeric(file6$align_identity)
f <- as.numeric(file7$align_identity)
g <- as.numeric(file8$align_identity)
h <- as.numeric(file9$align_identity)
i <- as.numeric(file10$align_identity)

###Trying to plot histo of alignment ID
hist(a,breaks=100,col=rgb(1,0,0,0.5),main="",xlab="align identity", xlim = c(70,100))
par(new=TRUE)
hist(b,breaks=50,col=rgb(0,1,0,0.5),main="",xlab="",ylab="", xlim = c(70,100), axes=FALSE)
par(new=TRUE)
hist(c,breaks=50,col=rgb(0,0,1,0.5),main="",xlab="",ylab="", xlim = c(70,100),axes=FALSE)
par(new=TRUE)
hist(d,breaks=50,col=rgb(0,0,1,0.5),main="",xlab="",ylab="", xlim = c(70,100),axes=FALSE)
legend("topleft", legend=c("ONT_raw","ONT_canu","ONT_ill","Ill"), fill=c("red","green","blue","gray"))
dev.off()

ggplot() + geom_density(aes(x=x), fill="red", data=data.frame(x=a), alpha=.5) + 
  geom_density(aes(x=x), fill="blue", data=data.frame(x=b), alpha=.5) +
  geom_density(aes(x=x), fill="green", data=data.frame(x=c), alpha=.5) +
  geom_density(aes(x=x), fill="yellow", data=data.frame(x=d), alpha=.5) +
  scale_colour_manual(values=c("red"="ONT_raw", "blue"="ONT_canu", "green"="ONT_ill", "yellow"="Ill"), name="Samples")
#scale_colour_manual(values=c("Bo_1H"="red", "Bo_2H"="blue", "Bo_3H"="green", "Bo_4H"="yellow", "Bo_5H"="purple", "Bo_6H"="orange"), name="Samples")
dev.off()


bmmisrate <- cbind(a,b,e,c,d,f,g,h,i)
boxplot.matrix(bmmisrate, use.cols = T,ylim=c(50,100),ylab="align Identity (%)",col = c('red', 'green', 'blue', 'gray','yellow','purple','orange','yellow','yellow'), names=c("ONT_raw","ONT_canu","ONT_canu2","ONT_ill","Ill",'y-raw','y_canu','y_canu_lsc','y_canu_lsc2'), main = 'Alignment Identity and corection')

a1 <- as.numeric(file1$mismatches)
a2 <- as.numeric(file2$mismatches)
b1 <- as.numeric(file1$insertion)
b2 <- as.numeric(file2$insertion)
c1 <- as.numeric(file1$deletion)
c2 <- as.numeric(file2$deletion)
d1 <- as.numeric(file1$perc_read_aligned)
d2 <- as.numeric(file2$perc_read_aligned)

bmmisrate <- cbind(a1,a2,b1,b2,c1,c2,d1,d2)
boxplot.matrix(bmmisrate, use.cols = T,ylim=c(0,100),ylab="Percentage",col = c('red','red', 'green','green', 'blue','blue', 'gray','gray'), names=c("mismatches","", "Insertion","", "deletion","", "aligned",""), main = 'Effect of canu correction')


a1 <- as.numeric(file1$mismatches)
a2 <- as.numeric(file3$mismatches)
b1 <- as.numeric(file1$insertion)
b2 <- as.numeric(file3$insertion)
c1 <- as.numeric(file1$deletion)
c2 <- as.numeric(file3$deletion)
d1 <- as.numeric(file1$perc_read_aligned)
d2 <- as.numeric(file3$perc_read_aligned)

bmmisrate <- cbind(a1,a2,b1,b2,c1,c2,d1,d2)
boxplot.matrix(bmmisrate, use.cols = T,ylim=c(0,100),ylab="Percentage",col = c('red','red', 'green','green', 'blue','blue', 'gray','gray'), names=c("mismatches","", "Insertion","", "deletion","", "aligned",""), main = 'Effect of LSC correction')

a1 <- as.numeric(file2$mismatches)
a2 <- as.numeric(file3$mismatches)
b1 <- as.numeric(file2$insertion)
b2 <- as.numeric(file3$insertion)
c1 <- as.numeric(file2$deletion)
c2 <- as.numeric(file3$deletion)
d1 <- as.numeric(file2$perc_read_aligned)
d2 <- as.numeric(file3$perc_read_aligned)

bmmisrate <- cbind(a1,a2,b1,b2,c1,c2,d1,d2)
boxplot.matrix(bmmisrate, use.cols = T,ylim=c(0,30),ylab="Percentage",col = c('red','red', 'green','green', 'blue','blue', 'gray','gray'), names=c("mis_canu","mis_LSC", "Ins_canu","Ins_LSC", "del_canu","del_LSC", "aligned_canu","aligned_LSC"), main = 'Comparing Canu and LSC correction')

bmmisrate <- cbind(a1,a2,b1,b2,c1,c2)
boxplot.matrix(bmmisrate, use.cols = T,ylim=c(0,30),ylab="Percentage",col = c('red','red', 'green','green', 'blue','blue'), names=c("mis_canu","mis_LSC", "Ins_canu","Ins_LSC", "del_canu","del_LSC"), main = 'Comparing Canu and LSC correction')

head(file5)
#Plot saturation
plot(log10(file5$reads),file5$genes_ONT, type="b", pch=18, col="red", xlab=expression('log'[10]*"(Total reads or Total bases)"),ylab="No. of genes", xlim = c(3,11), ylim = c(0,12000))
lines(log10(file5$reads.1),file5$genes._illumina, pch=18, col="blue", type="b", lty=2) 
lines(log10(file5$bases),file5$genes_ONT, pch=18, col="green", type="b", lty=2)
lines(log10(file5$bases.1),file5$genes._illumina, pch=18, col="black", type="b", lty=2)
legend("topleft", legend=c("ONT-reads", "Illumina-reads", "ONT-bases", "Illumina-bases"),
       col=c("red", "blue", "green", "black"), lty=c(2,2,2,2), cex=0.8)

