# Set working directory
##For combined transcripts
setwd("/home/abayega/R/tests/test27/all_combined/")
introns_1 <- read.table("coverage.1H.bed",header=F, sep = '\t')
introns_2 <- read.table("coverage.2H.bed",header=F, sep = '\t')
introns_3 <- read.table("coverage.3H.bed",header=F, sep = '\t')
introns_4 <- read.table("coverage.4H.bed",header=F, sep = '\t')
introns_5 <- read.table("coverage.5H.bed",header=F, sep = '\t')
introns_6 <- read.table("coverage.6H.bed",header=F, sep = '\t')
introns_female <- read.table("coverage.female.bed",header=F, sep = '\t')
introns_male <- read.table("coverage.male.bed",header=F, sep = '\t')

reads1=1247035
reads2=1703243
reads3=1287945
reads4=2338587
reads5=3729976
reads6=1357637
reads.female=1620587
reads.male=1673016

head(introns_1)
introns_1$V17 <- 10000*introns_1$V13/reads1
introns_2$V17 <- 10000*introns_2$V13/reads2
introns_3$V17 <- 10000*introns_3$V13/reads3
introns_4$V17 <- 10000*introns_4$V13/reads4
introns_5$V17 <- 10000*introns_5$V13/reads5
introns_6$V17 <- 10000*introns_6$V13/reads6
introns_female$V17 <- 10000*introns_female$V13/reads.female
introns_male$V17 <- 10000*introns_male$V13/reads.male
n <- 0.04
intron_cov <- cbind(introns_1$V17[which(introns_1$V17 < n)],introns_2$V17[which(introns_2$V17 < n)],introns_3$V17[which(introns_3$V17 < n)],
                    introns_4$V17[which(introns_4$V17 < n)],introns_5$V17[which(introns_5$V17 < n)],introns_6$V17[which(introns_6$V17 < n)],
                    introns_female$V17[which(introns_female$V17 < n)],introns_male$V17[which(introns_male$V17 < n)])
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0.005,n))


##This for zygotic_early
setwd("/home/abayega/R/tests/test27/zygotic_early/")
introns_1 <- read.table("coverage.1H.bed.zygotic_early",header=F, sep = '\t')
introns_2 <- read.table("coverage.2H.bed.zygotic_early",header=F, sep = '\t')
introns_3 <- read.table("coverage.3H.bed.zygotic_early",header=F, sep = '\t')
introns_4 <- read.table("coverage.4H.bed.zygotic_early",header=F, sep = '\t')
introns_5 <- read.table("coverage.5H.bed.zygotic_early",header=F, sep = '\t')
introns_6 <- read.table("coverage.6H.bed.zygotic_early",header=F, sep = '\t')
introns_female <- read.table("coverage.female.bed.zygotic_early",header=F, sep = '\t')
introns_male <- read.table("coverage.male.bed.zygotic_early",header=F, sep = '\t')

reads1=1247035
reads2=1703243
reads3=1287945
reads4=2338587
reads5=3729976
reads6=1357637
reads.female=1620587
reads.male=1673016

head(introns_1)
introns_1$V17 <- 10000*introns_1$V13/reads1
introns_2$V17 <- 10000*introns_2$V13/reads2
introns_3$V17 <- 10000*introns_3$V13/reads3
introns_4$V17 <- 10000*introns_4$V13/reads4
introns_5$V17 <- 10000*introns_5$V13/reads5
introns_6$V17 <- 10000*introns_6$V13/reads6
introns_female$V17 <- 10000*introns_female$V13/reads.female
introns_male$V17 <- 10000*introns_male$V13/reads.male
n <- 0.04
intron_cov <- cbind(introns_1$V17[which(introns_1$V17 < n)],introns_2$V17[which(introns_2$V17 < n)],introns_3$V17[which(introns_3$V17 < n)],
                    introns_4$V17[which(introns_4$V17 < n)],introns_5$V17[which(introns_5$V17 < n)],introns_6$V17[which(introns_6$V17 < n)],
                    introns_female$V17[which(introns_female$V17 < n)],introns_male$V17[which(introns_male$V17 < n)])
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0,n), names=c("1H","2H","3H","4H","5H","6H","adult_female","adult_male"), 
               main = 'intron coverage among zygotic genes', ylab='intron coverage (RPG10K)')


#Let us try to plot gene expressions
introns_1 <- read.table("coverage.1H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_2 <- read.table("coverage.2H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_3 <- read.table("coverage.3H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_4 <- read.table("coverage.4H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_5 <- read.table("coverage.5H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_6 <- read.table("coverage.6H.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_female <- read.table("coverage.female.bed.zygotic_early_gene_expn",header=F, sep = '\t')
introns_male <- read.table("coverage.male.bed.zygotic_early_gene_expn",header=F, sep = '\t')

head(introns_1)
intron_cov <- cbind(introns_1$V25,introns_1$V26,introns_1$V27,introns_1$V28,introns_1$V29,introns_1$V30,introns_1$V31,introns_1$V32)
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0,200000), names=c("1H","2H","3H","4H","5H","6H","adult_female","adult_male"),
               main = 'absolute gene expression among zygotic genes', ylab='transcripts/embryo')


##This for matenral degraded
setwd("/home/abayega/R/tests/test27/maternal_degraded/")
introns_1 <- read.table("coverage.1H.bed.maternal_degraded",header=F, sep = '\t')
introns_2 <- read.table("coverage.2H.bed.maternal_degraded",header=F, sep = '\t')
introns_3 <- read.table("coverage.3H.bed.maternal_degraded",header=F, sep = '\t')
introns_4 <- read.table("coverage.4H.bed.maternal_degraded",header=F, sep = '\t')
introns_5 <- read.table("coverage.5H.bed.maternal_degraded",header=F, sep = '\t')
introns_6 <- read.table("coverage.6H.bed.maternal_degraded",header=F, sep = '\t')
introns_female <- read.table("coverage.female.bed.maternal_degraded",header=F, sep = '\t')
introns_male <- read.table("coverage.male.bed.maternal_degraded",header=F, sep = '\t')

reads1=1247035
reads2=1703243
reads3=1287945
reads4=2338587
reads5=3729976
reads6=1357637
reads.female=1620587
reads.male=1673016

head(introns_1)
introns_1$V18 <- 10000*introns_1$V13/reads1
introns_2$V18 <- 10000*introns_2$V13/reads2
introns_3$V18 <- 10000*introns_3$V13/reads3
introns_4$V18 <- 10000*introns_4$V13/reads4
introns_5$V18 <- 10000*introns_5$V13/reads5
introns_6$V18 <- 10000*introns_6$V13/reads6
introns_female$V18 <- 10000*introns_female$V13/reads.female
introns_male$V18 <- 10000*introns_male$V13/reads.male
n <- 100
m = 30
intron_cov <- cbind(introns_1$V18[which(introns_1$V18 < n)],introns_2$V18[which(introns_2$V18 < n)],introns_3$V18[which(introns_3$V18 < n)],
                    introns_4$V18[which(introns_4$V18 < n)],introns_5$V18[which(introns_5$V18 < n)],introns_6$V18[which(introns_6$V18 < n)],
                    introns_female$V18[which(introns_female$V18 < n)],introns_male$V18[which(introns_male$V18 < n)])
intron_cov <- cbind(introns_1$V18[which(introns_1$V18 < n & introns_1$V13 == 0)],introns_2$V18[which(introns_2$V18 < n & introns_2$V13 <= m)],
                    introns_3$V18[which(introns_3$V18 < n & introns_3$V13 <= m)],introns_4$V18[which(introns_4$V18 < n & introns_4$V13 <= m)],
                    introns_5$V18[which(introns_5$V18 < n & introns_5$V13 <= m)],introns_6$V18[which(introns_6$V18 < n & introns_6$V13 <= m)],
                    introns_female$V18[which(introns_female$V18 < n & introns_female$V13 == 0)],introns_male$V18[which(introns_male$V18 < n & introns_male$V13 == 0)])
intron_cov2 <- cbind(introns_1$V13,introns_1$V18,introns_2$V13,introns_2$V18,introns_3$V13,introns_3$V18,introns_4$V13,introns_4$V18,introns_5$V13,introns_5$V18,
                     introns_6$V13,introns_6$V18,introns_female$V13,introns_female$V18introns_male$V13,introns_male$V18)
colnames(intron_cov2) <- c('y1','V2','V3','V4','V5','V6','V7','V8','V9','V10','V11','V12','V13','V14')
intron_cov3 <- subset(data.frame(intron_cov2), y1 <= 5 & V11 <= 5 & V13 <= 5, select = c(V2,V4,V6,V8,V10,V12,V14))
boxplot(intron_cov3)
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0,0.6), names=c("1H","2H","3H","4H","5H","6H","adult_female","adult_male"), 
               main = 'intron coverage among maternal degraded genes', ylab='intron coverage (RPG10K)')


#Let us try to plot gene expressions
introns_1 <- read.table("coverage.1H.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_2 <- read.table("coverage.2H.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_3 <- read.table("coverage.3H.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_4 <- read.table("coverage.4H.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_5 <- read.table("coverage.5H.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_6 <- read.table("coverage.6H.bedmaternal_degraded_gene_expn",header=F, sep = '\t')
introns_female <- read.table("coverage.female.bed.maternal_degraded_gene_expn",header=F, sep = '\t')
introns_male <- read.table("coverage.male.bed.maternal_degraded_gene_expn",header=F, sep = '\t')

head(introns_1)
intron_cov <- cbind(introns_1$V25,introns_1$V26,introns_1$V27,introns_1$V28,introns_1$V29,introns_1$V30,introns_1$V31,introns_1$V32)
intron_cov <- cbind(introns_1$V17,introns_1$V18,introns_1$V19,introns_1$V20,introns_1$V21,introns_1$V22,introns_1$V23,introns_1$V24)
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0,500000), names=c("1H","2H","3H","4H","5H","6H","adult_female","adult_male"),
               main = 'absolute gene expression among BUSCO maternal-degraded genes', ylab='transcripts/embryo')
