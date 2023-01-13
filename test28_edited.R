# Set working directory
setwd("/home/abayega/R/tests/test28")
exonsf_all <- read.table("exons_mandalorion_editted",header=F, row=1, sep = '\t')
intronsf_all <- read.table("introns_mandalorion_editted",header=F, row=1, sep = '\t')

exonsf_md <- read.table("exons_mandalorion_editted.maternal_degraded_genes",header=F, row=1, sep = '\t')
intronsf_md <- read.table("introns_mandalorion_editted.maternal_degraded_genes",header=T, row=1, sep = '\t')

exonsf <- read.table("exons_mandalorion_editted.zygotic_early_genes",header=F, row=1, sep = '\t')
intronsf <- read.table("introns_mandalorion_editted.zygotic_early_genes",header=F, row=1, sep = '\t')

expZyg <- read.table('allzygoticEarlygenes_expandedList',header=F)
expZyg <- read.table('allzygoticEarlygenes',header=F)

colnames_forAllExonsIntrons <- c('scaffold',	'start',	'end',	'length',	'RPG10K_1H',	'RPG10K_2H',	'RPG10K_3H',	'RPG10K_4H',	'RPG10K_5H',
                       'RPG10K_6H',	'RPG10K_fem',	'RPG10K_mal',	'std_1H',	'std_2H',	'std_3H',	'std_4H',	'std_5H',	'std_6H',	'std_fem',	'std_mal')

colnames(exonsf) <- colnames_forAllExonsIntrons
colnames(intronsf) <- colnames_forAllExonsIntrons
n <- 5
intron_cov <- cbind(intronsf$RPG10K_1H,intronsf$RPG10K_2H,intronsf$RPG10K_3H,intronsf$RPG10K_4H,intronsf$RPG10K_5H,intronsf$RPG10K_6H,intronsf$RPG10K_fem,intronsf$RPG10K_mal)
boxplot.matrix(intron_cov, use.cols = TRUE, ylim = c(0.005,n))

exon_cov <- cbind(exonsf$RPG10K_1H,exonsf$RPG10K_2H,exonsf$RPG10K_3H,exonsf$RPG10K_4H,exonsf$RPG10K_5H,exonsf$RPG10K_6H,exonsf$RPG10K_fem,exonsf$RPG10K_mal)
boxplot.matrix(exon_cov, use.cols = TRUE, ylim = c(0.005,n))

library(ggplot2)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

#ploting gene expression among exons and intron
#160120
#Trying to expand zygotic genes
colnames(exonsf_all) <- colnames_forAllExonsIntrons
colnames(intronsf_all) <- colnames_forAllExonsIntrons

exonsf_exp <- subset(exonsf_all, rownames(exonsf_all) %in% expZyg$V1)
intronsf_exp <- subset(intronsf_all, rownames(intronsf_all) %in% expZyg$V1)

time.points <- c(rep(c("1H","2H","3H","4H","5H","6H","female","male"), each=(nrow(exonsf_exp)+nrow(intronsf_exp))))

ftr <- c( rep(rep(c("exon","intron"), c(nrow(exonsf_exp),nrow(intronsf_exp))), times = 8))
length(ftr)
head(ftr)
ftr_mat <- c(exonsf_exp$RPG10K_1H,intronsf_exp$RPG10K_1H,
             exonsf_exp$RPG10K_2H,intronsf_exp$RPG10K_2H,
             exonsf_exp$RPG10K_3H,intronsf_exp$RPG10K_3H,
             exonsf_exp$RPG10K_4H,intronsf_exp$RPG10K_4H,
             exonsf_exp$RPG10K_5H,intronsf_exp$RPG10K_5H,
             exonsf_exp$RPG10K_6H,intronsf_exp$RPG10K_6H,
             exonsf_exp$RPG10K_fem,intronsf_exp$RPG10K_fem,
             exonsf_exp$RPG10K_mal,intronsf_exp$RPG10K_mal)

cov.data=data.frame(time.points,ftr,ftr_mat)
cov.data$ftr <- factor(cov.data$ftr, levels = c('intron','exon'),ordered = TRUE)
head(cov.data)
cov.data2 <- cov.data[cov.data$ftr_mat<2,]
ggplot(cov.data, aes(x=time.points, y=ftr_mat, fill=ftr)) + 
  geom_boxplot() +
  ggtitle("Assessing gene expression in exons and introns among zygotic early genes") +
  xlab("timepoints") +
  ylab("RPG10K") +
  ylim(c(0,1)) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00")) + #999999=blue,#56B4E9=grey,#E69F00=yellow not sure though
  scale_fill_manual(values=c("#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999")) +
  my.theme

write.csv(cov.data)

###What I can say is that using the whole expanded list of zygotic genes returns a much bigger list of geenes but the intron-exon coverage does not make a lot of sense
#This works beautifully.
##Now let us see how introns coverage looks across all genes
colnames(exonsf_all) <- c('scaffold',	'start',	'end',	'length',	'RPG10K_1H',	'RPG10K_2H',	'RPG10K_3H',	'RPG10K_4H',	'RPG10K_5H',
                      'RPG10K_6H',	'RPG10K_fem',	'RPG10K_mal',	'std_1H',	'std_2H',	'std_3H',	'std_4H',	'std_5H',	'std_6H',	'std_fem',	'std_mal')
colnames(intronsf_all) <- c('scaffold',	'start',	'end',	'length',	'RPG10K_1H',	'RPG10K_2H',	'RPG10K_3H',	'RPG10K_4H',	'RPG10K_5H',
                        'RPG10K_6H',	'RPG10K_fem',	'RPG10K_mal',	'std_1H',	'std_2H',	'std_3H',	'std_4H',	'std_5H',	'std_6H',	'std_fem',	'std_mal')

time.points <- c(rep(c("1H","2H","3H","4H","5H","6H","female","male"), each=(nrow(exonsf_all)+nrow(intronsf_all))))

ftr <- c( rep(rep(c("exon","intron"), c(nrow(exonsf_all),nrow(intronsf_all))), times = 8))

ftr_mat <- c(exonsf_all$RPG10K_1H,intronsf_all$RPG10K_1H,
             exonsf_all$RPG10K_2H,intronsf_all$RPG10K_2H,
             exonsf_all$RPG10K_3H,intronsf_all$RPG10K_3H,
             exonsf_all$RPG10K_4H,intronsf_all$RPG10K_4H,
             exonsf_all$RPG10K_5H,intronsf_all$RPG10K_5H,
             exonsf_all$RPG10K_6H,intronsf_all$RPG10K_6H,
             exonsf_all$RPG10K_fem,intronsf_all$RPG10K_fem,
             exonsf_all$RPG10K_mal,intronsf_all$RPG10K_mal)

cov.data=data.frame(time.points,ftr,ftr_mat)

ggplot(cov.data, aes(x=time.points, y=ftr_mat, fill=ftr)) + 
  geom_boxplot() +
  ggtitle("Assessing gene expression in exons and introns among ALL genes") +
  xlab("timepoints") +
  ylab("RPG10K") +
  ylim(c(0,2)) +
  my.theme

#Unfortunately, we don't see a trend among intron coverage. So, let us check maternal degraded genes.
colnames(exonsf_md) <- c('scaffold',	'start',	'end',	'length',	'RPG10K_1H',	'RPG10K_2H',	'RPG10K_3H',	'RPG10K_4H',	'RPG10K_5H',
                          'RPG10K_6H',	'RPG10K_fem',	'RPG10K_mal',	'std_1H',	'std_2H',	'std_3H',	'std_4H',	'std_5H',	'std_6H',	'std_fem',	'std_mal')
colnames(intronsf_md) <- c('scaffold',	'start',	'end',	'length',	'RPG10K_1H',	'RPG10K_2H',	'RPG10K_3H',	'RPG10K_4H',	'RPG10K_5H',
                            'RPG10K_6H',	'RPG10K_fem',	'RPG10K_mal',	'std_1H',	'std_2H',	'std_3H',	'std_4H',	'std_5H',	'std_6H',	'std_fem',	'std_mal')

time.points <- c(rep(c("1H","2H","3H","4H","5H","6H","female","male"), each=(nrow(exonsf_md)+nrow(intronsf_md))))

ftr <- c( rep(rep(c("exon","intron"), c(nrow(exonsf_md),nrow(intronsf_md))), times = 8))

ftr_mat <- c(exonsf_md$RPG10K_1H,intronsf_md$RPG10K_1H,
             exonsf_md$RPG10K_2H,intronsf_md$RPG10K_2H,
             exonsf_md$RPG10K_3H,intronsf_md$RPG10K_3H,
             exonsf_md$RPG10K_4H,intronsf_md$RPG10K_4H,
             exonsf_md$RPG10K_5H,intronsf_md$RPG10K_5H,
             exonsf_md$RPG10K_6H,intronsf_md$RPG10K_6H,
             exonsf_md$RPG10K_fem,intronsf_md$RPG10K_fem,
             exonsf_md$RPG10K_mal,intronsf_md$RPG10K_mal)

cov.data=data.frame(time.points,ftr,ftr_mat)

ggplot(cov.data, aes(x=time.points, y=ftr_mat, fill=ftr)) + 
  geom_boxplot() +
  ggtitle("Assessing gene expression in exons and introns among maternal degraded genes") +
  xlab("timepoints") +
  ylab("RPG10K") +
  ylim(c(0,2)) +
  my.theme
###



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
