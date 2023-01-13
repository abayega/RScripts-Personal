setwd('/home/abayega/R/tests/test25')
raw_count <- read.table("AllvAll_combined",header=T,row.names=1, sep = '\t')
head(raw_count)

#number of un expressed genes at 1-hour
length(raw_count$RPG10K_1H[which(raw_count$RPG10K_1H < 0.01)]) #6928
#number of un expressed genes at 2-hour
length(raw_count$RPG10K_2H[which(raw_count$RPG10K_2H < 0.01)]) #7429
#number of un expressed genes at 3-hour
length(raw_count$RPG10K_3H[which(raw_count$RPG10K_3H < 0.01)]) #7258
#number of un expressed genes at 4-hour
length(raw_count$RPG10K_4H[which(raw_count$RPG10K_4H < 0.01)]) #7040
#number of un expressed genes at 5-hour
length(raw_count$RPG10K_5H[which(raw_count$RPG10K_5H < 0.01)]) #6861
#number of un expressed genes at 6-hour
length(raw_count$RPG10K_6H[which(raw_count$RPG10K_6H < 0.01)]) #6877

#number of expressed genes at 1-hour
length(raw_count$RPG10K_1H[which(raw_count$RPG10K_1H >= 0.01)]) #11043
#number of expressed genes at 2-hour
length(raw_count$RPG10K_2H[which(raw_count$RPG10K_2H >= 0.01)]) #10542
#number of expressed genes at 3-hour
length(raw_count$RPG10K_3H[which(raw_count$RPG10K_3H >= 0.01)]) #10713
#number of expressed genes at 4-hour
length(raw_count$RPG10K_4H[which(raw_count$RPG10K_4H >= 0.01)]) #10931
#number of expressed genes at 5-hour
length(raw_count$RPG10K_5H[which(raw_count$RPG10K_5H >= 0.01)]) #11110
#number of expressed genes at 6-hour
length(raw_count$RPG10K_6H[which(raw_count$RPG10K_6H >= 0.01)]) #11094

head(rownames(raw_count)[which(raw_count$RPG10K_1H >= 0.01 && raw_count$RPG10K_2H < 0.01)])

#I am trying to clean up the table
x1_6x <- subset(raw_count, select = c(RPG10K_1H,RPG10K_2H,RPG10K_3H,RPG10K_4H,RPG10K_5H,RPG10K_6H))
x1_6x[x1_6x < 0.01] = 0
raw_count2 <- cbind(raw_count[,1:15], x1_6x, raw_count[,22:(ncol(raw_count))])

#row_sub = apply(x1_6x, 1, function(row) all(row !=0 ))
#x1_6x <- x1_6x[row_sub,]

raw_count2 <- subset(raw_count2, RPG10K_1H >= 0.01 | RPG10K_2H >= 0.01 | RPG10K_3H >= 0.01 | RPG10K_4H >= 0.01 | RPG10K_5H >= 0.01 | RPG10K_6H >= 0.01)
raw_count3 <- subset(raw_count2, RPG10K_1H < 181.8971	& RPG10K_2H < 365.4984	& RPG10K_3H < 161.0578 & RPG10K_4H < 175.6243 & RPG10K_5H < 175.4735 & RPG10K_6H < 206.1444)


#row_sub = apply(a, 1, function(row) all(row !=0 ))
#a <- a[row_sub,]
#raw_count_clean[raw_count_clean < 0.2] = 0

#Assuming that RNA degradation is not specific but random, let us see the expression of genes accross timepoints
k = cbind(log10(raw_count3$TPE_1H[raw_count3$TPE_1H>0]), log10(raw_count3$TPE_2H[raw_count3$TPE_2H>0]), log10(raw_count3$TPE_3H[raw_count3$TPE_3H>0]), 
          log10(raw_count3$TPE_4H[raw_count3$TPE_4H>0]), log10(raw_count3$TPE_5H[raw_count3$TPE_5H>0]), log10(raw_count3$TPE_6H[raw_count3$TPE_6H>0]))
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green', 'blue', 'gray','yellow','purple'), names=c("1H","2H","3H","4H","5H","6H"),
               main = 'Comparison of expression levels accross timepoints', ylab="log10 (transcripts per embryo)")
na <- (raw_count3$TPE_1H[raw_count3$TPE_1H>0])
nn <- (raw_count3$TPE_2H[raw_count3$TPE_2H>0])
vioplot(na,nn, col = c('red', 'green'), names=c("1H","2H"), ylim = c(0,1000000))
#So the story is this; globally, the genes at 1 hpo have a high expression. They are degraded at 2 hpo and this seems global. If this is a stochastic process then the most highly expressed
#genes should have a higher effect. So, we use GFOLD to determine maternal transcripts and see if they have a higher gene expression

#Try ggplot
library(ggplot2)
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization
# Basic violin plot

expn <- c(log10(raw_count3$TPE_1H[raw_count3$TPE_1H>0]), log10(raw_count3$TPE_2H[raw_count3$TPE_2H>0]), log10(raw_count3$TPE_3H[raw_count3$TPE_3H>0]), 
              log10(raw_count3$TPE_4H[raw_count3$TPE_4H>0]), log10(raw_count3$TPE_5H[raw_count3$TPE_5H>0]), log10(raw_count3$TPE_6H[raw_count3$TPE_6H>0]))
fac <- c(rep('1H',length(log10(raw_count3$TPE_1H[raw_count3$TPE_1H>0]))), rep('2H', length(log10(raw_count3$TPE_2H[raw_count3$TPE_2H>0]))),
         rep('3H', length(log10(raw_count3$TPE_3H[raw_count3$TPE_3H>0]))), rep('4H', length(log10(raw_count3$TPE_4H[raw_count3$TPE_4H>0]))),
         rep('5H', length(log10(raw_count3$TPE_5H[raw_count3$TPE_5H>0]))), rep('6H', length(log10(raw_count3$TPE_6H[raw_count3$TPE_6H>0])))
         )
expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('expn', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)

p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "time-points") +
  labs(fill = 'timepoints') +
  scale_y_continuous(expression(log[10]*" (Transcripts per embryo)")) + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

# violin plot with median points
p + stat_summary(fun.y=median, geom="point", size=2, color="red")
p + geom_boxplot(width=0.1)
p + stat_summary(fun.data="mean_sdl", mult=1, geom="crossbar", width=0.2 )

#Let us try to determine if maternal degraded transcripts have a higher expression compared to other genes at 1 hpo
maternal_degraded_genes <- subset(raw_count3,RPG10K_1H != 0 & X1v2 =="down",select = c(TPE_1H))
not_maternal_degraded_genes <- subset(raw_count3,RPG10K_1H != 0 & X1v2 !="down",select = c(TPE_1H))

par(mfrow=c(1,1))
k = cbind(log10(maternal_degraded_genes$TPE_1H), head(log10(not_maternal_degraded_genes$TPE_1H), 1497))
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green'), names=c("maternal degraded genes","other genes"), 
               main = 'comparison of expression levels of maternal degraded genes versus other genes', ylab="log10 (transcripts per embryo)")
#So, indeed maternal genes have a higher expression. at 1 hpo. We now determine if indeed they were cleared

#140120
#Trying to see if the differences are statistically significant
length(log10(maternal_degraded_genes$TPE_1H))
length(head(log10(not_maternal_degraded_genes$TPE_1H), 1497))

#wilcox.test(log10(maternal_degraded_genes$TPE_1H),head(log10(not_maternal_degraded_genes$TPE_1H), 1497),paired = FALSE,conf.level = 0.99)
#t.test(log10(maternal_degraded_genes$TPE_1H),head(log10(not_maternal_degraded_genes$TPE_1H), 1497),paired = FALSE,conf.level = 0.99)

#wilcox.test(maternal_degraded_genes$TPE_1H,head(not_maternal_degraded_genes$TPE_1H, 1497),paired = FALSE,conf.level = 0.99)

#Data is not normally distributed so we shall use the wilcox.test
shapiro.test(head(not_maternal_degraded_genes$TPE_1H, 1497))
wilcox.test(maternal_degraded_genes$TPE_1H,head(not_maternal_degraded_genes$TPE_1H, 1497),paired = FALSE,conf.level = 0.99)


#using ggplot2
expn <- c(log10(maternal_degraded_genes$TPE_1H), log10(not_maternal_degraded_genes$TPE_1H))
fac <- c(rep('maternal_degraded_genes', length(log10(maternal_degraded_genes$TPE_1H))), rep('others', length(log10(not_maternal_degraded_genes$TPE_1H)))
         )
expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('expn', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)
p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  labs(x = "time-points") +
  labs(fill = '') +
  scale_y_continuous(expression(log[10]*" (Transcripts per embryo)")) + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#rownames(raw_count2)[which(raw_count2$TPE_1H == max(raw_count2$TPE_1H))]

#Let us confirm that the expression levels of maternal transcripts at 2 hpo are down
maternal_degraded_genes <- subset(raw_count3,RPG10K_1H != 0 & X1v2 =="down",select = c(TPE_2H))
not_maternal_degraded_genes <- subset(raw_count3,RPG10K_2H != 0 & RPG10K_1H != 0 & X1v2 !="down",select = c(TPE_2H))
k = cbind(log10(maternal_degraded_genes$TPE_2H), head(log10(not_maternal_degraded_genes$TPE_2H), 1497))
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green'), names=c("maternal degraded genes (2H)","other genes (2H)"), 
               main = 'comparison of expression levels of maternal degraded genes versus other genes at 2 hpo', ylab="log10 (transcripts per embryo)")
#And yes, at 2 hpo the maternal transcripts are degraded to levels similar to other genes. Now, the question is how to the genes get expressed to the levels we see at 3 hpo? 
#Are the maternal genes just remade again? (this would be so stupid). Is it new genes?

#140120
#Trying to see if the differences are statistically significant
length(log10(maternal_degraded_genes$TPE_2H))
length(head(log10(not_maternal_degraded_genes$TPE_2H), 1497))

#wilcox.test(log10(maternal_degraded_genes$TPE_2H),head(log10(not_maternal_degraded_genes$TPE_2H), 1497),paired = FALSE,conf.level = 0.99)
#t.test(log10(maternal_degraded_genes$TPE_2H),head(log10(not_maternal_degraded_genes$TPE_2H), 1497),paired = FALSE,conf.level = 0.99)

#wilcox.test(maternal_degraded_genes$TPE_1H,head(not_maternal_degraded_genes$TPE_1H, 1497),paired = FALSE,conf.level = 0.99)
#t.test(maternal_degraded_genes$TPE_1H,head(not_maternal_degraded_genes$TPE_1H, 1497),paired = FALSE,conf.level = 0.99)

#Data is normally distributed so we shall use the t.test
shapiro.test(head(not_maternal_degraded_genes$TPE_2H, 1497))
wilcox.test(maternal_degraded_genes$TPE_2H,head(not_maternal_degraded_genes$TPE_2H, 1497),paired = FALSE,conf.level = 0.99)

#using ggplot2
expn <- c(log10(maternal_degraded_genes$TPE_2H), log10(not_maternal_degraded_genes$TPE_2H))
fac <- c(rep('maternal_degraded_genes', length(log10(maternal_degraded_genes$TPE_2H))), rep('others', length(log10(not_maternal_degraded_genes$TPE_2H))))

expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('expn', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)
p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  labs(x = "time-points") +
  labs(fill = '') +
  scale_y_continuous(expression(log[10]*" (Transcripts per embryo)")) + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#First, let us see the expression levels of maternal degraded transcripts at 3 hpo to see how they look, my gues is that they are not that high
maternal_degraded_genes <- subset(raw_count3,RPG10K_3H != 0 & RPG10K_1H != 0 & X1v2 =="down",select = c(TPE_3H))
not_maternal_degraded_genes <- subset(raw_count3,RPG10K_3H != 0 & RPG10K_1H != 0 & X1v2 !="down",select = c(TPE_3H))
k = cbind(log10(maternal_degraded_genes$TPE_3H), log10(not_maternal_degraded_genes$TPE_3H))
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green'), names=c("maternal degraded genes (3H)","other genes (3H)"), 
               main = 'xcomparison of expression levels of maternal degraded genes versus other genes at 3 hpo', ylab="log10 (transcripts per embryo)")
#And I was right, the maternal genes that were formerly degraded at 2 hpo are not synthesised at levels above othe genes, they are very similar


#To see if it is new genes (zygotic), we find genes how many were absent at 1 and 2 hpo but present at 3 hpo and their expression levels
raw_count4 <- subset(raw_count3, RPG10K_1H < 0.01 & RPG10K_2H < 0.01 & RPG10K_3H >= 0.01)
raw_count4 <- subset(raw_count4, RPG10K_1H < 181.8971	& RPG10K_2H < 365.4984	& RPG10K_3H < 161.0578, select = c(TPE_3H,gene_length,gene_gc,transcript_length,transcript_gc))
nrow(raw_count4)
#So, only 203 genes meet this criteria. Now we compare their expression compared to the other samples, my gues is that the embryo needs them so it made enough of them although they have to be short genes
raw_count5 <- subset(raw_count3, RPG10K_1H >= 0.01 | RPG10K_2H >= 0.01 & RPG10K_3H >= 0.01)
raw_count5 <- subset(raw_count5, RPG10K_1H < 181.8971	& RPG10K_2H < 365.4984	& RPG10K_3H < 161.0578, select = c(TPE_3H,gene_length,gene_gc,transcript_length,transcript_gc))

k = cbind(log10(raw_count4$TPE_3H),log10(raw_count5$TPE_3H))
par(mgp=c(2.4,1,0),mar=c(4,4,4,2)+0.1)
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green'), names=c("zygotic genes (3H)","other genes (3H)"), cex.lab=1.5, cex.names = 1.5, cex.axis = 1.5,
               main = 'comparison of expression levels of zygotic genes versus other genes at 3 hpo', ylab="log10 (transcripts per embryo)")
#So, the zygotic genes actually have a much lower expression coz probably they are spatially expressed, and it's just starting to express them. But how long are these genes compared to the rest?

#using ggplot2
expn <- c(log10(raw_count4$TPE_3H), log10(raw_count5$TPE_3H))
fac <- c(rep('maternal_degraded_genes', length(log10(raw_count4$TPE_3H))), rep('others', length(log10(raw_count5$TPE_3H))))

expn_mat <- data.frame(expn, fac)
colnames(expn_mat) <- c('expn', 'fac')
expn_mat$fac <- as.factor(expn_mat$fac)
head(expn_mat)
p <- ggplot(expn_mat, aes(x=fac, y=expn, fill=fac)) + #, color=fac)
  geom_violin() +
  geom_boxplot(width=0.1, fill="white") +
  labs(x = "time-points") +
  labs(fill = '') +
  scale_y_continuous(expression(log[10]*" (Transcripts per embryo)")) + 
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=16, face="bold"),
    axis.text.y = element_text(size=16, face="bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size=14, face="bold")
  )
p

k = cbind(raw_count4$gene_length,raw_count4$transcript_length,raw_count5$gene_length,raw_count5$transcript_length)
k = cbind(log(raw_count4$gene_length),log(raw_count4$transcript_length),log(raw_count5$gene_length),log(raw_count5$transcript_length))
boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green', 'blue', 'gray'), names=c('zygote gene length','zygote transcript length','gene length (other genes)','transcript length (other genes)'),
               main = 'comparison of gene lengths of zygotic genes versus other genes', ylab="log( length (bp) )")

k = cbind(raw_count4$gene_gc,raw_count4$transcript_gc,raw_count5$gene_gc,raw_count5$transcript_gc)
boxplot.matrix(k,use.cols = TRUE, col = c('yellow','purple','orange','green'), names=c('zygote genes','zygote transcripts','other genes','other transcripts)'),
               main = 'comparison of gene GC of zygotic genes versus other genes', ylab=" GC content (%)")
#Okay, this is a bit bizzarre, there doesn't seem to be any visible differences in the length of genes, transcripts, or GC content between zygotic and other genes

#So, this is the idea, it seems to me that the embryo degrades maternal transcripts and then restarts expression of the same genes but this time the expression is not dominated by maternal transcripts
#it now controls who dominates the expression, and probably the area.
#so, probably it is important to determine the most highly expressed genes and see if these are enriched in anything,



#col = c('red', 'green', 'blue', 'gray','yellow','purple','orange','yellow','yellow'), names=c("ONT_raw","ONT_canu","ONT_canu2","ONT_ill","Ill",'y-raw','y_canu','y_canu_lsc','y_canu_lsc2')
#k = cbind(log10(sort(raw_count3$TPE_1H[raw_count3$TPE_1H>0], decreasing = T)[1:200],log10(sort(raw_count3$TPE_1H[raw_count3$TPE_1H>0], decreasing = T)[201:length(raw_count3$TPE_1H)])))
#boxplot.matrix(k,use.cols = TRUE, col = c('red', 'green'), names=c("200 most highly expressed genes","2H",
          
high_expn <- subset(raw_count3, select = c(TPE_3H))
head(high_expn)
#high_expn <- data.frame(high_expn)
#class(high_expn)
#high_expn <- high_expn[order(high_expn$TPE_1H,decreasing = T),]
#high_expn <- high_expn[with(high_expn, order(TPE_1H)), ]

write.table(high_expn, file = '3_hour_TPE', quote = F, sep = '\t')


#Now I will try gene set enrichment using 


##ill and ONT foldchange corr
setwd('/home/abayega/R/tests/test25')
ont<-read.table("ont_5v6.diff", sep = '\t', header = T,row.names = 1)
ill<-read.table("ill_Bo.E.5Hv6H.diff", sep = '\t', header = T,row.names = 1)

plot(ill$GFOLD.0.01.,ont$GFOLD.0.01., cex=0.25,pch = 20)
abline(a=0,b=1, col = 'red')

dfm <- cbind(ill$GFOLD.0.01.,ont$GFOLD.0.01.)
colnames(dfm) <- c('ill','ont')
dfm <- data.frame(dfm)

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(ont ~ ill, df1);
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
                     r2 = format(cor(df1$ill,df1$ont, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = dfm, aes(x = ill, y = ont)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of illumina and ONT Gfold") +
  scale_x_continuous("illumina Gfold") +
  scale_y_continuous(expression("ONT Gfold")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 0.7, y = 3, label = lm_eqn(dfm), parse = T, color = "red",size=6) +
  annotate(geom = "text", x = 0.7, y = 4, label = lm_eqn2(dfm), parse = T, color = "red",size=6) +
  my.theme

k <- cbind(ill$X1stRPKM,ill$X2ndRPKM,ont$X1stRPKM,ont$X2ndRPKM)
colnames(k) <- c('ill_5h_fpkm','ill_6h_fpkm','ONT_5h_fpkm','ONT_6h_fpkm')

k <- data.matrix(k)
boxplot.matrix(k, use.cols = T, ylim=c(0,10))

#Scatter of FPKMs
k <- data.frame(k)
plot(k$ill_5h_fpkm,k$ill_6h_fpkm, cex=0.25,pch = 20, xlim = c(0,1000), ylim = c(0,1000))
abline(a=0,b=1, col = 'red')
cor(k$ill_5h_fpkm,k$ill_6h_fpkm)


plot(k$ONT_5h_fpkm,k$ONT_6h_fpkm, cex=0.25,pch = 20, xlim = c(0,1000), ylim = c(0,1000))
abline(a=0,b=1, col = 'red')


?maPlot
install.packages("edgeR")

source("https://bioconductor.org/biocLite.R")
biocLite("edgeR")

library("edgeR")
maPlot(k$ill_5h_fpkm,k$ill_6h_fpkm, smooth.scatter=T, lowess=T)
maPlot(k$ONT_5h_fpkm,k$ONT_6h_fpkm, smooth.scatter=T, lowess=T)


install.packages("vioplot")
library("vioplot")


##13/01/2020
#Trying to find maternal genes reorganization

abs_exprn_1H = subset(raw_count, RPG10K_1H >= 0.01,select=c(TPE_1H))
abs_exprn_2H = subset(raw_count, RPG10K_2H >= 0.01,select=c(TPE_2H))
abs_exprn_3H = subset(raw_count, RPG10K_3H >= 0.01,select=c(TPE_3H))

first_3_hrs <- subset(raw_count,select=c(RPG10K_1H,RPG10K_2H,RPG10K_3H,TPE_1H,TPE_2H,TPE_3H))

write.table(first_3_hrs, file="first_3_hrs.txt", quote = F, sep = '\t')

raw_count <- read.table("first_3_hrs.edited2.txt",header=F,row.names=1, sep = '\t')
head(raw_count)

colnames(raw_count) <- c('RPG10K_1H','RPG10K_2H','RPG10K_3H','TPE_1H','TPE_2H','TPE_3H')
raw_count <- raw_count[73:nrow(raw_count),]


#Now let us the percentage destabilisation
raw_count$perc_destn <- 100*(raw_count$TPE_1H-raw_count$TPE_2H)/raw_count$TPE_1H

#Now let us get the relative expression of 3 hrs AEL to 1 hr AEL
raw_count$rel3to1 <- (raw_count$TPE_3H/raw_count$TPE_1H)

#Try to remove inf and NaN
raw_count$remove <- (raw_count$rel3to1 + raw_count$perc_destn)
raw_count_clean <- raw_count[!is.na(raw_count$remove),]

dim(raw_count_clean_sorted)
#Plotting 
dd <- raw_count_clean
raw_count_clean_sorted <- dd[with(dd, order(-perc_destn)), ]

par(mar = c(5,7,4,2) + 0.1)
plot(seq(1,length(raw_count_clean_sorted$rel3to1),by=1),raw_count_clean_sorted$rel3to1,type="p",cex.axis = 1.5,cex.lab = 1.5,
     main = "Assessing reorganization of maternal transcripts in developing B. oleae embryo",
     ylab="Transcripts per embryo at\n 3 hrs AEL relative to 1 hr AEL", xlab="Genes sorted: Most downregulated to least downregulated at 2 hrs AEL)")
#points(seq(1,length(raw_count_clean_sorted$rel3to1),by=1),raw_count_clean_sorted$rel3to1, cex = .5, col = "dark red")
abline(lm(raw_count_clean_sorted$rel3to1~seq(1,length(raw_count_clean_sorted$rel3to1))) , col="blue")

plot(seq(1,length(raw_count_clean_sorted$rel3to1),by=1),raw_count_clean_sorted$rel3to1,type="p",cex.axis = 1.5,cex.lab = 1.5,xaxt="none", yaxt="none", xlab="",ylab="",
     main = "Assessing reorganization of maternal transcripts in developing B. oleae embryo")
axis(1, font = 2,cex.axis=1.5)
axis(2, font = 2,cex.axis=1.5)
mtext(side=1, line=2, "Genes sorted: Most downregulated to least downregulated at 2 hrs AEL", font=2,cex=1.5) #col='blue'
mtext(side=2, line=3, "Transcripts per embryo at\n 3 hrs AEL relative to 1 hr AEL", font=2, cex=1.5) #col='blue'
abline(lm(raw_count_clean_sorted$rel3to1~seq(1,length(raw_count_clean_sorted$rel3to1))) , col="blue")

plot(raw_count_clean_sorted$rel3to1[(length(raw_count_clean_sorted$rel3to1)-1000):length(raw_count_clean_sorted$rel3to1)]) #, yrange, type="n", xlab="Age (days)", ylab="Circumference (mm)" ) 
lines(raw_count_clean_sorted$rel3to1, fitted(fit), col="blue")
barplot(raw_count_clean_sorted$TPE_1H, ylim = c(0,10000000))
barplot(log(raw_count_clean_sorted$TPE_3H, 2), ylim = c(0,20))

