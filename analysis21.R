rm(list = ls())
setwd("~/R/tests/test21")
a<-read.table('ERCC_expected_and_RPG10k.txt', header = T, row.names = 1, sep = '\t')
a<-read.table('rsem_gene_counts.txt', header = T, row.names = 1, sep = '\t')

head(a)

row_sub = apply(a, 1, function(row) all(row !=0 ))
a <- a[row_sub,]

#erccs <- a
#erccs <- a[1:53,]

tail(a)
tail(erccs)

plot(log(a$X1H_exp),log(a$X1H_TPM),xlab="log10(number_of_ERCC_molecules)",ylab="log10(number_of_sequenced_reads)")
glm(a$X1H_TPM ~ offset(log(a$X1H_exp)), family=poisson(link=log))
new_y <- log(a$X1H_exp) - 13.99
lines(log(a$X1H_exp),new_y,col="red")

glm(a$X2H_TPM ~ offset(log(a$X2H_exp)), family=poisson(link=log))
glm(a$X3H_TPM ~ offset(log(a$X3H_exp)), family=poisson(link=log))
glm(a$X4H_TPM ~ offset(log(a$X4H_exp)), family=poisson(link=log))
glm(a$X5H_TPM ~ offset(log(a$X5H_exp)), family=poisson(link=log))
glm(a$X6H_TPM ~ offset(log(a$X6H_exp)), family=poisson(link=log))
glm(a$fem.1 ~ offset(log(a$fem)), family=poisson(link=log))
glm(a$mal.1 ~ offset(log(a$mal)), family=poisson(link=log))

#Getting beta for Illumina rsem quantified TPMs
tpm_5H <- c(a$X5H_TPM_rsem)
tpm_6H <- c(a$X6H_TPM_rsem)
ercc_exp <- c(a$X5H_exp)

tpm_5H_df <- as.data.frame(cbind(ercc_exp, tpm_5H))
tpm_6H_df <- as.data.frame(cbind(ercc_exp, tpm_6H))

row_sub = apply(tpm_5H_df, 1, function(row) all(row !=0 )) #Remove rows with zeros
tpm_5H_df <- tpm_5H_df[row_sub,]
row_sub = apply(tpm_6H_df, 1, function(row) all(row !=0 )) #Remove rows with zeros
tpm_6H_df <- tpm_6H_df[row_sub,]

#5H beta
glm(tpm_5H_df$tpm_5H ~ offset(log(tpm_5H_df$ercc_exp)), family=poisson(link=log))
#6H beta
glm(tpm_6H_df$tpm_6H ~ offset(log(tpm_6H_df$ercc_exp)), family=poisson(link=log))

#Plot correlations per sample
plot(log2(a$X1H_exp),log2(a$X1H_TPM), col=rgb(1,0,0,0.5), pch=15,xlab="Log2(Expected No. molecules)",ylab="Log2(Pooled sample ERCC TPM)")
par(new=TRUE)
plot(log2(a$X2H_exp),log2(a$X2H_TPM), col=rgb(0,1,0,0.5), pch=16, axes = False )
par(new=TRUE)
plot(log2(a$X3H_exp),log2(a$X3H_TPM), col=rgb(0,0,1,0.5), pch=17, axes = False )
par(new=TRUE)
plot(log2(a$X4H_exp),log2(a$X4H_TPM), col=rgb(0.5,0.5,0.5,0.5), pch=18, axes = False )
par(new=TRUE)
plot(log2(a$X5H_exp),log2(a$X5H_TPM), col=rgb(1,0.9,0,0.5), pch=19, axes = False )
par(new=TRUE)
plot(log2(a$X6H_exp),log2(a$X6H_TPM), col=rgb(1,0,1,0.5), pch=20, axes = False )
legend("topleft", legend=c("Bo_1H","Bo_2H","Bo_3H","Bo_4H","Bo_5H","Bo_6H"), pch = c(15,16,17,18,19,20), col=c("red","green","blue","gray","yellow","purple"))
dev.off()

#Plot all ERCC data with natural log
plot(log(all_ercc_exp),log(all_ercc_tpm),xlab="ln(Expected No. molecules)",ylab="ln(Pooled sample ERCC TPM)")
glm(all_ercc_tpm ~ offset(log(all_ercc_exp)), family=poisson(link=log))
new_y <- log(all_ercc_exp) - 13.54
lines(log(all_ercc_exp),new_y,col="red")

#Get Pearson correlation
cor(all_ercc_exp,all_ercc_tpm, method=c("pearson"))
cor(log(all_ercc_exp),log(all_ercc_tpm), method=c("pearson"))

p <- log(all_ercc_tpm)
q <- log(all_ercc_exp)
number<-match(c("-Inf"),log(all_ercc_tpm))
negative_number<--1*number
p <- p[negative_number]
q <- q[negative_number]
cor(p,q,method=c("pearson"))

lm(all_ercc_tpm ~ all_ercc_exp)
all_ercc_tpm = 1.153 + 1.225e-06(all_ercc_exp)

#Try to plot TPM vs time
plot(c(1,2,3,4,5,6),a[1,9:14], lty=c(1), lwd=c(3), col=c("green"))
lines(c(1,2,3,4,5,6),a[1,9:14], lty=c(1), lwd=c(3), col=c("red"))


f <- read.table("Bo_E_1H_C010_10_pass_edited.trimmed_stranded_edited_readlengths.txt", header = F)
head(f)
boxplot(f$V1)
head(a)



#Restarting new analysis of ERCC to make them absolute counts
library(ggplot2)
###This is old code, we are doing new code below this one
plot(log(a$X1H_exp),log(a$X1H_TPM),xlab="log10(Expected_number_of_ERCC_molecules)",ylab="log10(Observed_normalised_number_of_ERCC)")

all_ercc_tpm <- c(a$X1H_TPM,a$X2H_TPM,a$X3H_TPM,a$X4H_TPM,a$X5H_TPM,a$X6H_TPM)
all_ercc_exp <- c(a$X1H_exp,a$X2H_exp,a$X3H_exp,a$X4H_exp,a$X5H_exp,a$X6H_exp)
all_ercc <- as.data.frame(cbind(all_ercc_tpm, all_ercc_exp))
head(all_ercc)

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(log2(all_ercc_tpm) ~ log2(all_ercc_exp), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = all_ercc, aes(x = log2(all_ercc_exp), y = log2(all_ercc_tpm))) +
  geom_point() +
  geom_smooth(method = "lm") +
  scale_x_continuous("log2(Expected_number_of_ERCC_molecules)") +
  scale_y_continuous("log2(Observed_ERCC)") +
  my.theme

ggplot(data = erccs, aes(x = log2(X1H_exp), y = log2(X1H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("log2(Expected No. off ERCC molecules)") +
  scale_y_continuous("log2(Observed ERCC)") + #, limits = c(1,6)) +
  annotate(geom = "text", x = 20, y = 7, label = lm_eqn(erccs), parse = T, color = "red") +
  my.theme

#This is the plot used for the PhD committee slides
ggplot(data = all_ercc, aes(x = log2(all_ercc_exp), y = log2(all_ercc_tpm))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("log2(Expected No. of ERCC molecules)") +
  scale_y_continuous("log2(Observed ERCC)") + #, limits = c(1,6)) +
  annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#Plotting ERCC length
ercc_ratios<-read.table('ERCC_ratios_full_edited', header = T, sep = '\t')
head(ercc_ratios)

ggplot(data = ercc_ratios, aes(x = ercc_length, y = len_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous("ratio Observed:Expected ERCC length") + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 10, y = 0.5, label = lm_eqn(ercc_ratios), parse = T, color = "red") +
  my.theme
### Old code ends here

#Ploting correlation of ONTs

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(log2(X6H_TPM) ~ log2(X6H_exp), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
#Ploting correlation of ONTs 5H
ggplot(data = a, aes(x = log2(X5H_exp), y = log2(X5H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation exp vs RPG10 Bo_E_5H ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(a), parse = T, color = "red") +
  my.theme
#Ploting correlation of ONTs 6H
ggplot(data = a, aes(x = log2(X6H_exp), y = log2(X6H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation exp vs RPG10 Bo_E_6H ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(a), parse = T, color = "red") +
  my.theme

#biases in 5H ONT
a$exp_ratio_5H <- log2(a$X5H_abs)/log2(a$X5H_exp)
ggplot(data = a, aes(x = ercc_gc, y = exp_ratio_5H)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with GC 5H ONT") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  my.theme

#Plotting ERCC length
ggplot(data = a, aes(x = ercc_length, y = exp_ratio_5H)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with length 5H ONT") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  my.theme


#Plot length and GC ratios but for expression ###This code is probably wrong. Use other code down the page or search for Correlation of abs ratios with length 6H ONT
all_ercc_abs <- c(a$X1H_abs,a$X2H_abs,a$X3H_abs,a$X4H_abs,a$X5H_abs,a$X6H_abs)
all_ercc_exp <- c(a$X1H_exp,a$X2H_exp,a$X3H_exp,a$X4H_exp,a$X5H_exp,a$X6H_exp)
all_ercc2 <- as.data.frame(cbind(all_ercc_abs, all_ercc_exp))
all_ercc2$ratio <- all_ercc_abs/all_ercc_exp
head(all_ercc2)
all_ercc2$len <- a$ercc_length
all_ercc2$GC <- a$ercc_gc
head(all_ercc2)

lm_eqn = function(df1){
  m = lm(log2(ratio) ~ GC, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = all_ercc2, aes(x = GC, y = log10(ratio))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous("log2(Observed:Expected)") + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

ggplot(data = all_ercc2, aes(x = len, y = log10(ratio))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous("log2(Observed:Expected)") + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme
###

#Printing read lengths of male contigs
male<-read.table('male_contigs_readlength', header = F, sep = '\t')
boxplot(male$V1, ylab = "scaffold_length")
hist(male$V1, breaks = 50, ylim = c(0,100))
sum(male$V1)


#For Illumina correlations and sensitivity and LLD (lower limit of detection)
#Plot correlations per sample
library(ggplot2)

#Doing correlations for 5 Hour timepoint
head(tpm_5H_df)
my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(log2(tpm_5H) ~ log2(ercc_exp), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = tpm_5H_df, aes(x = log2(ercc_exp), y = log2(tpm_5H))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("correlation exp vs tpm Bo_E_5H") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" TPM")) +
  annotate(geom = "text", x = 12, y = 9, label = lm_eqn(tpm_5H_df), parse = T, color = "red") +
  my.theme

#Doing correlations for 6 Hour timepoint
head(tpm_6H_df)

lm_eqn = function(df1){
  m = lm(log2(tpm_6H) ~ log2(ercc_exp), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = tpm_6H_df, aes(x = log2(ercc_exp), y = log2(tpm_6H))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("correlation exp vs tpm Bo_E_6H") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" TPM")) +
  annotate(geom = "text", x = 12, y = 9, label = lm_eqn(tpm_6H_df), parse = T, color = "red") +
  my.theme

#Plot length and GC ratios for 5H Illumina timepoint
bo_5H <-read.table('rsem_gene_counts.txt', header = T, row.names = 1, sep = '\t')
head(bo_5H)
bo_5H$exp_ratio <- log2(bo_5H$X5H_abs)/log2(bo_5H$X5H_exp)

bo_5H <- as.data.frame(bo_5H)

row_sub = apply(bo_5H, 1, function(row) all(row !=0 ))
bo_5H <- bo_5H[row_sub,]

lm_eqn = function(df1){
  m = lm(exp_ratio ~ ercc_gc, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = bo_5H, aes(x = ercc_gc, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ratios with GC") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

ggplot(data = bo_5H, aes(x = ercc_length, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ratios with length") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#Plot length and GC ratios for 6H Illumina timepoint
bo_6H <-read.table('rsem_gene_counts.txt', header = T, row.names = 1, sep = '\t')
head(bo_6H)
bo_6H$exp_ratio <- log2(bo_6H$X6H_abs)/log2(bo_6H$X5H_exp)

bo_6H <- as.data.frame(bo_6H)

row_sub = apply(bo_6H, 1, function(row) all(row !=0 ))
bo_6H <- bo_6H[row_sub,]

lm_eqn = function(df1){
  m = lm(exp_ratio ~ ercc_gc, df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = bo_6H, aes(x = ercc_gc, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ratios with GC 6H") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

ggplot(data = bo_6H, aes(x = ercc_length, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ratios with length 6H") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

bo_counts <-read.table('rsem_gene_read_counts', header = T, sep = '\t')
head(bo_counts)
bo_counts$count_ratio <- (bo_counts$X5H_count)/(bo_counts$X6H_count)

row_sub = apply(bo_counts, 1, function(row) all(row !=0 ))
bo_counts <- bo_counts[row_sub,]

ggplot(data = bo_counts, aes(x = read_len, y = count_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of count ratios with length 5-6H") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression("Counts ratio (5H/6H)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme
max()

###We are tyring to see correlation of ont and illumina quantification
gene_exp <- read.table('ont_illumina_correlation', header = T, row.names = 1, sep = '\t')
head(gene_exp)
plot(gene_exp$X5H_RPG10,gene_exp$X5H_TPM)

add_1 <- apply()

gene_exp <- gene_exp + 1

#ONT RPG10 and Illumina TPM 5H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                        r2 = format(cor(log(gene_exp$X5H_RPG10),log(gene_exp$X5H_TPM), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X5H_RPG10), y = log(X5H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT RPG10 and Illumina TPM Bo_5H") +
  scale_x_continuous("ln(RPG10 +1)") +
  scale_y_continuous(expression("ln(TPM +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 0.7, y = 11, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme

#ONT RPG10 and Illumina TPM 6H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)",  
                   list(           
                     r2 = format(cor(log(gene_exp$X6H_RPG10),log(gene_exp$X6H_TPM), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X6H_RPG10), y = log(X6H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT RPG10 and Illumina TPM Bo_6H") +
  scale_x_continuous("ln(RPG10 +1)") +
  scale_y_continuous(expression("ln(TPM +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 0.7, y = 11, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme

#ONT abs/emb and Illumina abs/emb 5H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(gene_exp$X5H_abs_emb),log(gene_exp$X5H_abs_emb.1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X5H_abs_emb), y = log(X5H_abs_emb.1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT abs and Illumina abs Bo_5H") +
  scale_x_continuous("ln(Transcripts per embryo (ONT) +1)") +
  scale_y_continuous(expression("ln(Transcripts per embryo (Ill) +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 2, y = 11, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme

#ONT abs/emb and Illumina abs/emb 6H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(gene_exp$X6H_abs_emb),log(gene_exp$X6H_abs_emb.1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X6H_abs_emb), y = log(X6H_abs_emb.1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT abs and Illumina abs Bo_6H") +
  scale_x_continuous("ln(Transcripts per embryo (ONT) +1)") +
  scale_y_continuous(expression("ln(Transcripts per embryo (Ill) +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 0.7, y = 11, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme

#See if Illumina has some length bias
plot(gene_exp$length,gene_exp$X5H_TPM, xlim = c(0,10000))

###Repeating ONT ratio bias assesment
a<-read.table('ERCC_expected_and_RPG10k.txt', header = T, row.names = 1, sep = '\t')

#We are going to plot ratios using absolute numbers (but not the per embryo abs)
a$abs_5H_ratio <- log2(a$X5H_abs)/log2(a$X5H_exp)
a$abs_6H_ratio <- log2(a$X6H_abs)/log2(a$X6H_exp)
#a$abs_5H_ratio2 <- log2(a$X5H_abs/a$X5H_exp) This is bad coz the 
a$abs_5H_ratio3 <- a$X5H_abs/a$X5H_exp
b <- a

row_sub = apply(b, 1, function(row) all(row !=0 ))
b <- b[row_sub,]

#GC 5H time point
ggplot(data = b, aes(x = ercc_gc, y = abs_5H_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with GC 5H ONT") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

ggplot(data = b, aes(x = ercc_length, y = abs_5H_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with length 5H ONT") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#GC 6H time point
ggplot(data = b, aes(x = ercc_gc, y = abs_6H_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with GC 6H ONT") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

ggplot(data = b, aes(x = ercc_length, y = abs_6H_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with length 6H ONT") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme


###Correlation of ERCC quantification between ONT and Illumina in absolute and then TPM
#5H
gene_exp <- read.table('ercc_ONT_illumina', header = T, row.names = 1, sep = '\t')
head(gene_exp)

lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(gene_exp$X5H_abs_emb),log(gene_exp$X5H_abs_emb.1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X5H_abs_emb), y = log(X5H_abs_emb.1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT abs_emb and Illumina abs_emb Bo_5H") +
  scale_x_continuous("ln(ONT Transcripts/embryo)") +
  scale_y_continuous(expression("ln(Illumina Transcripts/embryo)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 9, y = 16, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme

#6H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(gene_exp$X6H_abs_emb),log(gene_exp$X6H_abs_emb.1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_exp, aes(x = log(X6H_abs_emb), y = log(X6H_abs_emb.1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT abs_emb and Illumina abs_emb Bo_6H") +
  scale_x_continuous("ln(ONT Transcripts/embryo)") +
  scale_y_continuous(expression("ln(Illumina Transcripts/embryo)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 9, y = 16, label = lm_eqn(gene_exp), parse = T, color = "red") +
  my.theme
