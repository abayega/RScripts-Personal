setwd("~/R/tests/test9")
a<-read.table('ERCC_expected_and_RPG10k.txt', header = T, row.names = 1, sep = '\t')
erccs <- a
head(a)
plot(log(a$X1H_exp),log(a$X1H_TPM),xlab="log10(number_of_ERCC_molecules)",ylab="log10(number_of_sequenced_reads)")
glm(a$X1H_TPM ~ offset(log(a$X1H_exp)), family=poisson(link=log))
new_y <- log(a$X1H_exp) - 13.99
lines(log(a$X1H_exp),new_y,col="red")

glm(a$X2H_TPM ~ offset(log(a$X2H_exp)), family=poisson(link=log))
glm(a$X3H_TPM ~ offset(log(a$X3H_exp)), family=poisson(link=log))
glm(a$X4H_TPM ~ offset(log(a$X4H_exp)), family=poisson(link=log))
glm(a$X5H_TPM ~ offset(log(a$X5H_exp)), family=poisson(link=log))
glm(a$X6H_TPM ~ offset(log(a$X6H_exp)), family=poisson(link=log))


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

#Plotting ERCC GC ratios
ggplot(data = ercc_ratios, aes(x = ercc_GC, y = GC_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  scale_x_continuous("ERCC GC (%)") +
  scale_y_continuous("ratio Observed:Expected ERCC") + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 10, y = 0.5, label = lm_eqn(ercc_ratios), parse = T, color = "red") +
  my.theme


#Plot length and GC ratios but for expression
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

#Printing read lengths of male contigs
male<-read.table('male_contigs_readlength', header = F, sep = '\t')
boxplot(male$V1, ylab = "scaffold_length")
hist(male$V1, breaks = 50, ylim = c(0,100))
sum(male$V1)
