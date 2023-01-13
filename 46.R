# Set working directory
rm(list = ls())
setwd("/home/abayega/R/tests/test46")
raw_count <- read.table("ERCC_C1H",header=T,row.names=1)

head(raw_count)
colnames(raw_count) = c('transcript_id','expected','B02C1HF','B05C1HF','B06C1HF','B09C1HF','B12C1HF','B03C1HM','B07C1HM','B08C1HM','B10C1HM','B11C1HM',
                        'B02C1HFr','B05C1HFr','B06C1HFr','B09C1HFr','B12C1HFr','B03C1HMr','B07C1HMr','B08C1HMr','B10C1HMr','B11C1HMr','SD','CV','total','tot_chop')

filtered <-raw_count[,12:ncol(raw_count)]
unfiltered <-raw_count[,1:12]

#For ERCCs, making correlation plot 
library(ggplot2)
my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  df1 = subset(df1, rownames(df1)!='ERCC-00116',select=c(expected, B02C1HF))
  #print(length(rownames(df1)))
  m = lm(log2(expected) ~ log2(B02C1HF), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
#Ploting correlation for un corrected and unchopped reads
unfil_total = subset(raw_count, SD < 300,select=c(expected, total))
ggplot(data = unfil_total, aes(x = log2(expected), y = log2(total))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs unfiltered") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" counts")) +
  annotate(geom = "text", x = 9, y = 11, label = lm_eqn(unfil_total,expected,total), parse = T, color = "red",size=6) +
  my.theme

#Ploting correlation for corrected and chopped reads
fil_total = subset(raw_count, SD < 300,select=c(expected, tot_chop))
ggplot(data = fil_total, aes(x = log2(expected), y = log2(tot_chop))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs filtered") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" counts")) +
  annotate(geom = "text", x = 9, y = 11, label = lm_eqn(fil_total), parse = T, color = "red",size=6) +
  my.theme

#Ploting correlation of unfiltered B02
unfil_B02 = subset(raw_count, B02C1HF > 1,select=c(expected, B02C1HF))
ggplot(data = unfil_B02, aes(x = log2(expected), y = log2(B02C1HF))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs C1H B02") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" counts")) +
  annotate(geom = "text", x = 14, y = 11, label = lm_eqn(unfil_B02), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of filtered B02
fil_B02 = subset(raw_count, B02C1HFr > 1,select=c(expected, B02C1HFr))
ggplot(data = fil_B02, aes(x = log2(expected), y = log2(B02C1HFr))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs C1H B02r") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" counts")) +
  annotate(geom = "text", x = 13, y = 11, label = lm_eqn(fil_B02), parse = T, color = "red",size=6) +
  my.theme

#Ploting correlation of ONTs 5H
ont_exprn_erccs_1H = subset(ont_exprn_erccs, X5H_TPM >= 0.01,select=c(X5H_exp, X5H_TPM))
ggplot(data = ont_exprn_erccs_1H, aes(x = log2(X5H_exp), y = log2(X5H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs RPG10K Bo_E_5H ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10K")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(ont_exprn_erccs_1H), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of ONTs 6H
ont_exprn_erccs_1H = subset(ont_exprn_erccs, X6H_TPM >= 0.01,select=c(X6H_exp, X6H_TPM))
ggplot(data = ont_exprn_erccs_1H, aes(x = log2(X6H_exp), y = log2(X6H_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs RPG10K Bo_E_6H ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10K")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(ont_exprn_erccs_1H), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of ONTs female
ont_exprn_erccs_1H = subset(ont_exprn_erccs, fem_TPM >= 0.01,select=c(fem, fem_TPM))
ggplot(data = ont_exprn_erccs_1H, aes(x = log2(fem), y = log2(fem_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs RPG10K female ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10K")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(ont_exprn_erccs_1H), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of ONTs male
ont_exprn_erccs_1H = subset(ont_exprn_erccs, mal_TPM >= 0.01,select=c(fem, mal_TPM))
ggplot(data = ont_exprn_erccs_1H, aes(x = log2(fem), y = log2(mal_TPM))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs RPG10K male ONT") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" RPG10K")) +
  annotate(geom = "text", x = 16, y = 6, label = lm_eqn(ont_exprn_erccs_1H), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of ONTs end here


#biases in male heads ONT
ont_exprn_erccs_male = subset(ont_exprn_erccs, mal_TPM >= 0.02,select=c(mal, mal_abs,ercc_length,ercc_gc))
ont_exprn_erccs_male = subset(ont_exprn_erccs_male, rownames(ont_exprn_erccs_male)!='ERCC-00116')
ont_exprn_erccs_male$exp_ratio <- log2(ont_exprn_erccs_male$mal_abs/ont_exprn_erccs_male$mal)

ggplot(data = ont_exprn_erccs_male, aes(x = ercc_gc, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with GC male ONT") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  my.theme

#Plotting ERCC length
ggplot(data = ont_exprn_erccs_male, aes(x = ercc_length, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of abs ratios with length male ONT") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Observed:Expected)")) + #, limits = c(1,6)) +
  my.theme


####Getting beta for Illumina rsem quantified TPMs
setwd("/home/abayega/R/tests/test24")
raw_count <- read.table("rsem_combined",header=T,row.names=1,sep = '\t')
head(raw_count)

raw_count_clean = subset(raw_count, rownames(raw_count)!='ERCC-00116',select=c(exp_5H, TPM_5H, TPM_6H,ercc_length, ercc_gc))
#Am going to set a TPM cut off of 0.2 and remove all ERCCs with < 0.2 TPM
raw_count_clean[raw_count_clean < 0.2] = 0

#5H beta
glm(raw_count_clean$TPM_5H[which(raw_count_clean$TPM_5H >=0.2)] ~ offset(log(raw_count_clean$exp_5H[which(raw_count_clean$TPM_5H >=0.2)])), family=poisson(link=log)) #we use 57 ERCCs
#6H beta
glm(raw_count_clean$TPM_6H[which(raw_count_clean$TPM_6H >=0.2)] ~ offset(log(raw_count_clean$exp_5H[which(raw_count_clean$TPM_6H >=0.2)])), family=poisson(link=log)) #we use 75 ERCCs


#Correlation of ERCC quantification between ONT and Illumina in absolute and then TPM
setwd('/home/abayega/R/tests/test24')
raw_count <- read.table('rsem_rpg10k_combined', header = T, row.names = 1, sep = '\t')
raw_count2 <- read.table("ercc_rsem_combined",header=T,row.names=1,sep = '\t')
head(gene_exp)

ercc_raw_count <- raw_count[1:92,]
ercc_raw_count$expected <- raw_count2$exp_5H

#For ERCCs, let us try to make plots for detection limits
#Ploting correlation of illumina
#Same theme as line 41
lm_eqn = function(df1){
  df1 = subset(df1, rownames(df1)!='ERCC-00116',select=c(expected, TPM_6H))
  print(length(rownames(df1)))
  m = lm(log2(TPM_6H) ~ log2(expected), df1);
  eq <- substitute(italic(y) == a + b %.% italic(x)*","
                   ~~italic(r)^2~"="~r2~~", p < "~~p, 
                   list(a = format(coef(m)[1], digits = 2),
                        b = format(coef(m)[2], digits = 2),              
                        r2 = format(summary(m)$r.squared, digits = 3),              
                        p = format(coef(summary(m))[8], digits = 3, scientific=T)             
                   ))
  as.character(as.expression(eq));                 
}
#Ploting correlation of ill 5H
ill_exprn_erccs_5H = subset(ercc_raw_count, TPM_5H >= 0.1,select=c(expected, TPM_5H))
ggplot(data = ill_exprn_erccs_5H, aes(x = log2(expected), y = log2(TPM_5H))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs TPM Bo_E_5H illumina") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" TPM")) +
  annotate(geom = "text", x = 17, y = 9, label = lm_eqn(ill_exprn_erccs_5H), parse = T, color = "red",size=6) +
  my.theme
#Ploting correlation of ill 6H
ill_exprn_erccs_6H = subset(ercc_raw_count, TPM_6H >= 0.1,select=c(expected, TPM_6H))
ggplot(data = ill_exprn_erccs_6H, aes(x = log2(expected), y = log2(TPM_6H))) +
  geom_point() +
  geom_smooth(method = "lm") +
  ggtitle("correlation expected vs TPM Bo_E_6H illumina") +
  scale_x_continuous(expression(log[2]*" Transcript molecules")) +
  scale_y_continuous(expression(log[2]*" TPM")) +
  annotate(geom = "text", x = 14, y = 10.5, label = lm_eqn(ill_exprn_erccs_6H), parse = T, color = "red",size=6) +
  my.theme

#Now to Spearman correlation between ONT and Illumina, first for ERCCs then genes
#5H
lm_eqn = function(df1){
  print(length(rownames(df1)))
  df1 = subset(df1, rownames(df1)!='ERCC-00116')
  print(length(rownames(df1)))
  df1 = subset(df1, TPM_5H >= 0.1)
  print(length(rownames(df1)))
  df1 = subset(df1, RPG10K_5H >= 0.01)
  print(length(rownames(df1)))
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(df1$TPE_5H_ont),log(df1$TPE_5H), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = ercc_raw_count, aes(x = log(TPE_5H_ont), y = log(TPE_5H))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC ONT abs_emb and Illumina abs_emb Bo_5H") +
  scale_x_continuous("log(ONT Transcripts/embryo)") +
  scale_y_continuous(expression("log(Illumina Transcripts/embryo)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 7, y = 13, label = lm_eqn(ercc_raw_count), parse = T, color = "red",size=6) +
  my.theme

#6H
lm_eqn = function(df1){
  print(length(rownames(df1)))
  df1 = subset(df1, rownames(df1)!='ERCC-00116')
  print(length(rownames(df1)))
  df1 = subset(df1, TPM_6H >= 0.1)
  print(length(rownames(df1)))
  df1 = subset(df1, RPG10K_6H >= 0.01)
  print(length(rownames(df1)))
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log(df1$TPE_6H_ont),log(df1$TPE_6H), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = ercc_raw_count, aes(x = log(TPE_6H_ont), y = log(TPE_6H))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC ONT abs_emb and Illumina abs_emb Bo_6H") +
  scale_x_continuous("log(ONT Transcripts/embryo)") +
  scale_y_continuous(expression("log(Illumina Transcripts/embryo)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 8, y = 14, label = lm_eqn(ercc_raw_count), parse = T, color = "red",size=6) +
  my.theme
##End

##Trying to plot length and GC bias in Illumina
#GC 5H time point
ercc_raw_count_5H = subset(ercc_raw_count, TPM_5H >= 0.1,select=c(expected, abs_5H,gene_length,gene_gc))
ercc_raw_count_5H = subset(ercc_raw_count_5H, rownames(ercc_raw_count_5H)!='ERCC-00116')
ercc_raw_count_5H$exp_ratio <- log2(ercc_raw_count_5H$expected/ercc_raw_count_5H$abs_5H)

ggplot(data = ercc_raw_count_5H, aes(x = gene_gc, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC abs ratios with GC 5H illumina") +
  scale_x_continuous("ERCC GC content (%)") +
  scale_y_continuous(expression(log[2]*" (Expected:Observed)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#length 5H time point
ggplot(data = ercc_raw_count_5H, aes(x = gene_length, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC abs ratios with length 5H illumina") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Expected:Observed)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#GC 6H time point
ercc_raw_count_6H = subset(ercc_raw_count, TPM_6H >= 0.1,select=c(expected, abs_6H,gene_length,gene_gc))
ercc_raw_count_6H = subset(ercc_raw_count_6H, rownames(ercc_raw_count_6H)!='ERCC-00116')
ercc_raw_count_6H$exp_ratio <- log2(ercc_raw_count_6H$expected/ercc_raw_count_6H$abs_6H)

ggplot(data = ercc_raw_count_6H, aes(x = gene_gc, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC abs ratios with length 6H ONT") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Expected:Observed)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

#length 6H time point
ggplot(data = ercc_raw_count_6H, aes(x = gene_length, y = exp_ratio)) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ERCC abs ratios with length 6H illumina") +
  scale_x_continuous("ERCC length (bp)") +
  scale_y_continuous(expression(log[2]*" (Expected:Observed)")) + #, limits = c(1,6)) +
  #annotate(geom = "text", x = 17, y = 7, label = lm_eqn(all_ercc), parse = T, color = "red") +
  my.theme

###Gene expression correlation btn ONT abs/emb and Illumina abs/emb 5H
gene_raw_count <- raw_count[93:nrow(raw_count),]
gene_raw_count_5_6H = subset(gene_raw_count, RPG10K_5H >= 0.01, select=c(TPE_5H, TPE_6H, TPE_5H_ont, TPE_6H_ont, gene_length, gene_gc, transcript_length, transcript_gc))

lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log((gene_raw_count_5_6H$TPE_5H)+1),log((gene_raw_count_5_6H$TPE_5H_ont)+1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_raw_count_5_6H, aes(x = log(TPE_5H+1), y = log(TPE_5H_ont+1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT TPE and Illumina TPE Bo_5H") +
  scale_x_continuous("log (Transcripts per embryo (ONT) +1)") +
  scale_y_continuous(expression("log (Transcripts per embryo (ill) +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 7, y = 17, label = lm_eqn(gene_raw_count_5_6H), parse = T, color = "red",size=6) +
  my.theme

#Gene expression correlation btn ONT abs/emb and Illumina abs/emb 6H
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log((gene_raw_count_5_6H$TPE_6H)+1),log((gene_raw_count_5_6H$TPE_6H_ont)+1), method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_raw_count_5_6H, aes(x = log(TPE_6H+1), y = log(TPE_6H_ont+1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of ONT TPE and Illumina TPE Bo_6H") +
  scale_x_continuous("log (Transcripts per embryo (ONT) +1)") +
  scale_y_continuous(expression("log (Transcripts per embryo (ill) +1)")) + #, limits = c(1,6)) +
  annotate(geom = "text", x = 7, y = 17, label = lm_eqn(gene_raw_count_5_6H), parse = T, color = "red",size=6) +
  my.theme
#End

##Now let us try to see if log expression has correlation with length for illumina and ONT
#I will not save the graphs but by changing a few variables u can plot the correlation of gene expression with gene length or GC and transcript length or GC
lm_eqn = function(df1){
  eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                   list(           
                     r2 = format(cor(log((gene_raw_count_5_6H$TPE_6H)+1), gene_raw_count_5_6H$transcript_gc, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = gene_raw_count_5_6H, aes(x = transcript_gc, y = log(TPE_6H+1))) +
  geom_point() +
  geom_smooth(method = "lm", formula = y~(x)) +
  ggtitle("Correlation of Illumina TPE and gene length Bo_6H") +
  scale_x_continuous("gene length (bp)" , limits = c(0,100)) +
  scale_y_continuous(expression("log (Transcripts per embryo (ill) +1)")) + #, limits = c(0,100000)) +
  annotate(geom = "text", x = 20000, y = 17, label = lm_eqn(gene_raw_count_5_6H), parse = T, color = "red",size=6) +
  my.theme
#End

#Now, let us try to plot the correlation of expression accross timepoints
gene_raw_count_ont <- gene_raw_count[,41:46]
head(gene_raw_count_ont)

plot_cor <- function(df){
  cols = colnames(df)
  n = 1
  for (col in cols){
    for (col2 in cols){
      lm_eqn = function(df1){
        eq <- substitute(italic(r) == r2*", (Spearman's rho)", 
                         list(           
                           r2 = format(cor(log((df$col)+1), log((df$col2)+1), method=c("spearman")), digits = 3)
                         ))
        as.character(as.expression(eq));
      }
      ggplot(data = df, aes_string(x = col, y = col2)) +
        geom_point() +
        geom_smooth(method = "lm", formula = y~(x)) +
        #ggtitle("Correlation of Illumina TPE and gene length Bo_6H") +
        #scale_x_continuous("gene length (bp)" , limits = c(0,100)) +
        #scale_y_continuous(expression("log (Transcripts per embryo (ill) +1)")) + #, limits = c(0,100000)) +
        annotate(geom = "text", x = 8, y = 12, label = lm_eqn(gene_raw_count_5_6H), parse = T, color = "red",size=6) +
        my.theme
      #print(cols)
    }
    n = n + 1
  }
}
plot_cor(gene_raw_count_ont) #This function refused to work, not sure why

par(mfrow=c(6,6), omi = rep(0,4), mar = c(2,2,2,0.1), mgp=c(0.8,0,0))
par(mfrow=c(6,6))
reg<-lm(log(TPE_2H_ont+1) ~ log(TPE_1H_ont+1), data = gene_raw_count_ont)
coeff=coefficients(reg)
corn = format(cor(log((gene_raw_count_ont$TPE_2H_ont)+1), log((gene_raw_count_ont$TPE_1H_ont)+1), method=c("spearman")), digits = 2)
eq = paste0("r = ",corn)
plot(log((gene_raw_count_ont$TPE_1H_ont)+1),log((gene_raw_count_ont$TPE_2H_ont)+1), main=eq, ylab='a', xlab='v', pch = 20, cex = 0.1)
abline(reg, col="blue")

par(mfrow=c(6,6), mar = c(2,2,2,0.1), mgp=c(1,0,0))
plot_cor <- function(df){
  x = 1
  y = 1 
  for (i in names(df)){
    for (j in names(df)){
      a = df[[i]]+1
      b = df[[j]]+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      eq = paste0("r = ",corn)
      #plot(log(a),log(b), main=eq, ylab=paste('log (',x,'H)', sep = ""), xlab=paste('log (',y,'H)', sep = ""), pch = 20, cex = 0.1,xlim = c(0,20), ylim = c(0,20))
      plot(log(a),log(b), main=eq, ylab=paste('log (',x,'H)', sep = ""), xlab=paste('log (',y,'H)', sep = ""), pch = 20, cex = 0.1)
      abline(reg, col="blue")
      abline(a=0,b=1, col="red")
      x = x + 1
    }
    y = y + 1
    x = 1
  }
}
plot_cor(gene_raw_count_ont)
par(mfrow=c(1,1))

dev.off()
#Save the graph as correlation of ONT TPE accross timepoints
#Shit, I almost spent the whole afternoon today 100918 just to derive this function lines 423 - 453,

#K, now I try to plot expression correlation as one graph
setwd("/home/abayega/R/tests/test24")
exp_cor <- read.table("expression_correlation",header=T,row.names=1, sep = '\t')
head(exp_cor)
install.packages("gplots")
library('gplots')

colnames(exp_cor) <- c('1H','2H','3H','4H','5H','6H')
heatmap.2(as.matrix(exp_cor), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
           key=TRUE, density.info = 'none', key.xlab = 'Spearman rho') 
#heatmap.2(as.matrix(exp_cor), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", key=FALSE,
#          trace='none',lwid = c(.01,.99),lhei = c(.01,.99),margins = c(5,15 ))





#plot(ont_exprn_erccs_male$ercc_length,ont_exprn_erccs_male$mal_abs, ylim = c(0,10000000))

#low_expn = apply(ont_exprn, 1, function(row) all(row !=0 ))
#ont_exprn_erccs <- a[row_sub,]
#erccs <- a[1:53,]

head(erccs)

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


