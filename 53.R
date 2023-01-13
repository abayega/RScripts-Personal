rm(list=ls())

setwd('/home/abayega/R/tests/test53')
total_reads <- read.table("number_of_reads.tsv",header=F, sep = '\t')
summary_file <- read.table("/media/abayega/Padlock_DT/reads/nanopore/covid19/B004/B004_twist1/sequencing_summary_PAH62448_851e1099.txt_edited",header=T, sep = '\t')
flagstats <- read.table("flagstats.tsv",header=T, sep = '\t')
#summary_file <- fread("/media/abayega/Padlock_DT/reads/nanopore/covid19/B004/B004_twist1/sequencing_summary_PAH62448_851e1099.txt_edited",header=T, sep = '\t')
raw_count1 <- read.table("all_coverage.bed",header=F, sep = '\t')
covseq_report2 <- read.csv("V4testSamples_test.csv", header=T)
covseq_report3 <- read.csv("covseq_report_edited.csv", header=T)
covseq_report5 <- read.csv("LATEST_REPORT.csv", header=T)
variants1 <- read.table("sars_variants_v3",header=T, sep = '\t')
variants1 <- read.table("sars_variants_v3r",header=T, sep = '\t')
v4_concentrations <- read.table("V4_PCR_concn_of_serial_dilution.txt", header = T)

library(ggplot2)
library(plyr)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14, face="bold"),
                  axis.title.x = element_text(size=14, face="bold"), #axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(size=14, face="bold"), #axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

# Function to calculate the mean and the standard deviation
# for each group
#+++++++++++++++++++++++++
# data : a data frame
# varname : the name of a column containing the variable to be summariezed
# groupnames : vector of column names to be used as
# grouping variables
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

head(total_reads_b)
##Find total number of reads generated per barcode, and pass/fail categories
total_reads_a = aggregate(summary_file$protocol, by = list(summary_file$protocol, summary_file$filter), FUN = length)
colnames(total_reads_a) <- c('protocol','filter','reads')
total_reads_a$reads2 <- total_reads_a$reads/1000000
ggplot(total_reads_a, aes(x=protocol, y=reads2, fill=filter)) + 
  geom_col() +
  ggtitle("Total reads distribution among protocols") +
  xlab("Protocol") +
  ylab("Total number of reads (millions)") +
  my.theme

##Find total number of reads generated per barcode, and pass/fail categories
total_reads_b = subset(summary_file, select = c("protocol",'filter','length'))
colsx = c('length')
total_reads_b = aggregate(x = total_reads_b[ , colnames(total_reads_b) %in% colsx], by = list(total_reads_b$protocol, total_reads_b$filter), FUN = sum)
colnames(total_reads_b) <- c('protocol','filter','bases')
total_reads_b$bases2 <- total_reads_b$bases/1000000000
ggplot(total_reads_b, aes(x=protocol, y=bases2, fill=filter)) + 
  geom_col() +
  ggtitle("Distribution of total bases among protocols") +
  xlab("Protocol") +
  ylab("Total number of bases (gigabases)") +
  my.theme
ggsave(filename = "Distribution of total bases among protocols.png", width = 25, height = 15, units = "cm")

##Find total number of pass reads (before pychopping) generated across protocols
#Do this on the cluster for memory issues
total_reads_x = aggregate(summary_file$protocol, by = list(summary_file$protocol, summary_file$dilution, summary_file$barcode, summary_file$filter), FUN = length)
total_reads_x = read.table('total_reads_x',header = T, sep = '\t')
colnames(total_reads_x) <- c('protocol','part_p_ml','barcode','filter','reads')
total_reads_x_pass <- subset(total_reads_x, filter == 'PASS' & protocol != 'unclassified')
total_reads_x_pass$part_p_ml <- as.numeric(total_reads_x_pass$part_p_ml, options(scipen = 999))
total_reads_x_pass$reads <- total_reads_x_pass$reads/10000
total_reads_x_pass_sum <- data_summary(total_reads_x_pass, varname="reads", groupnames=c("protocol", "part_p_ml"))

#Turn your 'treatment' column into a character vector
total_reads_x_pass_sum$part_p_ml <- as.character(total_reads_x_pass_sum$part_p_ml )
#Then turn it back into a factor with the levels in the correct order
total_reads_x_pass_sum$part_p_ml  <- factor(total_reads_x_pass_sum$part_p_ml , levels=sort(unique(total_reads_x_pass_sum$part_p_ml), decreasing = T),ordered = F)

ggplot(total_reads_x_pass_sum, aes(x=protocol, y=reads, group=part_p_ml, color=part_p_ml)) + 
  geom_errorbar(aes(ymin=reads-sd, ymax=reads+sd), width=.1, position = position_dodge(width=0.6)) +
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  labs(y = "Total number of pass reads (x10,000)") +
  labs(x = "Protocol") +
  #ylim(0,100) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

##I want to have a global look at how many pychopped pass reads were generated per barcode
colnames(total_reads) <- c('protocol', 'dilution','pass_reads')
#Turn your 'treatment' column into a character vector
total_reads$dilution <- as.character(total_reads$dilution )
#Then turn it back into a factor with the levels in the correct order
total_reads$dilution  <- factor(total_reads$dilution , levels=unique(total_reads$dilution),ordered = TRUE)

#Summarize the data :
total_reads$pass_reads <- total_reads$pass_reads/10000
df3 <- data_summary(total_reads, varname="pass_reads", 
                    groupnames=c("protocol", "dilution"))

ggplot(df3, aes(x=protocol, y=pass_reads, group=dilution, color=dilution)) + 
  geom_errorbar(aes(ymin=pass_reads-sd, ymax=pass_reads+sd), width=.1, position = position_dodge(width=0.6)) +
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  labs(y = "Total number of pass reads (x10,000)") +
  labs(x = "") +
  #ylim(0,100) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

##Now, get number of SARS-CoV2 mapped reads across protocols
#Turn your 'treatment' column into a character vector
flagstats$concn <- as.character(flagstats$concn )
#Then turn it back into a factor with the levels in the correct order
flagstats$concn  <- factor(flagstats$concn , levels=unique(flagstats$concn),ordered = TRUE)
df3 <- data_summary(flagstats, varname="mapping_rate", 
                    groupnames=c("Protocol", "concn"))

ggplot(df3, aes(x=Protocol, y=mapping_rate, group=concn, color=concn)) + 
  geom_errorbar(aes(ymin=mapping_rate-sd, ymax=mapping_rate+sd), width=.1, position = position_dodge(width=0.6)) +
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  labs(y = "Read mapping rate (%)") +
  labs(x = "Protocol") +
  #ylim(0,100) +
  my.theme

##Try plotting read length distribution
p <- ggplot(summary_file, aes(x=protocol, y=length, color=filter)) +
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(y = "read length (bp)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

##Plotting all_coverage.bed
colnames(raw_count1) <- c('protocol', 'dilution','coverage')
#Turn your 'treatment' column into a character vector
raw_count1$dilution <- as.character(raw_count1$dilution )
#Then turn it back into a factor with the levels in the correct order
raw_count1$dilution  <- factor(raw_count1$dilution , levels=unique(raw_count1$dilution),ordered = TRUE)

#Summarize the data :
df3 <- data_summary(raw_count1, varname="coverage", 
                    groupnames=c("protocol", "dilution"))
head(df3)

#Try ggplot line plots http://www.sthda.com/english/wiki/ggplot2-line-plot-quick-start-guide-r-software-and-data-visualization
# coverage plotting

head(raw_count1)

df1 <- df3[which(df3$protocol == "artic_v3"),]

# Standard deviation of the mean
ggplot(df1, aes(x=dilution, y=coverage, group=protocol, color=protocol)) + 
  geom_errorbar(aes(ymin=coverage-sd, ymax=coverage+sd), width=.1) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()+
  labs(y = "genome breadth of coverage (%)") +
  labs(x = "concentration (viral particle copies/ml)") +
  ylim(0,100) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

# Use position_dodge to move overlapped errorbars horizontally
ggplot(df3, aes(x=dilution, y=coverage, group=protocol, color=protocol)) + 
  geom_errorbar(aes(ymin=coverage-sd, ymax=coverage+sd), width=.1, 
                position=position_dodge(0.05)) +
  geom_line() + geom_point()+
  scale_color_brewer(palette="Paired")+theme_minimal()


###Okay, now I want to plot some summary stats from CovSeq samples that have been processed in the lab

#Plot a violin of concentrations from all samples that have concentration >= 0
dim(covseq_report2) #The whole dataset has 68373 samples

#Now subset only samples with conce >= 0
covseq_report2_sub <- subset(covseq_report2, covseq_report2$Concentration >= 0, select = c(Concentration))
dim(covseq_report2_sub) #We have 66336 samples left
summary(covseq_report2_sub$Concentration)

#violin
p <- ggplot(covseq_report2_sub, aes(x="All samples", y=Concentration)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(y = "concentration (ng/ul)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#Make a violin plot that shows concentration of PCR products for each lineage
p <- ggplot(covseq_report3, aes(x=PANGOLIN, y=Concentration, fill=Concentration)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Lineage') +
  labs(y = "concentration (ng/ul)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#Attempt to make concentration for each category in lineages you want to assay
levels(factor(covseq_report3$PANGOLIN))
lineages <- c( "A.2.5","AY.10","AY.12","AY.20","AY.24","AY.3","AY.4","AY.5","B.1.1.7","B.1.351","B.1.525","B.1.617.2","B.1.621","B.1.621.1","P.1")
lin_assay <- c( "AY.10","AY.3","AY.4","B.1.1.7","B.1.351","B.1.617.2","B.1.621","P.1")

covseq_report3_lin <- subset(covseq_report3, covseq_report3$Concentration > 0 & PANGOLIN %in% lin_assay, select = c(Concentration, PANGOLIN, grade, avgCov ))
covseq_report3_lin$grade <- factor(covseq_report3_lin$grade, levels=c("low","mid","high"),ordered = TRUE)

p <- ggplot(covseq_report3_lin, aes(x=PANGOLIN, y=Concentration, fill=grade), color=grade) +
  geom_boxplot(width=0.5) + theme_minimal() +
  scale_fill_manual(values = c("yellow", "red", "green")) +
  labs(x = "") +
  labs(fill = 'grade') +
  labs(y = "concentration (ng/ul)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

#Now, make a violin of lineage coverage
p <- ggplot(covseq_report3, aes(x=PANGOLIN, y=avgCov, fill=avgCov)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'Lineage') +
  labs(y = "average coverage (reads/base)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

##Now, try to see if concentration is correlated with average coverage
#Subset the dataframe to get equal numbers of samples with concentration and average coverage
covseq_report3_sub = subset(covseq_report3, covseq_report3$Concentration > 0,select = c(avgCov,Concentration))

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45, face="bold"),
                  axis.title.y = element_text(vjust=.2, face="bold"),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(Concentration ~ avgCov, df1);
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
                     r2 = format(cor(covseq_report3_sub$avgCov,covseq_report3_sub$Concentration, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = covseq_report3_sub, aes(x = avgCov, y = Concentration)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  ggtitle("correlation of coverage with concentration") +
  scale_x_continuous("average coverage (reads/base)") +
  scale_y_continuous("concentration (ng/ul)") +
  geom_smooth(method = "loess") +
  annotate(geom = "text", x = 300, y = 210, label = lm_eqn(covseq_report3_sub), parse = T, color = "red") +
  annotate(geom = "text", x = 175, y = 230, label = lm_eqn2(covseq_report3_sub), parse = T, color = "red") +
  my.theme

#Some summary stats I did with Pandas in Python, see if I can do it in R
aggregate(covseq_report3$grade, by = list(covseq_report3$PANGOLIN, covseq_report3$grade), FUN = length)
lin_cats = as.table(table(covseq_report3$PANGOLIN,covseq_report3$grade) )

write.csv(lin_cats, file = "lineage_categorisation_on_concentraiton", quote = F)


covseq_report5_sub <- subset(covseq_report5, GISAID_ID!="NA", select = c(Sample.Name, GISAID_ID))
write.csv(covseq_report5_sub, file = "GISAID_submitted_samples.csv", quote = F)
length(which(covseq_report5$GISAID_ID!="NA" & covseq_report5$PANGOLIN=="B.1.1.7"))

#Make a grouped boxplots for the concentration of the variants

variants1$Variant <- factor(variants1$Variant, levels=c("WT","SA","UK","BRZ","IND"),ordered = TRUE)
variants1$Dilution <- factor(variants1$Dilution, levels=c("_1e2","_1e3","_1e4","_1e5","_1e6"),ordered = TRUE)
variants1$ID <- factor(variants1$ID, levels=c("D","E","F"),ordered = TRUE)

#p <- ggplot(variants1, aes(x=Variant, y=Concentration, fill=ID), color=ID) +
p <- ggplot(variants1, aes(x=Variant, y=Concentration, fill=Dilution), color=Dilution) +
  geom_boxplot(width=0.5) + theme_minimal() +
  #scale_fill_manual(values = c("yellow", "red", "green")) +
  labs(x = "") +
  labs(fill = 'Dilution') +
  labs(y = "concentration (ng/ul)") + 
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p

##Plotting the concentrations of the serial dilutions I did for cell cultured variants
v4_concentrations_sub <- v4_concentrations[1:(nrow(v4_concentrations)-2),]
v4_concentrations_sub <- subset(v4_concentrations_sub, Dilution != "1e-2")
v4_concentrations_sub$Sample <- factor(v4_concentrations_sub$Sample, levels=c("WT","UK","SA","BRZ","IND"),ordered = TRUE)

p <- ggplot(v4_concentrations_sub, aes(x=Sample, y=Concentration, fill=Dilution), color=Dilution) +
  geom_boxplot(width=0.5) + theme_minimal() +
  #scale_fill_manual(values = c("yellow", "red", "green")) +
  labs(x = "") +
  labs(fill = 'Dilution') +
  labs(y = "concentration (ng/ul)") + 
  ggtitle("Post PCR concentration of serial diluted variants") +
  #ylim(0,3000) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p
