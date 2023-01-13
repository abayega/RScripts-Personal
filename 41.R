setwd('/home/abayega/R/tests/test41')
raw_count1 <- read.table("Bo_E_1H_C010_10_pass_edited.trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count2 <- read.table("Bo_E_2H_C010_09_pass_edited.trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count3 <- read.table("Bo_E_3H_C010_08_pass.trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count4 <- read.table("Bo_E_4H_C010_06_pass.trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count5 <- read.table("Bo_E_5H_combined_pass_trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count6 <- read.table("Bo_E_6H_C010_07_pass.trimmed_stranded_edited_readlengths_GC.txt",header=F, sep = '\t')
raw_count_Fem <- read.table("Bo_Heads_female_C010_12_01_pass.both_adapters_stranded_readlengths_GC.txt",header=F, sep = '\t')
raw_count_Mal <- read.table("Bo_Heads_male_C010_12_02_pass.both_adapters_stranded_readlengths_GC.txt",header=F, sep = '\t')
totRNA <- read.table("total_RNA_extracted_from_capitata_single_embryos",header=T, sep = '\t')

raw_countNON <- raw_count_Mal
colnames(raw_countNON) <- c('readname','length', 'GC')
head(raw_countNON)

summary(raw_countNON$length)
summary(raw_countNON$GC)


readlengths <- matrix(raw_count6$V2, raw_count2$V2, raw_count3$V2, raw_count4$V2, raw_count5$V2,
                      raw_count1$V2, raw_count_Fem$V2, raw_count_Mal$V2)
No.exons <- matrix(raw_countNON$No.of.exons, raw_countZygo$No.of.exons)
intronLen <- matrix(raw_countNON$totalIntronLength, raw_countZygo$totalIntronLength)
exonLen <- matrix(raw_countNON$totalExonLength, raw_countZygo$totalExonLength)

#Try ggplot violin plots
library(ggplot2)
#http://www.sthda.com/english/wiki/ggplot2-violin-plot-quick-start-guide-r-software-and-data-visualization

# Gene length

readlengths <- c(raw_count1$V2, raw_count2$V2, raw_count3$V2, raw_count4$V2, raw_count5$V2,raw_count6$V2, raw_count_Fem$V2, raw_count_Mal$V2)
conditions <- c(rep('1H',length(raw_count1$V2)), rep('2H', length(raw_count2$V2)), rep('3H', length(raw_count3$V2)), rep('4H', length(raw_count4$V2)),
         rep('5H',length(raw_count5$V2)), rep('6H', length(raw_count6$V2)), rep('female', length(raw_count_Fem$V2)), rep('male', length(raw_count_Mal$V2)))
expn_mat <- data.frame(readlengths, conditions)

colnames(expn_mat) <- c('readLength', 'conditions')
#Turn your 'treatment' column into a character vector
expn_mat$conditions <- as.character(expn_mat$conditions)
#Then turn it back into a factor with the levels in the correct order
expn_mat$conditions <- factor(expn_mat$conditions, levels=unique(expn_mat$conditions),ordered = TRUE)

head(expn_mat)

p <- ggplot(expn_mat, aes(x=conditions, y=readlengths, fill=conditions)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'timepoint') +
  labs(y = "Read length (bases)") + 
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

# GC

readlengths <- c(raw_count1$V3, raw_count2$V3, raw_count3$V3, raw_count4$V3, raw_count5$V3,raw_count6$V3, raw_count_Fem$V3, raw_count_Mal$V3)
conditions <- c(rep('1H',length(raw_count1$V3)), rep('2H', length(raw_count2$V3)), rep('3H', length(raw_count3$V3)), rep('4H', length(raw_count4$V3)),
                rep('5H',length(raw_count5$V3)), rep('6H', length(raw_count6$V3)), rep('female', length(raw_count_Fem$V3)), rep('male', length(raw_count_Mal$V3)))
expn_mat <- data.frame(readlengths, conditions)

colnames(expn_mat) <- c('readLength', 'conditions')
#Turn your 'treatment' column into a character vector
expn_mat$conditions <- as.character(expn_mat$conditions)
#Then turn it back into a factor with the levels in the correct order
expn_mat$conditions <- factor(expn_mat$conditions, levels=unique(expn_mat$conditions),ordered = TRUE)

head(expn_mat)

p <- ggplot(expn_mat, aes(x=conditions, y=readlengths, fill=conditions)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "") +
  labs(fill = 'timepoint') +
  labs(y = "Read length (bases)") + 
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

#total_RNA_extracted_from_capitata_single_embryos

readlengths <- c(totRNA$X1H, totRNA$X2H, totRNA$X3H, totRNA$X4H, totRNA$X5H, totRNA$X6H, totRNA$X7H, totRNA$X8H, totRNA$X9H, totRNA$X10H
                 , totRNA$X11H, totRNA$X12H, totRNA$X13H, totRNA$X14H, totRNA$X15H)
conditions <- c(rep('1H',length(totRNA$X1H)), rep('2H', length(totRNA$X2H)), rep('3H', length(totRNA$X3H)), rep('4H', length(totRNA$X4H)),
                rep('5H',length(totRNA$X5H)), rep('6H', length(totRNA$X6H)), rep('7H', length(totRNA$X7H)), rep('8H', length(totRNA$X8H)),
                rep('9H',length(totRNA$X9H)), rep('10H', length(totRNA$X10H)), rep('11H', length(totRNA$X11H)), rep('12H', length(totRNA$X12H)),
                rep('13H',length(totRNA$X13H)), rep('14H', length(totRNA$X14H)), rep('15H', length(totRNA$X15H)))
expn_mat <- data.frame(readlengths, conditions)

colnames(expn_mat) <- c('Total RNA', 'conditions')

#Turn your 'treatment' column into a character vector
expn_mat$conditions <- as.character(expn_mat$conditions)
#Then turn it back into a factor with the levels in the correct order
expn_mat$conditions <- factor(expn_mat$conditions, levels=unique(expn_mat$conditions),ordered = TRUE)

head(expn_mat)

p <- ggplot(expn_mat, aes(x=conditions, y=readlengths, fill=conditions)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  geom_boxplot(width=0.1) + theme_minimal() +
  labs(x = "timepoint") +
  labs(fill = 'timepoint') +
  labs(y = "Total RNA (ng/embryo)") + 
  #ylim(0,100) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
p


