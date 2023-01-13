rm(list=ls())

setwd('/home/abayega/R/tests/test53/B004_11_2')
output_workspace = '/home/banthony/R/tests/test53/B004_11_2'
summary_file <- read.table("sequencing_summary_PAI76911_dc3ee531.txt_edited",header=T, sep = '\t')
consensus_variants <- read.table("consensus_variant.txt", header=T, sep = '\t')
consensus_variants <- read.table("consensus_variant2.txt", header=T, sep = '\t')

#Remove Ebbs Twist1 samples
summary_file <- subset(summary_file, protocol != "ebb")

#library("devtools")
library(gridExtra)
library(grid)
library(ggplot2)
library(lattice)
library(plyr)
library("plot3D")
library(plotly)
library("gg3D")

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14, face="bold"),
                  axis.title.x = element_text(size=14, face="bold"), #axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(size=14, face="bold"), #axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

summary_file <- raw_count1
head(total_reads_b)

##Find total number of reads generated per protocol, and pass/fail categories
total_reads_a = aggregate(summary_file$protocol, by = list(summary_file$protocol, summary_file$pass_filter), FUN = length)
colnames(total_reads_a) <- c('protocol','filter','reads')
total_reads_a$reads2 <- total_reads_a$reads/1000000
ggplot(total_reads_a, aes(x=protocol, y=reads2, fill=filter)) + 
  geom_col() +
  ggtitle("Total reads distribution among protocols") +
  xlab("Protocol") +
  ylab("Total number of reads (millions)") +
  my.theme

##Find total number of reads generated per barcode, and pass/fail categories
total_reads_b = subset(summary_file, select = c("protocol",'pass_filter','length'))
colsx = c('length')
total_reads_b = aggregate(x = total_reads_b[ ,colnames(total_reads_b) %in% colsx], by = list(total_reads_b$protocol, total_reads_b$pass_filter), FUN = sum)
colnames(total_reads_b) <- c('protocol','filter','bases')
total_reads_b$bases2 <- total_reads_b$bases/1000000000
total_reads_b1 <- subset(total_reads_b, protocol != "ebb")
ggplot(total_reads_b1, aes(x=protocol, y=bases2, fill=filter)) + 
  geom_col() +
  ggtitle("Distribution of total bases among protocols") +
  xlab("Protocol") +
  ylab("Total number of bases (gigabases)") +
  my.theme
ggsave(filename = "Distribution of total bases among protocols.png", width = 25, height = 15, units = "cm")

#Now, let's focus on pass reads
summary_file_pass <- subset(summary_file, pass_filter == "PASS")
ggplot(summary_file_pass, aes(x=protocol, y=length, color=filter)) +
  geom_violin(fill='#A4A4A4', color="darkred") +
  labs(x = "") +
  labs(y = "read length (bp)") + 
  ylim(0,1500) +
  theme(
    #plot.title = element_text(color="red", size=14, face="bold.italic"),
    #axis.text.x = element_text(size=14,angle=90,hjust=.5,vjust=.5,face="plain"),
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")

##Plotting read lengths
summary_file_pass_V3 <- subset(summary_file_pass, protocol == "artic_v3")
artic_v3 <- ggplot(summary_file_pass_V3, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=10) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(400,600) +
  my.theme 
ggsave(filename = "Artic V3 pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")

summary_file_pass_V4.1 <- subset(summary_file_pass, protocol == "artic_v4.1")
artic_v4.1 <- ggplot(summary_file_pass_V4.1, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(400,600) +
  my.theme 
ggsave(filename = "Artic V4.1 pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")  

summary_file_pass_midnight <- subset(summary_file_pass, protocol == "midnight")
midnight <- ggplot(summary_file_pass_midnight, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(100,1500) +
  my.theme 
ggsave(filename = "Midnight pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")  

summary_file_pass_snap <- subset(summary_file_pass, protocol == "snap")
snap <- ggplot(summary_file_pass_snap, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(100,900) +
  my.theme 
ggsave(filename = "SNAP pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")  

summary_file_pass_ebb <- subset(summary_file_pass, protocol == "entebbe")
entebbe <- ggplot(summary_file_pass_ebb, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(100,2200) +
  my.theme 
ggsave(filename = "Entebbe pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")

png(filename = paste0(output_workspace,"Pass read length distribution across protocols.png"),width = 30,height=20, units = "cm",res = 1200)
grid.arrange(artic_v3,artic_v4.1,midnight,snap,nrow=2,ncol=2)
dev.off()

summary_file_pass_unclass <- subset(summary_file_pass, protocol == "unclassified" & pass_filter == "PASS")
unclass_pass <- ggplot(summary_file_pass_unclass, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(100,2200) +
  my.theme 
ggsave(filename = "Unclassified pass read length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")

summary_file_fail <- subset(summary_file, pass_filter == "FAIL")
fail <- ggplot(summary_file_fail, aes(x=length)) + 
  geom_histogram(color="darkblue", fill="lightblue",binwidth=20) +
  geom_vline(aes(xintercept=mean(length)),color="blue", linetype="dashed", size=1) +
  labs(x = "read length (bp)") +
  #labs(y = "count") + 
  xlim(100,2200) +
  my.theme 
ggsave(filename = "All FAIL reads length distribution.png", path=output_workspace, width = 20, height = 15, units = "cm")


##Let us try to find total number of pass reads per variant per protocol
total_reads_var = aggregate(summary_file_pass$protocol, by = list(summary_file_pass$protocol, summary_file_pass$variant), FUN = length)
colnames(total_reads_var) <- c('protocol','variant','total_pass_reads')

#Remove unclassified and negative ctrl reads
total_reads_var <- subset(total_reads_var, protocol != "unclassified" & variant != 'NegativeCtrl')

#divide reads by 1,000,000
total_reads_var$total_pass_reads2 <- total_reads_var$total_pass_reads/1000000

ggplot(total_reads_var, aes(x=protocol, y=total_pass_reads2, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Total number of pass reads accross protocols") +
  labs(y = "Total number of pass reads (millions)") +
  labs(x = "Protocol") +
  #ylim(0,100) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )

##I processed the B004-11-2 results with the Artic pipeline and will try to plot the summaries here
#Starting with plotting post-PCR amplicon concentrations
pcr_conc <- read.table("PCR_concentrations.txt", header=T, sep = '\t')
pcr_conc$Variant[pcr_conc$Variant == "WT"] <- "B.1"
pcr_conc$Variant[pcr_conc$Variant == "UK"] <- "B.1.1.7"
pcr_conc$Variant[pcr_conc$Variant == "SA"] <- "B.1.351"
pcr_conc$Variant[pcr_conc$Variant == "BRZ"] <- "P.1"
pcr_conc$Variant[pcr_conc$Variant == "IND"] <- "B.1.617.2"
pcr_conc$Variant[pcr_conc$Variant == "OMI"] <- "BA.1"

pcr_conc$Level[pcr_conc$Level == "High"] <- "1e-3"
pcr_conc$Level[pcr_conc$Level == "Medium"] <- "1e-5"
pcr_conc$Level[pcr_conc$Level == "Low"] <- "1e-6"
colnames(pcr_conc) <- c('Protocol','Variant','Dilution','PCR_Concentration')

levels_protocols = c('Artic_v3','Artic_v4.1','SNAP','Qiaseq','Midnight','Entebbe')
pcr_conc$Protocol <- factor(pcr_conc$Protocol , levels=levels_protocols,ordered = T)

levels_lev = c('1e-3','1e-5','1e-6')
pcr_conc$Dilution <- factor(pcr_conc$Dilution , levels=levels_lev,ordered = T)

px <- ggplot(pcr_conc, aes(x=Variant, y=PCR_Concentration, color=Dilution)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Concentration (ng/ul)") +
  #labs(y = "Average length of gaps (as %ge of genome)") +
  labs(x = "Variants") +
  #ylim(0,10) +
  my_theme2
px + facet_grid(rows = vars(Protocol))
ggsave(filename = "Post-PCR concentration of amplicons.png", width = 25, height = 15, units = "cm")

#Then plot number of pychopped reads. For SNAP, we are plotting the number of aligned reads as we had a huge number of unalignabe reads in the Medium and Low samples
pychop_reads <- read.table("number_of_pchop_reads.txt", header=F, sep = '\t')
colnames(pychop_reads) <- c('protocol','variant','dilution','reads')
pychop_reads$variant[pychop_reads$variant == "WT"] <- "B.1"
pychop_reads$variant[pychop_reads$variant == "UK"] <- "B.1.1.7"
pychop_reads$variant[pychop_reads$variant == "SA"] <- "B.1.351"
pychop_reads$variant[pychop_reads$variant == "BRZ"] <- "P.1"
pychop_reads$variant[pychop_reads$variant == "IND"] <- "B.1.617.2"

pychop_reads$dilution[pychop_reads$dilution == "High"] <- "1e-3"
pychop_reads$dilution[pychop_reads$dilution == "Medium"] <- "1e-5"
pychop_reads$dilution[pychop_reads$dilution == "Low"] <- "1e-6"

pychop_reads$reads2 <- pychop_reads$reads/10000

levels_protocols = c('artic_v3','artic_v4.1','snap','midnight','entebbe')
pychop_reads$protocol <- factor(pychop_reads$protocol , levels=levels_protocols,ordered = T)

#levels_dilution = c('1e-3','1e-5',"1e-6")
#pychop_reads$dilution <- factor(pychop_reads$dilution , levels=levels_dilution,ordered = T)

p1 <- ggplot(pychop_reads, aes(x=variant, y=reads2, group=dilution, color=dilution)) +
  geom_point(position = position_dodge(width=0.6))+
  #scale_color_brewer(palette="Paired")+theme_minimal()+
  #ggtitle("Total number of Ns in consensus ") +
  labs(y = "Total number of reads (x10,000)") +
  labs(x = "Variant") +
  my_theme2
p1 + facet_grid(rows = vars(protocol))
ggsave(filename = "Total number of pychopped reads1.png", width = 25, height = 15, units = "cm")

# Stacked
ggplot(pychop_reads, aes(x = protocol, y = reads2, fill = variant)) + 
  geom_bar(stat = "identity", color = "black")
  #scale_fill_manual(values = c("#DADAEB", "#9E9AC8", "#6A51A3"))
# Grouped
ggplot(pychop_reads, aes(fill=variant, y=reads2, x=protocol)) + 
  geom_bar(position="dodge", stat="identity")

#You can attempt a grouped and stacked bar by following the example at https://stackoverflow.com/questions/46597278/how-to-plot-a-stacked-and-grouped-bar-chart-in-ggplot

#Total number of Ns
#Turn it into a factor with the levels in the correct order
consensus_variants$No_reads <- as.character(consensus_variants$No_reads )
levels_reads = c("100","1000","5000","10000","20000","40000","80000","100000","200000","400000","800000","1000000","1500000","2000000","max")
levels_reads = c("100","1000","5000","10000","20000","40000","80000","100000","200000","400000")
consensus_variants$No_reads <- factor(consensus_variants$No_reads , levels=levels_reads,ordered = T)

levels_reads = c('WT','UK','SA','BRZ','IND','OM')
consensus_variants$variant <- factor(consensus_variants$variant , levels=levels_reads,ordered = T)

levels_reads = c('artic_v3','artic_v4.1','entebbe','midnight','snap')
consensus_variants$protocol <- factor(consensus_variants$protocol , levels=levels_reads,ordered = T)

levels_reads = c('Yes','No')
consensus_variants$correct_lineage2 <- factor(consensus_variants$correct_lineage2 , levels=levels_reads,ordered = T)

consensus_variants$consensus_Ns2 <- consensus_variants$consensus_Ns/10000
ggplot(consensus_variants, aes(x=No_reads, y=consensus_Ns2, group=protocol, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Total number of Ns in consensus ") +
  labs(y = "Total number of Ns (x10,000)") +
  labs(x = "Number of reads sampled") +
  scale_y_continuous() +
  #xlim(100,800000) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Total number of Ns in consensus2.png", width = 25, height = 15, units = "cm")

#Number of gaps
ggplot(consensus_variants, aes(x=No_reads, y=No_gaps, group=protocol, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Total number of gaps in consensus ") +
  labs(y = "Total number of gaps (Ns>=5)") +
  labs(x = "Number of reads sampled") +
  scale_y_continuous() +
  #xlim(100,800000) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Total number of gaps in consensus.png", width = 25, height = 15, units = "cm")

#Size of gaps
ggplot(consensus_variants, aes(x=No_reads, y=ave_len_gap, group=protocol, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Average length of gaps in consensus ") +
  labs(y = "Average length of gaps (bp)") +
  labs(x = "Number of reads sampled") +
  scale_y_continuous() +
  ylim(0,3000) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Average length of gaps in consensus.png", width = 25, height = 15, units = "cm")

#Percentage of reads per amplicon, we plot the coefficient of variation
ggplot(consensus_variants, aes(x=No_reads, y=CoV_per_per_ampn, group=protocol, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Coefficient of variation in percentage of reads per amplicon") +
  labs(y = "Coefficient of variation") +
  labs(x = "Number of reads sampled") +
  scale_y_continuous() +
  #ylim(0,3000) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Coefficient of variation in percentage of reads per amplicon.png", width = 25, height = 15, units = "cm")

#Correct lineage detected
ggplot(consensus_variants, aes(x=No_reads, y=correct_lineage2, group=protocol, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) + theme_minimal() +
  scale_color_brewer(palette="Paired")+theme_minimal()+
  ggtitle("Detection of correct lineage") +
  labs(y = "Correct lineage") +
  labs(x = "Number of reads sampled") +
  #scale_y_continuous() +
  #ylim(0,3000) +
  theme(
    axis.text.x = element_text(size=14),
    axis.text.y = element_text(size=14),
    axis.title.x = element_text(size=14, face="bold"),
    axis.title.y = element_text(size=14, face="bold")
  )
ggsave(filename = "Detection of correct lineage.png", width = 25, height = 15, units = "cm")

##Now, in order to represent the data in the well, I will use facets disregard the failed code for 3D representation
#Total number of Ns
#Try 1
qplot(x=0, y=0, z=0, geom="blank") + 
  theme_void() +
  axes_3D()

data(iris)
ggplot(iris, aes(x=Petal.Width, y=Sepal.Width, z=Petal.Length, color=Species)) + 
  theme_void() +
  axes_3D() +
  stat_3D()

ggplot(consensus_variants1, aes(x=No_reads, y=consensus_Ns, z=variant, color=protocol)) + 
  theme_void() +
  axes_3D() +
  stat_3D()

#Try 2
consensus_variants1$protocol[which(consensus_variants1$protocol == "BRZ")] <- 1
consensus_variants1$protocol[which(consensus_variants1$protocol == "IND")] <- 2
consensus_variants1$protocol[which(consensus_variants1$protocol == "SA")] <- 3
consensus_variants1$protocol[which(consensus_variants1$protocol == "UK")] <- 4
consensus_variants1$protocol[which(consensus_variants1$protocol == "WT")] <- 5
consensus_variants1$correct_lineage2 <- as.factor(consensus_variants1$correct_lineage2)
consensus_variants1$variant <- as.numeric(consensus_variants1$variant)
x <- reads <- consensus_variants1$No_reads
y <- Ns <- consensus_variants1$consensus_Ns
z <- variant <- consensus_variants1$variant
points3D(x, y, z, bty = "b2",pch = 19, cex = 0.5, ticktype = "detailed", theta = 15, phi = 20)

#Try 3
data(mtcars)
mtcars$am[which(mtcars$am == 0)] <- 'Automatic'
mtcars$am[which(mtcars$am == 1)] <- 'Manual'
mtcars$am <- as.factor(mtcars$am)
fig <- plot_ly(mtcars, x = ~wt, y = ~hp, z = ~qsec, color = ~am, colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Weight'),
                                   yaxis = list(title = 'Gross horsepower'),
                                   zaxis = list(title = '1/4 mile time')))
fig

fig <- plot_ly(consensus_variants1, x = ~variant, y = ~No_reads, z = ~consensus_Ns, color = ~correct_lineage2, colors = c('#BF382A', '#0C4B8E'))
fig <- fig %>% add_markers()
fig <- fig %>% layout(scene = list(xaxis = list(title = 'Variant'),
                                   yaxis = list(title = 'Reads sampled'),
                                   zaxis = list(title = 'Total Ns')))
fig

#Try 4
#https://ggplot2.tidyverse.org/reference/facet_grid.html
consensus_variants <- read.table("consensus_variant2.txt", header=T, sep = '\t')
consensus_variants1 <- consensus_variants[which(consensus_variants$No_reads <= 100000),]
consensus_variants1$variant[consensus_variants1$variant == "WT"] <- "B.1"
consensus_variants1$variant[consensus_variants1$variant == "UK"] <- "B.1.1.7"
consensus_variants1$variant[consensus_variants1$variant == "SA"] <- "B.1.351"
consensus_variants1$variant[consensus_variants1$variant == "BRZ"] <- "P.1"
consensus_variants1$variant[consensus_variants1$variant == "IND"] <- "B.1.617.2"
consensus_variants1$correct_lineage2 <- as.factor(consensus_variants1$correct_lineage2)

levels_reads = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
consensus_variants1$No_reads <- factor(consensus_variants1$No_reads , levels=levels_reads,ordered = T)

#levels_variants = c('WT','UK','SA','BRZ','IND')
#levels_variants = c('B.1','B.1.1.7','B.1.351','P.1','B.1.617.2')
#consensus_variants1$variant <- factor(consensus_variants1$variant , levels=levels_variants,ordered = T)

levels_protocols = c('artic_v3','artic_v4.1','snap','midnight','entebbe')
consensus_variants1$protocol <- factor(consensus_variants1$protocol , levels=levels_protocols,ordered = T)

levels_lineage = c('Yes','No')
consensus_variants1$correct_lineage2 <- factor(consensus_variants1$correct_lineage2 , levels=levels_lineage,ordered = T)

my_theme2 = theme(
  axis.text.x = element_text(size=14),
  axis.text.y = element_text(size=14),
  axis.title.x = element_text(size=14, face="bold"),
  axis.title.y = element_text(size=14, face="bold"),
  strip.text.x = element_text(size=12, color="black",
                              face="bold"),
  strip.text.y = element_text(size=12, color="black",
                              face="bold"),
  strip.background = element_rect(colour="black", fill="white", size=1.5, linetype="solid")
)

consensus_variants1$consensus_Ns2 <- consensus_variants1$consensus_Ns/10000
levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000")
#levels_reads2 = c("20000","40000","80000","100000","200000")
consensus_variants2 <- consensus_variants1[consensus_variants1$No_reads %in% levels_reads2,]
p1 <- ggplot(consensus_variants2, aes(x=No_reads, y=consensus_Ns2, group=variant, color=variant)) +
  geom_point(position = position_dodge(width=0.6))+
  #scale_color_brewer(palette="Paired")+theme_minimal()+
  #ggtitle("Total number of Ns in consensus ") +
  labs(y = "Total number of Ns (x10,000)") +
  labs(x = "Number of reads sampled") +
  my_theme2
p1 + facet_grid(rows = vars(protocol))
ggsave(filename = "Total number of Ns in consensus1.png", width = 25, height = 15, units = "cm")

#Total number of gaps, leave this one out, it is meaningless
#Number of gaps
consensus_variants1$No_gaps2 <- (100*consensus_variants1$No_gaps/29903)
p2 <- ggplot(consensus_variants1, aes(x=No_reads, y=No_gaps, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Total number of gaps (Ns>=5)") +
  labs(x = "Number of reads sampled") +
  my_theme2
p2 + facet_grid(rows = vars(protocol))
ggsave(filename = "Total number of gaps in consensus2.png", width = 25, height = 15, units = "cm")

#Size of gaps, change between ave_len_gap and ave_len_gap2 to plot raw or percentage values
consensus_variants1$ave_len_gap2 <- (100*consensus_variants1$ave_len_gap/29903)
levels_reads2 = c("5000","10000","20000","40000","80000","100000","200000","400000")
consensus_variants3 <- consensus_variants1[which(consensus_variants1$No_reads %in% levels_reads2),]
p3 <- ggplot(consensus_variants3, aes(x=No_reads, y=ave_len_gap, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Average length of gaps (bp)") +
  #labs(y = "Average length of gaps (as %ge of genome)") +
  labs(x = "Number of reads sampled") +
  #ylim(0,10) +
  my_theme2
p3 + facet_grid(rows = vars(protocol))
ggsave(filename = "Average length of gaps in consensus raw.png", width = 25, height = 15, units = "cm")

#Percentage of reads per amplicon, we plot the coefficient of variation
p4 <- ggplot(consensus_variants1, aes(x=No_reads, y=CoV_per_per_ampn, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Coefficient of variation") +
  labs(x = "Number of reads sampled") +
  my_theme2
p4 + facet_grid(rows = vars(protocol))
ggsave(filename = "Coefficient of variation in percentage of reads per amplicon2.png", width = 25, height = 15, units = "cm")

#Correct lineage detected
p5 <- ggplot(consensus_variants1, aes(x=No_reads, y=correct_lineage2, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Correct lineage called?") +
  labs(x = "Number of reads sampled") +
  my_theme2
p5 + facet_grid(rows = vars(protocol))
ggsave(filename = "Correct lineage detected.png", width = 25, height = 15, units = "cm")

####
## Now, let us do the plots for the High Medium Low analysis
hi_med_low <- read.table("combined_results_for_high_med_low_analysis.txt", header = T, sep = '\t')

hi_med_low$variant[hi_med_low$variant == "WT"] <- "B.1"
hi_med_low$variant[hi_med_low$variant == "UK"] <- "B.1.1.7"
hi_med_low$variant[hi_med_low$variant == "SA"] <- "B.1.351"
hi_med_low$variant[hi_med_low$variant == "BRZ"] <- "P.1"
hi_med_low$variant[hi_med_low$variant == "IND"] <- "B.1.617.2"
hi_med_low$variant[hi_med_low$variant == "OMI"] <- "BA.1"
#hi_med_low$variant[hi_med_low$variant == "OMI"] <- "BA.1"

hi_med_low$dilution[hi_med_low$dilution == "High"] <- "1e-3"
hi_med_low$dilution[hi_med_low$dilution == "Medium"] <- "1e-5"
hi_med_low$dilution[hi_med_low$dilution == "Low"] <- "1e-6"

hi_med_low$correct_lineage <- 'No'
hi_med_low$correct_lineage[hi_med_low$lineage == hi_med_low$variant] <- "Yes"

hi_med_low1 <- hi_med_low[which(hi_med_low$No_reads %in% c(100,1000,5000,10000,20000,40000,80000,100000,200000)),]
levels_reads = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
hi_med_low1$No_reads <- factor(hi_med_low1$No_reads , levels=levels_reads,ordered = T)

levels_protocols = c('artic_v3','artic_v4.1','Qiaseq','snap','midnight','entebbe')
hi_med_low1$protocol <- factor(hi_med_low1$protocol , levels=levels_protocols,ordered = T)

levels_lineage = c('Yes','No')
hi_med_low1$correct_lineage <- factor(hi_med_low1$correct_lineage , levels=levels_lineage,ordered = T)

hi_med_low1$consensus_Ns2 <- hi_med_low1$consensus_Ns/10000
#levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000")
levels_reads2 = c("20000","40000","80000","100000","200000")
#consensus_variants2 <- consensus_variants1[consensus_variants1$No_reads %in% levels_reads2,]
hi_med_low1_h <- subset(hi_med_low1, dilution == "1e-6" & hi_med_low1$No_reads %in% levels_reads2)
p1 <- ggplot(hi_med_low1_h, aes(x=No_reads, y=consensus_Ns2, group=variant, color=variant)) +
  geom_point(position = position_dodge(width=0.6))+
  labs(y = "Total number of Ns (x10,000)") +
  labs(x = "Number of reads sampled") +
  my_theme2
p1 + facet_grid(rows = vars(protocol))
ggsave(filename = "Total number of Ns in consensus2 Low.png", width = 25, height = 15, units = "cm")

#Size of gaps, change between ave_len_gap and ave_len_gap2 to plot raw or percentage values
levels_reads2 = c("5000","10000","20000","40000","80000","100000","200000") #,"400000")
hi_med_low1_i <- subset(hi_med_low1, dilution == "1e-6" & hi_med_low1$No_reads %in% levels_reads2)
p3 <- ggplot(hi_med_low1_i, aes(x=No_reads, y=ave_len_gap, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Average length of gaps (bp)") +
  #labs(y = "Average length of gaps (as %ge of genome)") +
  labs(x = "Number of reads sampled") +
  #ylim(0,10) +
  my_theme2
p3 + facet_grid(rows = vars(protocol), scales="free_y")
ggsave(filename = "Average length of gaps in consensus raw Low2.png", width = 25, height = 15, units = "cm")

#Correct lineage detected
levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
hi_med_low1_j <- subset(hi_med_low1, dilution == "1e-6" & hi_med_low1$No_reads %in% levels_reads2)
p5 <- ggplot(hi_med_low1_j, aes(x=No_reads, y=correct_lineage, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "Correct lineage called?") +
  labs(x = "Number of reads sampled") +
  my_theme2
p5 + facet_grid(rows = vars(protocol))
ggsave(filename = "Correct lineage detected Low.png", width = 25, height = 15, units = "cm")

#Confidence of correct lineage detected
levels_reads2 = c("100","1000","5000","10000","20000","40000","80000","100000","200000")
hi_med_low1_k <- subset(hi_med_low1, dilution == "1e-3" & hi_med_low1$No_reads %in% levels_reads2)
p6 <- ggplot(hi_med_low1_k, aes(x=No_reads, y=ambiguity, group=variant, color=variant)) + 
  geom_point(position = position_dodge(width=0.6)) +
  labs(y = "ambiguity") +
  labs(x = "Number of reads sampled") +
  ylim(0.7,1) +
  my_theme2
p6 + facet_grid(rows = vars(protocol))
ggsave(filename = "Confidence of correct lineage detected High.png", width = 25, height = 15, units = "cm")

