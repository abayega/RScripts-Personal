rm(list=ls())

setwd('/home/abayega/R/tests/test54')
summary_file <- read.table("B002_05_5AB_summary",header=T, sep = '\t')
illumina_exp <- read.table("GCF_000347755.3_Ccap_2.1_embryos_4-8Hrs_rsem.isoforms.results",header=T, sep = '\t')
ont_panhandle <- read.table("Ccap21_B002_05_5A_pass_liqa",header=T, sep = '\t')
ont_regular <- read.table("Ccap21_B002_05_5B_pass_liqa",header=T, sep = '\t')

library(ggplot2)
library(plyr)

head(summary_file)

my.theme <- theme(axis.text = element_text(colour="black", size=12),
                  text = element_text(size=10),
                  title = element_text(size=14),
                  axis.title.x=  element_text(vjust=-0.3),
                  axis.title.y = element_text(vjust=.1),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

#Aggregate results just as in pandas
total_reads = aggregate(summary_file$Primer.type, by = list(summary_file$Primer.type), FUN = length)
total_reads = data.frame(total_reads)
colnames(total_reads) <- c("Primer_type", "total_reads")

fail_pass <- aggregate(summary_file$Primer.type, by = list(summary_file$Quality, summary_file$Primer.type), FUN = length)
colnames(fail_pass) <- c('Quality','Primer_type','Total_reads')

barcoded <- aggregate(summary_file$Primer.type, by = list(summary_file$Primer.type, summary_file$Quality, summary_file$Barcode), FUN = length)
colnames(barcoded) <- c('Primer_type','Quality','Barcode','Total_reads')

head(barcoded)

#Plot total reads
ggplot(total_reads, aes(x=Primer_type, y=total_reads)) + #, fill=condition)) + 
  geom_bar(stat="identity") + #geom_bar(position="fill", stat="identity") +
  ggtitle("Total reads from the 2 experiments") +
  xlab("Primer type") +
  ylab("Total number of reads") +
  #ylim(c(0,1)) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00")) + #999999=blue,#56B4E9=grey,#E69F00=yellow not sure though
  #scale_fill_manual(values=c("#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999")) +
  my.theme

#Plot pass fail
ggplot(fail_pass, aes(x=Primer_type, y=Total_reads, fill=Quality)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Percentage of pass and fail reads from the 2 experiments") +
  xlab("Primer type") +
  ylab("%ge of total reads") +
  #ylim(c(0,1)) +
  my.theme

#Plot all-reads barcoded
barcoded$Barcode = factor(barcoded$Barcode , levels = c('unclassified','B01'), ordered = F)
ggplot(barcoded, aes(x=Primer_type, y=Total_reads, fill=Barcode)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Percentage of barcoded reads from the 2 experiments") +
  xlab("Primer type") +
  ylab("%ge of total reads") +
  #ylim(c(0,1)) +
  my.theme

#Plot pass-reads barcoded
barcoded_pass <- subset(barcoded, Quality == "PASS")
ggplot(barcoded_pass, aes(x=Primer_type, y=Total_reads, fill=Barcode)) + 
  geom_bar(position="fill", stat="identity") +
  ggtitle("Percentage of barcoded Pass reads from the 2 experiments") +
  xlab("Primer type") +
  ylab("%ge of total reads") +
  #ylim(c(0,1)) +
  my.theme
#Get percentage of unclassified reads
> 100*562563/(4308950+562563)
[1] 11.54801
> 100*170977/(170977+11123654)
[1] 1.51379

#Plot pass read lengths
p <- ggplot(summary_file, aes(x=Primer.type, y=read.length)) + #, color=fac)
  geom_violin(fill='#A4A4A4', color="darkred") +
  #geom_boxplot(width=0.1) + theme_minimal() +
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
p


GCF_000347755.3_Ccap_2.1_genomic_B002_05_5B.all.pychop.cutadapt_minimap2.sort <- c(0.0,0.102484705332,0.179333628237,0.241079849432,0.274799061871,0.301012828324,0.325152529358,0.345202523933,0.381431898041,0.425965880163,0.487605685527,0.515474113827,0.534568452505,0.553043075459,0.574587064843,0.584677789555,0.604992780416,0.617478904626,0.631949371103,0.632423025882,0.637505946767,0.644573209919,0.647619623914,0.644917496432,0.646574244865,0.646695266791,0.651784447431,0.654799562652,0.650363483094,0.654463622478,0.662916378022,0.652790181366,0.663813609542,0.671957550516,0.66924081694,0.673219934398,0.681528715582,0.686709705958,0.68521988432,0.699239648452,0.70452288149,0.708339245324,0.706968359026,0.704744059493,0.712070059176,0.717491006819,0.715053875623,0.722315191173,0.726782569504,0.728483136221,0.731487818517,0.730617712602,0.720842062214,0.72646332201,0.736526921119,0.738722008463,0.743892565915,0.740750169013,0.748003138224,0.744898299851,0.759068298098,0.769766219025,0.789442714897,0.793077545842,0.799216278701,0.801536561141,0.809021141278,0.813544857403,0.816651782361,0.836071628287,0.840413811523,0.846293807851,0.851053308072,0.845884837205,0.847779456319,0.855948436313,0.879839833741,0.896401058316,0.899699949087,0.912117215995,0.916263260247,0.929640356222,0.929275203859,0.93592723661,0.946418585629,0.946948578201,0.96107684475,0.968444576131,0.970823282949,0.976169113535,0.982983899911,0.983758022919,0.989064208391,0.993640089139,0.997385509085,1.0,0.979305250682,0.938237086126,0.831138941517,0.56406859022)
GCF_000347755.3_Ccap_2.1_genomic_B002_05_5A.all.porechop.pychop.cutadapt_minimap2.sort <- c(0.0,0.135674834857,0.256664971363,0.347867338046,0.42174885914,0.475774510766,0.528241672937,0.576865155963,0.639253313029,0.705328507171,0.775318761655,0.820160032677,0.847608085507,0.877168823469,0.909635370969,0.919408646444,0.94049584155,0.954775055048,0.969087976333,0.970651422816,0.971081692805,0.974764764259,0.984110743961,0.983419734323,0.986697320925,0.984626274824,0.994684083679,0.99560807361,0.990844370155,0.992876751442,1.0,0.995657643885,0.997297428592,0.99858823856,0.995788509412,0.988786212325,0.992499025944,0.990477550118,0.986645767839,0.996866167197,0.98988369822,0.990203922198,0.985365863331,0.979433292785,0.98248285612,0.981372481954,0.97942337873,0.979939900998,0.980756819135,0.97932622099,0.979426352946,0.975978244598,0.970526505722,0.967415475245,0.971141177135,0.967708931275,0.969133580986,0.957943587044,0.960509344493,0.952684180836,0.95283784869,0.952180546839,0.96045481719,0.951148493708,0.950432698933,0.944668667323,0.943773428151,0.939746338987,0.933884158232,0.936478666441,0.928643588729,0.925663423779,0.92057652213,0.913831990475,0.904919254979,0.901949004084,0.898322442744,0.892722984448,0.890510167359,0.885947719222,0.877246153099,0.873080267164,0.868207509104,0.866775919553,0.859791467766,0.849776289348,0.846966646145,0.840731696919,0.829649766177,0.82904699163,0.825321289739,0.814714242234,0.801958818998,0.796473372335,0.784364345489,0.768800270455,0.744385918473,0.699266062504,0.587412719138,0.338582825486)


pdf("B002_05_5AB_geneBody_coverage.geneBodyCoverage.curves.pdf")
png("B002_05_5AB_geneBody_coverage.geneBodyCoverage.curves.png")
x=1:100
icolor = colorRampPalette(c("#7fc97f","#beaed4","#fdc086","#ffff99","#386cb0","#f0027f"))(2)
plot(x,GCF_000347755.3_Ccap_2.1_genomic_B002_05_5B.all.pychop.cutadapt_minimap2.sort,type='l',xlab="Gene body percentile (5'->3')", ylab="Coverage",lwd=0.8,col=icolor[1])
lines(x,GCF_000347755.3_Ccap_2.1_genomic_B002_05_5A.all.porechop.pychop.cutadapt_minimap2.sort,type='l',col=icolor[2])
legend(60,0.2,fill=icolor[1:2], legend=c('Regular primers','Panhandle primers'))
dev.off()

##Some summary stats
summary(summary_file$read.length[which(summary_file$Quality == 'PASS' & summary_file$Primer.type == 'Panhandle primer')])
summary(summary_file$read.length[which(summary_file$Quality == 'PASS' & summary_file$Primer.type == 'Regular primer')])
hist(summary_file$read.length[which(summary_file$Quality == 'PASS' & summary_file$Primer.type == 'Panhandle primer' & summary_file$read.length < 1500)],
     xlab="read length (bp)", col="darkmagenta", ylim=c(0,1300000), main="Histogram of Panhandle primers")
hist(summary_file$read.length[which(summary_file$Quality == 'PASS' & summary_file$Primer.type == 'Regular primer' & summary_file$read.length < 1500)],
     xlab="read length (bp)", col="blue", ylim=c(0,1300000), main="Histogram of Regular primers")


###Now, let us try to plot correlation of isoform expression
head(ont_panhandle)

#Get intersection of the 3 experiments
illumina_isoforms <- illumina_exp$transcript_id
panhandle_isoforms <- ont_panhandle$IsoformName
regular_isoforms <- ont_regular$IsoformName

intersect_all <- intersect(intersect(illumina_isoforms,panhandle_isoforms),regular_isoforms)
length(regular_isoforms)

illumina_exp2 <- subset(illumina_exp, (illumina_exp$transcript_id %in% intersect_all), select=c(transcript_id,TPM))
ont_panhandle2 <- subset(ont_panhandle, (ont_panhandle$IsoformName %in% intersect_all), select=c(IsoformName,RelativeAbundance))
ont_regular2 <- subset(ont_regular, (ont_regular$IsoformName %in% intersect_all), select=c(IsoformName,RelativeAbundance))

#order all of them
illumina_exp3 <- illumina_exp2[order(illumina_exp2$transcript_id),]
ont_panhandle3 <- ont_panhandle2[order(ont_panhandle2$IsoformName),]
ont_regular3 <- ont_regular2[order(ont_regular2$IsoformName),]

head(final_df)
dim(ont_regular3)

final_df <- cbind.data.frame(illumina_exp3$transcript_id, illumina_exp3$TPM, ont_panhandle3$RelativeAbundance, ont_regular3$RelativeAbundance)
colnames(final_df) <- c('transcript_id','Illumina_TPM','Panhandle_RA','Regular_RA')

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45, face="bold"),
                  axis.title.y = element_text(vjust=.2, face="bold"),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

lm_eqn = function(df1){
  m = lm(Illumina_TPM ~ Regular_RA, df1);
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
                     r2 = format(cor(final_df$Illumina_TPM,final_df$Regular_RA, method=c("spearman")), digits = 3)
                   ))
  as.character(as.expression(eq));                 
}

ggplot(data = final_df, aes(x = Panhandle_RA, y = Illumina_TPM)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  ggtitle("correlation of Illumina TPM and Panhandle primer relative abundance") +
  scale_x_continuous("Panhandle RA") +
  scale_y_continuous("Illumina TPM") +
  geom_smooth(method = "loess") +
  annotate(geom = "text", x = 0.2, y = 900, label = lm_eqn(final_df), parse = T, color = "red") +
  annotate(geom = "text", x = 0.2, y = 800, label = lm_eqn2(final_df), parse = T, color = "red") +
  my.theme

ggplot(data = final_df, aes(x = Regular_RA, y = Illumina_TPM)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  ggtitle("correlation of Illumina TPM and Regular primer relative abundance") +
  scale_x_continuous("Regular RA") +
  scale_y_continuous("Illumina TPM") +
  geom_smooth(method = "loess") +
  annotate(geom = "text", x = 0.2, y = 900, label = lm_eqn(final_df), parse = T, color = "red") +
  annotate(geom = "text", x = 0.2, y = 800, label = lm_eqn2(final_df), parse = T, color = "red") +
  my.theme

ggplot(data = final_df, aes(x = Regular_RA, y = Panhandle_RA)) +
  geom_point() +
  #geom_smooth(method = "lm") +
  ggtitle("correlation of Panhandle and Regular primer relative abundance") +
  scale_x_continuous("Regular RA") +
  scale_y_continuous("Panhandle RA") +
  geom_smooth(method = "loess") +
  annotate(geom = "text", x = 0.2, y = 0.9, label = lm_eqn(final_df), parse = T, color = "red") +
  annotate(geom = "text", x = 0.2, y = 0.8, label = lm_eqn2(final_df), parse = T, color = "red") +
  my.theme
