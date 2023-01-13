rm(list = ls())

outdir <- "/home/abayega/R/tests/test8"
outprefix <- "out"
print(outdir)

pngpdfdir <- sprintf("%s/pngpdf", outdir)

if (! file.exists(pngpdfdir)){
    dir.create(pngpdfdir)
}

setwd("/home/abayega/R/tests/test8")


#Bo_1H<-read.table('Bo_E_1H_sequencing_summary.txt', header = T)
#Bo_2H<-read.table('Bo_E_2H_sequencing_summary.txt', header = T)
#Bo_3H<-read.table('Bo_E_3H_sequencing_summary.txt', header = T)
#Bo_4H<-read.table('Bo_E_4H_sequencing_summary.txt', header = T)
#Bo_5H<-read.table('Bo_E_5H_sequencing_summary.txt', header = T)
#Bo_6H<-read.table('Bo_E_6H_sequencing_summary.txt', header = T)

Bo_1H<-read.table('Bo_E_1H_sequence_lengths.txt', header = T)
Bo_2H<-read.table('Bo_E_2H_sequence_lengths.txt', header = T)
Bo_3H<-read.table('Bo_E_3H_sequence_lengths.txt', header = T)
Bo_4H<-read.table('Bo_E_4H_sequence_lengths.txt', header = T)
Bo_5H<-read.table('Bo_E_5H_sequence_lengths.txt', header = T)
Bo_6H<-read.table('Bo_E_6H_sequence_lengths.txt', header = T)

Bo_1H<-read.table('Bo_E_1H_C010_10_pass_edited_readlengths.txt', header = F)
Bo_2H<-read.table('Bo_E_2H_C010_09_pass_edited_readlengths.txt', header = F)
Bo_3H<-read.table('Bo_E_3H_C010_08_pass_edited_readlengths.txt', header = F)
Bo_4H<-read.table('Bo_E_4H_C010_06_pass_edited_readlengths.txt', header = F)
Bo_5H<-read.table('Bo_E_5H_C010_05C_pass_edited_readlengths.txt', header = F)
Bo_6H<-read.table('Bo_E_6H_C010_07_pass_edited_readlengths.txt', header = F)

#Bo_1H<-read.table('x1', header = T)
#Bo_2H<-read.table('x2', header = T)
#Bo_3H<-read.table('x3', header = T)
#Bo_4H<-read.table('x4', header = T)
#Bo_5H<-read.table('x5', header = T)
#Bo_6H<-read.table('x6', header = T)

library("ggplot2") #, lib="/home/banthony/Rpackages")

pngpath <- sprintf('%s/%s_read_length_disn.png', pngpdfdir, outprefix)
print(sprintf("Generating %s", pngpath))
png(pngpath)

ggplot() + geom_density(aes(x=x), fill="red", data=data.frame(x=Bo_1H$sequence_length_template), alpha=.5) + 
  geom_density(aes(x=x), fill="blue", data=data.frame(x=Bo_2H$sequence_length_template), alpha=.5) +
  geom_density(aes(x=x), fill="green", data=data.frame(x=Bo_3H$sequence_length_template), alpha=.5) +
  geom_density(aes(x=x), fill="yellow", data=data.frame(x=Bo_4H$sequence_length_template), alpha=.5) +
  geom_density(aes(x=x), fill="purple", data=data.frame(x=Bo_5H$sequence_length_template), alpha=.5) +
  geom_density(aes(x=x), fill="orange", data=data.frame(x=Bo_6H$sequence_length_template), alpha=.5)
dev.off()

ggplot() + geom_density(aes(x=x), fill="red", data=data.frame(x=Bo_1H$sequence_length_template), alpha=.5)

a=rnorm(100,mean=10,sd=1)
b=rnorm(100,mean=13,sd=1)
hist(a,xlim=c(5,18),ylim=c(0,30),breaks=10,col=rgb(1,1,0,0.7),main="",xlab="number")
par(new=TRUE)
hist(b,xlim=c(5,18),ylim=c(0,30),breaks=10,col=rgb(0,1,1,0.4),main="",xlab="",ylab="")

b <- Bo_2H$sequence_length_template
head(as.numeric(b))
max(as.numeric(b))
b = as.numeric(b)
length(b[which(b>0)])
summary(b)

a <- as.numeric(Bo_1H$sequence_length_template)
b <- as.numeric(Bo_2H$sequence_length_template)
c <- as.numeric(Bo_3H$sequence_length_template)
d <- as.numeric(Bo_4H$sequence_length_template)
e <- as.numeric(Bo_5H$sequence_length_template)
f <- as.numeric(Bo_6H$sequence_length_template)

hist(a,breaks=30,col=rgb(1,0,0,0.5),main="",xlab="number")
par(new=TRUE)
hist(b,breaks=30,col=rgb(0,1,0,0.5),main="",xlab="",ylab="", axes=FALSE)
par(new=TRUE)
hist(c,breaks=30,col=rgb(0,0,1,0.5),main="",xlab="",ylab="", axes=FALSE)
par(new=TRUE)
hist(d,breaks=30,col=rgb(0.5,0.5,0.5,0.5),main="",xlab="",ylab="", axes=FALSE)
par(new=TRUE)
hist(e,breaks=30,col=rgb(1,0.9,0,0.5),main="",xlab="",ylab="", axes=FALSE)
par(new=TRUE)
hist(f,breaks=30,col=rgb(1,0,1,0.5),main="",xlab="",ylab="", axes=FALSE)
legend("top", legend=c("Bo_1H","Bo_2H","Bo_3H","Bo_4H","Bo_5H","Bo_6H"), fill=c("red","green","blue","gray","yellow","purple"))
dev.off()

#rgb(1,1,0,0.7)
#rgb(0,1,1,0.4)
#red: rgb(1,0,0,1)
#green: rgb(0,1,0,1)
#blue: rgb(0,0,1,1)
#gray: rgb(0.5,0.5,0.5,1)
#yellow: rgb(1,0.9,0,1)
#purple: rgb(1,0,1,1)
#orange: rgb(1,0.5,0,1)

#Make box plot of read distribution
bmmisrate <- cbind(a,b,c,d,e,f)
boxplot.matrix(bmmisrate, use.cols = T, ylim=c(0,10000),col = c('red', 'green', 'blue', 'gray', 'yellow','purple'), names=c("Bo_1H","Bo_2H","Bo_3H","Bo_4H","Bo_5H","Bo_6H"), main = 'Read length distribution')


k <- read.table('col2_template_name_length', header = F)
k <- k
summary(k$V1)

a <- as.numeric(Bo_1H$sequence_length_template)
summary(a)
Bo_1H<-read.table('x1', header = T)
a$x<-as.numeric(Bo_1H$sequence_length_template)
a$bins=findInterval(a$x,seq(0,400000,200))
v<-tapply(a$x,a$bins,FUN=sum)

summary(f)
length(a)
length(b)
length(c)
length(d)
length(e)
length(f)

cols <- sample(c("red","green","pink"),100,TRUE)
plot(rnorm(100),rnorm(100),col=cols,pch=16,cex=4)
plot(rnorm(100),rnorm(100),col=addTrans(cols,200),pch=16,cex=4)

library(ggridges)
install.packages("ggridges")
library(ggplot2)
data=data.frame(x=Bo_2H$sequence_length_template)

a<-read.table("together_dimention_sum_value_B2H")
head(a)
a_corrected<-a$V1*100
data <- data.frame(x=a_corrected,y=a$V2)
ggplot(data, aes(x, y, height = c(0,200000000))) + geom_ridgeline()


head(b)
