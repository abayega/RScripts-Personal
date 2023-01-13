rm(list = ls())
bwam <- read.table('Bo.E.5H_all_rnaseq_combined_pass.no_adapters.readlengths', header = F, sep = '\t')
graphm <- read.table('Bo.E.5H_all_rnaseq_combined_pass.single_adapter.readlengths', header = F, sep = '\t')
lastm <- read.table('Bo.E.5H_all_rnaseq_combined_pass.trimmed_stranded.readlengths', header = F, sep = '\t')
lastn <- read.table('Bo.E.5H_all_rnaseq_combined_fail.readlengths', header = F, sep = '\t')

f <- cbind(bwam$V1,graphm$V1,lastm$V1,lastn$V1)
boxplot.matrix(f, use.cols = T, col = c('red', 'green', 'blue','yellow'), names=c('no adapter','Single adapter','Full length','fail'), ylim=c(0,20000),main = "read lengths dist'n")
summary(bwam)
summary(graphm)
summary(lastm)
summary(lastn)

sum(bwam)
sum(graphm)
sum(lastm)
sum(lastn)
sum(sum(bwam),sum(graphm),sum(lastm),sum(lastn))

length(bwam$V1)
length(graphm$V1)
length(lastm$V1)
length(lastn$V1)

