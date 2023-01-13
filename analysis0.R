size = read.table('lengths.txt', header = F, sep = '\t')

boxplot(size$V1)
summary(size$V1)
sum(size$V1<5)
head(size$V1<5)

which(size$V1<5)
which(size$V1==1)

hist(size$align_identity, breaks = 50, xlab = 'Percentage alignment identity')
summary(size$align_identity)
plot(size$align_identity,size$read_identity)
plot(size$align_identity,size$read_length)
