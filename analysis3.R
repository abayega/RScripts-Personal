size = read.table('C004_FailPass_2DnTemplate_readlength.txt')
size2 = read.table('merged.filtered_subreads.readlengths.txt')

f <- cbind(size$V1,size2$V1)
boxplot.matrix(f, use.cols = T, col = c('red', 'green'), names=c('ONT','PacBio'), main = "read lengths dist'n")


lnmisrate <- cbind(size3$V1,size2$V1,size$V1)
boxplot.matrix(lnmisrate, use.cols = T, col = c('red', 'green', 'blue'), names=c('Medfly1','Medfly2','Medfly3'), main = 'read lengths')

print(length(size$V1), length(size2$V1), length(size3$V1))

head(size)
hist(size$V1[which(size$V1>10000)], breaks = 20, main = 'combined_C004_02-04-06_passfail', xlab = 'sizes')
hist(size$V1, breaks = 50, xlim = c(0,800000), ylim = c(0,100), xlab = 'read lengths (bp)')
sum(as.numeric(size$V1))
k <- size$V1[which(size$V1>1000)]
sum(as.numeric(k))

summary(size$V1)

x = c('1000', '5000', '10000', '20000', '30000', '50000', '60000')
y = c('1864246', '14809707', '22008097', '51735451', '50362626', '93617830', '28008067')

plot(x,y, xlab = 'bins of sizes', ylab = 'total bases(bp)', main = 'Bins of bases and Total bases')

boxplot(size2$V1)
fail <- read.table('summaryFail.txt')
hist(fail$V1, breaks = seq(0,300000, 1000), xlim = c(0,20000), ylim = c(0,500))
max(fail)
pass <- read.table('summaryPass.txt')
hist(pass$V1, breaks = seq(0,300000, 1000), xlim = c(0,120000), ylim = c(0,10000))
max(pass)

#This code subsets a dataframe and returns the rows that do not have a zero. Two codes achieve the same thing
e <- matrix(c(1,2,0,1,2,3), nrow=3, ncol=2)
head(e)
m <- function(c){
  n <- match(0,c)
  return(is.na(n))
}
m(c)
data_matrix <- e
s <- subset(data_matrix,apply(data_matrix,1,m))
dim(s)

s <- subset(data_matrix,apply(data_matrix,1,function(b) length(b[b != 0])==200)) #This is one line code that performs same as code above
dim(s)

s <- data_matrix[as.vector(apply(data_matrix,1,function(i) length(which(i>0))==200)),] #This is one line code that performs same as code above