
data(fruitfly)
colnames(fruitfly) ## check if arrays are arranged in the default order
gnames <- rownames(fruitfly)
assay <- rep(c("A", "B", "C"), each = 12)
time.grp <- rep(c(1:12), 3)
size <- rep(3, nrow(fruitfly))


out1 <- mb.long(fruitfly, times=12, reps=size, rep.grp = assay, time.grp = time.grp)
summary(out1)
plotProfile(out1, type="b", gnames=gnames, legloc=c(2,15), pch=c("A","B","C"), xlab="Hour")

setwd("~/R/tests/test9")
datax <- read.table("Bo_E_all_combined.txt", header = T, sep = "\t", row.names = 1)
data2 <- cbind(datax$A_Bo_E_1H,datax$A_Bo_E_2H,datax$A_Bo_E_3H,datax$A_Bo_E_4H,datax$A_Bo_E_5H,datax$A_Bo_E_6H)
dim(data2)

#library(timecourse)
gnames <- rownames(datax)
assay <- rep(c("A", "B", "C"), each = 6)
time.grp <- rep(c(1:6), 3)
size <- rep(3, nrow(datax))


out1 <- mb.long(datax, times=6, reps=size, rep.grp = assay, time.grp = time.grp)
summary(out1)
plotProfile(out1, type="b", gnames=gnames, legloc=c(2,15), pch=c("A","B","C"), xlab="Hour")

gnames <- rownames(data2)
assay <- rep(c("A"), each = 6)
time.grp <- rep(c(1:6), 1)
size <- rep(1, nrow(data2))


out1 <- mb.long(data2, times=6, reps=size, rep.grp = assay, time.grp = time.grp)
summary(out1)
plotProfile(out1, type="b", gnames=gnames, legloc=c(2,15), pch=c("A","B","C"), xlab="Hour")
