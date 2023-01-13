#https://primer-monitor.neb.com/lineages

rm(list = ls())

setwd('/home/abayega/R/tests/test51_new/')

fiocruz_F = read.table('fiocruz_F', header = F, sep = '\t')
fiocruz_R = read.table('fiocruz_R', header = F, sep = '\t')

artic_V3_F = read.table('artic_V3_F', header = F, sep = '\t')
artic_V3_R = read.table('artic_V3_R', header = F, sep = '\t')

artic_V4_F = read.table('artic_V4_F', header = F, sep = '\t')
artic_V4_R = read.table('artic_V4_R', header = F, sep = '\t')

nether_F = read.table('nether_F', header = F, sep = '\t')
nether_R = read.table('nether_R', header = F, sep = '\t')

ebbs_F = read.table('ebbs_F', header = F, sep = '\t')
ebbs_R = read.table('ebbs_R', header = F, sep = '\t')

midnight_F = read.table('midnight_F', header = F, sep = '\t')
midnight_R = read.table('midnight_R', header = F, sep = '\t')

#Load ALL mutations overlapping Artic V3 primers
overlaps_1 = read.table('overlapping_variants', header = F, sep = '\t')

#Load variants
B.1.1.519 = read.table('B.1.1.519', header = F, sep = '\t')
B.1.1.7 = read.table('B.1.1.7', header = F, sep = '\t')
B.1.351 = read.table('B.1.351', header = F, sep = '\t')
B.1.617 = read.table('B.1.617', header = F, sep = '\t')
B.1.617.2 = read.table('B.1.617.2', header = F, sep = '\t')
P.1 = read.table('P.1', header = F, sep = '\t')
A.2.5.2 = read.table('A.2.5.2', header = F, sep = '\t')
B.1.1.529 = read.table('B.1.1.529', header = F, sep = '\t')
  
head(artic_F)

#Make empty sketch
plot(1, type="n", xlab="genome pos", ylab="primers", xlim=c(0, 30000), ylim=c(0, 3))

points.default(fiocruz_F$V1, fiocruz_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(fiocruz_R$V1, fiocruz_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(artic_F$V1, artic_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(artic_R$V1, artic_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(artic_V4_F$V1, artic_V4_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(artic_V4_R$V1, artic_V4_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(nether_F$V1, nether_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(nether_R$V1, nether_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(ebbs_F$V1, ebbs_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(ebbs_R$V1, ebbs_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(midnight_F$V1, midnight_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(midnight_R$V1, midnight_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

#Plot variants
points.default(A.2.5.2$V1, A.2.5.2$V2, type="p", pch=20, col="red")
points.default(P.1$V1, P.1$V2, type="p", pch=20, col="red")
points.default(B.1.617.2$V1, B.1.617.2$V2, type="p", pch=20, col="red")
points.default(B.1.617$V1, B.1.617$V2, type="p", pch=20, col="red")
points.default(B.1.351$V1, B.1.351$V2, type="p", pch=20, col="red")
points.default(B.1.1.519$V1, B.1.1.519$V2, type="p", pch=20, col="red")
points.default(B.1.1.7$V1, B.1.1.7$V2, type="p", pch=20, col="red")
points.default(B.1.1.529$V1, B.1.1.529$V2, type="p", pch=20, col="red")

#Plot overlappings
points.default(overlaps_1$V1, overlaps_1$V2, type="p", pch=20, col="green")
