

rm(list = ls())

setwd('/home/abayega/R/tests/test51')

fiocruz_F = read.table('fiocruz_F', header = F, sep = '\t')
fiocruz_R = read.table('fiocruz_R', header = F, sep = '\t')

artic_F = read.table('artic_F', header = F, sep = '\t')
artic_R = read.table('artic_R', header = F, sep = '\t')

nether_F = read.table('nether_F', header = F, sep = '\t')
nether_R = read.table('nether_R', header = F, sep = '\t')

ebbs_F = read.table('ebbs_F', header = F, sep = '\t')
ebbs_R = read.table('ebbs_R', header = F, sep = '\t')

midnight_F = read.table('midnight_F', header = F, sep = '\t')
midnight_R = read.table('midnight_R', header = F, sep = '\t')

head(artic_F)

plot(1, type="n", xlab="genome pos", ylab="primers", xlim=c(0, 30000), ylim=c(0, 3))

points.default(fiocruz_F$V1, fiocruz_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(fiocruz_R$V1, fiocruz_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(artic_F$V1, artic_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(artic_R$V1, artic_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(nether_F$V1, nether_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(nether_R$V1, nether_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(ebbs_F$V1, ebbs_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(ebbs_R$V1, ebbs_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

points.default(midnight_F$V1, midnight_F$V2, type="p", pch=1, col="black", bg=NA, cex=1)
points.default(midnight_R$V1, midnight_R$V2, type="p", pch=23, col="blue", bg=NA, cex=1)

