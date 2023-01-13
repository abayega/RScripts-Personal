rm(list = ls())
setwd('/home/abayega/R/tests/test43')

boleae_6H <- read.table("Boleae_6H_transcripts_per_embryo",header=T,sep = '\t')
boleae_zygotic <- read.table("potential_early_zygotic_genes_full-exuastic_list.txt",header=F, sep = '\t')

boleae_zygotic_tpe <- subset(boleae_6H, as.character(boleae_6H$gene) %in% as.character(boleae_zygotic$V1))
dim(boleae_zygotic_tpe)

sum(boleae_zygotic_tpe$X6H_abs_emb)

#S0, the zygotic genes only produced 26744825 molecules per embryo at 6 hours AEL