
rm(list = ls())
setwd("/home/abayega/R/tests/test32")
library('gplots') #has the heatmap.2 function

Rcount1 <- read.table("zGalaxy304-Run1-TEannot.gff3_edited.CountSummary",header=F, sep = '\t')
Rlength1 <- read.table("zGalaxy304-Run1-TEannot.gff3_edited.LengthSummary",header=F, sep = '\t')

###counts which = number of copies
Rcount1$freq <- 1
#Find which class is more prevalent in out annotated TEs
aggregate(Rcount1$freq, by=list(Category=Rcount1$V7), FUN=sum)
#Category     x
#1     Class_I  8926
#2    Class_II  27412
#3    Class_NA   481
#4 Class_noCat 16472

#Find which order is more prevalent in out annotated TEs
aggregate(Rcount1$freq, by=list(Category=Rcount1$V8), FUN=sum)
#Category     x
#1               DIRS   204
#2           Helitron 21299
#3      Helitron|LARD    70
#4       Helitron|LTR    26
#5       Helitron|TIR   201
#6      Helitron|TRIM    21
#7               LARD  2860
#8               LINE  1823
#9                LTR  2059
#10          LTR|DIRS     5
#11           LTR|TIR    30
#12          Maverick    84
#13              MITE   653
#14             noCat 16808
#15               PLE    48
#16          PLE|LARD     9
#17 PotentialHostGene   293
#18              SINE   145
#19         SINE|LARD     4
#20               SSR   188
#21               TIR  4844
#22              TRIM  1617

#Find which family is more prevalent in out annotated TEs
aggregate(Rcount1$freq, by=list(Category=Rcount1$V9), FUN=sum)
#       Category     x
#1           complete  2144
#2         incomplete 19636
#3             LINE_?    92
#4             LINE_I   965
#5        LINE_Jockey   480
#6            LINE_R2    20
#7           LINE_RTE    73
#8       LINE_unknown   193
#9              LTR_?     3
#10       LTR_Bel-Pao   631
#11         LTR_Copia    71
#12         LTR_Gypsy  1310
#13       LTR_unknown    44
#14     noCat_unknown 16808
#15         TIR_CACTA    71
#16           TIR_hAT   468
#17        TIR_Merlin    10
#18          TIR_MuDR    14
#19             TIR_P   142
#20 TIR_PIF-Harbinger   103
#21      TIR_PiggyBac    83
#22   TIR_Tc1-Mariner  3724
#23       TIR_Transib    59
#24       TIR_unknown   170

##Now the length
Rcount1 <- read.table("zGalaxy304-Run1-TEannot.gff3_edited.CountSummary",header=F, sep = '\t')
Rlength1 <- read.table("zGalaxy304-Run1-TEannot.gff3_edited.LengthSummary",header=F, sep = '\t')

###counts which = number of copies
Rcount1$freq <- 1
#Find which class is more prevalent in out annotated TEs
##Now the length
#Find which class is more longer in out annotated TEs
aggregate(Rlength1$V11, by=list(Category=Rlength1$V7), FUN=sum)
#Category       x
#1     Class_I 4348379
#2    Class_II 11885280
#3    Class_NA  236993
#4 Class_noCat 8602884

#Find which order is longer in out annotated TEs
aggregate(Rlength1$V11, by=list(Category=Rlength1$V8), FUN=sum)
#        Category       x
#1               DIRS  118767
#2           Helitron 9015399
#3      Helitron|LARD   31011
#4       Helitron|LTR   13782
#5       Helitron|TIR   96851
#6      Helitron|TRIM    8275
#7               LARD 1314240
#8               LINE  921510
#9                LTR 1083297
#10          LTR|DIRS    2725
#11           LTR|TIR   13927
#12          Maverick   39446
#13              MITE  257292
#14             noCat 8769198
#15               PLE   21847
#16          PLE|LARD    3614
#17 PotentialHostGene  156806
#18              SINE   55106
#19         SINE|LARD    1116
#20               SSR   80187
#21               TIR 2311284
#22              TRIM  757856

#Find which family is longer in out annotated TEs
aggregate(Rlength1$V11, by=list(Category=Rlength1$V9), FUN=sum)
#           Category       x
#1           complete  909237
#2         incomplete 8341328
#3             LINE_?   47504
#4             LINE_I  496798
#5        LINE_Jockey  234758
#6            LINE_R2    6764
#7           LINE_RTE   35465
#8       LINE_unknown  100221
#9              LTR_?    2948
#10       LTR_Bel-Pao  353390
#11         LTR_Copia   31762
#12         LTR_Gypsy  675050
#13       LTR_unknown   20147
#14     noCat_unknown 8769198
#15         TIR_CACTA   34914
#16           TIR_hAT  233231
#17        TIR_Merlin    4436
#18          TIR_MuDR    6148
#19             TIR_P   68739
#20 TIR_PIF-Harbinger   48803
#21      TIR_PiggyBac   31513
#22   TIR_Tc1-Mariner 1770587
#23       TIR_Transib   28046
#24       TIR_unknown   84867
#25               NA  2737682   