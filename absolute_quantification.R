rm(list = ls())
setwd("~/R/tests/test49_rep/absolute_quantification/")
expresn_table <- read.table('../C0-15H_counts_matrix.tsv.tpm_edited.tsv', header = T, row.names = 1, sep = '\t')
ercc_table <- read.table('ERCC_expression.txt', header = T, row.names = 1, sep = '\t')
correct_column_order <- read.table('correct_column_order')
#samples_filtered <- read.table('samples_filtered', header = T, row.names = 1, sep = '\t')
#delete <- read.table('delete1', header = T, row.names = 1, sep = '\t')

#correct the order of the columns
col_order <- c(correct_column_order$V1)
expresn_table <- expresn_table[, col_order[2:length(col_order)]]
ercc_table <- ercc_table[, c('transcript_id',	'Molecules_per_2ul_of_1740_dilution_used_per_embryo', col_order[2:length(col_order)])]

head(ercc_table)
'''
expresn_table1 = expresn_table[,1:4]
ercc_table1 = ercc_table[,2:6]

ercc_table3 = subset(ercc_table1, "B02C1HF" != 0 & Molecules_per_2ul_of_1740_dilution_used_per_embryo > 5000, select = c("Molecules_per_2ul_of_1740_dilution_used_per_embryo","B02C1HF"))
coef = glm(delete$B02C1HF ~ offset(log(delete$Molecules_per_2ul_of_1740_dilution_used_per_embryo)), family=poisson(link=log))$coefficients[[1]]
head(delete$B02C1HF *(2.71^(-coef) ))
'''
#Make a function that takes the ERCC table containing the expected No of molecules and the whole expression matrix
#and calculates the conversion factor for each sample then computes the absolute quantification for that sample, combines
#all columns and returns a data table
make_abs_expn <- function(ercc, whole, itemx){
  new_abs = ''
  column_names = c()
  excluded_samples = c()
  for (i in colnames(ercc)) {
    if (nchar(i) <= 12 & i %in% colnames(itemx)) { #remove bad samples
      new_ercc = subset(ercc, i != 0 & Molecules_per_2ul_of_1740_dilution_used_per_embryo > 5000, select = c("Molecules_per_2ul_of_1740_dilution_used_per_embryo",i))
      column_names = c(column_names,i)
      #coef = glm(round(new_ercc[[i]], digits = 0) ~ offset(log(Molecules_per_2ul_of_1740_dilution_used_per_embryo)), family=poisson(link=log), data = new_ercc)$coefficients[[1]]
      coef = glm(new_ercc[[i]] ~ offset(log(Molecules_per_2ul_of_1740_dilution_used_per_embryo)), family=poisson(link=log), data = new_ercc)$coefficients[[1]]
      abs = (whole[[i]] * (2.71^(-coef)))
      new_abs = cbind(new_abs,abs)
    }
    else{
      excluded_samples = c(excluded_samples,i)
    }
  }
  new_abs <- new_abs[,2:length(colnames(new_abs))]
  colnames(new_abs) <- column_names
  rownames(new_abs) <- rownames(whole)
  write.table(new_abs, file = 'C0-15H_counts_matrix.TPE.tsv', quote = FALSE, sep = '\t')
  return(excluded_samples)
}
excluded_samples_all = make_abs_expn(ercc_table,expresn_table,expresn_table)

#read the expression profile normalized to transcripts per embryo
expresn_table_tpe = read.table('C0-15H_counts_matrix.TPE.tsv', header = T, row.names = 1, sep = '\t')
expresn_table_tpe_ercc = expresn_table_tpe[1:92,]
expresn_table_tpe_genes= expresn_table_tpe[93:nrow(expresn_table_tpe),]

#create sex and timepoint factors
sample_type1 = gsub('B[0-9]+','',colnames(expresn_table_tpe))
sample_type = gsub('HF','H',sample_type1)
sample_type = gsub('HM','H',sample_type)
sample_type = gsub('HX','H',sample_type)

sample_sex = gsub('C[0-9]+H','',sample_type1)

total_transcripts = colSums(expresn_table_tpe_genes)

#We now try ggplot boxplots
library(ggplot2)

my.theme <- theme(axis.text = element_text(colour="black", size=15),
                  text = element_text(size=16),
                  title = element_text(size=19),
                  axis.title.x=  element_text(vjust=-0.45),
                  axis.title.y = element_text(vjust=.2),
                  axis.ticks = element_line(colour="black"),
                  axis.line = element_line())

total_transcripts_df = data.frame(sample_type,sample_sex,total_transcripts)
total_transcripts_df$sample_sex <- factor(total_transcripts_df$sample_sex, levels = unique(sample_sex),ordered = TRUE)
total_transcripts_df$sample_type <- factor(total_transcripts_df$sample_type, levels = unique(sample_type),ordered = TRUE)
head(total_transcripts_df)

ggplot(total_transcripts_df, aes(x=sample_type, y=total_transcripts, fill=sample_sex)) + 
  geom_boxplot() +
  ggtitle("Plotting total expressed transcripts per embryo") +
  xlab("timepoints") +
  ylab("Total transcripts per embryo (molecules)") +
  #ylim(c(0,1)) +
  #scale_fill_manual(values=c("#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00", "#999999", "#E69F00")) + #999999=blue,#56B4E9=grey,#E69F00=yellow not sure though
  #scale_fill_manual(values=c("#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999", "#56B4E9", "#999999")) +
  my.theme
ggsave(filename = "Plotting total expressed transcripts per embryo.png", width = 25, height = 15, units = "cm")

#Try plotting same graph for ERCCs  ''
expresn_table_tpe_ercc1 <- cbind(expresn_table_tpe_ercc,ercc_table$Molecules_per_2ul_of_1740_dilution_used_per_embryo)
expresn_table_tpe_ercc1 <- expresn_table_tpe_ercc1[expresn_table_tpe_ercc1$`ercc_table$Molecules_per_2ul_of_1740_dilution_used_per_embryo` > 5000,]
total_transcripts_ercc1 = colSums(expresn_table_tpe_ercc1)

total_transcripts_ercc_df = data.frame(sample_type,sample_sex,total_transcripts_ercc1[1:160])
total_transcripts_ercc_df$sample_sex <- factor(total_transcripts_ercc_df$sample_sex, levels = unique(sample_sex),ordered = TRUE)
total_transcripts_ercc_df$sample_type <- factor(total_transcripts_ercc_df$sample_type, levels = unique(sample_type),ordered = TRUE)

#ggplot(total_transcripts_ercc_df, aes(x=sample_type, y=total_transcripts_ercc, fill=sample_sex)) + 
ggplot(total_transcripts_ercc_df, aes(x=seq(1,length(sample_sex)), y=total_transcripts_ercc)) + 
  geom_bar(stat="identity") + #geom_boxplot() +
  ggtitle("Plotting total expressed ERCCs per embryo") +
  xlab("embryos") +
  ylab("Total transcripts (molecules)") +
  my.theme

#Let me try plotting un-normalized ERCCs
ercc_table_ed = ercc_table[3:ncol(ercc_table)]
ercc_table_ed = subset(ercc_table_ed, select = colnames(expresn_table_tpe))

barplot(height = colSums(ercc_table_ed), names.arg = seq(1,ncol(ercc_table_ed)),
        main = "Barplot of unnormalized ERCC TPM per embryo", col = "lightblue", ylab = "TPM", xlab = "embryos", cex.lab=1.5, cex.names = 1.5, cex.axis = 1.25)

###I will draw the correlations of the ERCCs expected to observed

plot(log(a$X1H_exp),log(a$X1H_TPM),xlab="log10(number_of_ERCC_molecules)",ylab="log10(number_of_sequenced_reads)")
glm(ercc_table$B02C1HF ~ offset(log(ercc_table$Molecules_per_2ul_.of_1740_dilution_used_per_embryo)), family=poisson(link=log))$coefficients[[1]]
new_y <- log(a$X1H_exp) - 13.99
lines(log(a$X1H_exp),new_y,col="red")








#Now, let us try to plot the correlation of expression accross timepoints using the absolute quantification results
bulk_sex <- read.table('../123', header = F, row.names = 1, sep = '\t')
bulk_sex <- read.table('../C1-15H_counts_matrix_tpm_edited_bulk_sex.tsv.averages', header = T, row.names = 1, sep = '\t')
bulk_sex <- read.table('../C1-15H_counts_matrix_tpm_edited_bulk.tsv.averages', header = T, row.names = 1, sep = '\t')


female_ids <- c('C1HF','C2HF','C3HF','C4HF','C5HF','C6HF','C7HF','C8HF','C9HF','C10HF','C11HF','C12HF','C13HF','C14HF','C15HF')
male_ids <- c('C1HM','C2HM','C3HM','C4HM','C5HM','C6HM','C7HM','C8HM','C9HM','C10HM','C11HM','C12HM','C13HM','C14HM','C15HM')
bulk_ids <- c('C1H','C2H','C3H','C4H','C5H','C6H','C7H','C8H','C9H','C10H','C11H','C12H','C13H','C14H','C15H')

C1H <- new_abs[,1:10]
expresn_table2 <- bulk_sex[91:length(rownames(bulk_sex)),]
to_take = names(expresn_table2) %in% male_ids
expresn_table2 <- expresn_table2[to_take]
  
par(mfrow=c(15,15), mar = c(2,2,2,0.1), mgp=c(1,0,0))
plot_cor <- function(df){
  x = 1
  y = 1 
  for (i in names(df)){
    for (j in names(df)){
      a = as.numeric(df[[i]])+1
      b = as.numeric(df[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      eq = paste0("r = ",corn)
      #plot(log(a),log(b), main=eq, ylab=paste('log (',x,'H)', sep = ""), xlab=paste('log (',y,'H)', sep = ""), pch = 20, cex = 0.1,xlim = c(0,20), ylim = c(0,20))
      #plot(log(a),log(b), main=eq, ylab=paste('log (',x,'H)', sep = ""), xlab=paste('log (',y,'H)', sep = ""), pch = 20, cex = 0.1)
      plot(log(a),log(b), main=eq, ylab=j, xlab=i, pch = 20, cex = 0.1)
      abline(reg, col="blue")
      abline(a=0,b=1, col="red")
      x = x + 1
    }
    y = y + 1
    x = 1
  }
}
plot_cor(as.data.frame(expresn_table2))
par(mfrow=c(1,1))

#dev.off()

#Let us try to get the correlation of the bulk samples we previously sequenced to the single embryo samples 'quasi-pooled'
C5_6H_bulk_sex <- read.table('C5_6H_counts_matrix_bulk.tsv.tpm.tsv_NewEd', header = T, row.names = 1, sep = '\t')
C1_15H_bulk_sex <- read.table('C1-15H_counts_matrix_tpm_edited_bulk.tsv.averages_NewEd', header = T, row.names = 1, sep = '\t')

C5_6H_bulk_sex2 <- C5_6H_bulk_sex[78:length(rownames(C5_6H_bulk_sex)),]
C1_15H_bulk_sex2 <- C1_15H_bulk_sex[78:length(rownames(C1_15H_bulk_sex)),]

colnames(C5_6H_bulk_sex2) = c('C5HB','C5HB2','C6HB')
to_take = c('C5HB','C6HB')
to_take = names(C5_6H_bulk_sex2) %in% to_take
C5_6H_bulk_sex2 = C5_6H_bulk_sex2[to_take]

par(mfrow=c(2,15), mar = c(2,2,2,0.1), mgp=c(1,0,0))
plot_cor <- function(df1, df2){
  x = 1
  y = 1 
  for (i in names(df1)){
    for (j in names(df2)){
      a = as.numeric(df1[[i]])+1
      b = as.numeric(df2[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      eq = paste0("r = ",corn)
      plot(log(a),log(b), main=eq, ylab=j, xlab=i, pch = 20, cex = 0.1)
      abline(reg, col="blue")
      abline(a=0,b=1, col="red")
      x = x + 1
    }
    y = y + 1
    x = 1
  }
}
plot_cor(as.data.frame(C5_6H_bulk_sex2),as.data.frame(C1_15H_bulk_sex2))


#Try to get a table with means of expression
bulk_sex <- read.table('../C1-15H_counts_matrix_tpm_edited_bulk_sex.tsv', header = F, row.names = 1, sep = '\t')
bulk_sex <- read.table('../123', header = F, row.names = 1, sep = '\t')
bulk_sex <- read.table('../C1-15H_counts_matrix.tsv.tpm_edited.tsv', header = T, row.names = 1, sep = '\t')

female_ids <- c('C1HF','C2HF','C3HF','C4HF','C5HF','C6HF','C7HF','C8HF','C9HF','C10HF','C11HF','C12HF','C13HF','C14HF','C15HF')
male_ids <- c('C1HM','C2HM','C3HM','C4HM','C5HM','C6HM','C7HM','C8HM','C9HM','C10HM','C11HM','C12HM','C13HM','C14HM','C15HM')


##This code refused to work. I wanted to use it to return a dataframe containing the means for each timepoints
'''
average_sex <- function(table){
  sex_table = ''
  x = 1
  for (i in seq(1,length(female_ids))){
    sample_x = female_ids[x]
    columns_to_take = which(table[1,]==sample_x)
    new_table = table[columns_to_take]
    #new_table = table[,columns_to_take[1]:columns_to_take[length(columns_to_take)]] #subset(table, table[1,]==sample_x)
    #new_table = data.frame(new_table)
    print(typeof(new_table))
    print(rowMeans(as.data.frame(new_table[2:3,])))
    new_table$sample_x = rowMeans(new_table[2:length(rownames(new_table)),2:ncol(new_table)])
    sex_table = cbind(sex_table,new_table$sample_x)
    x = x + 1
  }
  return(head(sex_table))
}
average_sex(bulk_sex)
'''
#Writing a function that can return for us a table of correlation coefficiences
return_correlation_table <- function(df){
  x = 1
  y = 1
  columns_x = ''
  for (i in names(df)){
    columns_x1 = c()
    for (j in names(df)){
      a = as.numeric(df[[i]])+1
      b = as.numeric(df[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      columns_x1 = c(columns_x1, corn)
      x = x + 1
    }
    y = y + 1
    x = 1
    columns_x = cbind(columns_x, columns_x1)
  }
  #print(columns_x)
  columns_x = columns_x[,2:ncol(columns_x)]
  colnames(columns_x) = colnames(df)
  rownames(columns_x) = colnames(df)
  write.table(columns_x, file="C1-15H_counts_matrix_tpm_edited2_bulk_sex.averages.correlation_coefficiencies.txt", quote = FALSE, sep='\t')
  #return(columns_x)
}
return_correlation_table(as.data.frame(expresn_table))

#Now, we shall draw a heatmap for the correlation coefficiences as this is easier to see
correlation_coefficiencies_bulk = read.table('correlation_coefficiencies_bulk.txt', header = T, row.names = 1, sep = '\t')
correlation_coefficiencies_female = read.table('correlation_coefficiencies_female.txt', header = T, row.names = 1, sep = '\t')
correlation_coefficiencies_male = read.table('correlation_coefficiencies_male.txt', header = T, row.names = 1, sep = '\t')

library('gplots')

heatmap.2(as.matrix(correlation_coefficiencies_male), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=TRUE, density.info = 'none', key.xlab = 'Spearman rho') 

return_correlation_table1 <- function(df1,df2){
  x = 1
  y = 1
  columns_x = ''
  for (i in names(df1)){
    columns_x1 = c()
    for (j in names(df2)){
      a = as.numeric(df1[[i]])+1
      b = as.numeric(df2[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      columns_x1 = c(columns_x1, corn)
      x = x + 1
    }
    y = y + 1
    x = 1
    columns_x = cbind(columns_x, columns_x1)
  }
  #print(columns_x)
  columns_x = columns_x[,2:ncol(columns_x)]
  colnames(columns_x) = colnames(df1)
  rownames(columns_x) = colnames(df2)
  write.table(columns_x, file="correlation_coefficiencies_bulk_vs_SE.txt", quote = FALSE, sep='\t')
  #return(columns_x)
}
corrn_x = return_correlation_table1(as.data.frame(C5_6H_bulk_sex2),as.data.frame(C1_15H_bulk_sex2))

correlation_coefficiencies_bulk_vs_SE = read.table('correlation_coefficiencies_bulk_vs_SE.txt', header = T, row.names = 1, sep = '\t')
heatmap.2(as.matrix(correlation_coefficiencies_bulk_vs_SE), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=TRUE, density.info = 'none', key.xlab = 'Spearman rho')


##I am trying to filter the embryos to remove 'bad samples'
##Now that we have the matrix of the single embryo expression correlation, let us remove embryos that have average correlation of less than 0.4
corrn_table = read.table('correlation_coefficiencies_complete.txt', header = T, row.names = 1, sep = '\t')
corrn_table = subset(corrn_table, select=(which(colMeans(corrn_table) >= 0.4)))

#Now, get the original expression profile and only take columns in that appear in the corrn_table
new_expressn_table = subset(expresn_table, select=colnames(corrn_table))
new_expressn_table = new_expressn_table[92:nrow(new_expressn_table),] #remove ERCCs

write.table(new_expressn_table,file = 'C1-15H_counts_matrix.tsv.tpm_edited2.tsv', quote=F, sep = '\t')

#Now, we use python to create averages of expression. We shall replot them here
new_bulk_averages = read.table('C1-15H_counts_matrix_tpm_edited2_bulk.tsv.averages', header = T, row.names = 1, sep = '\t')

return_correlation_table <- function(df){
  x = 1
  y = 1
  columns_x = ''
  for (i in names(df)){
    columns_x1 = c()
    for (j in names(df)){
      a = as.numeric(df[[i]])+1
      b = as.numeric(df[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      columns_x1 = c(columns_x1, corn)
      x = x + 1
    }
    y = y + 1
    x = 1
    columns_x = cbind(columns_x, columns_x1)
  }
  columns_x = columns_x[,2:ncol(columns_x)]
  colnames(columns_x) = colnames(df)
  rownames(columns_x) = colnames(df)
  write.table(columns_x, file="C1-15H_counts_matrix_tpm_edited2_bulk.averages.correlation_coefficiencies.txt", quote = FALSE, sep='\t')
  #return(columns_x)
}
return_correlation_table(as.data.frame(new_bulk_averages))

matrix = read.table('C1-15H_counts_matrix_tpm_edited2_bulk.averages.correlation_coefficiencies.txt', header = T, row.names = 1, sep = '\t')

heatmap.2(as.matrix(matrix), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none',
          key=TRUE, breaks=seq(0, 1, length.out=101), density.info = 'none', key.xlab = 'Spearman rho')

#Now, let us do sex plots
new_bulk_averages = read.table('C1-15H_counts_matrix_tpm_edited2_bulk_sex.tsv.averages', header = T, row.names = 1, sep = '\t')

return_correlation_table <- function(df){
  x = 1
  y = 1
  columns_x = ''
  for (i in names(df)){
    columns_x1 = c()
    for (j in names(df)){
      a = as.numeric(df[[i]])+1
      b = as.numeric(df[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      columns_x1 = c(columns_x1, corn)
      x = x + 1
    }
    y = y + 1
    x = 1
    columns_x = cbind(columns_x, columns_x1)
  }
  columns_x = columns_x[,2:ncol(columns_x)]
  colnames(columns_x) = colnames(df)
  rownames(columns_x) = colnames(df)
  write.table(columns_x, file="C1-15H_counts_matrix_tpm_edited2_bulk_sex.averages.correlation_coefficiencies.txt", quote = FALSE, sep='\t')
  #return(columns_x)
}
return_correlation_table(as.data.frame(new_bulk_averages))

matrix = read.table('C1-15H_counts_matrix_tpm_edited2_bulk_sex.averages.correlation_coefficiencies.txt', header = T, row.names = 1, sep = '\t')

#Female
female_matrix = subset(matrix, rownames(matrix) %in% female_ids, select = colnames(matrix) %in% female_ids)

heatmap.2(as.matrix(female_matrix), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none',
          key=TRUE, breaks=seq(0, 1, length.out=101), density.info = 'none', key.xlab = 'Spearman rho')

#male
male_matrix = subset(matrix, rownames(matrix) %in% male_ids, select = colnames(matrix) %in% male_ids)

heatmap.2(as.matrix(male_matrix), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none',
          key=TRUE, breaks=seq(0, 1, length.out=101), density.info = 'none', key.xlab = 'Spearman rho')


#Now I will try to do the correlation of the 5 and 6 hour timepoints with the new filtered bulk
new_bulk_averages = read.table('C1-15H_counts_matrix_tpm_edited2_bulk.tsv.averages', header = T, row.names = 1, sep = '\t')
new_bulk_averages = subset(new_bulk_averages, rownames(new_bulk_averages) %in% rownames(C5_6H_bulk_sex2))

return_correlation_table1 <- function(df1,df2){
  x = 1
  y = 1
  columns_x = ''
  for (i in names(df1)){
    columns_x1 = c()
    for (j in names(df2)){
      a = as.numeric(df1[[i]])+1
      b = as.numeric(df2[[j]])+1
      reg<-lm(log(b) ~ log(a))
      coeff=coefficients(reg)
      corn = format(cor(log(b), log(a), method=c("spearman")), digits = 2)
      columns_x1 = c(columns_x1, corn)
      x = x + 1
    }
    y = y + 1
    x = 1
    columns_x = cbind(columns_x, columns_x1)
  }
  #print(columns_x)
  columns_x = columns_x[,2:ncol(columns_x)]
  colnames(columns_x) = colnames(df1)
  rownames(columns_x) = colnames(df2)
  write.table(columns_x, file="correlation_coefficiencies_bulk_vs_SE_filtered.txt", quote = FALSE, sep='\t')
  #return(columns_x)
}
return_correlation_table1(as.data.frame(C5_6H_bulk_sex2),as.data.frame(new_bulk_averages))

correlation_coefficiencies_bulk_vs_SE_fil = read.table('correlation_coefficiencies_bulk_vs_SE_filtered.txt', header = T, row.names = 1, sep = '\t')
heatmap.2(as.matrix(correlation_coefficiencies_bulk_vs_SE_fil), Rowv=FALSE, Colv=FALSE, dendrogram='none', notecol="black", trace='none', 
          key=TRUE, breaks=seq(0, 1, length.out=101), density.info = 'none', key.xlab = 'Spearman rho')
