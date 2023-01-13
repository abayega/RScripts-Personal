#https://www.analyticsvidhya.com/blog/2016/03/practical-guide-principal-component-analysis-python/
rm(list = ls())
setwd("/home/abayega/R/tests/test50")
library('gplots') #has the heatmap.2 function
library('ggplot2')
library('Matrix')
library('irlba')
library('gridExtra')

micro_abund <- read.csv("phyla_dataset_d3.csv",header=T, row.names = 1)

metadata_x = micro_abund[,1178:1183]
micro_abund_x = micro_abund[,1:1177]

micro_abund_x1 = micro_abund_x[,(which(colSums(micro_abund_x) != 0))]
micro_abund_x1 = log(micro_abund_x1 + 1)
#prin_comp <- prcomp(log(micro_abund_x1), scale. = T, center = T)
#prin_comp <- prcomp(micro_abund_x1, scale. = T, center = T)

#Data is so sparse, we need to center it
#https://stats.stackexchange.com/questions/35185/dimensionality-reduction-svd-or-pca-on-a-large-sparse-matrix

m <- 500000
n <- 100
i <- unlist(lapply(1:m, function(i) rep(i, sample(25:50,1))))
j <- sample(1:n, length(i), replace=TRUE)
x <- sparseMatrix(i, j, x=runif(length(i)))

x = as.matrix(micro_abund_x1)
m = nrow(x)
n = ncol(x)
n_comp <- 10
system.time({
  xt.x <- crossprod(x)
  x.means <- colMeans(x)
  xt.x <- (xt.x - m * tcrossprod(x.means)) / (m-1)
  svd.0 <- irlba(xt.x, nu=0, nv=n_comp, tol=1e-10)
})

system.time(pca <- prcomp(x, center=TRUE))

''' #Seems to produce same result
get.col <- function(i) x[,i] # Emulates reading a column
system.time({
  xt.x <- matrix(numeric(), n, n)
  x.means <- rep(numeric(), n)
  for (i in 1:n) {
    i.col <- get.col(i)
    x.means[i] <- mean(i.col)
    xt.x[i,i] <- sum(i.col * i.col)
    if (i < n) {
      for (j in (i+1):n) {
        j.col <- get.col(j)
        xt.x[i,j] <- xt.x[j,i] <- sum(j.col * i.col)
      }    
    }
  }
  xt.x <- (xt.x - m * outer(x.means, x.means, `*`)) / (m-1)
  svd.0 <- svd(xt.x / m)
}
)
system.time(pca <- prcomp(x, center=TRUE))
'''

prin_comp <- pca

new_prin_comp = as.data.frame(prin_comp$x)
new_prin_comp$sample_id = metadata_x$uc_cd
#new_prin_comp$sample_id = metadata_x$col_site
#new_prin_comp$sample_id = metadata_x$stool_biopsy

#new_prin_comp$sex = edit_colnames(new_prin_comp)[(nrow(new_prin_comp)+1):(nrow(new_prin_comp)*2)]
#metadata_x$

##get some details about the results
summary(prin_comp)
names(prin_comp)
prin_comp$x[1:5,1:5]
# Eigenvalues
eig <- (prin_comp$sdev)^2
# Variances in percentage
variance <- eig*100/sum(eig)
# Cumulative variances
cumvar <- cumsum(variance)
eig.cumvar <- data.frame(eig = eig, variance = variance,
                         cumvariance = cumvar)

#ggplot(as.data.frame(prin_comp$x), aes(x=PC1,y=PC2)) + geom_point(size=4) +
#ggplot(as.data.frame(new_prin_comp), aes(x=PC1,y=PC2)) + geom_point(size=4) +
#ggplot(new_prin_comp, aes(x=PC1,y=PC2, label=sample_id, color=sex)) +
ggplot(new_prin_comp, aes(x=PC1,y=PC2, color=sample_id)) + geom_point(size=4) +
  #geom_label(aes(fill = sample_id), colour = "white", fontface = "bold")+
  theme_bw(base_size=32) +
  labs(x=paste0("PC1: ",round(variance[1],1),"%"),
       y=paste0("PC2: ",round(variance[2],1),"%")) +
  theme(legend.position="top")

biplot(prin_comp, scale = 0, xlabs=rep(".", dim(prin_comp$x)[1]), main="PCA of 1100 top differentially expressed; log then z-score")


new_prin_comp_x = as.data.frame(prin_comp$x[,1:5])
new_prin_comp_x$sample_id = metadata_x$uc_cd

plot_cor <- function(df){
  plot_list <- list()
  x = 1
  y = 1 
  z = 1
  for (i in names(df[,1:5])){
    for (j in names(df[,1:5])){
      g = (ggplot(df, aes_string(x=i,y=j, color="sample_id")) + geom_point() +
        #theme_bw(base_size=32) +
        labs(x=paste0(i,": ",round(variance[x],1),"%"),
             y=paste0(j,": ",round(variance[y],1),"%")) +
        theme(legend.position="top"))
      #print(
      #  ggplot(new_prin_comp, aes(x=i,y=j, color=sample_id)) + geom_point(size=4) +
      #    theme_bw(base_size=32) +
      #    labs(x=paste0(i,": ",round(variance[x],1),"%"),
      #         y=paste0(j,": ",round(variance[y],1),"%")) +
      #    theme(legend.position="top")
      #  )
      #print(g)
      plot_list[[z]] <- g
      y = y + 1
      z = z + 1
    }
    x = x + 1
    y = 1
  }
  grid.arrange(grobs=plot_list,ncol=5)
}
plot_cor(new_prin_comp_x)

dev.off()
#1





























