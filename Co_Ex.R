library(ggplot2)
library(reshape2)
library(qgraph)
co_ex <- function() { 
  data = read.delim(file.choose(), header=TRUE, sep="\t")
  
  #selecting pearson correlation or cross-correlation matrix  
  n <- readline(prompt="\n \t Enter an integer \n
       \t Pearson Correlation Matrix: Enter 1 \n
       \t Cross-Correlation Matrix  : Enter 2 \n")
  
  if (n == 1) corr(data)
  if (n == 2) xcorr(data) 
}

#calculate correlation matrix
corr <- function(data){
  
  #creating color vector
  c <- vector( ,length=ncol(data))
  for(i in 1:ncol(data)){
    if ((max(data[,i]) > 2) & (min(data[,i]) > 0.5)) c[i] <- "green"
    else if ((max(data[,i]) < 2) & (min(data[,i]) < 0.5)) c[i] <- "red"
    else if ((max(data[,i]) > 2) & (min(data[,i]) < 0.5)) c[i] <- "blue"
    else c[i]  <- "black"
  } 
  cor_matrix <- cor(data)
  
  cor_matrix[,which(c == "black")] <- 0
  cor_matrix[which(c == "black"),] <- 0
  
  #asking for cut-off value
  x <- readline(prompt="Enter cut-off value(0-1): ")
  new_cor_matrix <- matrix(, ncol(data), ncol(data))
  new_cor_matrix <- cor_matrix
  
  #draw network 
  for (i in seq_len(nrow(new_cor_matrix))){
    for(j in seq_len(ncol(new_cor_matrix))){
      if (abs(new_cor_matrix[i,j]) < x) new_cor_matrix[i,j] <- 0
      else new_cor_matrix[i,j] = new_cor_matrix[i,j]
    }
  }
  new_cor_matrix[lower.tri(new_cor_matrix,diag = T)] <- 0
  newplot(cor_matrix,new_cor_matrix,c)
  
}

xcorr <- function(data)  {
  
  #create an empty cross-correlation matrix  
  cor_matrix <- matrix(, ncol(data),ncol(data)) 
  
  #creating color vector
  c <- vector( ,length=ncol(data))
  for(i in 1:ncol(data)){
    if ((max(data[,i]) > 2) & (min(data[,i]) > 0.5)) c[i] <- "green"
    else if ((max(data[,i]) < 2) & (min(data[,i]) < 0.5)) c[i] <- "red"
    else if ((max(data[,i]) > 2) & (min(data[,i]) < 0.5)) c[i] <- "blue"
    else c[i]  <- "black"
  } 
  
  #creating a cross-correlation matrix
  for (i in seq_len(ncol(data))){
    for (j in seq_len(ncol(data))){
      if (abs(min(ccf(data[,i],data[,j],plot=F)$acf)) < max(ccf(data[,i],data[,j],plot=F)$acf))
        cor_matrix[i,j] <- max(ccf(data[,i],data[,j],plot=F)$acf)
      if (abs(min(ccf(data[,i],data[,j],plot=F)$acf)) > max(ccf(data[,i],data[,j],plot=F)$acf))
        cor_matrix[i,j] <- min(ccf(data[,i],data[,j],plot=F)$acf)
    }
  }
  
  #naming the rows and columns of the new matrix
  rownames(cor_matrix) <- colnames(data)
  colnames(cor_matrix) <- colnames(data)
  
  cor_matrix[,which(c == "black")] <- 0
  cor_matrix[which(c == "black"),] <- 0
  new_cor_matrix <- matrix(, ncol(data), ncol(data))
  new_cor_matrix <- cor_matrix
  
  #asking for cut-off value
  x <- readline(prompt="Enter cut-off value(0-1): ")
  
  #draw network 
  for (i in seq_len(nrow(new_cor_matrix))){
    for(j in seq_len(ncol(new_cor_matrix))){
      if (abs(new_cor_matrix[i,j]) < x) new_cor_matrix[i,j] <- 0
      else new_cor_matrix[i,j] = new_cor_matrix[i,j]
    }
  }  
  new_cor_matrix[lower.tri(new_cor_matrix, diag = T) ] <- 0
  newplot(cor_matrix,new_cor_matrix,c)
}

#function for plotting graphs and heatmaps 
newplot <- function(cor_matrix,new_cor_matrix,c){
  a <- readline(prompt="\n \t Enter an integer \n
       \t For Correlation Heatmap   : Enter 1 \n
       \t For Co-Expression Network : Enter 2 \n
       \t For Both                  : Enter 3 \n")
  if(a == 1){
    print(qplot(x=Var1, y=Var2, data=melt(cor_matrix), fill=value, 
                geom="tile",main = "Correlation Heatmap") + scale_fill_gradient2(limits=c(-1, 1)))
  }
  if(a == 2){
    qgraph(new_cor_matrix, border.color = c, directed = FALSE,esize =10, cut = NULL, 
           maximum = 1, minimum = 0.5, edge.labels = F,border.width = 4)
    title("Co-expression network", line = 3)
  }
  if(a == 3){
    print(qplot(x=Var1, y=Var2, data=melt(cor_matrix), fill=value, 
                geom="tile",main = "Correlation Heatmap") + scale_fill_gradient2(limits=c(-1, 1)))
    readline("Press Enter for Co-Expression Network") 
    qgraph(new_cor_matrix, border.color = c, directed = FALSE,esize =10, cut = NULL, 
           maximum = 1, minimum = 0.5, edge.labels = F,border.width = 4) 
    title("Co-expression network", line = 3)
  }
}





