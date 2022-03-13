
# install.packages("signnet")
# install.packages("ggraph")
# install.packages("igraph")
# install.packages("zoo")
# install.packages("devtools")
# devtools::install_github("furkangursoy/signed.backbones")



library(igraph)
library(signnet)
library(ggraph)
library(devtools)
library(zoo)
library(signed.backbones)



# Function reciprocity calculations
reciprocity_info <- function(g){
  g_n <- delete_edges(g, E(g)[E(g)$sign == 1])
  g_p <- delete_edges(g, E(g)[E(g)$sign == -1])
  
  res <- c(length(E(g)), length(E(g)) - reciprocity(g) * length(E(g)), reciprocity(g_p) * length(E(g_p)), reciprocity(g_n) * length(E(g_n)), reciprocity(g) * length(E(g)) - (reciprocity(g_p) * length(E(g_p)) + reciprocity(g_n) * length(E(g_n))))  
  return(res)
}


#Function for structural balance calculations
structuralbalance_info <- function(g){
  
  ppp = count_signed_triangles(g)[1]
  ppn = count_signed_triangles(g)[2]
  pnn = count_signed_triangles(g)[3]
  nnn = count_signed_triangles(g)[4]
  
  res <- c(ppp[[1]], ppn[[1]], pnn[[1]], nnn[[1]] )
  return(res)
}



data <- read.csv('input/migration.csv', header=FALSE, stringsAsFactors = FALSE, encoding = 'UTF-8')
data[data==""] <- NA
data <- zoo::na.locf(data)
names(data) <- c('year', 'source', 'target', 'value')
data = data[data$value != 0, ]
rownames(data) <- NULL
head(data)



rec_data <- data.frame(matrix(ncol = 7, nrow = 0))
sb_data  <- data.frame(matrix(ncol = 6, nrow = 0))

for (year in 2008:2020){
  for (sigma in seq(5, 100,5)){
    net  <- data[data$year == year, c('source', 'target', 'value')]
    backbone_edges <- signed.backbones::extract(net, directed = TRUE, significance_threshold = paste(sigma,'pc'), vigor_threshold = 0)
    names(backbone_edges)[3] <- 'sign'
    backbone  <- graph_from_data_frame( backbone_edges, directed = TRUE,  vertices = NULL)
    backbone_undirected <- as.undirected(backbone, mode ="collapse", edge.attr.comb = "mean")
    backbone_undirected <- delete_edges(backbone_undirected, E(backbone_undirected)[(E(backbone_undirected)$sign == 0)])
    
    rec_data <- rbind(rec_data, c(sigma, year, reciprocity_info(backbone)))
    sb_data  <- rbind(sb_data,  c(sigma, year, structuralbalance_info(backbone_undirected)))
  }
}

names(rec_data) <- c('sigma', 'year', 'edges', 'nonrec', 'recpos', 'recneg', 'posneg')
names(sb_data)  <- c('sigma', 'year', 'ppp', 'ppn', 'pnn', 'nnn')


rec_data = round(rec_data)
rec_data$reciprocated <- (rec_data$recpos + rec_data$recneg)/(rec_data$recpos + rec_data$recneg + rec_data$nonrec + rec_data$posneg)
rec_data$conflict     <- (rec_data$posneg)/(rec_data$recpos + rec_data$recneg + rec_data$nonrec + rec_data$posneg)
write.csv(rec_data, 'output/reciprocity.csv')


sb_data$sb  <- (sb_data$ppp + sb_data$pnn)/(sb_data$ppp+sb_data$ppn+sb_data$pnn+sb_data$nnn)
sb_data$wsb <- (sb_data$ppp + sb_data$pnn + sb_data$nnn)/(sb_data$ppp+sb_data$ppn+sb_data$pnn+sb_data$nnn)
write.csv(sb_data, 'output/structural_balance.csv')
