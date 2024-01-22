library(RTransferEntropy)
library(igraph)
library(future)
library(ggplot2)
library(broom)
plan(multisession)

setwd("C:/Users/twokr/PVDiffuse")
df <- read.csv("ts_pv_tract.csv")
neighbor <- read.csv("neighbors.csv")
centroids <- read.csv("centroids.csv")
#df <- df[, colSums(df) >= 30]
# all
te_all <- matrix(999,ncol(df),ncol(df))
pval_all <- matrix(999,ncol(df),ncol(df))
for (i in 1:ncol(df)){
  int_neigh <- neighbor["nbr_FID"][neighbor["src_FID"] == i - 1,] + 1
  for (j in 1:length(int_neigh)){
      int <- transfer_entropy(df[13:252,i],df[13:252,int_neigh[j]], quiet = TRUE)$coef
      te_all[i,int_neigh[j]] <- int[1,2]
      te_all[int_neigh[j],i] <- int[2,2]
      pval_all[i,int_neigh[j]] <- int[1,4]
      pval_all[int_neigh[j],i] <- int[2,4]
  }
}

#2000-2009
te_first <- matrix(999,ncol(df),ncol(df))
pval_first <- matrix(999,ncol(df),ncol(df))
for (i in 1:ncol(df)){
  int_neigh <- neighbor["nbr_FID"][neighbor["src_FID"] == i - 1,] + 1
  for (j in 1:length(int_neigh)){
    int <- transfer_entropy(df[13:132,i],df[13:132,int_neigh[j]], quiet = TRUE, nboot = 100)$coef
    te_first[i,int_neigh[j]] <- int[1,2]
    te_first[int_neigh[j],i] <- int[2,2]
    pval_first[i,int_neigh[j]] <- int[1,4]
    pval_first[int_neigh[j],i] <- int[2,4]
  }
}

#2010-2019
te_sec <- matrix(999,ncol(df),ncol(df))
pval_sec <- matrix(999,ncol(df),ncol(df))
for (i in 1:ncol(df)){
  int_neigh <- neighbor["nbr_FID"][neighbor["src_FID"] == i - 1,] + 1
  for (j in 1:length(int_neigh)){
    int <- transfer_entropy(df[133:252,i],df[133:252,int_neigh[j]], quiet = TRUE, nboot = 100)$coef
    te_sec[i,int_neigh[j]] <- int[1,2]
    te_sec[int_neigh[j],i] <- int[2,2]
    pval_sec[i,int_neigh[j]] <- int[1,4]
    pval_sec[int_neigh[j],i] <- int[2,4]
  }
}

vmShp<- readOGR(dsn = "C:/Users/twokr/PVDiffuse/tl_2021_50_tract.shp", layer = "tl_2021_50_tract")

str(vmShp,max.level = 2)

names(vmShp@data)

vmShp_tidied <- tidy(vmShp)
vmShp_tidied["events"] <- 0
vmShp_tidied["events2"] <- 0
vmShp_tidied["events3"] <- 0
for (i in 1:nrow(vmShp_tidied)){
  
  vmShp_tidied[i,"events"] <- as.numeric(colSums(df[13:252,])[as.numeric(vmShp_tidied[i,"id"]) + 1])
  vmShp_tidied[i,"events2"] <- as.numeric(colSums(df[13:132,])[as.numeric(vmShp_tidied[i,"id"]) + 1])
  vmShp_tidied[i,"events3"] <- as.numeric(colSums(df[133:252,])[as.numeric(vmShp_tidied[i,"id"]) + 1])
  
}

g <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group),
               fill = "white", color = "gray") +
  coord_fixed()
g

g1 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = events),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g1

g2 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = events2),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g2

g3 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = events3),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g3


load("dat.RData")
te_all_adj <- (pval_all < 0.001)*1
graph <- graph_from_adjacency_matrix(te_all_adj, mode = "directed")
edge_list <- get.edgelist(graph)
from <- centroids[edge_list[,1],]
rownames(from) <- seq_len(nrow(from))
to <- centroids[edge_list[,2],]
rownames(to) <- seq_len(nrow(to))
fl <- cbind(from,to)
colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

nodes_to_remove <- which(colSums(df[13:252,]) < 1)
te_all_filtered <- te_all_adj[-nodes_to_remove,-nodes_to_remove]
centroids_filtered <- centroids[-nodes_to_remove,]
graph_filtered <- graph_from_adjacency_matrix(te_all_filtered, mode = "directed")
edge_list_filtered <- get.edgelist(graph_filtered)
from_filtered <- centroids_filtered[edge_list_filtered[,1],]
rownames(from_filtered) <- seq_len(nrow(from_filtered))
to_filtered <- centroids_filtered[edge_list_filtered[,2],]
rownames(to_filtered) <- seq_len(nrow(to_filtered))
fl_filtered <- cbind(from_filtered,to_filtered)
colnames(fl_filtered) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

gn <- g1 + 
  geom_curve(data=fl_filtered,aes(x=orig.lon,y=orig.lat,xend=dest.lon,yend=dest.lat),
             linewidth=0.1,col='red',alpha=1, curvature = 0.1,
             arrow = arrow(length = unit(0.005, "npc"), type = "closed")) +
  theme_minimal()
gn
ggsave("all.png",gn)


te_first_adj <- (pval_first < 0.001)*1
graph <- graph_from_adjacency_matrix(te_first_adj, mode = "directed")
edge_list <- get.edgelist(graph)
from <- centroids[edge_list[,1],]
rownames(from) <- seq_len(nrow(from))
to <- centroids[edge_list[,2],]
rownames(to) <- seq_len(nrow(to))
fl <- cbind(from,to)
colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

nodes_to_remove <- which(colSums(df[13:132,]) < 1)
te_first_filtered <- te_first_adj[-nodes_to_remove,-nodes_to_remove]
centroids_filtered <- centroids[-nodes_to_remove,]
graph_filtered <- graph_from_adjacency_matrix(te_first_filtered, mode = "directed")
edge_list_filtered <- get.edgelist(graph_filtered)
from_filtered <- centroids_filtered[edge_list_filtered[,1],]
rownames(from_filtered) <- seq_len(nrow(from_filtered))
to_filtered <- centroids_filtered[edge_list_filtered[,2],]
rownames(to_filtered) <- seq_len(nrow(to_filtered))
fl_filtered <- cbind(from_filtered,to_filtered)
colnames(fl_filtered) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

gn1 <- g2 + 
  geom_curve(data=fl_filtered,aes(x=orig.lon,y=orig.lat,xend=dest.lon,yend=dest.lat),
             linewidth=0.1,col='red',alpha=1, curvature = 0.1,
             arrow = arrow(length = unit(0.005, "npc"), type = "closed")) +
  theme_minimal()
gn1
ggsave("first.png",gn1)


te_sec_adj <- (pval_sec < 0.05)*1
graph <- graph_from_adjacency_matrix(te_sec_adj, mode = "directed")
edge_list <- get.edgelist(graph)
from <- centroids[edge_list[,1],]
rownames(from) <- seq_len(nrow(from))
to <- centroids[edge_list[,2],]
rownames(to) <- seq_len(nrow(to))
fl <- cbind(from,to)
colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

nodes_to_remove <- which(colSums(df[133:252,]) < 1)
te_sec_filtered <- te_sec_adj[-nodes_to_remove,-nodes_to_remove]
centroids_filtered <- centroids[-nodes_to_remove,]
graph_filtered <- graph_from_adjacency_matrix(te_sec_filtered, mode = "directed")
edge_list_filtered <- get.edgelist(graph_filtered)
from_filtered <- centroids_filtered[edge_list_filtered[,1],]
rownames(from_filtered) <- seq_len(nrow(from_filtered))
to_filtered <- centroids_filtered[edge_list_filtered[,2],]
rownames(to_filtered) <- seq_len(nrow(to_filtered))
fl_filtered <- cbind(from_filtered,to_filtered)
colnames(fl_filtered) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")

gn2 <- g3 + 
  geom_curve(data=fl_filtered,aes(x=orig.lon,y=orig.lat,xend=dest.lon,yend=dest.lat),
             linewidth=0.1,col='red',alpha=1, curvature = 0.1,
             arrow = arrow(length = unit(0.005, "npc"), type = "closed")) +
  theme_minimal()
gn2
ggsave("sec.png",gn2)

