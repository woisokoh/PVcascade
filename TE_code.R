library(RTransferEntropy)
library(igraph)
library(future)
library(ggplot2)
library(broom)
library(rgdal)
library(sf)
library(patchwork)
plan(multisession)

setwd("C:/Users/twokr/PVDiffuse")
df_vm <- read.csv("ts_pv_tract_vm.csv")
neighbor_vm <- read.csv("neighbors_vm.csv")
centroids_vm <- read.csv("centroids_vm.csv")

df_sea <- read.csv("ts_pv_tract_sea.csv")
neighbor_sea <- read.csv("neighbors_sea.csv")
centroids_sea <- read.csv("centroids_sea.csv")

neighbor_count_vm <- matrix(0,nrow(centroids_vm),1)
neighbor_count_sea <- matrix(0,nrow(centroids_sea),1)
for (i in 1:nrow(centroids_vm)){
  neighbor_count_vm[i,1] <- sum(neighbor_vm$src_FID == i - 1)
}
for (i in 1:nrow(centroids_sea)){
  neighbor_count_sea[i,1] <- sum(neighbor_sea$src_FID == i - 1)
}

neighbor_vm <- neighbor_vm[neighbor_vm["src_FID"]<neighbor_vm["nbr_FID"],]
neighbor_sea <- neighbor_sea[neighbor_sea["src_FID"]<neighbor_sea["nbr_FID"],]

te_cal <- function(df,neighbor){
  te_all <- matrix(999,ncol(df),ncol(df))
  pval_all <- matrix(999,ncol(df),ncol(df))
  for (i in 1:nrow(neighbor)){
   from <- neighbor["nbr_FID"][i,] + 1
   to <- neighbor["src_FID"][i,] + 1
   int <- transfer_entropy(df[,from],df[,to], quiet = TRUE,nboot = 100)$coef
   te_all[from,to] <- int[1,2]
   te_all[to,from] <- int[2,2]
   pval_all[from,to] <- int[1,4]
   pval_all[to,from] <- int[2,4]
   if (i == 0){
     print("TRANSFER ENTROPY CALCULATION")
     print("RUNNING ... ... 0 %")
   }
   if (i == floor(nrow(neighbor)/4)){print("RUNNING ... ... 25 %")}
   if (i == floor(nrow(neighbor)/2)){print("RUNNING ... ... 50 %")}
   if (i == floor(nrow(neighbor)/4*3)){print("RUNNING ... ... 75 %")}
  }
  print("DONE")
  return(list(te_all,pval_all))
}
plot_causal_flows <- function(df,pval,te,centroids,shp_tidied,th = 1, sig = 0.001){
  nodes_to_remove <- as.numeric(which(colSums(df) < th))
  signi_adj <- (pval < sig)*1
  graph <- graph_from_adjacency_matrix(signi_adj, mode = "directed")
  edge_list <- get.edgelist(graph)
  from <- centroids[edge_list[,1],]
  rownames(from) <- seq_len(nrow(from))
  to <- centroids[edge_list[,2],]
  rownames(to) <- seq_len(nrow(to))
  fl <- cbind(from,to)
  colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")
  
  from_trans <- st_as_sf(fl, coords = c("orig.lon", "orig.lat"), crs = 4326)
  from_trans <- st_transform(from_trans, 3857)$geometry
  from_trans <- st_coordinates(from_trans)
  to_trans <- st_as_sf(fl, coords = c("dest.lon", "dest.lat"), crs = 4326)
  to_trans <- st_transform(to_trans, 3857)$geometry
  to_trans <- st_coordinates(to_trans)
  fl$orig.lon <- from_trans[,1]
  fl$orig.lat <- from_trans[,2]
  fl$dest.lon <- to_trans[,1]
  fl$dest.lat <- to_trans[,2]
  
  idx <- matrix(0,nrow(edge_list),1)
  idx2 <- matrix(0,nrow(edge_list),1)
  for (i in 1:nrow(from)){
    fl[i,7] <- te[edge_list[i,1],edge_list[i,2]]
    if (sum((as.numeric(nodes_to_remove) == edge_list[i,1])*1) +
        sum((as.numeric(nodes_to_remove) == edge_list[i,2])*1) == 0){
      idx[i,1] <- 1
    }
    idx2[i,1] <- nrow(fl[fl["dest"] == fl[i,"orig"] & fl["orig"] == fl[i,"dest"],])
  }
  colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat","te")
  fl_filtered <- fl[idx == 1 & idx2 == 0,]
  med_te <- median(fl_filtered$te)
  fl_thick <- fl[fl["te"] > med_te & idx == 1 & idx2 == 0,]
  fl_filtered_two <- fl[idx == 1 & idx2 == 1,]
  fl_thick_two <- fl[idx == 1 & fl["te"] > med_te & idx2 == 1,]
  
  shorten_factor <- 0.9
  fl_filtered$short_xend <- fl_filtered$orig.lon + shorten_factor * (fl_filtered$dest.lon - fl_filtered$orig.lon)
  fl_filtered$short_yend <- fl_filtered$orig.lat + shorten_factor * (fl_filtered$dest.lat - fl_filtered$orig.lat)
  fl_filtered$short_x <- (1 - shorten_factor) * fl_filtered$short_xend + shorten_factor * fl_filtered$orig.lon
  fl_filtered$short_y <- (1 - shorten_factor) * fl_filtered$short_yend + shorten_factor * fl_filtered$orig.lat
  
  fl_thick$short_xend <- fl_thick$orig.lon + shorten_factor * (fl_thick$dest.lon - fl_thick$orig.lon)
  fl_thick$short_yend <- fl_thick$orig.lat + shorten_factor * (fl_thick$dest.lat - fl_thick$orig.lat)
  fl_thick$short_x <- (1 - shorten_factor) * fl_thick$short_xend + shorten_factor * fl_thick$orig.lon
  fl_thick$short_y <- (1 - shorten_factor) * fl_thick$short_yend + shorten_factor * fl_thick$orig.lat
  
  fl_filtered_two$short_xend <- fl_filtered_two$orig.lon + shorten_factor * (fl_filtered_two$dest.lon - fl_filtered_two$orig.lon)
  fl_filtered_two$short_yend <- fl_filtered_two$orig.lat + shorten_factor * (fl_filtered_two$dest.lat - fl_filtered_two$orig.lat)
  fl_filtered_two$short_x <- (1 - shorten_factor) * fl_filtered_two$short_xend + shorten_factor * fl_filtered_two$orig.lon
  fl_filtered_two$short_y <- (1 - shorten_factor) * fl_filtered_two$short_yend + shorten_factor * fl_filtered_two$orig.lat
  
  fl_thick_two$short_xend <- fl_thick_two$orig.lon + shorten_factor * (fl_thick_two$dest.lon - fl_thick_two$orig.lon)
  fl_thick_two$short_yend <- fl_thick_two$orig.lat + shorten_factor * (fl_thick_two$dest.lat - fl_thick_two$orig.lat)
  fl_thick_two$short_x <- (1 - shorten_factor) * fl_thick_two$short_xend + shorten_factor * fl_thick_two$orig.lon
  fl_thick_two$short_y <- (1 - shorten_factor) * fl_thick_two$short_yend + shorten_factor * fl_thick_two$orig.lat
  
  g <- ggplot() +
    geom_polygon(data = shp_tidied, 
                 aes(x = long, y = lat, group = group, fill = event),
                 color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    coord_fixed()
  
  gn <- g +
    geom_curve(data=fl_filtered,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.1, col='red',alpha=1, curvature = 0,
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
    geom_curve(data=fl_thick,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.8, col='red',alpha=1, curvature = 0,
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
    geom_curve(data=fl_filtered_two,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.1, col='yellow',alpha=1, curvature = 0) +
    geom_curve(data=fl_thick_two,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.8, col='yellow',alpha=1, curvature = 0) +
    theme_minimal()+
    theme(axis.text.x = element_blank(),   # Remove x-axis tick labels
          axis.text.y = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  red_count <- nrow(fl_filtered) + nrow(fl_thick)
  yell_count <- nrow(fl_filtered_two) + nrow(fl_thick_two)
  print(paste0("One-way causal flows: ", red_count))
  print(paste0("Two-way causal flows: ", yell_count))
  print(paste0("One-way causal flow percentage: ", red_count / (red_count + yell_count)))
  return(gn)
}

out_vm = te_cal(df_vm[13:252,],neighbor_vm) # may take some computation time.
te_vm = out_vm[[1]]
pval_vm = out_vm[[2]]

out_sea = te_cal(df_sea[13:252,],neighbor_sea) # may take some computation time.
te_sea = out_sea[[1]]
pval_sea = out_sea[[2]]

vmShp<- readOGR(dsn = "C:/Users/twokr/PVDiffuse/tl_2019_50_tract_vm_census.shp", layer = "tl_2019_50_tract_vm_census")
vmShp <- spTransform(vmShp,CRS("+init=epsg:3857"))
vmShp_tidied <- tidy(vmShp)
evSum_vm <- as.numeric(colSums(df_vm[13:252,]))
vmShp_tidied["event"] <- evSum_vm[as.numeric(vmShp_tidied$id) + 1]

seaShp<- readOGR(dsn = "C:/Users/twokr/PVDiffuse/tl_2019_53_tract_seat_census.shp", layer = "tl_2019_53_tract_seat_census")
seaShp <- spTransform(seaShp,CRS("+init=epsg:3857"))
seaShp_tidied <- tidy(seaShp)
evSum_sea <- as.numeric(colSums(df_sea[13:252,]))
seaShp_tidied["event"] <- evSum_sea[as.numeric(seaShp_tidied$id) + 1]

gn_vm <- plot_causal_flows(df_vm[13:252,],pval_vm,te_vm,centroids_vm,vmShp_tidied,th = 1, sig = 0.01)
gn_vm <- gn_vm + ggspatial::annotation_scale(location = "br", width_hint = 0.2,
                                    style = "ticks")
gn_vm
ggsave("causal_vm.png",gn_vm)

gn_sea <- plot_causal_flows(df_sea[13:252,],pval_sea,te_sea,centroids_sea,seaShp_tidied,th = 1, sig = 0.01)
gn_sea <- gn_sea+ ggspatial::annotation_scale(location = "bl", width_hint = 0.2,
                                    style = "ticks")
gn_sea
ggsave("causal_sea.png",gn_sea)

(gn_sea | gn_vm) + plot_layout(ncol = 2)


##
nodes_to_remove <- as.numeric(which(colSums(df_vm) < 1))
signi_adj <- (pval_vm < 0.01)*1
graph <- graph_from_adjacency_matrix(signi_adj, mode = "directed")
ind <- degree(graph, v = V(graph), 
       mode = "in", 
       loops = TRUE, normalized = FALSE)
outd <- degree(graph, v = V(graph), 
               mode = "out", 
               loops = TRUE, normalized = FALSE)
indp <- ind/neighbor_count_vm
outdp <- outd/neighbor_count_vm
vmShp_tidied["ind"] <- ind[as.numeric(vmShp_tidied$id) + 1]
vmShp_tidied["outd"] <- outd[as.numeric(vmShp_tidied$id) + 1]
vmShp_tidied["indp"] <- indp[as.numeric(vmShp_tidied$id) + 1]
vmShp_tidied["outdp"] <- outdp[as.numeric(vmShp_tidied$id) + 1]


g1 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = ind),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()

g2 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = outd),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g3 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = indp),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()

g4 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = outdp),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()


g1
g2
g3
g4


Hub <- hub.score(graph)$vector
Authority <- authority.score(graph)$vector
vmShp_tidied["hub"] <- Hub[as.numeric(vmShp_tidied$id) + 1]
vmShp_tidied["authority"] <- Authority[as.numeric(vmShp_tidied$id) + 1]
g5 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = hub),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g6 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = authority),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g5
g6



closeness <- closeness(graph)
vmShp_tidied["closeness"] <- closeness[as.numeric(vmShp_tidied$id) + 1]
g7 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = closeness),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g7

btw <- betweenness(graph,normalized = TRUE)
vmShp_tidied["btw"] <- btw[as.numeric(vmShp_tidied$id) + 1]
g8 <- ggplot() +
  geom_polygon(data = vmShp_tidied, 
               aes(x = long, y = lat, group = group, fill = btw),
               color = "gray") +
  scale_fill_gradient(low = "lightblue", high = "darkblue") +
  coord_fixed()
g8
