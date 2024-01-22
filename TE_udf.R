tract_cal <- function(pv,census){
  for (i in 1:nrow(pv)){
    for (j in 1:nrow(census)){
      if (length(st_within(st_point(c(pv$Long[i],pv$Lat[i])), census$geometry[j])[[1]]) != 0){
        if (st_within(st_point(c(pv$Long[i],pv$Lat[i])), census$geometry[j])[[1]] == 1){
          pv[i,"tract"] = j
          break
        }
      }
    }
  }
  return(pv)
}

te_cal <- function(df,neighbor){
  te_all <- matrix(999,nrow(df),nrow(df))
  pval_all <- matrix(999,nrow(df),nrow(df))
  print("TRANSFER ENTROPY CALCULATION")
  print("RUNNING ... ... 0 %")
  k <- 0
  for (i in 1:nrow(df)){
    from <- i
    for (j in 1:length(neighbor[[i]])){
      if (neighbor[[i]][j] > i){
        to <- neighbor[[i]][j]
        int <- transfer_entropy(df[from,],df[to,], quiet = TRUE,nboot = 100)$coef
        te_all[from,to] <- int[1,2]
        te_all[to,from] <- int[2,2]
        pval_all[from,to] <- int[1,4]
        pval_all[to,from] <- int[2,4]
      }
    }
    if (i == floor(nrow(df)/4 & k == 0)){
      print("RUNNING ... ... 25 %")
      k <- k + 1}
    if (i == floor(nrow(df)/2) & k == 1){
      print("RUNNING ... ... 50 %")
      k <- k + 1}
    if (i == floor(nrow(df)/4*3) & k == 2){
      print("RUNNING ... ... 75 %")
      k <- k + 1}
  }
  print("DONE")
  return(list(te_all,pval_all))
}

plot_causal_flows <- function(df,pval,te,centroids,census,th = 1, sig = 0.001){
  nodes_to_remove <- as.numeric(which(rowSums(df) < th))
  signi_adj <- (pval < sig)*1
  graph <- graph_from_adjacency_matrix(signi_adj, mode = "directed")
  edge_list <- get.edgelist(graph)
  from <- centroids[edge_list[,1],]
  from_xy <- st_coordinates(from)
  to <- centroids[edge_list[,2],]
  to_xy <- st_coordinates(to)
  fl <- data.frame(from$ID,from_xy,to$ID,to_xy)
  colnames(fl) <- c("orig","orig.lon","orig.lat","dest","dest.lon","dest.lat")
  
  idx <- matrix(0,nrow(edge_list),1)
  idx2 <- matrix(0,nrow(edge_list),1)
  for (i in 1:nrow(from)){
    fl[i,7] <- te[edge_list[i,1],edge_list[i,2]]
    if ((sum((as.numeric(nodes_to_remove) == edge_list[i,1])*1) +
        sum((as.numeric(nodes_to_remove) == edge_list[i,2])*1) == 0) &
        fl[i,7] > 0){
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
  
  events <- rowSums(df)
  census["event"] <- events
  
  g <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = event),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()
  gn <- g +
    geom_curve(data=fl_filtered,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.1, col='red',alpha=1, curvature = 0,
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
    geom_curve(data=fl_thick,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 1, col='red',alpha=1, curvature = 0,
               arrow = arrow(length = unit(0.01, "npc"), type = "closed")) +
    geom_curve(data=fl_filtered_two,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 0.1, col='yellow',alpha=1, curvature = 0) +
    geom_curve(data=fl_thick_two,aes(x=short_x,y=short_y,xend=short_xend,yend=short_yend),
               linewidth = 1, col='yellow',alpha=1, curvature = 0) +
    theme_minimal()+
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank())
  red_count <- nrow(fl_filtered)
  yell_count <- nrow(fl_filtered_two) / 2
  print(paste0("One-way causal flows: ", red_count))
  print(paste0("Two-way causal flows: ", yell_count))
  print(paste0("One-way causal flow percentage: ", red_count / (red_count + yell_count)))
  return(gn)
}

time_tract_df <- function(pv,census){
  pv <- pv[!is.na(pv$tract),]
  df <- matrix(0,nrow(census),length(time))
  for (i in 1:nrow(pv)){
    tract <- pv$tract[i]
    ym <- which(lubridate::year(time) == pv$Year[i] & lubridate::month(time) == pv$Month[i])
    df[tract,ym] <- 1
  }
  return(df)
}

plot_net_prop <- function(df,te,pval,neighbor,census,filename){
  #nodes_to_remove <- as.numeric(which(rowSums(df) < 1))
  signi_adj <- (pval < 0.01 & te > 0)*1
  graph <- graph_from_adjacency_matrix(signi_adj, mode = "directed")
  ind <- degree(graph, v = V(graph), 
                mode = "in", 
                loops = TRUE, normalized = FALSE)
  outd <- degree(graph, v = V(graph), 
                 mode = "out", 
                 loops = TRUE, normalized = FALSE)
  neighbor_count <- matrix(0,length(neighbor),1)
  for (i in 1:length(neighbor)){
    neighbor_count[i] <- length(neighbor[[i]])
  }
  indp <- ind/neighbor_count
  outdp <- outd/neighbor_count
  Hub <- hub.score(graph)$vector
  Authority <- authority.score(graph)$vector
  closeness <- closeness(graph)
  btw <- betweenness(graph,normalized = TRUE)
  
  census$Indegree <- ind
  census$Outdegree <- outd
  census$Indegree_per <- indp
  census$Outdegree_per <- outdp
  census$Hub_score <- Hub
  census$Authority_score <- Authority
  census$Closeness_centrality <- closeness
  census$Betweenness_centrality <- btw
  
  g1 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Indegree),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
      ) + 
    labs(title = "Indegree")

  g2 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Outdegree),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Outdegree")
  
  g3 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Hub_score),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Hub score")
  
  g4 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Authority_score),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Authority score")
  
  g5 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Indegree_per),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+
    theme(
    legend.title = element_blank(),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank()
  ) + 
    labs(title = "Indegree %")
  
  
  g6 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Outdegree_per),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Outdegree %")
  
  
  g7 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Closeness_centrality),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal()+ 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Closeness centrality")
  
  g8 <- ggplot() +
    geom_sf(data = census, 
            aes(geometry = geometry, fill = Betweenness_centrality),
            color = "gray") +
    scale_fill_gradient(low = "lightblue", high = "darkblue") +
    theme_minimal() + 
    theme(
      legend.title = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank()
    ) + 
    labs(title = "Betweenness centrality")
  
  ggsave(filename,gridExtra::grid.arrange(g1,g2,g3,g4,g5,g6,g7,g8, ncol = 4))
}

sea_adjust <- function(neighbor){
  river_north <- c(33,48,49,50,56,53,55,41)
  river_south <- c(58,59,61,62,69,68,63,64,65)
  river_west <- c(104,116,122)
  river_east <- c(97,117)
  for (i in 1:8){
    for (j in 1:9){
      sf_poly1 <-st_sf(geometry = st_sfc(sf_data$geometry[river_north[i]]))
      sf_poly2 <-st_sf(geometry = st_sfc(sf_data$geometry[river_south[j]]))
      if (as.numeric(st_distance(sf_poly1,sf_poly2)) <= 1000){
        neighbor[[river_north[i]]] <- c(neighbor[[river_north[i]]],river_south[j])
        neighbor[[river_south[j]]] <- c(neighbor[[river_south[j]]],river_north[i])
      }
    }
  }
  
  for (i in 1:3){
    for (j in 1:2){
      sf_poly1 <-st_sf(geometry = st_sfc(sf_data$geometry[river_west[i]]))
      sf_poly2 <-st_sf(geometry = st_sfc(sf_data$geometry[river_east[j]]))
      if (as.numeric(st_distance(sf_poly1,sf_poly2)) <= 1000){
        neighbor[[river_west[i]]] <- c(neighbor[[river_west[i]]],river_east[j])
        neighbor[[river_east[j]]] <- c(neighbor[[river_east[j]]],river_west[i])
      }
    }
  }
  return(neighbor)
}

clust_detect<- function(target,graph){
  int <- as.numeric(unique(neighbors(graph, target, mode = "all"))) 
  clust1 <- c(target,int)
  i <- 0
  while(i < 100){
    target <- int
    neigh <- c()
    for (j in 1:length(target)){
      int <- as.numeric(unique(neighbors(graph, target[j], mode = "all")))
      neigh <- c(neigh,int)
      clust1 <- c(clust1,int)
    }
    int <- unique(neigh)
    i <- i + 1
  }
  clust1 <- unique(clust1)
}