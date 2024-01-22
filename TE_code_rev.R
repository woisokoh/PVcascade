library(RTransferEntropy)
library(igraph)
library(future)
library(ggplot2)
library(broom)
library(sf)
library(ineq)
library(plyr)
library(patchwork)
library(tidyverse)
library(fs)
library(network)
library(ergm)
library(ergm.count)
library(ergmito)
plan(multisession)

# Set your directory and load pre-defined functions
setwd("C:/Users/twokr/PVDiffuse")
dump("TE_udf", file="TE_udf");source("TE_udf.R")

# Load census data
load("VT_acs.RData")
load("Seattle_ACS.RData")

# Load pv adoption data
pv_vm <- read.csv("PV_heat_vm.csv")
pv_sea <- read.csv("PV_heat_sea.csv")
pv_sea$Lat <- pv_sea$Latitude
pv_sea$Long <- pv_sea$Longitude
pv_sea <- pv_sea %>%
  dplyr::mutate(CompletedDate = as.Date(CompletedDate),
                Year = lubridate::year(CompletedDate),
                Month = lubridate::month(CompletedDate))
time <- seq(as.Date("1999-01-01"), as.Date("2019-12-01"), by = "1 month")

# Locate adoption events into each tract
pv_vm <- tract_cal(pv_vm,VT_acs_2019)
pv_sea <- tract_cal(pv_sea,seattle_tracts)

# Create matrix at monthly and tract levels
df_vm <- time_tract_df(pv_vm,VT_acs_2019)
df_sea <- time_tract_df(pv_sea, seattle_tracts)

# Find neighbors which share borders 
# For Seattle, tracts with the shortest distance < 1 km were also defined as
# neighbors as well if the tract faces water bodies.
poly <- VT_acs_2019$geometry
sf_data <- st_sf(ID = c(1:nrow(VT_acs_2019)), geometry = st_sfc(poly))
neighbor_vm  <- st_touches(sf_data)
centroid_vm <- st_centroid(sf_data)

poly <- seattle_tracts$geometry
sf_data <- st_sf(ID = c(1:nrow(seattle_tracts)), geometry = st_sfc(poly))
neighbor_sea  <- st_touches(sf_data)
centroid_sea <- st_centroid(sf_data)
neighbor_sea <- sea_adjust(neighbor_sea)

# Calculate transfer entropy
set.seed(12345)
out_vm = te_cal(df_vm[,13:252],neighbor_vm) # may take some computation time.
te_vm = out_vm[[1]]
pval_vm = out_vm[[2]]

set.seed(12345)
out_sea = te_cal(df_sea[,13:252],neighbor_sea) # may take some computation time.
te_sea = out_sea[[1]]
pval_sea = out_sea[[2]]

# Plot causal flows
gn_vm <- plot_causal_flows(df_vm[,13:252],pval_vm,te_vm,centroid_vm,VT_acs_2019,th = 14, sig = 0.01)+ 
  ggspatial::annotation_scale(location = "br", width_hint = 0.2,style = "ticks")
gn_vm
ggsave("causal_vm.png",gn_vm)

gn_sea <- plot_causal_flows(df_sea[,13:252],pval_sea,te_sea,centroid_sea,seattle_tracts,th = 14, sig = 0.01)+ 
  ggspatial::annotation_scale(location = "br", width_hint = 0.2,style = "ticks")
gn_sea
ggsave("causal_sea.png",gn_sea)

# Plot transfer entropy distributions
nodes_bye <- rowSums(df_vm) < 14 
te_vm2 <- te_vm
te_vm2[nodes_bye,] <- 0
te_vm2[,nodes_bye] <- 0
te_vm_dist <- te_vm2[pval_vm < 0.01 & te_vm2 > 0]

nodes_bye <- rowSums(df_sea) < 14 
te_sea2 <- te_sea
te_sea2[nodes_bye,] <- 0
te_sea2[,nodes_bye] <- 0
te_sea_dist <- te_sea2[pval_sea < 0.01 & te_sea2 > 0]

tot_dist <- data.frame(x = c(te_vm_dist,te_sea_dist),
                       group = rep(c("Vermont", "Seattle"), 
                                   c(length(te_vm_dist), length(te_sea_dist)))
                       )
mu <- ddply(tot_dist, "group", summarise, grp.mean=median(x))

ted <- ggplot(tot_dist, aes(x=x, fill=group)) + 
  geom_histogram(aes(y=..density..), alpha=0.5, 
                 position="identity")+
  geom_vline(data=mu, aes(xintercept=grp.mean, color=group),
             linetype="dashed",linewidth = 1)+
  theme_minimal()+
  labs(title="Transfer Entropy Density",x="TE", y = "Density")
ted
ggsave("te_dist.png",ted, width = 5,height = 2)

#graph_vm <- graph_from_adjacency_matrix(te_vm * ((pval_vm < 0.01 & te_vm > 0)*1), mode = "directed",weighted = TRUE)
int <- (pval_vm < 0.01 & te_vm > 0)*1
nodes_bye <- rowSums(df_vm) < 14 
int[nodes_bye,] <- 0
int[,nodes_bye] <- 0
graph_vm <- graph_from_adjacency_matrix(int, mode = "directed")
transitivity(graph_vm, type = "global")
#ineq(page_rank(graph_vm)$vector, type = "Gini")
ineq(degree(graph_vm, mode = "total"), type = "Gini")
ineq(degree(graph_vm, mode = "total")[degree(graph_vm, mode = "total")>0], type = "Gini")
#Lasym(page_rank(graph_vm)$vector)
Lasym(degree(graph_vm, mode = "total"))
Lasym(degree(graph_vm, mode = "total")[degree(graph_vm, mode = "total")>0])
#REAT::hoover(page_rank(graph_vm)$vector)
REAT::hoover(degree(graph_vm, mode = "total"))
REAT::hoover(degree(graph_vm, mode = "total")[degree(graph_vm, mode = "total")>0])

#graph_sea <- graph_from_adjacency_matrix(te_sea * ((pval_sea < 0.01 & te_sea > 0)*1), mode = "directed",weighted = TRUE)
int <- (pval_sea < 0.01 & te_sea > 0)*1
nodes_bye <- rowSums(df_sea) < 14 
nodes_bye <- rowSums(df_sea) < 14 
int[nodes_bye,] <- 0
int[,nodes_bye] <- 0
graph_sea <- graph_from_adjacency_matrix(int, mode = "directed")
transitivity(graph_sea, type = "global")
#ineq(page_rank(graph_sea)$vector, type = "Gini")
ineq(degree(graph_sea, mode = "total"), type = "Gini")
ineq(degree(graph_sea, mode = "total")[degree(graph_sea, mode = "total")>0], type = "Gini")
#Lasym(page_rank(graph_sea)$vector)
Lasym(degree(graph_sea, mode = "total"))
Lasym(degree(graph_sea, mode = "total")[degree(graph_sea, mode = "total")>0])
#REAT::hoover(page_rank(graph_sea)$vector)
REAT::hoover(degree(graph_sea, mode = "total"))
REAT::hoover(degree(graph_sea, mode = "total")[degree(graph_sea, mode = "total")>0])

plot_net_prop(df_vm,te_vm,pval_vm,neighbor_vm,VT_acs_2019,"gg_vm.png")
plot_net_prop(df_sea,te_sea,pval_sea,neighbor_sea,seattle_tracts,"gg_sea.png")

dat <- VT_acs_2019 |> group_by(id) |>
  summarize(
    Area = Area,
    PopDensity = PopDensity,
    HomeOwn = HomeOwn,
    SingleFamily = SingleFamily,
    Edu = Edu,
    Income = Income,
    White = White,
    Wealth = Wealth,
    Gini = Gini
  )

mat <- (pval_vm < 0.01 & te_vm > 0)*1
del <- which(colSums(mat) > 0 | rowSums(mat) > 0)
mat2 <- mat[del,del]
dat2 <- dat[del,]

net <- network(mat2, directed = TRUE)
net
net %e% "flow" <- mat2
net %v% "PopDensity" <- dat2$PopDensity
net %v% "HomeOwn" <- dat2$HomeOwn
net %v% "SingleFamily" <- dat2$SingleFamily
net %v% "Edu" <- dat2$Edu
net %v% "Income" <- dat2$Income
net %v% "White" <- dat2$White
net %v% "Wealth" <- dat2$Wealth
net %v% "Gini" <- dat2$Gini

summary(net~triadcensus)

f0 <- ergm(net ~ edges)
summary(f0)

f1 <- ergm(net ~ edges + mutual + transitiveties + cyclicalties)
summary(f1)

f2 <- ergm(net ~ edges + mutual + transitiveties + cyclicalties + 
             diff("PopDensity") + diff("Edu") + 
             diff("White"))
summary(f2)

f3 <- ergm(net ~ edges + mutual + transitiveties + cyclicalties + 
             nodeicov("PopDensity") + nodeocov("PopDensity") +
             nodeicov("HomeOwn") + nodeocov("HomeOwn") +
             nodeicov("SingleFamily") + nodeocov("SingleFamily") +
             nodeicov("Edu") + nodeocov("Edu") +
             nodeicov("Income") + nodeocov("Income") +
             nodeicov("White") + nodeocov("White") +
             nodeicov("Wealth") + nodeocov("Wealth") +
             nodeicov("Gini") + nodeocov("Gini")
)
summary(f3)


poly <- seattle_tracts$geometry
sf_data <- st_sf(ID = c(1:nrow(seattle_tracts)), geometry = st_sfc(poly))
neighbor_vm  <- st_touches(sf_data)
centroid_vm <- st_centroid(sf_data)


seattle_tracts$X <- st_coordinates(centroid_sea)[,1]
seattle_tracts$Y <- st_coordinates(centroid_sea)[,2]
seattle_tracts$id <- c(1:135)
ggplot() +
  geom_sf(data = sf_data,
            aes(geometry = geometry),
            color = "gray") +
  geom_text(data = seattle_tracts, aes(x = X,y = Y, label = id),size = 3) +
  theme_minimal()

poly <- VT_acs_2019$geometry
sf_data <- st_sf(ID = c(1:nrow(VT_acs_2019)), geometry = st_sfc(poly))
neighbor_vm  <- st_touches(sf_data)
centroid_vm <- st_centroid(sf_data)

VT_acs_2019$X <- st_coordinates(centroid_vm)[,1]
VT_acs_2019$Y <- st_coordinates(centroid_vm)[,2]
VT_acs_2019$id <- c(1:248)
ggplot() +
  geom_sf(data = sf_data,
          aes(geometry = geometry),
          color = "gray") +
  geom_text(data = VT_acs_2019, aes(x = X,y = Y, label = id),size = 3) +
  theme_minimal()

dat <- VT_acs_2019 |> group_by(id) |>
  summarize(
    Area = Area,
    PopDensity = PopDensity,
    HomeOwn = HomeOwn,
    SingleFamily = SingleFamily,
    Edu = Edu,
    Income = Income,
    White = White,
    Wealth = Wealth,
    Gini = Gini
  )

eds <- get.edgelist(graph_vm)
graph_vm2 <- delete_edges(graph_vm, edges = c(which(eds[,1] == 122 & eds[,2] == 135)))
c_vm <- clusters(graph_vm2)
nn_vm <- c(104, unlist(sapply(1:10, function(degree) ego(graph_vm2, order = degree, nodes = 104, mindist = 1))))
clust1_vm <- unique(nn_vm)
clust1_vm
nn_vm <- c(135, unlist(sapply(1:30, function(degree) ego(graph_vm2, order = degree, nodes = 135, mindist = 1))))
clust2_vm <- unique(nn_vm)


mat <- (pval_vm < 0.01 & te_vm > 0)*1
nodes_bye <- rowSums(df_vm) < 14 
mat[nodes_bye,] <- 0
mat[,nodes_bye] <- 0
mat[122,135] <- 0
graph_vm2 <- graph_from_adjacency_matrix(mat,mode = "directed")
clust1_vm <- clust_detect(104,graph_vm2)
clust2_vm <- clust_detect(122,graph_vm2)
mat0 <- (pval_vm < 0.01 & te_vm > 0)*1
nodes_bye <- rowSums(df_vm) < 14 
mat0[nodes_bye,] <- 0
mat0[,nodes_bye] <- 0
mat0[122,135] <- 1
del <- which(colSums(mat0) > 0 | rowSums(mat0) > 0)
mat0 <- mat[del,del]
dat0 <- dat[del,]

sort_clust1_vm <- sort(clust1_vm)
sort_clust2_vm <- sort(clust2_vm)
mat1 <- mat[sort_clust1_vm,sort_clust1_vm]
dat1 <- dat[sort_clust1_vm,]
mat2 <- mat[sort_clust2_vm,sort_clust2_vm]
dat2 <- dat[sort_clust2_vm,]

net0 <- network(mat0, directed = TRUE)
net0
net0 %e% "flow" <- mat0
net0 %v% "PopDensity" <- dat0$PopDensity
net0 %v% "Pop" <- dat0$PopDensity * dat0$Area
net0 %v% "HomeOwn" <- dat0$HomeOwn
net0 %v% "SingleFamily" <- dat0$SingleFamily
net0 %v% "Edu" <- dat0$Edu
net0 %v% "Income" <- dat0$Income
net0 %v% "White" <- dat0$White
net0 %v% "Wealth" <- dat0$Wealth
net0 %v% "Gini" <- dat0$Gini

net1 <- network(mat1, directed = TRUE)
net1
net1 %e% "flow" <- mat1
net1 %v% "PopDensity" <- dat1$PopDensity
net1 %v% "Pop" <- dat1$PopDensity * dat1$Area
net1 %v% "HomeOwn" <- dat1$HomeOwn
net1 %v% "SingleFamily" <- dat1$SingleFamily
net1 %v% "Edu" <- dat1$Edu
net1 %v% "Income" <- dat1$Income
net1 %v% "White" <- dat1$White
net1 %v% "Wealth" <- dat1$Wealth
net1 %v% "Gini" <- dat1$Gini

net2 <- network(mat2, directed = TRUE)
net2
net2 %e% "flow" <- mat2
net2 %v% "PopDensity" <- dat2$PopDensity
net2 %v% "Pop" <- dat2$PopDensity * dat2$Area
net2 %v% "HomeOwn" <- dat2$HomeOwn
net2 %v% "SingleFamily" <- dat2$SingleFamily
net2 %v% "Edu" <- dat2$Edu
net2 %v% "Income" <- dat2$Income
net2 %v% "White" <- dat2$White
net2 %v% "Wealth" <- dat2$Wealth
net2 %v% "Gini" <- dat2$Gini

f0_1 <- ergm(net1 ~ edges)
summary(f0_1)
f0_2 <- ergm(net2 ~ edges)
summary(f0_2)

f1_1 <- ergm(net1 ~ edges + mutual + transitiveties + cyclicalties)
summary(f1_1)

f1_2 <- ergm(net2 ~ edges + mutual + transitiveties + cyclicalties)
summary(f1_2)

f2_1 <- ergm(net1 ~ edges + diff("Pop") + diff("PopDensity") + diff("HomeOwn") + 
               diff("SingleFamily") +diff("Edu") + diff("Income") + 
               diff("White")+diff("Wealth")+diff("Gini"))
summary(f2_1)

f2_2 <- ergm(net2 ~ edges + diff("Pop") + diff("PopDensity") + diff("HomeOwn") + 
               diff("SingleFamily") +diff("Edu") + diff("Income") + 
               diff("White")+diff("Wealth")+diff("Gini"))
summary(f2_2)

f3_1 <- ergm(net1 ~ edges + absdiff("Pop") + absdiff("PopDensity") + absdiff("HomeOwn") + 
               absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("White")+absdiff("Wealth")+absdiff("Gini"))
summary(f3_1)

f3_2 <- ergm(net2 ~ edges + absdiff("Pop") + absdiff("PopDensity") + absdiff("HomeOwn") + 
               absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("White")+absdiff("Wealth")+absdiff("Gini"))
summary(f3_2)

f4_1 <- ergm(net1 ~ edges + nodeocov("Pop") + nodeocov("PopDensity") + nodeocov("HomeOwn") + 
               nodeocov("SingleFamily") +nodeocov("Edu") + nodeocov("Income") + 
               nodeocov("White")+nodeocov("Wealth")+nodeocov("Gini"))
summary(f4_1)

f4_2 <- ergm(net2 ~ edges + nodeocov("Pop") + nodeocov("PopDensity") + nodeocov("HomeOwn") + 
               nodeocov("SingleFamily") +nodeocov("Edu") + nodeocov("Income") + 
               nodeocov("White")+nodeocov("Wealth")+nodeocov("Gini"))
summary(f4_2)

f5_1 <- ergm(net1 ~ edges + nodeicov("Pop") + nodeicov("PopDensity") + nodeicov("HomeOwn") + 
               nodeicov("SingleFamily") +nodeicov("Edu") + nodeicov("Income") + 
               nodeicov("White")+nodeicov("Wealth")+nodeicov("Gini"))
summary(f5_1)

f5_2 <- ergm(net2 ~ edges + nodeicov("Pop") + nodeicov("PopDensity") + nodeicov("HomeOwn") + 
               nodeicov("SingleFamily") +nodeicov("Edu") + nodeicov("Income") + 
               nodeicov("White")+nodeicov("Wealth")+nodeicov("Gini"))
summary(f5_2)

f6_1 <- ergm(net1 ~ edges + nodecov("Pop") + nodecov("PopDensity") + nodecov("HomeOwn") + 
               nodecov("SingleFamily") +nodecov("Edu") + nodecov("Income") + 
               nodecov("White")+nodecov("Wealth")+nodecov("Gini"))
summary(f6_1)

f6_2 <- ergm(net2 ~ edges + nodecov("Pop") + nodecov("PopDensity") + nodecov("HomeOwn") + 
               nodecov("SingleFamily") +nodecov("Edu") + nodecov("Income") + 
               nodecov("White")+nodecov("Wealth")+nodecov("Gini"))
summary(f6_2)

f7_1 <- ergm(net1 ~ edges + mutual + transitiveties + cyclicalties + nodecov("Pop") +
              nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_1)

f7_2 <- ergm(net2 ~ edges + mutual + transitiveties + cyclicalties + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_2)

f7_0 <- ergm(net0 ~ edges + mutual + transitiveties + cyclicalties + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_0)


## seattle
mat <- (pval_sea < 0.01 & te_sea > 0)*1
nodes_bye <- rowSums(df_sea) < 14 
mat[nodes_bye,] <- 0
mat[,nodes_bye] <- 0
graph_sea2 <- graph_from_adjacency_matrix(mat,mode = "directed")
clust1_sea <- clust_detect(33,graph_sea2)
clust2_sea <- clust_detect(40,graph_sea2)
clust3_sea <- clust_detect(107,graph_sea2)
clust4_sea <- clust_detect(113,graph_sea2)
#mat0 <- (pval_sea < 0.01 & te_sea > 0)*1
#mat0[nodes_bye,] <- 0
#mat0[,nodes_bye] <- 0
#del <- which(colSums(mat0) > 0 | rowSums(mat0) > 0)
#mat0 <- mat[del,del]
#dat0 <- dat[del,]

dat <- seattle_tracts |> group_by(id) |>
  summarize(
    Area = Area,
    PopDensity = PopDensity,
    HomeOwn = HomeOwn,
    SingleFamily = SingleFamily,
    Edu = Edu,
    Income = Income,
    White = White,
    Wealth = Poverty,
    Gini = Gini
  )

sort_clust1_sea <- sort(clust1_sea)
sort_clust2_sea <- sort(clust2_sea)
sort_clust3_sea <- sort(clust3_sea)
sort_clust4_sea <- sort(clust4_sea)
mat1 <- mat[sort_clust1_sea,sort_clust1_sea]
dat1 <- dat[sort_clust1_sea,]
mat2 <- mat[sort_clust2_sea,sort_clust2_sea]
dat2 <- dat[sort_clust2_sea,]
mat3 <- mat[sort_clust3_sea,sort_clust3_sea]
dat3 <- dat[sort_clust3_sea,]
mat4 <- mat[sort_clust4_sea,sort_clust4_sea]
dat4 <- dat[sort_clust4_sea,]

net0 <- network(mat, directed = TRUE)
net0
net0 %e% "flow" <- mat
net0 %v% "PopDensity" <- dat$PopDensity
net0 %v% "Pop" <- dat$PopDensity * dat$Area
net0 %v% "HomeOwn" <- dat$HomeOwn
net0 %v% "SingleFamily" <- dat$SingleFamily
net0 %v% "Edu" <- dat$Edu
net0 %v% "Income" <- dat$Income
net0 %v% "White" <- dat$White
net0 %v% "Wealth" <- dat$Wealth
net0 %v% "Gini" <- dat$Gini

net1 <- network(mat1, directed = TRUE)
net1
net1 %e% "flow" <- mat1
net1 %v% "PopDensity" <- dat1$PopDensity
net1 %v% "Pop" <- dat1$PopDensity * dat1$Area
net1 %v% "HomeOwn" <- dat1$HomeOwn
net1 %v% "SingleFamily" <- dat1$SingleFamily
net1 %v% "Edu" <- dat1$Edu
net1 %v% "Income" <- dat1$Income
net1 %v% "White" <- dat1$White
net1 %v% "Wealth" <- dat1$Wealth
net1 %v% "Gini" <- dat1$Gini

net2 <- network(mat2, directed = TRUE)
net2
net2 %e% "flow" <- mat2
net2 %v% "PopDensity" <- dat2$PopDensity
net2 %v% "Pop" <- dat2$PopDensity * dat2$Area
net2 %v% "HomeOwn" <- dat2$HomeOwn
net2 %v% "SingleFamily" <- dat2$SingleFamily
net2 %v% "Edu" <- dat2$Edu
net2 %v% "Income" <- dat2$Income
net2 %v% "White" <- dat2$White
net2 %v% "Wealth" <- dat2$Wealth
net2 %v% "Gini" <- dat2$Gini

net3 <- network(mat3, directed = TRUE)
net3
net3 %e% "flow" <- mat3
net3 %v% "PopDensity" <- dat3$PopDensity
net3 %v% "Pop" <- dat3$PopDensity * dat3$Area
net3 %v% "HomeOwn" <- dat3$HomeOwn
net3 %v% "SingleFamily" <- dat3$SingleFamily
net3 %v% "Edu" <- dat3$Edu
net3 %v% "Income" <- dat3$Income
net3 %v% "White" <- dat3$White
net3 %v% "Wealth" <- dat3$Wealth
net3 %v% "Gini" <- dat3$Gini

net4 <- network(mat4, directed = TRUE)
net4
net4 %e% "flow" <- mat4
net4 %v% "PopDensity" <- dat4$PopDensity
net4 %v% "Pop" <- dat4$PopDensity * dat4$Area
net4 %v% "HomeOwn" <- dat4$HomeOwn
net4 %v% "SingleFamily" <- dat4$SingleFamily
net4 %v% "Edu" <- dat4$Edu
net4 %v% "Income" <- dat4$Income
net4 %v% "White" <- dat4$White
net4 %v% "Wealth" <- dat4$Wealth
net4 %v% "Gini" <- dat4$Gini

f7_1s <- ergm(net1 ~ edges + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_1s)

f7_2 <- ergm(net2 ~ edges + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_2)

f7_3 <- ergm(net3 ~ edges + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_3)

f7_4 <- ergm(net4 ~ edges +nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_4)

f7_0 <- ergm(net0 ~ edges + mutual + transitiveties + cyclicalties + nodecov("Pop") +
               nodecov("SingleFamily") + nodecov("Income") + nodecov("Wealth") +
               absdiff("Pop") + absdiff("SingleFamily") +absdiff("Edu") + absdiff("Income") + 
               absdiff("Wealth"))
summary(f7_0)