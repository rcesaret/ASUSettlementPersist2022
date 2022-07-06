#### SpatialNets.R
#### Rudolf Cesaretti, 7/4/2022

#### "SpatialNets" 
#### 
#### 
#### 
#### 

pak <- c("rgdal", "sp", "sf", "GISTools", "lwgeom", "tidyverse", "tidyr",  
         "data.table", "zoo", "igraph", "raster", "terra", "stars", "ggplot2", 
         "cowplot", "deldir", "cccd", "ggraph", "geosphere", "statnet")
# Install packages not yet installed
ip <- pak %in% rownames(installed.packages())
if (any(ip == FALSE)) {
  install.packages(pak[!ip])
}
# load packages
invisible(lapply(pak, library, character.only = TRUE))
rm(pak,ip)

#library(scales)
#library(pracma)
#library(modelsummary)

###############################################################
########################  SpatialNets  ########################
###############################################################

SpatialNets <- function(PointsList,
                        Coords,
                        CDmatsList,
                        RelNeighbor = F,
                        Gabriel = F,
                        BetaSkel = F,
                        Delaunay = F,
                        KNN = F,
                        MaxDist = F,
                        RelNeighbor = F,
                        KNN_k = 6,
                        MaxDist_t,
                        NodeStats = T,
                        PlotResults = F,
                        ggBasemap
                        
                        )

for (i in length(PointsList)){

nodes <- Pts_List[[i]]@data
xy <- data.frame(x = nodes$East, y = nodes$North)
locs <- st_as_sf(Pts_List[[i]])
CD.mat <- CD.mats[[i]]/3600 #convert to hrs




mergelist <- list()
i=1

if (RelNeighbor == T){
  rng1 <- cccd::rng(CD.mat)
  mergelist[[i]] <- rng1
  i=i+1
}
if (Gabriel == T){
  gg1 <- cccd::gg(CD.mat)
  mergelist[[i]] <- gg1
  i=i+1
}
if (BetaSkel == T){
  beta_s <- cccd::gg(CD.mat, r = 1.5)
  mergelist[[i]] <- beta_s
  i=i+1
}
if (Delaunay == T){
  deltri <- deldir(nodes[, "East"], nodes[, "North"])
  mergelist[[i]] <- deltri
  i=i+1
}
if (KNN == T){
  knn <- nng(CD.mat, k = KNN_k)
  el1 <- as.data.frame(get.edgelist(knn))
  colnames(el1) <- c("from", "to")
  knng <- graph_from_data_frame(el1)
  mergelist[[i]] <- knng
  i=i+1
}
if (MaxDist == T){
  mdnet <- network(event2dichot(CD.mat, method = "absolute",thresh = 1,leq = TRUE),directed = FALSE)
  mergelist[[i]] <- mdnet
  i=i+1
}

l=length(mergelist)

#The union of two or more graphs are created. The graphs may have identical or overlapping vertex sets.
igraph::union(..., byname = "auto")

#convert to matrices
#add together
#option for binary or weights or both?

#determine distances among all settlements using shortest path distances weighted graph 
#https://igraph.org/r/doc/distances.html 
#https://book.archnetworks.net/exploratory#WalksPathsDistance

#create weighted and unweighted undirected graphs
#output revised distance matrix
#output network
#add centrality metrics to points or df and output
#weighted binary network edgeweights --> weighted centrality?

if (NodeStats == T){
  
}


if (plot ==T){
  
  
  
  
}

}











#### Relative Neighborhood Graph
#rng1 <- rng(nodes[, c("East", "North")])
rng1 <- cccd::rng(CD.mat)

ggraph(rng1, layout = "kk") +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()

ggraph(rng1,
       layout = "manual",
       x = nodes[, c("East")],
       y = nodes[, c("North")]) +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()

#### Gabriel Graphs

gg1 <- cccd::gg(CD.mat)

ggraph(gg1, layout = "stress") +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()

ggraph(gg1,
       layout = "manual",
       x = nodes[, c("East")],
       y = nodes[, c("North")]) +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()

#### Beta Skeletons

beta_s <- cccd::gg(CD.mat, r = 1.5)

ggraph(beta_s,
       layout = "manual",
       x = nodes[, c("East")],
       y = nodes[, c("North")]) +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()


#### mapping igraphs in ggplot2 ####

nodes$ID = c(1:nrow(nodes))
list.edge.attributes(gg1)
# Extract edgelist from network object
edgelist <- get.edgelist(gg1)
# Create dataframe of beginning and ending points of edges
edges <- as.data.frame(matrix(NA, nrow(edgelist), 4))
colnames(edges) <- c("X1", "Y1", "X2", "Y2")
for (i in seq_len(nrow(edgelist))) {
  
  edges[i, ] <- c(nodes[which(nodes$ID == edgelist[i, 1]), "East"],
                  nodes[which(nodes$ID == edgelist[i, 1]), "North"],
                  nodes[which(nodes$ID == edgelist[i, 2]), "East"],
                  nodes[which(nodes$ID == edgelist[i, 2]), "North"])
  }
##ggp1 = 
ggplot() +  geom_stars(data = Hillshade.s)+
  scale_fill_gradientn(colours = c("turquoise", "black", "white"), 
                       values = scales::rescale(c(-9999, -161, 254)), guide="none")+
  geom_segment(data = edges,aes(x = X1,y = Y1,xend = X2,yend = Y2),col = "black",size = 1) +
  geom_point(data = nodes[, c("East", "North")],aes(East, North),alpha = 0.8,col = "black",
    fill = "white",shape = 21,size = 1.5,show.legend = FALSE) +theme_void()

###############

#### Delaunay Triangulation 

dt1 <- deldir(nodes[, "East"], nodes[, "North"])
plot(dt1)
# Extract Voronoi polygons for plotting
Voronoi_mapdat <- as.data.frame(dt1$dirsgs)
# Extract network for plotting
network_mapdat <- as.data.frame(dt1$delsgs)

ggplot() +  geom_stars(data = Hillshade.s)+
  scale_fill_gradientn(colours = c("turquoise", "black", "white"), 
                       values = scales::rescale(c(-9999, -161, 254)), guide="none")+
  geom_segment(data = Voronoi_mapdat,aes(x = x1,y = y1,xend = x2,yend = y2),col = "black",size = 1) +
  geom_segment(data = network_mapdat,aes(x = x1,y = y1,xend = x2,yend = y2),col = "red",size = 1) +
  geom_point(data = nodes[, c("East", "North")],aes(East, North),alpha = 0.8,col = "black",
             fill = "white",shape = 21,size = 3,show.legend = FALSE) +theme_void()

#### K-nearest Neighbors

# Calculate k=6 nearest neighbor graph
nn <- nng(CD.mat, k = 6)
el1 <- as.data.frame(get.edgelist(nn))
colnames(el1) <- c("from", "to")
g <- graph_from_data_frame(el1)

# Plot both graphs
ggraph(g, layout = "manual",
       x = nodes[, "East"], y = nodes[, "North"]) +
  geom_edge_link(width = 1) +
  geom_node_point(size = 2) +
  labs(edge_color = "K") +
  theme_graph()


#### Maximum Distance Networks


# Note we use the leq=TRUE argument here as we want nodes less than
# the threshold to count.
mdnet <- network(event2dichot(CD.mat, method = "absolute",thresh = 1,leq = TRUE),directed = FALSE)

# Plot network
ggraph(mdnet,
       layout = "manual",
       x = nodes[, "East"],
       y = nodes[, "North"]) +
  geom_edge_link() +
  geom_node_point(size = 2) +
  theme_graph()


















  ggnewscale::new_scale_fill()+
  geom_stars(data = lcp_dens.s, alpha=0.5)+
  scale_fill_gradientn(colours = c(NA, "black", "black"), 
                       values = scales::rescale(c(0, 50, 273)), guide="none")+
  geom_sf(data = CatchLims.sf, aes(geometry = geometry), color="black", size=2, alpha = 0) +
  geom_sf(data = Catch.sf, aes(geometry = geometry),color="black", size=0.35, alpha = 0) +
  geom_sf(data = Sites.sf, aes(geometry = geometry),color="black", size=0.1, fill="red") +
  labs(title= "Least Cost Path Density", subtitle = "Period 9, CL (c.AD 100-550)")+
  coord_sf()+
  theme_void()+
  theme(plot.title = element_text(hjust = 0.5, face="bold", size=16),
        plot.subtitle = element_text(hjust = 0.5, face="bold", size=14),
        plot.background = element_rect(fill = "white", colour = NA))















