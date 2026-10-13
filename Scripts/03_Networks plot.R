#Networks plot
#Author: Laura Paltrinieri
#17/04/2026
#Modified: 07/09/2026

#Network and ternary plots
#This script builds unipartite overlap networks for each niche dimension (Temporal, Spatial, Trophic)
#and hunting-pressure treatment (Non-hunted, Hunted, Overhunted), plus the combined overlap network
#(with module hulls where modularity is significant) and ternary plots showing the relative
#contribution of each dimension to the combined overlap between primate pairs.
#All figures are exported as PDF file.

##Libraries
library(igraph)
library(dplyr)
library(tidygraph)
library(ggraph)
library(scales)
library(ggplot2)
library(ggpubr)
library(ggforce)
library(ggtern)
rm(list = ls())

#### Upload data ####
#Upload centrality data for each dimension to set node size range and plot structure
#and extract global min/max strength for consistent node sizing across plots
TR_centr <- read.csv("Results/TR_centrality.csv")
T_centr <- read.csv("Results/T_centrality.csv")
S_centr <- read.csv("Results/S_centrality.csv")
centr <- rbind(S_centr,T_centr,TR_centr)

global_limits <- range(centr$Strength_norm, na.rm = TRUE)
global_limits

#### Trophic ####
#Upload Trophic overlap matrices
TR_OH <- read.csv("Results/TR_Overhunted.csv", row.names = 1)
TR_HU <- read.csv("Results/TR_Hunted.csv", row.names = 1)
TR_NH <- read.csv("Results/TR_Non_hunted.csv", row.names = 1)

#Split trophic centrality data by hunting-pressure treatment
df_NH <- TR_centr %>%
  filter(Hunting_Pressure == "Non_hunted")
df_HU <- TR_centr %>%
  filter(Hunting_Pressure == "Hunted")
df_OH <- TR_centr %>%
  filter(Hunting_Pressure == "Overhunted")

#Plot the unipartite networks
##Non-hunted
#Convert overlap matrix to R matrix format and build the graph
TR_NH <- as.matrix(TR_NH)
graph_obj <- graph_from_adjacency_matrix(TR_NH, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges as dataframe
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Create nodes dataframe from graph vertex names
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength and Louvain module membership (pre-computed in "02_Metrics" script)
membership_NH <- readRDS("Results/TR_membership_NH.rds")
nodes_df <- nodes_df %>%
  left_join(df_NH %>% select(Genus, Strength_norm), by = c("name" = "Genus")) %>%
  left_join(membership_NH, by = c("name" = "Genus"))

#Combine nodes and edges 
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)

#Create a manual networks layout
#Build an unweighted graph, used only as reference for the fixed manual layout
g_layout <- graph_from_adjacency_matrix(
  TR_NH,
  mode = "undirected",
  weighted = FALSE,
  diag = FALSE)

#Set fixed node coordinates, kept consistent across all plots
fixed_coords <- data.frame(
  name = c("Ateles",
           "Sapajus",
           "Cebus",
           "Callicebus",
           "Saguinus",
           "Lagothrix",
           "Alouatta",
           "Saimiri",
           "Pithecia",
           "Cacajao",
           "Aotus"),
  
  x = c(3.8,
        -0.4,
        2.9252476,
        1.4599457,
        -0.3799586,
        3.3,
        3.8506908,
        -0.9401148,
        1.4,
        1.3,
        1.1066913),
  
  y = c(-0.72443006,
        -1.4,
        1.58043183,
        -0.3,
        1.20965325,
        -1.8,
        0.8,
        -0.32135395,
        -2.57215214,
        -1.5,
        2.18047752),
  
  stringsAsFactors = FALSE)

#Attach fixed coordinates to nodes
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name") %>%
  filter(!is.na(x))

#Extract node data with flipped y-coordinate, used later to draw module hulls
hull_df <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble() %>%
  mutate(y_plot = -y)

#Plot
NH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y)  +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) + 
  scale_edge_color_gradientn(
    colors = c("#D8CCB9", "#FFE6B4","#F3D97C", "#EABE1F"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_classic() +
  ggtitle("Non_hunted - Trophic")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE)) +
  geom_mark_hull(data = hull_df,aes(x = x, y = y_plot, group = module, colour = module), 
                 concavity = 10, expand = unit(3, "mm"), alpha = 0, color = "black") # Modules
NH

##Hunted
#Convert overlap matrix to R matrix format and build the graph
TR_HU <- as.matrix(TR_HU)
graph_obj <- graph_from_adjacency_matrix(TR_HU, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

#Create nodes dataframe from graph vertex names
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength
nodes_df <- nodes_df %>%
  left_join(df_HU %>% select(Genus, Strength_norm), by = c("name" = "Genus")) 

#Combine nodes and edges
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
HU <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) +
  scale_edge_color_gradientn(
    colors = c("#D8CCB9", "#FFE6B4","#F3D97C", "#EABE1F"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_classic() +
  ggtitle("Hunted - Trophic")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))
HU

#Overhunted
#Convert overlap matrix to R matrix format and build the graph
TR_OH <- as.matrix(TR_OH)
graph_obj <- graph_from_adjacency_matrix(TR_OH, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges 
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Create nodes dataframe from graph vertex names
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength
nodes_df <- nodes_df %>%
  left_join(df_OH %>% select(Genus, Strength_norm), by = c("name" = "Genus")) 

#Combine nodes and edges
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
OH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) + 
  scale_edge_color_gradientn(
    colors = c("#D8CCB9", "#FFE6B4","#F3D97C", "#EABE1F"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Overhunted - Trophic")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))
OH

Trophic1 <- ggarrange(NH, HU, OH,
                      nrow = 1, common.legend = T, legend = "bottom")
Trophic1

#### Temporal ####
#Upload Temporal matrices
T_OH <- read.csv("Results/T_Overhunted.csv", row.names = 1)
T_HU <- read.csv("Results/T_Hunted.csv", row.names = 1)
T_NH <- read.csv("Results/T_Non_hunted.csv", row.names = 1)

#Split trophic centrality data by hunting-pressure treatment
df_NH <- T_centr %>%
  filter(Hunting_Pressure == "Non_hunted")
df_HU <- T_centr %>%
  filter(Hunting_Pressure == "Hunted")
df_OH <- T_centr %>%
  filter(Hunting_Pressure == "Overhunted")

#Plot the unipartite networks
#Overhunted
#Convert overlap matrix to R matrix format and build the graph
T_OH <- as.matrix(T_OH)
graph_obj <- graph_from_adjacency_matrix(T_OH, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges 
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength 
nodes_df <- nodes_df %>%
  left_join(df_OH %>% select(Genus, Strength_norm), by = c("name" = "Genus"))

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
OH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) + 
  scale_edge_color_gradientn(
    colors = c("#CCD7E2", "#DCDDF9", "#9FBCF0","#5155E1"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Overhunted - Temporal")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

OH

##Hunted
#Convert overlap matrix to R matrix format and build the graph
T_HU <- as.matrix(T_HU)
graph_obj <- graph_from_adjacency_matrix(T_HU, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength
nodes_df <- nodes_df %>%
  left_join(df_HU %>% select(Genus, Strength_norm), by = c("name" = "Genus")) 

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
HU <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) + 
  scale_edge_color_gradientn(
    colors = c("#CCD7E2", "#DCDDF9", "#9FBCF0","#5155E1"),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    limits = c(0, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Hunted - Temporal")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

HU

##Non-hunted
#Convert overlap matrix to R matrix format and build the graph
T_NH <- as.matrix(T_NH)
graph_obj <- graph_from_adjacency_matrix(T_NH, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength 
nodes_df <- nodes_df %>%
  left_join(df_NH %>% select(Genus, Strength_norm), by = c("name" = "Genus"))

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

#Plot
NH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.8) + 
  scale_edge_color_gradientn(
    colors = c("#CCD7E2", "#DCDDF9", "#9FBCF0","#5155E1"),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    limits = c(0, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0, 0.25, 0.5, 0.75, 1),
    limits = c(0, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Non_hunted - Temporal")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

NH

Temporal1 <- ggarrange(NH, HU, OH,
                       nrow = 1, common.legend = T, legend = "bottom")
Temporal1

#### Spatial ####
#Upload Spatial matrices
S_OH <- read.csv("Results/S_Overhunted.csv", row.names = 1)
S_HU <- read.csv("Results/S_Hunted.csv", row.names = 1)
S_NH <- read.csv("Results/S_Non_hunted.csv", row.names = 1)

df_NH <- S_centr %>%
  filter(Hunting_Pressure == "Non_hunted")
df_HU <- S_centr %>%
  filter(Hunting_Pressure == "Hunted")
df_OH <- S_centr %>%
  filter(Hunting_Pressure == "Overhunted")

#Plot the unipartite networks
#Overhunted
#Convert overlap matrix to R matrix format and build the graph
S_OH <- as.matrix(S_OH)
graph_obj <- graph_from_adjacency_matrix(S_OH, mode = "undirected", weighted = TRUE, diag = FALSE)

#Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength
nodes_df <- nodes_df %>%
  left_join(df_OH %>% select(Genus, Strength_norm), by = c("name" = "Genus")) 

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
OH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.9) + 
  scale_edge_color_gradientn(
    colors = c("#E9D9D5", "#FCCFC5","#FAAF8C", "#F46444"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Overhunted - Spatial")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

OH

##Hunted
#Convert overlap matrix to R matrix format and build the graph
S_HU <- as.matrix(S_HU)
graph_obj <- graph_from_adjacency_matrix(S_HU, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength
nodes_df <- nodes_df %>%
  left_join(df_HU %>% select(Genus, Strength_norm), by = c("name" = "Genus"))

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

# Plot
HU <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.9) + 
  scale_edge_color_gradientn(
    colors = c("#E9D9D5", "#FCCFC5","#FAAF8C", "#F46444"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Hunted - Spatial")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

HU 


##Non-hunted
#Convert overlap matrix to R matrix format and build the graph
S_NH <- as.matrix(S_NH)
graph_obj <- graph_from_adjacency_matrix(S_NH, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges 
edge_df <- igraph::as_data_frame(graph_obj, what = "edges") 

#Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add node strength and Louvain module membership (pre-computed in "02_Metrics" script)
membership_NH <- readRDS("Results/S_membership_NH.rds")
nodes_df <- nodes_df %>%
  left_join(df_NH %>% select(Genus, Strength_norm), by = c("name" = "Genus")) %>%
  left_join(membership_NH, by = c("name" = "Genus"))

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

#Extract node data with flipped y-coordinate, used later to draw module hulls
hull_df <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble() %>%
  mutate(y_plot = -y)

# Plot
NH <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = weight), alpha = 0.9) + 
  scale_edge_color_gradientn(
    colors = c("#E9D9D5", "#FCCFC5","#FAAF8C", "#F46444"),
    values = scales::rescale(c(0.001, 0.25, 0.5, 0.75, 1)),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  scale_edge_width_continuous(
    breaks = c(0.001, 0.25, 0.5, 0.75, 1),
    limits = c(0.001, 1),
    name = "Niche overlap") +
  geom_node_point(aes(size = Strength_norm), color = "grey30") +
  scale_size_continuous(
    range = c(2, 10),
    limits = global_limits,
    breaks = round(seq(global_limits[1], global_limits[2], length.out = 5), 2)) +
  geom_node_text(aes(label = name), repel = TRUE, 
                 size = 6) +
  theme_classic() +
  ggtitle("Non_hunted - Spatial")+ theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))+
  geom_mark_hull(data = hull_df,aes(x = x, y = y_plot, group = module, colour = module), 
                 concavity = 10, expand = unit(3, "mm"), alpha = 0, color = "black") # Modules

NH

Spatial1 <- ggarrange(NH, HU, OH,
                      nrow = 1, common.legend = T, legend = "bottom")
Spatial1

ggarrange(Temporal1, Spatial1, Trophic1,
          nrow = 3)


#### Combined ####
#Upload matrices of combined overlap and relative contribution of each dimension per treatment
combined_objects <- readRDS("Results/combined_objects.rds")
list2env(combined_objects, envir = .GlobalEnv)

#Create the ternary plot structure
res <- 300
Spatial <- seq(0, 1, length.out = res)
Temporal <- seq(0, 1, length.out = res)
grid <- expand.grid(Spatial = Spatial, Temporal = Temporal)
grid$Trophic <- 1 - grid$Spatial - grid$Temporal
grid <- grid[grid$Trophic >= 0, ]

#Select vertex colors
my_colors <- c("#D4310C", "#729CE9", "#EEC944" )
rgb_S <- col2rgb(my_colors[1])/255
rgb_T <- col2rgb(my_colors[2])/255
rgb_Tr <- col2rgb(my_colors[3])/255

#Blend the three vertex colors according to each grid point's ternary coordinates
grid$col <- rgb(
  red   = grid$Spatial*rgb_S[1] + grid$Temporal*rgb_T[1] + grid$Trophic*rgb_Tr[1],
  green = grid$Spatial*rgb_S[2] + grid$Temporal*rgb_T[2] + grid$Trophic*rgb_Tr[2],
  blue  = grid$Spatial*rgb_S[3] + grid$Temporal*rgb_T[3] + grid$Trophic*rgb_Tr[3])

##Non-hunted
#Create a dataframe with all unique primate pairs
I_pairs <- expand.grid(
  Prim1 = rownames(C_S_NH),
  Prim2 = colnames(C_S_NH),
  stringsAsFactors = FALSE)
I_pairs <- I_pairs[I_pairs$Prim1 < I_pairs$Prim2, ]

#Extract Spatial, Temporal, Trophic contribution and total overlap for each pair
I_pairs$Spatial <- mapply(function(x,y) C_S_NH[x,y], I_pairs$Prim1, I_pairs$Prim2)
I_pairs$Temporal <- mapply(function(x,y) C_T_NH[x,y], I_pairs$Prim1, I_pairs$Prim2)
I_pairs$Trophic <- mapply(function(x,y) C_TR_NH[x,y], I_pairs$Prim1, I_pairs$Prim2)
I_pairs$Total_overlap <- mapply(function(x,y) C_NH[x,y], I_pairs$Prim1, I_pairs$Prim2)

#Assign a blended color to each pair based on its Spatial/Temporal/Trophic contribution
I_pairs$col <- mapply(function(S,T,Tr){
  rgb(
    red   = S*rgb_S[1] + T*rgb_T[1] + Tr*rgb_Tr[1],
    green = S*rgb_S[2] + T*rgb_T[2] + Tr*rgb_Tr[2],
    blue  = S*rgb_S[3] + T*rgb_T[3] + Tr*rgb_Tr[3]
  )
}, I_pairs$Spatial, I_pairs$Temporal, I_pairs$Trophic)

# Create the graph
C_NH <- as.matrix(C_NH)
graph_obj <- graph_from_adjacency_matrix(C_NH, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

# Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Add a column of strength and modularity membership in the nodes
C_membership_NH <- readRDS("Results/C_membership_NH.rds")
nodes_df <- nodes_df %>%
  left_join(C_membership_NH, by = c("name" = "Genus"))

#Combine edges and nodes
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)

#Attach fixed coordinates to nodes
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

#Extract nodes and edges as a dataframe
nodes <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble()

edges_named <- graph_tbl %>%
  activate(edges) %>%
  as_tibble()

#Match edge colors to the corresponding pair in I_pairs, so networks and ternary plot share the same colors
#Add primate names to each edge based on node index
edges_named$from_name <- nodes$name[edges_named$from]
edges_named$to_name   <- nodes$name[edges_named$to]

#Standardize primate pair order (alphabetical) to match with I_pairs
edges_named <- edges_named %>%
  mutate(
    Prim1 = pmin(from_name, to_name),
    Prim2 = pmax(from_name, to_name))

#Standardize primate pair order (alphabetical) in I_pairs as well
I_pairs2 <- I_pairs %>%
  mutate(
    Prim1 = pmin(Prim1, Prim2),
    Prim2 = pmax(Prim1, Prim2))

#Match each edge with its corresponding blended color
edges_colored <- edges_named %>%
  left_join(I_pairs2, by = c("Prim1", "Prim2"))

#Assign the blended color to the graph edges
graph_tbl <- graph_tbl %>%
  activate(edges) %>%
  mutate(col = edges_colored$col)

#Extract node data with flipped y-coordinate, used later to draw module hulls
hull_df <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble() %>%
  mutate(y_plot = -y)

# Network Plot
C_NH_plot <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = col), alpha = 0.8) +
  scale_edge_color_identity() +
  scale_edge_width_continuous(breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    limits = c(0.001, 3),
    name = "Niche overlap") +
  geom_node_point(size = 5, color = "grey30") +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_classic() +
  ggtitle("Non_hunted - Combined") +
  theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE)) +
  geom_mark_hull(
    data = hull_df,
    aes(x = x, y = y_plot, group = module, colour = module),
    concavity = 10,
    expand = unit(3, "mm"),
    alpha = 0,
    color = "black")

C_NH_plot

#Ternary plot
I_tern_points <- ggtern(I_pairs, aes(x = Spatial, y = Temporal, z = Trophic)) +
  geom_point(aes(fill = col, size = Total_overlap), 
             shape = 21,    # cerchio pieno con bordo
             color = "black", # bordo nero
             alpha = 0.9) +
  scale_fill_identity()+  scale_color_identity() +
  scale_size_continuous(
    breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    range = c(1, 9),       # dimensione minima e massima dei punti
    limits = c(0.1, 3)) + # usa i colori così come sono
  scale_T_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_L_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_R_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  theme_bw() +
  theme_showgrid() +
  theme(
    axis.title=element_blank(),
    tern.axis.title.T = element_text(color=my_colors[2], size=12),
    tern.axis.title.L = element_text(color=my_colors[1], size=12),
    tern.axis.title.R = element_text(color=my_colors[3], size=12),
    tern.axis.text.T = element_text(color=my_colors[2], size=9),
    tern.axis.text.L = element_text(color=my_colors[1], size=9),
    tern.axis.text.R = element_text(color=my_colors[3], size=9),
    panel.ontop = FALSE,
    panel.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom")

I_tern_points <- I_tern_points +
  ggtitle("Non_hunted")

I_tern_points

##Hunted
#Create a dataframe with all unique primate pairs
H_pairs <- expand.grid(
  Prim1 = rownames(C_S_HU),
  Prim2 = colnames(C_S_HU),
  stringsAsFactors = FALSE
)
H_pairs <- H_pairs[H_pairs$Prim1 < H_pairs$Prim2, ]

#Extract Spatial, Temporal, Trophic contribution and total overlap for each pair
H_pairs$Spatial <- mapply(function(x,y) C_S_HU[x,y], H_pairs$Prim1, H_pairs$Prim2)
H_pairs$Temporal <- mapply(function(x,y) C_T_HU[x,y], H_pairs$Prim1, H_pairs$Prim2)
H_pairs$Trophic <- mapply(function(x,y) C_TR_HU[x,y], H_pairs$Prim1, H_pairs$Prim2)
H_pairs$Total_overlap <- mapply(function(x,y) C_HU[x,y], H_pairs$Prim1, H_pairs$Prim2)

#Assign a blended color to each pair based on its Spatial/Temporal/Trophic contribution
H_pairs$col <- mapply(function(S,T,Tr){
  rgb(
    red   = S*rgb_S[1] + T*rgb_T[1] + Tr*rgb_Tr[1],
    green = S*rgb_S[2] + T*rgb_T[2] + Tr*rgb_Tr[2],
    blue  = S*rgb_S[3] + T*rgb_T[3] + Tr*rgb_Tr[3]
  )
}, H_pairs$Spatial, H_pairs$Temporal, H_pairs$Trophic)

##Hunted network
# Create the graph
C_HU <- as.matrix(C_HU)
graph_obj <- graph_from_adjacency_matrix(C_HU, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges (links) with color_val
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

# Create df with nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Combine nodes and edges
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

nodes <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble()

edges_named <- graph_tbl %>%
  activate(edges) %>%
  as_tibble()

#Match edge colors to the corresponding pair in I_pairs, so networks and ternary plot share the same colors
edges_named$from_name <- nodes$name[edges_named$from]
edges_named$to_name   <- nodes$name[edges_named$to]

edges_named <- edges_named %>%
  mutate(
    Prim1 = pmin(from_name, to_name),
    Prim2 = pmax(from_name, to_name))

H_pairs2 <- H_pairs %>%
  mutate(
    Prim1 = pmin(Prim1, Prim2),
    Prim2 = pmax(Prim1, Prim2))


edges_colored <- edges_named %>%
  left_join(H_pairs2, by = c("Prim1", "Prim2"))


graph_tbl <- graph_tbl %>%
  activate(edges) %>%
  mutate(col = edges_colored$col)

# Network Plot
C_HU_plot <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = col), alpha = 0.8) +
  scale_edge_color_identity() +
  scale_edge_width_continuous(breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    limits = c(0.001, 3),
    name = "Niche overlap") +
  geom_node_point(size = 5, color = "grey30") +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_classic() +
  ggtitle("Hunted - Combined") +
  theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

C_HU_plot


H_tern_points <- ggtern(H_pairs, aes(x = Spatial, y = Temporal, z = Trophic)) +
  geom_point(aes(fill = col, size = Total_overlap), 
             shape = 21,    # cerchio pieno con bordo
             color = "black", # bordo nero
             alpha = 0.8) +
  scale_fill_identity()+
  scale_color_identity() +
  scale_size_continuous(
    breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    range = c(1, 9),       # dimensione minima e massima dei punti
    limits = c(0.1, 3)) +
  scale_T_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_L_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_R_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  theme_bw() +
  theme_showgrid() +
  theme(
    tern.panel.grid.minor = element_line( size=0.7),
    tern.panel.grid.major = element_line(size=0.7),
    axis.title=element_blank(),
    tern.axis.title.T = element_text(color=my_colors[2], size=12),
    tern.axis.title.L = element_text(color=my_colors[1], size=12),
    tern.axis.title.R = element_text(color=my_colors[3], size=12),
    tern.axis.text.T = element_text(color=my_colors[2], size=9),
    tern.axis.text.L = element_text(color=my_colors[1], size=9),
    tern.axis.text.R = element_text(color=my_colors[3], size=9),
    panel.ontop = FALSE,
    panel.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom")

H_tern_points <- H_tern_points +
  ggtitle("Hunted")
H_tern_points

##Overhunted
#Create a dataframe with all unique primate pairs
O_pairs <- expand.grid(
  Prim1 = rownames(C_S_OH),
  Prim2 = colnames(C_S_OH),
  stringsAsFactors = FALSE)
O_pairs <- O_pairs[O_pairs$Prim1 < O_pairs$Prim2, ]

#Extract Spatial, Temporal, Trophic contribution and total overlap for each pair
O_pairs$Spatial <- mapply(function(x,y) C_S_OH[x,y], O_pairs$Prim1, O_pairs$Prim2)
O_pairs$Temporal <- mapply(function(x,y) C_T_OH[x,y], O_pairs$Prim1, O_pairs$Prim2)
O_pairs$Trophic <- mapply(function(x,y) C_TR_OH[x,y], O_pairs$Prim1, O_pairs$Prim2)
O_pairs$Total_overlap <- mapply(function(x,y) C_OH[x,y], O_pairs$Prim1, O_pairs$Prim2)

#Assign a blended color to each pair based on its Spatial/Temporal/Trophic contribution
O_pairs$col <- mapply(function(S,T,Tr){
  rgb(
    red   = S*rgb_S[1] + T*rgb_T[1] + Tr*rgb_Tr[1],
    green = S*rgb_S[2] + T*rgb_T[2] + Tr*rgb_Tr[2],
    blue  = S*rgb_S[3] + T*rgb_T[3] + Tr*rgb_Tr[3]
  )
}, O_pairs$Spatial, O_pairs$Temporal, O_pairs$Trophic)


##Overunted network
# Create the graph
C_OH <- as.matrix(C_OH)
graph_obj <- graph_from_adjacency_matrix(C_OH, mode = "undirected", weighted = TRUE, diag = FALSE)

# Extract the edges
edge_df <- igraph::as_data_frame(graph_obj, what = "edges")

# Nodes
nodes_df <- data.frame(name = V(graph_obj)$name, stringsAsFactors = FALSE)

#Combine nodes and edges
graph_tbl <- tbl_graph(nodes = nodes_df, edges = edge_df, directed = FALSE)
graph_tbl <- graph_tbl %>%
  activate(nodes) %>%
  left_join(fixed_coords, by = "name")

nodes <- graph_tbl %>%
  activate(nodes) %>%
  as_tibble()

edges_named <- graph_tbl %>%
  activate(edges) %>%
  as_tibble()

#Match edge colors to the corresponding pair in I_pairs, so networks and ternary plot share the same colors
edges_named$from_name <- nodes$name[edges_named$from]
edges_named$to_name   <- nodes$name[edges_named$to]

edges_named <- edges_named %>%
  mutate(
    Prim1 = pmin(from_name, to_name),
    Prim2 = pmax(from_name, to_name))

O_pairs2 <- O_pairs %>%
  mutate(
    Prim1 = pmin(Prim1, Prim2),
    Prim2 = pmax(Prim1, Prim2))

edges_colored <- edges_named %>%
  left_join(O_pairs2, by = c("Prim1", "Prim2"))

graph_tbl <- graph_tbl %>%
  activate(edges) %>%
  mutate(col = edges_colored$col)

# Network Plot
C_OH_plot <- ggraph(graph_tbl, layout = "manual", x = x, y = -y) +
  geom_edge_link(aes(width = weight, color = col), alpha = 0.8) +
  scale_edge_color_identity() +
  scale_edge_width_continuous(breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    limits = c(0.001, 3),
    name = "Niche overlap") +
  geom_node_point(size = 5, color = "grey30") +
  geom_node_text(aes(label = name), repel = TRUE, size = 6) +
  theme_classic() +
  ggtitle("Overhunted - Combined") +
  theme_void() + 
  coord_fixed(
    xlim = range(graph_tbl %>% activate(nodes) %>% pull(x), na.rm = TRUE),
    ylim = range(-(graph_tbl %>% activate(nodes) %>% pull(y)), na.rm = TRUE))

C_OH_plot


#Ternary plots
O_tern_points <- ggtern(O_pairs, aes(x = Spatial, y = Temporal, z = Trophic)) +
  geom_point(aes(fill = col, size = Total_overlap), 
             shape = 21,    # cerchio pieno con bordo
             color = "black", # bordo nero
             alpha = 0.8) +
  scale_fill_identity()+
  scale_size_continuous(
    breaks = c(0.1, 0.5, 1, 1.5, 2, 3),
    range = c(1, 9),       # dimensione minima e massima dei punti
    limits = c(0.1, 3)) +
  scale_T_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_L_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  scale_R_continuous(breaks = seq(0,1,by=0.1), labels = seq(0,1,by=0.1)) +
  theme_bw() +
  theme_showgrid() +
  theme(
    tern.panel.grid.minor = element_line( size=0.7),
    tern.panel.grid.major = element_line( size=0.7),
    axis.title=element_blank(),
    tern.axis.title.T = element_text(color=my_colors[2], size=12),
    tern.axis.title.L = element_text(color=my_colors[1], size=12),
    tern.axis.title.R = element_text(color=my_colors[3], size=12),
    tern.axis.text.T = element_text(color=my_colors[2], size=9),
    tern.axis.text.L = element_text(color=my_colors[1], size=9),
    tern.axis.text.R = element_text(color=my_colors[3], size=9),
    panel.ontop = FALSE,
    panel.background = element_rect(fill = NA, color = NA),
    legend.position = "bottom")

O_tern_points <- O_tern_points +
  ggtitle("Overhunted")
O_tern_points

tern <- grid.arrange(I_tern_points, H_tern_points, O_tern_points, ncol = 3)
tern

net <- ggarrange(C_NH_plot, C_HU_plot, C_OH_plot, ncol = 3, common.legend = T)
net

comb <- ggarrange(net, tern, nrow = 2)
comb
ggarrange(Spatial1, Trophic1,
          nrow = 2 , legend = "none")

all <- ggarrange(Temporal1, Spatial1, Trophic1,
nrow = 3, common.legend = T)
all

ggsave("Figures/Spatial_networks.pdf", Spatial1, width = 14, height = 6)
ggsave("Figures/Trophic_networks.pdf", Trophic1, width = 14, height = 6)
ggsave("Figures/Temporal_networks.pdf", Temporal1, width = 14, height = 6)
ggsave("Figures/Ternary_plots.pdf", tern, width = 14, height = 6)
ggsave("Figures/Combined_networks.pdf", net, width = 14, height = 6)
ggsave("Figures/All_networks.pdf", all, width = 14, height = 14)
ggsave("Figures/comb_networks.pdf", comb, width = 14, height = 11)

