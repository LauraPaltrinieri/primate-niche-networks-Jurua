#Network metrics
#Author: Laura Paltrinieri
#17/04/2026
#Modified: 07/09/2026

#This script calculates network-level and node-level metrics (strength centrality, modularity,
#and nestedness) for each niche dimension (Temporal, Spatial, Trophic) and the combined overlap
#network, across hunting-pressure treatments (Non-hunted, Hunted, Overhunted). All metrics are
#tested against null models to assess statistical significance.
#Centrality data, module membership (when significative), and summary plots are exported for use in the network
#and ternary plot script.

# Load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(openxlsx)
library(officer)
library(spaa)
library(tibble)
library(igraph)
library(bipartite)
library(ggpubr)
library(tidytext)
library(ggh4x)
#clean_env
rm(list = ls())

##### Temporal #####
##Upload the matrices
T_OH <- read.csv("Results/T_Overhunted.csv", row.names = 1)
T_HU <- read.csv("Results/T_Hunted.csv", row.names = 1)
T_NH <- read.csv("Results/T_Non_hunted.csv", row.names = 1)

#Select the same order of species in each matrix based on body size
species_order <- c("Ateles", "Lagothrix", "Alouatta",  "Cacajao",
                   "Sapajus", "Cebus", "Pithecia", "Callicebus",
                   "Aotus", "Saimiri", "Saguinus")

reorder_matrix <- function(mat, order) {
  species_present <- order[order %in% rownames(mat)]
  mat[species_present, species_present]}

T_OH <- reorder_matrix(T_OH, species_order)
T_HU <- reorder_matrix(T_HU, species_order)
T_NH <- reorder_matrix(T_NH, species_order)

## Strength centrality
set.seed(123)
n_sim <- 1000

calc_strength_null <- function(mat){
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  diag(mat) <- 0
  
  g <- graph_from_adjacency_matrix(
    mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  N <- vcount(g)
  
  
  # strength osservata
  s_obs <- igraph::strength(g, mode = "all", weights = E(g)$weight)
  s_norm <- s_obs / (N - 1)
  
  # matrice null
  s_null <- matrix(NA, nrow = vcount(g), ncol = n_sim)
  rownames(s_null) <- names(s_obs)
  
  for (i in 1:n_sim) {
    
    g_rand <- g
    
    # stesso null model della modularità: shuffle pesi
    E(g_rand)$weight <- sample(E(g)$weight)
    
    s_null[, i] <- igraph::strength(g_rand, mode = "all")[names(s_obs)]
  }
  
  # statistiche null
  s_mean <- rowMeans(s_null)
  s_sd   <- apply(s_null, 1, sd)
  
  # z-score
  s_z <- (s_obs - s_mean) / s_sd
  
  data.frame(
    Genus = names(s_obs),
    Strength_obs = s_obs,
    Strength_norm = s_norm,
    Strength_mean_null = s_mean,
    Strength_sd_null = s_sd,
    Strength_z = s_z
  )
}

str_OH <- calc_strength_null(T_OH)
str_HU <- calc_strength_null(T_HU)
str_NH <- calc_strength_null(T_NH)

str_OH$Hunting_Pressure <- "Overhunted"
str_HU$Hunting_Pressure <- "Hunted"
str_NH$Hunting_Pressure <- "Non_hunted"

T_strength_null <- rbind(str_NH, str_HU, str_OH)
T_strength_null$Dimension <- "Temporal"

T_strength_null$Hunting_Pressure <- factor(
  T_strength_null$Hunting_Pressure,
  levels = c("Non_hunted", "Hunted", "Overhunted"),
  labels = c("Non_hunted", "Hunted", "Overhunted")
)

T_strength_null$Significance <- ifelse(T_strength_null$Strength_z > 1.96, "Positive",
                                       ifelse(T_strength_null$Strength_z < -1.96, "Negative", "n.s."))

# tabella riassuntiva solo dei significativi
T_strength_null[T_strength_null$Significance != "n.s.",
                c("Genus", "Hunting_Pressure", "Strength_z", "Significance")]

T_centr_plot <- ggplot(T_strength_null, aes(
  x = Strength_norm,
  y = reorder_within(Genus, Strength_norm, Hunting_Pressure)
)) +
  
  geom_segment(aes(
    x = 0, xend = Strength_norm,
    y = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    yend = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    color = Strength_z
  )) +
  
  # punti (testa del lollipop)
  geom_point(aes(color = Strength_z), size = 2) +
  
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "gold", "#F4A582", "#B2182B"),
    values = scales::rescale(c(-6, -2, 0, 1, 3)),
    limits = c(-6, 3)
  )+
  
  # linea zero (non necessaria ma utile)
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Hunting_Pressure, scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  scale_y_reordered() +
  
  labs(
    x = "Strength (normalized)",
    y = NULL,
    color = "z-score\n(null model)"
  ) +
  
  theme_bw()

T_centr_plot

#Modularity
###Overhunted
T_OH <- as.matrix(T_OH)
set.seed(123)
n_sim <- 1000

g_over <- igraph::graph_from_adjacency_matrix(T_OH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_OH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_OH <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_OH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}


Z_score_OH <- (Q_real_OH - mean(Q_null_OH)) / sd(Q_null_OH)

Q_real_OH #0
Z_score_OH #-1.79

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0
mod_T_OH <- c(obs = Q_real_OH, Z = Z_score_OH, ratio = Q_ratio_OH)

##Nestedness
mat <- as.matrix(T_OH)
diag(mat) <- 0

nested_real_OH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_OH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_OH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_OH <- (nested_real_OH - mean(nested_null_OH)) / sd(nested_null_OH)

nested_real_OH #65.87
Z_score_OH #1.52

N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH #0.89
nest_T_OH <- c(obs = nested_real_OH, Z = Z_score_OH, ratio = N_ratio_OH)

##Nested contribution
contrib_OH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_OH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_OH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_OH$contribution[k] <- (nested_real_OH - mean(null_k)) / sd(null_k)
}

contrib_OH[order(-contrib_OH$contribution), ]

null_df_OH <- data.frame(node = rep(names(null_list_OH), sapply(null_list_OH, length)),
                         WNODA_null = unlist(null_list_OH))
lev_OH <- contrib_OH$node[order(contrib_OH$contribution)]
null_df_OH$node <- factor(null_df_OH$node, levels = lev_OH)
contrib_OH$node_f <- factor(contrib_OH$node, levels = lev_OH)
x_lab_OH <- max(null_df_OH$WNODA_null, nested_real_OH) + 1.5

plot_contrib_T_OH <- ggplot(subset(contrib_OH, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#6F74E7") +
  geom_point(size = 3, colour = "#6F74E7") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Temporal-OH") +
  theme_bw()

plot_contrib_T_OH

###Hunted
T_HU <- as.matrix(T_HU)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(T_HU, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_HU <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_HU <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_HU[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_HU <- (Q_real_HU - mean(Q_null_HU)) / sd(Q_null_HU)

Q_real_HU #0
Z_score_HU #-1.65

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0.005
mod_T_HU <- c(obs = Q_real_HU, Z = Z_score_HU, ratio = Q_ratio_HU)

##Nestedness
mat <- as.matrix(T_HU)
diag(mat) <- 0

nested_real_HU <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_HU <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_HU[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_HU <- (nested_real_HU - mean(nested_null_HU)) / sd(nested_null_HU)

nested_real_HU #73.81 
Z_score_HU #3.05

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU #1.02
nest_T_HU <- c(obs = nested_real_HU, Z = Z_score_HU, ratio = N_ratio_HU)

##Nested contribution
contrib_HU <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_HU <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_HU[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_HU$contribution[k] <- (nested_real_HU - mean(null_k)) / sd(null_k)
}

contrib_HU[order(-contrib_HU$contribution), ]

null_df_HU <- data.frame(node = rep(names(null_list_HU), sapply(null_list_HU, length)),
                         WNODA_null = unlist(null_list_HU))
lev_HU <- contrib_HU$node[order(contrib_HU$contribution)]
null_df_HU$node <- factor(null_df_HU$node, levels = lev_HU)
contrib_HU$node_f <- factor(contrib_HU$node, levels = lev_HU)
x_lab_HU <- max(null_df_HU$WNODA_null, nested_real_HU) + 1.5

plot_contrib_T_HU <- ggplot(subset(contrib_HU, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#6F74E7") +
  geom_point(size = 3, colour = "#6F74E7") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Temporal-HU") +
  theme_bw()

plot_contrib_T_HU

###Non-hunted
T_NH <- as.matrix(T_NH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(T_NH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_NH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_NH <- numeric(n_sim)

for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_NH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_NH <- (Q_real_NH - mean(Q_null_NH)) / sd(Q_null_NH)

Q_real_NH #0
Z_score_NH #-3.02

Q_ratio_NH <- Q_real_NH / max(Q_null_NH)
Q_ratio_NH #0
mod_T_NH <- c(obs = Q_real_NH, Z = Z_score_NH, ratio = Q_ratio_NH)


##Nestedness
mat <- as.matrix(T_NH)
diag(mat) <- 0

nested_real_NH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_NH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_NH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_NH <- (nested_real_NH - mean(nested_null_NH)) / sd(nested_null_NH)

nested_real_NH #71.82
Z_score_NH #3.61

N_ratio_NH <- nested_real_NH / max(nested_null_NH)
N_ratio_NH #1.01
nest_T_NH <- c(obs = nested_real_NH, Z = Z_score_NH, ratio = N_ratio_NH)

##Nested contribution
contrib_NH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_NH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_NH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_NH$contribution[k] <- (nested_real_NH - mean(null_k)) / sd(null_k)
}

contrib_NH[order(-contrib_NH$contribution), ]

null_df_NH <- data.frame(node = rep(names(null_list_NH), sapply(null_list_NH, length)),
                         WNODA_null = unlist(null_list_NH))
lev_NH <- contrib_NH$node[order(contrib_NH$contribution)]
null_df_NH$node <- factor(null_df_NH$node, levels = lev_NH)
contrib_NH$node_f <- factor(contrib_NH$node, levels = lev_NH)
x_lab_NH <- max(null_df_NH$WNODA_null, nested_real_NH) + 1.5

plot_contrib_T_NH <- ggplot(subset(contrib_NH, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#6F74E7") +
  geom_point(size = 3, colour = "#6F74E7") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Temporal-NH") +
  theme_bw()

plot_contrib_T_NH

## Modularity plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(Q_null_NH),
                mean(Q_null_HU),
                mean(Q_null_OH)),
  sd_null = c(sd(Q_null_NH),
              sd(Q_null_HU),
              sd(Q_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  Q_real = c(Q_real_NH, Q_real_HU, Q_real_OH),
  colore = c("#6F74E7","#6F74E7",  "#6F74E7")  # Non_hunted, Hunted, Overhunted
)

# Plot
mod <- ggplot() +
  # Punto null model con barre di errore
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  # Punto modularità reale
  geom_point(data = df_real, aes(x = rete, y = Q_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Modularity")+
  theme_bw(base_size = 14)

#Nestedness plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(nested_null_NH),
                mean(nested_null_HU),
                mean(nested_null_OH)),
  sd_null = c(sd(nested_null_NH),
              sd(nested_null_HU),
              sd(nested_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  nested_real = c(nested_real_NH, nested_real_HU, nested_real_OH),
  colore = c("#6F74E7","#6F74E7", "#6F74E7")  # Non_hunted, Hunted, Overhunted
)

nest <- ggplot() +
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  geom_point(data = df_real, aes(x = rete, y = nested_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Nestedness (WNODA)") +
  theme_bw(base_size = 14)
nest
mod
T_nest_contrib <- rbind(
  data.frame(Genus = contrib_NH$node, Nest_contrib = contrib_NH$contribution, Hunting_Pressure = "Non_hunted"),
  data.frame(Genus = contrib_HU$node, Nest_contrib = contrib_HU$contribution, Hunting_Pressure = "Hunted")
)
T_nest_contrib$Dimension <- "Temporal"

write.csv(T_nest_contrib, "Results/T_nested_contribution.csv", row.names = FALSE)
#Combined plot
T_met <- ggarrange(mod, nest, nrow = 1)
T_met

Nest_Con <- ggarrange(plot_contrib_T_NH, plot_contrib_T_HU, plot_contrib_T_OH)
Nest_Con

##### Spatial #####
##Calculate Modularity
S_OH <- read.csv("Results/S_Overhunted.csv", row.names = 1)
S_HU <- read.csv("Results/S_Hunted.csv", row.names = 1)
S_NH <- read.csv("Results/S_Non_hunted.csv", row.names = 1)


#Select the order of the species
species_order <- c("Alouatta", "Aotus", "Ateles", "Cacajao", "Callicebus", 
                   "Cebus", "Lagothrix", "Pithecia", "Saguinus", "Saimiri",
                   "Sapajus")


#Create a function to reorder matrices
reorder_matrix <- function(mat, order) {
  species_present <- order[order %in% rownames(mat)]
  mat[species_present, species_present]
}

#Reorder primates in matrices
S_OH <- reorder_matrix(S_OH, species_order)
S_HU <- reorder_matrix(S_HU, species_order)
S_NH <- reorder_matrix(S_NH, species_order)

set.seed(123)
n_sim <- 1000

calc_strength_null <- function(mat){
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  diag(mat) <- 0
  
  g <- graph_from_adjacency_matrix(
    mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  N <- vcount(g)
  
  
  # strength osservata
  s_obs <- igraph::strength(g, mode = "all", weights = E(g)$weight)
  s_norm <- s_obs / (N - 1)
  
  # matrice null
  s_null <- matrix(NA, nrow = vcount(g), ncol = n_sim)
  rownames(s_null) <- names(s_obs)
  
  for (i in 1:n_sim) {
    
    g_rand <- g
    
    # stesso null model della modularità: shuffle pesi
    E(g_rand)$weight <- sample(E(g)$weight)
    
    s_null[, i] <- igraph::strength(g_rand, mode = "all")[names(s_obs)]
  }
  
  # statistiche null
  s_mean <- rowMeans(s_null)
  s_sd   <- apply(s_null, 1, sd)
  
  # z-score
  s_z <- (s_obs - s_mean) / s_sd
  
  data.frame(
    Genus = names(s_obs),
    Strength_obs = s_obs,
    Strength_norm = s_norm,
    Strength_mean_null = s_mean,
    Strength_sd_null = s_sd,
    Strength_z = s_z
  )
}

str_NH <- calc_strength_null(S_NH)
str_HU <- calc_strength_null(S_HU)
str_OH <- calc_strength_null(S_OH)

str_NH$Hunting_Pressure <- "Non_hunted"
str_HU$Hunting_Pressure <- "Hunted"
str_OH$Hunting_Pressure <- "Overhunted"

S_strength_null <- rbind(str_NH, str_HU, str_OH)

S_strength_null$Dimension <- "Spatial"

S_strength_null$Hunting_Pressure <- factor(
  S_strength_null$Hunting_Pressure,
  levels = c("Non_hunted", "Hunted", "Overhunted"),
  labels = c("Non_hunted", "Hunted", "Overhunted")
)
S_strength_null$Significance <- ifelse(S_strength_null$Strength_z > 1.96, "Positive",
                                       ifelse(S_strength_null$Strength_z < -1.96, "Negative", "n.s."))

# tabella riassuntiva solo dei significativi
S_strength_null[S_strength_null$Significance != "n.s.",
                c("Genus", "Hunting_Pressure", "Strength_z", "Significance")]

S_centr_plot <- ggplot(S_strength_null, aes(
  x = Strength_norm,
  y = reorder_within(Genus, Strength_norm, Hunting_Pressure)
)) +
  
  geom_segment(aes(
    x = 0, xend = Strength_norm,
    y = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    yend = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    color = Strength_z
  )) +
  
  # punti (testa del lollipop)
  geom_point(aes(color = Strength_z), size = 2) +
  
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "gold", "#F4A582", "#B2182B"),
    values = scales::rescale(c(-6, -2, 0, 1, 3)),
    limits = c(-6, 3)
  )+
  
  # linea zero (non necessaria ma utile)
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Hunting_Pressure, scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  scale_y_reordered() +
  
  labs(
    x = "Strength (normalized)",
    y = NULL,
    fill = "z-score\n(null model)"
  ) +
  
  theme_bw()

S_centr_plot

##Now I will calculate the Modularity in each treatment area and I will plot it

###Overhunted
S_OH <- as.matrix(S_OH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(S_OH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_OH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_OH <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_OH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_OH <- (Q_real_OH - mean(Q_null_OH)) / sd(Q_null_OH)

Q_real_OH #0.01
Z_score_OH #-1.61

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0.06
mod_S_OH <- c(obs = Q_real_OH, Z = Z_score_OH, ratio = Q_ratio_OH)


##Nestedness
mat <- as.matrix(S_OH)
diag(mat) <- 0

nested_real_OH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_OH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_OH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_OH <- (nested_real_OH - mean(nested_null_OH)) / sd(nested_null_OH)

nested_real_OH #65.08
Z_score_OH #1.43


N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH
nest_S_OH <- c(obs = nested_real_OH, Z = Z_score_OH, ratio = N_ratio_OH)

##Nested contribution
contrib_OH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_OH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_OH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_OH$contribution[k] <- (nested_real_OH - mean(null_k)) / sd(null_k)
}

contrib_OH[order(-contrib_OH$contribution), ]

null_df_OH <- data.frame(node = rep(names(null_list_OH), sapply(null_list_OH, length)),
                         WNODA_null = unlist(null_list_OH))
lev_OH <- contrib_OH$node[order(contrib_OH$contribution)]
null_df_OH$node <- factor(null_df_OH$node, levels = lev_OH)
contrib_OH$node_f <- factor(contrib_OH$node, levels = lev_OH)
x_lab_OH <- max(null_df_OH$WNODA_null, nested_real_OH) + 1.5

plot_contrib_S_OH <- ggplot(subset(contrib_OH, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#F46444") +
  geom_point(size = 3, colour = "#F46444") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Spatial-OH") +
  theme_bw()

plot_contrib_S_OH

###Hunted
S_HU <- as.matrix(S_HU)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(S_HU, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_HU <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_HU <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_HU[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_HU <- (Q_real_HU - mean(Q_null_HU)) / sd(Q_null_HU)

Q_real_HU #0.04
Z_score_HU #-2.08

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0.12
mod_S_HU <- c(obs = Q_real_HU, Z = Z_score_HU, ratio = Q_ratio_HU)

##Nestedness
mat <- as.matrix(S_HU)
diag(mat) <- 0

nested_real_HU <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_HU <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_HU[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_HU <- (nested_real_HU - mean(nested_null_HU)) / sd(nested_null_HU)

nested_real_HU #61.9
Z_score_HU #1.14

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU #0.89
nest_S_HU <- c(obs = nested_real_HU, Z = Z_score_HU, ratio = N_ratio_HU)

###Non_hunted
S_NH <- as.matrix(S_NH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(S_NH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_NH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_NH <- numeric(n_sim)

for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_NH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_NH <- (Q_real_NH - mean(Q_null_NH)) / sd(Q_null_NH)

Q_real_NH #0.17
Z_score_NH #2.58
membership(comm_real_over)

g_S_NH    <- g_over
comm_S_NH <- comm_real_over

saveRDS(list(g = g_S_NH,  comm = comm_S_NH),  "Results/S_modularity_NH.rds")

#Save the membership for the figure
membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(membership_NH, "Results/S_membership_NH.rds")

Q_ratio_NH <- Q_real_NH / max(Q_null_NH)
Q_ratio_NH #0.93
mod_S_NH <- c(obs = Q_real_NH, Z = Z_score_NH, ratio = Q_ratio_NH)

##Nestedness
mat <- as.matrix(S_NH)
diag(mat) <- 0

nested_real_NH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_NH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_NH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_NH <- (nested_real_NH - mean(nested_null_NH)) / sd(nested_null_NH)

nested_real_NH #55.7
Z_score_NH #-0.42

N_ratio_NH <- nested_real_NH / max(nested_null_NH)
N_ratio_NH #0.84
nest_S_NH <- c(obs = nested_real_NH, Z = Z_score_NH, ratio = N_ratio_NH)

## Modularity plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(Q_null_NH),
                mean(Q_null_HU),
                mean(Q_null_OH)),
  sd_null = c(sd(Q_null_NH),
              sd(Q_null_HU),
              sd(Q_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  Q_real = c(Q_real_NH, Q_real_HU, Q_real_OH),
  colore = c("#F46444","#F46444",  "#F46444")  # Non_hunted, Hunted, Overhunted
)

# Plot
mod <- ggplot() +
  # Punto null model con barre di errore
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  # Punto modularità reale
  geom_point(data = df_real, aes(x = rete, y = Q_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Modularity")+
  theme_bw(base_size = 14)


#Nestedness plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(nested_null_NH),
                mean(nested_null_HU),
                mean(nested_null_OH)),
  sd_null = c(sd(nested_null_NH),
              sd(nested_null_HU),
              sd(nested_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  nested_real = c(nested_real_NH, nested_real_HU, nested_real_OH),
  colore = c("#F46444","#F46444", "#F46444")  # Non_hunted, Hunted, Overhunted
)

nest <- ggplot() +
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  geom_point(data = df_real, aes(x = rete, y = nested_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Nestedness (WNODA)") +
  theme_bw(base_size = 14)
nest
#Combined plot
S_met <- ggarrange(mod, nest, nrow = 1)
S_met


##### Trophic #####
#Upload Trophic matrices
TR_OH <- read.csv("Results/TR_Overhunted.csv", row.names = 1)
TR_HU <- read.csv("Results/TR_Hunted.csv", row.names = 1)
TR_NH <- read.csv("Results/TR_Non_hunted.csv", row.names = 1)

# Select the order of the species
species_order <- c("Alouatta", "Aotus", "Ateles", "Cacajao", "Callicebus", 
                   "Cebus", "Lagothrix", "Pithecia", "Saguinus", "Saimiri",
                   "Sapajus")

#Create a function to reorder a matrices
reorder_matrix <- function(mat, order) {
  species_present <- order[order %in% rownames(mat)]
  mat[species_present, species_present]
}

#Reorder primates in matrices
TR_OH <- reorder_matrix(TR_OH, species_order)
TR_HU <- reorder_matrix(TR_HU, species_order)
TR_NH <- reorder_matrix(TR_NH, species_order)

#Strength centrality
set.seed(123)
n_sim <- 1000

calc_strength_null <- function(mat){
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  diag(mat) <- 0
  
  g <- graph_from_adjacency_matrix(
    mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  N <- vcount(g)
  
  
  # strength osservata
  s_obs <- igraph::strength(g, mode = "all", weights = E(g)$weight)
  s_norm <- s_obs / (N - 1)
  
  # matrice null
  s_null <- matrix(NA, nrow = vcount(g), ncol = n_sim)
  rownames(s_null) <- names(s_obs)
  
  for (i in 1:n_sim) {
    
    g_rand <- g
    
    # stesso null model della modularità: shuffle pesi
    E(g_rand)$weight <- sample(E(g)$weight)
    
    s_null[, i] <- igraph::strength(g_rand, mode = "all")[names(s_obs)]
  }
  
  # statistiche null
  s_mean <- rowMeans(s_null)
  s_sd   <- apply(s_null, 1, sd)
  
  # z-score
  s_z <- (s_obs - s_mean) / s_sd
  
  data.frame(
    Genus = names(s_obs),
    Strength_obs = s_obs,
    Strength_norm = s_norm,
    Strength_mean_null = s_mean,
    Strength_sd_null = s_sd,
    Strength_z = s_z
  )
}

str_NH <- calc_strength_null(TR_NH)
str_HU <- calc_strength_null(TR_HU)
str_OH <- calc_strength_null(TR_OH)

str_NH$Hunting_Pressure <- "Non_hunted"
str_HU$Hunting_Pressure <- "Hunted"
str_OH$Hunting_Pressure <- "Overhunted"

TR_strength_null <- rbind(str_NH, str_HU, str_OH)

TR_strength_null$Dimension <- "Trophic"

TR_strength_null$Hunting_Pressure <- factor(
  TR_strength_null$Hunting_Pressure,
  levels = c("Non_hunted", "Hunted", "Overhunted"),
  labels = c("Non_hunted", "Hunted", "Overhunted")
)

TR_strength_null$Significance <- ifelse(TR_strength_null$Strength_z > 1.96, "Positive",
                                        ifelse(TR_strength_null$Strength_z < -1.96, "Negative", "n.s."))

# tabella riassuntiva solo dei significativi
TR_strength_null[TR_strength_null$Significance != "n.s.",
                 c("Genus", "Hunting_Pressure", "Strength_z", "Significance")]


TR_centr_plot <- ggplot(TR_strength_null, aes(
  x = Strength_norm,
  y = reorder_within(Genus, Strength_norm, Hunting_Pressure)
)) +
  
  geom_segment(aes(
    x = 0, xend = Strength_norm,
    y = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    yend = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    color = Strength_z
  )) +
  
  # punti (testa del lollipop)
  geom_point(aes(color = Strength_z), size = 2) +
  
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "gold", "#F4A582", "#B2182B"),
    values = scales::rescale(c(-6, -2, 0, 1, 3)),
    limits = c(-6, 3)
  )+
  
  # linea zero (non necessaria ma utile)
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Hunting_Pressure, scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2)
  ) +
  
  scale_y_reordered() +
  
  labs(
    x = "Strength (normalized)",
    y = NULL,
    col = "z-score\n(null model)"
  ) +
  
  theme_bw()

TR_centr_plot

###Overhunted
TR_OH <- as.matrix(TR_OH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(TR_OH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_OH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_OH <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_OH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_OH <- (Q_real_OH - mean(Q_null_OH)) / sd(Q_null_OH)

Q_real_OH #0.03 (0)
Z_score_OH #-0.46 (-1.35)

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0.22 (0)
mod_TR_OH <- c(obs = Q_real_OH, Z = Z_score_OH, ratio = Q_ratio_OH)

##Nestedness
mat <- as.matrix(TR_OH)
diag(mat) <- 0

nested_real_OH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_OH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_OH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_OH <- (nested_real_OH - mean(nested_null_OH)) / sd(nested_null_OH)

nested_real_OH #60.3 (65.9)
Z_score_OH #0.61 (1.75)

N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH #0.88 (0.93)
nest_TR_OH <- c(obs = nested_real_OH, Z = Z_score_OH, ratio = N_ratio_OH)

##Nested contribution
contrib_OH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_OH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_OH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_OH$contribution[k] <- (nested_real_OH - mean(null_k)) / sd(null_k)
}

contrib_OH[order(-contrib_OH$contribution), ]

null_df_OH <- data.frame(node = rep(names(null_list_OH), sapply(null_list_OH, length)),
                         WNODA_null = unlist(null_list_OH))
lev_OH <- contrib_OH$node[order(contrib_OH$contribution)]
null_df_OH$node <- factor(null_df_OH$node, levels = lev_OH)
contrib_OH$node_f <- factor(contrib_OH$node, levels = lev_OH)
x_lab_OH <- max(null_df_OH$WNODA_null, nested_real_OH) + 1.5

plot_contrib_TR_OH <- ggplot(subset(contrib_OH, !is.na(contribution)),
                             aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#EABD20") +
  geom_point(size = 3, colour = "#EABD20") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Trophic-OH") +
  theme_bw()

plot_contrib_TR_OH

###Hunted
TR_HU <- as.matrix(TR_HU)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(TR_HU, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_HU <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_HU <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_HU[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_HU <- (Q_real_HU - mean(Q_null_HU)) / sd(Q_null_HU)

Q_real_HU #0.14 (0.14)
Z_score_HU #1.36 (1.02)

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0.62 (0.58)
mod_TR_HU <- c(obs = Q_real_HU, Z = Z_score_HU, ratio = Q_ratio_HU)

##Nestedness
mat <- as.matrix(TR_HU)
diag(mat) <- 0

nested_real_HU <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_HU <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_HU[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_HU <- (nested_real_HU - mean(nested_null_HU)) / sd(nested_null_HU)

nested_real_HU #46.42 (46.8)
Z_score_HU #-1.53 (-1.05)

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU #0.68 (0.7)
nest_TR_HU <- c(obs = nested_real_HU, Z = Z_score_HU, ratio = N_ratio_HU)

###Non_hunted
TR_NH <- as.matrix(TR_NH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(TR_NH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_NH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_NH <- numeric(n_sim)

for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_NH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_NH <- (Q_real_NH - mean(Q_null_NH)) / sd(Q_null_NH)

Q_real_NH #0.19 (0.2)
Z_score_NH #4.07 (4.85)
membership(comm_real_over)

g_TR_NH    <- g_over
comm_TR_NH <- comm_real_over

saveRDS(list(g = g_TR_NH, comm = comm_TR_NH), "Results/TR_modularity_NH.rds")

#Save the membership for the figure
membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(membership_NH, "Results/TR_membership_NH.rds")

Q_ratio_Non_hunted <- Q_real_NH/ max(Q_null_NH)
Q_ratio_Non_hunted #1.07 (1.09)
mod_TR_NH <- c(obs = Q_real_NH, Z = Z_score_NH, ratio = Q_ratio_Non_hunted)


##Nestedness
mat <- as.matrix(TR_NH)
diag(mat) <- 0

nested_real_NH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_NH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_NH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_NH <- (nested_real_NH - mean(nested_null_NH)) / sd(nested_null_NH)

nested_real_NH #52.4 (53.7)
Z_score_NH #.1.46 (-0.95)

nested_ratio_NH <- nested_real_NH/ max(nested_null_NH)
nested_ratio_NH #0.83 (0.83)
nest_TR_NH <- c(obs = nested_real_NH, Z = Z_score_NH, ratio = nested_ratio_NH)

## Modularity plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(Q_null_NH),
                mean(Q_null_HU),
                mean(Q_null_OH)),
  sd_null = c(sd(Q_null_NH),
              sd(Q_null_HU),
              sd(Q_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  Q_real = c(Q_real_NH, Q_real_HU, Q_real_OH),
  colore = c("#EABD20","#EABD20",  "#EABD20")  # Non_hunted, Hunted, Overhunted
)

# Plot
mod <- ggplot() +
  # Punto null model con barre di errore
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  # Punto modularità reale
  geom_point(data = df_real, aes(x = rete, y = Q_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Modularity")+
  theme_bw(base_size = 14)

mod

#Nestedness plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(nested_null_NH),
                mean(nested_null_HU),
                mean(nested_null_OH)),
  sd_null = c(sd(nested_null_NH),
              sd(nested_null_HU),
              sd(nested_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  nested_real = c(nested_real_NH, nested_real_HU, nested_real_OH),
  colore = c("#EABD20","#EABD20", "#EABD20")  # Non_hunted, Hunted, Overhunted
)

nest <- ggplot() +
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  geom_point(data = df_real, aes(x = rete, y = nested_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Nestedness (WNODA)") +
  theme_bw(base_size = 14)

#Combined plot
Tr_met <- ggarrange(mod, nest, nrow = 1)
Tr_met

##### Combined #####
C_OH <- as.matrix(read.csv("Results/Tot_Overhunted.csv", row.names = 1))
C_HU <- as.matrix(read.csv("Results/Tot_Hunted.csv", row.names = 1))
C_NH <- as.matrix(read.csv("Results/Tot_Non_hunted.csv", row.names = 1))

# Select the order of the species
species_order <- c("Alouatta", "Aotus", "Ateles", "Cacajao", "Callicebus", 
                   "Cebus", "Lagothrix", "Pithecia", "Saguinus", "Saimiri",
                   "Sapajus")

#Create a function to reorder a matrices
reorder_matrix <- function(mat, order) {
  species_present <- order[order %in% rownames(mat)]
  mat[species_present, species_present]
}

#Reorder primates in matrices
C_OH <- reorder_matrix(C_OH, species_order)
C_HU <- reorder_matrix(C_HU, species_order)
C_NH <- reorder_matrix(C_NH, species_order)

#Strength centrality
set.seed(123)
n_sim <- 1000

calc_strength_null <- function(mat){
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "numeric"
  diag(mat) <- 0
  
  g <- graph_from_adjacency_matrix(
    mat,
    mode = "undirected",
    weighted = TRUE,
    diag = FALSE
  )
  N <- vcount(g)
  
  
  # strength osservata
  s_obs <- igraph::strength(g, mode = "all", weights = E(g)$weight)
  s_norm <- s_obs / (N - 1)
  
  # matrice null
  s_null <- matrix(NA, nrow = vcount(g), ncol = n_sim)
  rownames(s_null) <- names(s_obs)
  
  for (i in 1:n_sim) {
    
    g_rand <- g
    
    # stesso null model della modularità: shuffle pesi
    E(g_rand)$weight <- sample(E(g)$weight)
    
    s_null[, i] <- igraph::strength(g_rand, mode = "all")[names(s_obs)]
  }
  
  # statistiche null
  s_mean <- rowMeans(s_null)
  s_sd   <- apply(s_null, 1, sd)
  
  # z-score
  s_z <- (s_obs - s_mean) / s_sd
  
  data.frame(
    Genus = names(s_obs),
    Strength_obs = s_obs,
    Strength_norm = s_norm,
    Strength_mean_null = s_mean,
    Strength_sd_null = s_sd,
    Strength_z = s_z
  )
}

str_NH <- calc_strength_null(C_NH)
str_HU <- calc_strength_null(C_HU)
str_OH <- calc_strength_null(C_OH)

str_NH$Hunting_Pressure <- "Non_hunted"
str_HU$Hunting_Pressure <- "Hunted"
str_OH$Hunting_Pressure <- "Overhunted"

C_strength_null <- rbind(str_NH, str_HU, str_OH)

C_strength_null$Hunting_Pressure <- factor(
  C_strength_null$Hunting_Pressure,
  levels = c("Non_hunted", "Hunted", "Overhunted"),
  labels = c("Non_hunted", "Hunted", "Overhunted")
)

C_strength_null$Significance <- ifelse(C_strength_null$Strength_z > 1.96, "Positive",
                                       ifelse(C_strength_null$Strength_z < -1.96, "Negative", "n.s."))

# tabella riassuntiva solo dei significativi
C_strength_null[C_strength_null$Significance != "n.s.",
                c("Genus", "Hunting_Pressure", "Strength_norm", "Strength_z", "Significance")]


C_centr_plot <- ggplot(C_strength_null, aes(
  x = Strength_norm,
  y = reorder_within(Genus, Strength_norm, Hunting_Pressure)
)) +
  
  geom_segment(aes(
    x = 0, xend = Strength_norm,
    y = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    yend = reorder_within(Genus, Strength_norm, Hunting_Pressure),
    color = Strength_z
  )) +
  
  # punti (testa del lollipop)
  geom_point(aes(color = Strength_z), size = 2) +
  
  scale_color_gradientn(
    colours = c("#2166AC", "#67A9CF", "gold", "#F4A582", "#B2182B"),
    values = scales::rescale(c(-6, -2, 0, 1, 3)),
    limits = c(-6, 3)
  )+
  
  # linea zero (non necessaria ma utile)
  geom_vline(xintercept = 0, linetype = "dashed") +
  
  facet_wrap(~Hunting_Pressure, scales = "free_y") +
  
  scale_x_continuous(
    limits = c(0, 3),
    breaks = seq(0, 3, 0.5)
  ) +
  
  scale_y_reordered() +
  
  labs(
    x = "Strength (normalized)",
    y = NULL,
    fill = "z-score\n(null model)"
  ) +
  
  theme_bw()

C_centr_plot


###Overhunted
C_OH <- as.matrix(C_OH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(C_OH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_OH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_OH <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_OH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_OH <- (Q_real_OH - mean(Q_null_OH)) / sd(Q_null_OH)

Q_real_OH # 0
Z_score_OH # -1.19 (-1.06)

Q_ratio_OH <- Q_real_OH/ max(Q_null_OH)
Q_ratio_OH #0
mod_C_OH <- c(obs = Q_real_OH, Z = Z_score_OH, ratio = Q_ratio_OH)

##Nestedness
mat <- as.matrix(C_OH)
diag(mat) <- 0

nested_real_OH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_OH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_OH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_OH <- (nested_real_OH - mean(nested_null_OH)) / sd(nested_null_OH)

nested_real_OH # 70.63 (67.46)
Z_score_OH # 2.41 (1.62)

nested_ratio_OH <- nested_real_OH/ max(nested_null_OH)
nested_ratio_OH #0.96 (0.91)
nest_C_OH <- c(obs = nested_real_OH, Z = Z_score_OH, ratio = nested_ratio_OH)

##Nested contribution
contrib_OH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_OH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_OH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_OH$contribution[k] <- (nested_real_OH - mean(null_k)) / sd(null_k)
}

contrib_OH[order(-contrib_OH$contribution), ]

plot_contrib_C_OH <- ggplot(subset(contrib_OH, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#BE956E") +
  geom_point(size = 3, colour = "#BE956E") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Combined-OH") +
  theme_bw()

plot_contrib_C_OH

###Hunted
C_HU <- as.matrix(C_HU)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(C_HU, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_HU <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_HU <- numeric(n_sim)
for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_HU[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_HU <- (Q_real_HU - mean(Q_null_HU)) / sd(Q_null_HU)

Q_real_HU # 0.06 (0.05)
Z_score_HU # 0.96 (0.19)

Q_ratio_HU <- Q_real_HU/ max(Q_null_HU)
Q_ratio_HU #0.48 (0.30)
mod_C_HU <- c(obs = Q_real_HU, Z = Z_score_HU, ratio = Q_ratio_HU)

##Nestedness
mat <- as.matrix(C_HU)
diag(mat) <- 0

nested_real_HU <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_HU <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_HU[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_HU <- (nested_real_HU - mean(nested_null_HU)) / sd(nested_null_HU)

nested_real_HU #61.9 (60)
Z_score_HU #1.32 (1.06)

nested_ratio_HU <- nested_real_HU/ max(nested_null_HU)
nested_ratio_HU #0.88 (0.84)
nest_C_HU <- c(obs = nested_real_HU, Z = Z_score_HU, ratio = nested_ratio_HU)

###Non_hunted
C_NH <- as.matrix(C_NH)
set.seed(123)
n_sim <- 1000

##Modularity
g_over <- igraph::graph_from_adjacency_matrix(C_NH, mode = "undirected", weighted = TRUE, diag = FALSE)
comm_real_over <- cluster_louvain(g_over, weights = E(g_over)$weight)

Q_real_NH <- modularity(g_over, membership(comm_real_over), weights = E(g_over)$weight)

Q_null_NH <- numeric(n_sim)

for (i in 1:n_sim) {
  g_rand <- g_over
  E(g_rand)$weight <- sample(E(g_over)$weight)
  comm_rand <- cluster_louvain(g_rand, weights = E(g_rand)$weight)
  Q_null_NH[i] <- modularity(g_rand, membership(comm_rand), weights = E(g_rand)$weight)
}

Z_score_NH <- (Q_real_NH - mean(Q_null_NH)) / sd(Q_null_NH)

Q_real_NH  # 0.08 (0.09)
Z_score_NH # 2.87 (2.95)
membership(comm_real_over)

g_C_NH    <- g_over
comm_C_NH <- comm_real_over

saveRDS(list(g = g_C_NH,  comm = comm_C_NH),  "Results/C_modularity_NH.rds")

#Save the membership for the figure
membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(membership_NH, "Results/C_membership_NH.rds")

Q_ratio_NH <- Q_real_NH/ max(Q_null_NH)
Q_ratio_NH #0.94 (0.73)
mod_C_NH <- c(obs = Q_real_NH, Z = Z_score_NH, ratio = Q_ratio_NH)





##Nestedness
mat <- as.matrix(C_NH)
diag(mat) <- 0

nested_real_NH <- nest.smdm(mat, weighted = TRUE, decreasing = "abund")$WNODAmatrix

nested_null_NH <- numeric(n_sim)

for(i in 1:n_sim){
  mat_rand <- mat
  mat_rand[upper.tri(mat_rand)] <- 0
  low <- mat[lower.tri(mat)]
  low[low > 0] <- sample(low[low > 0])
  mat_rand[lower.tri(mat_rand)] <- low
  mat_rand <- mat_rand + t(mat_rand)
  nested_null_NH[i] <- nest.smdm(mat_rand, weighted = TRUE, decreasing = "abund")$WNODAmatrix
}

Z_score_NH <- (nested_real_NH - mean(nested_null_NH)) / sd(nested_null_NH)

nested_real_NH # 62 (60)
Z_score_NH # 1.37 (0.82)

nested_ratio_NH <- nested_real_NH/ max(nested_null_NH)
nested_ratio_NH #0.92 (0.92)
nest_C_NH <- c(obs = nested_real_NH, Z = Z_score_NH, ratio = nested_ratio_NH)

##Nested contribution
contrib_NH <- data.frame(node = rownames(mat), contribution = NA_real_)
null_list_NH <- list()

for(k in 1:nrow(mat)){
  others <- setdiff(1:nrow(mat), k)
  if(sum(mat[k, others] > 0) == 0) next
  
  null_k <- numeric(n_sim)
  for(i in 1:n_sim){
    m <- mat
    w <- mat[k, others]
    w[w > 0] <- sample(w[w > 0])
    m[k, others] <- w
    m[others, k] <- w
    null_k[i] <- nest.smdm(m, weighted = TRUE, decreasing = "abund")$WNODAmatrix
  }
  null_list_NH[[rownames(mat)[k]]] <- null_k
  if(sd(null_k) > 0) contrib_NH$contribution[k] <- (nested_real_NH - mean(null_k)) / sd(null_k)
}

contrib_NH[order(-contrib_NH$contribution), ]

null_df_NH <- data.frame(node = rep(names(null_list_NH), sapply(null_list_NH, length)),
                         WNODA_null = unlist(null_list_NH))
lev_NH <- contrib_NH$node[order(contrib_NH$contribution)]
null_df_NH$node <- factor(null_df_NH$node, levels = lev_NH)
contrib_NH$node_f <- factor(contrib_NH$node, levels = lev_NH)
x_lab_NH <- max(null_df_NH$WNODA_null, nested_real_NH) + 1.5

plot_contrib_C_NH <- ggplot(subset(contrib_NH, !is.na(contribution)),
                            aes(x = contribution, y = reorder(node, contribution))) +
  geom_vline(xintercept = 0, colour = "grey50") +
  geom_vline(xintercept = c(-1.96, 1.96), linetype = "dashed", colour = "grey50") +
  geom_segment(aes(x = 0, xend = contribution, yend = reorder(node, contribution)),
               colour = "#BE956E") +
  geom_point(size = 3, colour = "#BE956E") +
  labs(x = "Nested contribution (z-score)", y = NULL, title = "Combined-NH") +
  theme_bw()

plot_contrib_C_NH

## Modularity plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(Q_null_NH),
                mean(Q_null_HU),
                mean(Q_null_OH)),
  sd_null = c(sd(Q_null_NH),
              sd(Q_null_HU),
              sd(Q_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  Q_real = c(Q_real_NH, Q_real_HU, Q_real_OH),
  colore = c("#BE956E","#BE956E",  "#BE956E")  # Non_hunted, Hunted, Overhunted
)

# Plot
mod <- ggplot() +
  # Punto null model con barre di errore
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  # Punto modularità reale
  geom_point(data = df_real, aes(x = rete, y = Q_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Modularity")+
  theme_bw(base_size = 14)

mod

#Nestedness plot
df_null_stats <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  mean_null = c(mean(nested_null_NH),
                mean(nested_null_HU),
                mean(nested_null_OH)),
  sd_null = c(sd(nested_null_NH),
              sd(nested_null_HU),
              sd(nested_null_OH))
)

df_real <- data.frame(
  rete = factor(c("Non_hunted", "Hunted", "Overhunted"),
                levels = c("Non_hunted", "Hunted", "Overhunted")),
  nested_real = c(nested_real_NH, nested_real_HU, nested_real_OH),
  colore = c("#BE956E","#BE956E", "#BE956E")  # Non_hunted, Hunted, Overhunted
)

nest <- ggplot() +
  geom_point(data = df_null_stats, aes(x = rete, y = mean_null), color = "gray40", size = 3) +
  geom_errorbar(data = df_null_stats, aes(x = rete,
                                          ymin = mean_null - sd_null,
                                          ymax = mean_null + sd_null),
                width = 0.05, color = "black") +
  geom_point(data = df_real, aes(x = rete, y = nested_real, color = colore), size = 4) +
  scale_color_identity() +
  labs(x = "", y = "Nestedness (WNODA)") +
  theme_bw(base_size = 14)

#Combined plot
C_met <- ggarrange(mod, nest, nrow = 1)
C_met


#### Summary table of network metrics ####
mod_all <- rbind(mod_T_NH, mod_T_HU, mod_T_OH,
                 mod_S_NH, mod_S_HU, mod_S_OH,
                 mod_TR_NH, mod_TR_HU, mod_TR_OH,
                 mod_C_NH, mod_C_HU, mod_C_OH)
nest_all <- rbind(nest_T_NH, nest_T_HU, nest_T_OH,
                  nest_S_NH, nest_S_HU, nest_S_OH,
                  nest_TR_NH, nest_TR_HU, nest_TR_OH,
                  nest_C_NH, nest_C_HU, nest_C_OH)

metrics_table <- data.frame(
  Network = rep(c("Temporal", "Spatial", "Trophic", "Combined"), each = 3),
  Hunting_Pressure = rep(c("Non-hunted", "Hunted", "Overhunted"), times = 4),
  `Observed modularity (z-score)` = sprintf("%.2f (%.2f)", mod_all[, "obs"], mod_all[, "Z"]),
  `Standardized modularity` = round(mod_all[, "ratio"], 2),
  `Observed nestedness (z-score)` = sprintf("%.2f (%.2f)", nest_all[, "obs"], nest_all[, "Z"]),
  `Standardized nestedness` = round(nest_all[, "ratio"], 2),
  check.names = FALSE,
  row.names = NULL
)
metrics_table

write.csv(metrics_table, "Results/Network_metrics_table.csv", row.names = FALSE)

library(flextable)

ft <- flextable(metrics_table)
ft <- autofit(ft)
ft <- merge_v(ft, j = "Network")
ft <- theme_booktabs(ft)

doc <- read_docx()
doc <- body_add_flextable(doc, ft)
print(doc, target = "Results/Network_metrics_table.docx")



## Nested contribution (all networks)
Nest_Con_all <- ggarrange(plot_contrib_T_NH, plot_contrib_T_HU, plot_contrib_T_OH,
                          plot_contrib_S_OH, plot_contrib_TR_OH, plot_contrib_C_NH,plot_contrib_C_OH,
                          ncol = 3, nrow = 3,
                          labels = "AUTO")
Nest_Con_all


#### Save figures and matrices ####
## Centrality
#Matrices
write.csv(T_strength_null, "Results/T_centrality.csv", row.names = TRUE)
write.csv(S_strength_null, "Results/S_centrality.csv", row.names = TRUE)
write.csv(TR_strength_null, "Results/TR_centrality.csv", row.names = TRUE)

#Plots
centr_plot <- ggarrange(T_centr_plot, S_centr_plot, TR_centr_plot, C_centr_plot, ncol = 1, common.legend = T)

ggsave("Figures/Centrality.jpg", centr_plot, width = 11, height = 14)
ggsave("Figures/Centrality.pdf", centr_plot, width = 11, height = 14)

## Modularity + Nestedness
#Temporal
ggsave("Figures/Metrics_Temporal.pdf", T_met, width = 6, height = 6)
ggsave("Figures/NestedContribution.pdf", Nest_Con, width = 12, height = 6)
ggsave("Figures/NestedContribution.jpg", Nest_Con, width = 12, height = 6)

#Spatial
ggsave("Figures/Metrics_Spatial.pdf", S_met, width = 6, height = 6)

#Trophic
ggsave("Figures/Metrics_Trophic.pdf", Tr_met, width = 6, height = 6)

#Combined
ggsave("Figures/Metrics_Combined.pdf", C_met, width = 6, height = 6)




ggsave("Figures/NestedContribution_all.png", Nest_Con_all,
       width = 14, height = 9, dpi = 300, bg = "white")
ggsave("Figures/NestedContribution_all.pdf", Nest_Con_all,
       width = 14, height = 9)
