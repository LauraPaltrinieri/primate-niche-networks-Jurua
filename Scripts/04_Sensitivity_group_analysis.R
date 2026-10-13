#Sensitivity analysis: group-level independent detections
#This script combines the Niche Overlap, Network Metrics, and Network/Ternary Plot scripts into
#a single pipeline, to test whether using the primate group (rather than the single individual)
#as the independent detection unit changes the overall pattern of results. All calculations follow
#the same methodology as the main analysis; the only difference is the independent detection unit.
#07/09/2026

# 01_Niche Overlap ####
#This script calculates pairwise niche overlap between primate genera along three dimensions
#(Temporal, Spatial, Trophic) for each hunting-pressure treatment (Non-hunted, Hunted, Overhunted),
#using activity pattern overlap (Temporal) and Pianka's index (Spatial, Trophic). It then combines
#the three dimensions into a single overlap index and calculates the relative contribution of each
#dimension to the combined overlap. All matrices are exported for use in the next scripts.

#Libraries
library(activity)   #solartime (Temporal)
library(overlap)    #temporal overlap estimation (Temporal)
library(circular)   #conversion to circular data (Temporal)
library(dplyr)      
library(tidyr)      
library(tibble)     
library(spaa)       #Pianka index (Spatial, Trophic)

#clean_env
rm(list = ls())

##### Upload the data #####
data <- read.csv("Data/Group_analysis/Primates_detections_g.csv")

#Filter for minimum 8 detections for primate genus in each site
data <- data %>%
  add_count(Vertebrate_ID, Hunting_Pressure, name = "n_detections") %>%
  filter(n_detections >= 10) 

#Create a function to reorder the matrices of this analysis with the same order
#Select the order
species_order <- c("Alouatta", "Aotus", "Ateles", "Cacajao", "Callicebus", 
                   "Cebus", "Lagothrix", "Pithecia", "Saguinus", "Saimiri",
                   "Sapajus")

#Create a function to reorder the matrices
reorder_matrix <- function(mat, order) {
  species_present <- order[order %in% rownames(mat)]
  mat[species_present, species_present]}


##### Temporal Overlap #####
#Temporal data
data_T <- data

#Set the DateTime format
data_T$DateTime <- as.POSIXct(data_T$DateTime, format = "%Y-%m-%d %H:%M:%S")

#Calculate the solartime: the sunset and sunrise exact times based on my study area and the period of the detection
solartime <- activity::solartime(data_T$DateTime, lat= data_T$Latitude , lon= data_T$Longitude, tz= 1  ,format ="%d/%m/%Y %H:%M:%S")
solartime

#Create a new column of the solar time, extracting it from "solartime"
data_T$solar <- solartime$solar
data_T$solar_input <- solartime$input

#Set the solar time as: circular
data_T$solar_circ <- circular(data_T$solar, units = "radians")
class(data_T$solar_circ) #has to be circular

#Create matrices of temporal overlap between every pair of primates in each site
overlap_matrices <- list()

hp <- unique(data_T$Hunting_Pressure)

#Loop
for (hunting_pressure in hp) {
  hp_data <- filter(data_T, Hunting_Pressure == hunting_pressure) #Select the treatment
  genera <- unique(hp_data$Vertebrate_ID) #Extract the unique primates genera in the Site
  n_gen <- length(genera) #Count of number of genera in the specific site to create a matrix with exact number of columns and rows
  
  mat <- matrix(NA, nrow = n_gen, ncol = n_gen, dimnames = list(genera, genera)) #Create an empty matrix with exact number of lines and columns
  
  for (i in seq_len(n_gen)) { #i = name of columns
    for (j in seq_len(n_gen)) { #j = name of lines
      g1 <- genera[i] #Assign the name of columns
      g2 <- genera[j] #Assign the name of lines
      
      if (g1 == g2) {
        mat[i, j] <- 1 #if the name of line and column is the same (same genus) so the value has to be 1
      } else { #if not
        n1 <- sum(hp_data$Vertebrate_ID == g1) #count the number of observation for each genus
        n2 <- sum(hp_data$Vertebrate_ID == g2)
        
        dhat_type <- ifelse(min(n1, n2) < 75, "Dhat1", "Dhat4") #use delta1 if the n of observations is <75, otherwise use delta4
        
        dhat_val <- overlapEst( #calculate the temporal overlap
          hp_data$solar_circ[hp_data$Vertebrate_ID == g1],
          hp_data$solar_circ[hp_data$Vertebrate_ID == g2],
          type = dhat_type #use delta1 or delta4 based of what selected before
        )
        
        mat[i, j] <- dhat_val #write the overlap value in the matrix
      }
    }
  }
  
  overlap_matrices[[hunting_pressure]] <- mat
}

#Extract each matrix for each treatment
T_OH <- as.data.frame(overlap_matrices[["Overhunted"]])
T_HU <- as.data.frame(overlap_matrices[["Hunted"]])
T_NH <- as.data.frame(overlap_matrices[["Non_hunted"]])

#Reorder temporal matrices based on the previous function
T_OH <- reorder_matrix(T_OH, species_order)
T_HU <- reorder_matrix(T_HU, species_order)
T_NH <- reorder_matrix(T_NH, species_order)

#Clean env.
rm(hunting_pressure, hp, overlap_matrices, solartime, data_T, hp_data, mat)

##### Spatial Overlap #####
#Spatial data
data_S <- data

#Count the occurrences of each primate genus in each Station_ID
data_S1 <- data_S %>%
  group_by(Tree_ID, Hunting_Pressure, Genus, Effort_H_Station, Focal_Genus, Sampled_Canopy) %>% 
  summarise(count = n()) %>%
  pivot_wider(names_from = Genus, values_from = count, values_fill = 0)

#Upload info about trees correctly sampled but with no primate detections 
#(and thus absent from the main detections dataset), needed to correctly 
#calculate foraging rates for each primate genus per plant genus
CT <- read.csv("Data/Missing_TreeID.csv")

missing_cols <- setdiff(names(data_S1), names(CT))
CT[missing_cols] <- 0
CT <- CT %>% select(all_of(names(data_S1)))

data_S1 <- rbind(data_S1, CT)

#Correct feeding interaction counts for the proportion of canopy actually sampled (standardized to 100% canopy coverage)
data_S2 <- data_S1 %>%
  mutate(across(all_of(c("Callicebus", "Sapajus", "Alouatta", "Aotus", 
                         "Ateles","Pithecia", "Cebus", "Saimiri", "Cacajao",
                         "Lagothrix", "Saguinus")),~ (.x / Sampled_Canopy) * 100))

#Sum corrected interactions per tree and calculate the Relative Abundance Index (RAI)  
#for each primate genus, standardized by the sampling effort
data_S3 <- data_S2 %>%
  group_by(Focal_Genus, Hunting_Pressure, Tree_ID, Effort_H_Station) %>%
  summarise(Callicebus = sum(Callicebus),
            Sapajus = sum(Sapajus),
            Alouatta = sum(Alouatta),
            Aotus = sum(Aotus),
            Ateles = sum(Ateles),
            Pithecia = sum(Pithecia),
            Cebus = sum(Cebus),
            Saimiri = sum(Saimiri),
            Cacajao = sum(Cacajao),
            Lagothrix = sum(Lagothrix),
            Saguinus = sum(Saguinus)) %>%
  mutate(RAI_Callicebus = (Callicebus/Effort_H_Station)*100,
         RAI_Sapajus = (Sapajus/Effort_H_Station)*100,
         RAI_Alouatta = (Alouatta/Effort_H_Station)*100,
         RAI_Aotus =(Aotus/Effort_H_Station)*100,
         RAI_Ateles = (Ateles/Effort_H_Station)*100,
         RAI_Pithecia = (Pithecia/Effort_H_Station)*100,
         RAI_Cebus = (Cebus/Effort_H_Station)*100,
         RAI_Saimiri = (Saimiri/Effort_H_Station)*100,
         RAI_Cacajao = (Cacajao/Effort_H_Station)*100,
         RAI_Lagothrix = (Lagothrix/Effort_H_Station)*100,
         RAI_Saguinus = (Saguinus/Effort_H_Station)*100) %>% 
  select(Hunting_Pressure,Focal_Genus, Tree_ID,
         RAI_Callicebus,RAI_Sapajus, RAI_Alouatta,
         RAI_Aotus, RAI_Ateles, RAI_Pithecia, RAI_Cebus,
         RAI_Saimiri,RAI_Cacajao, RAI_Lagothrix,RAI_Saguinus)

data_S <- data_S3
data_S <- data_S %>% ungroup()

#Create a different rai matrix for each treatment:
#Overhunted
S_OH <- data_S %>% 
  filter(Hunting_Pressure == "Overhunted") %>%
  select(Tree_ID, starts_with("RAI_")) %>%
  column_to_rownames(var = "Tree_ID") %>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Hunted
S_HU <- data_S %>% 
  filter(Hunting_Pressure == "Hunted") %>%
  select(Tree_ID, starts_with("RAI_")) %>%
  column_to_rownames(var = "Tree_ID")%>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Non_hunted
S_NH <- data_S %>% 
  filter(Hunting_Pressure == "Non_hunted") %>%
  select(Tree_ID, starts_with("RAI_")) %>%
  column_to_rownames(var = "Tree_ID")%>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Calculate the Spatial-overlap with Pianka index
S_OH <- as.data.frame(as.matrix(as.dist(niche.overlap(S_OH, method = "pianka"))))
diag(S_OH) <- 1

S_HU <- as.data.frame(as.matrix(as.dist(niche.overlap(S_HU, method = "pianka"))))
diag(S_HU) <- 1

S_NH <- as.data.frame(as.matrix(as.dist(niche.overlap(S_NH, method = "pianka"))))
diag(S_NH) <- 1

#Reorder Spatial matrices
S_OH <- reorder_matrix(S_OH, species_order)
S_HU <- reorder_matrix(S_HU, species_order)
S_NH <- reorder_matrix(S_NH, species_order)

#Clean env.
rm(list = c("CT", "data_S", "data_S1", "data_S2", "data_S3",ls(pattern = "^RAI_")))

##### Trophic Overlap #####
#Trophic data
data_TR <- data

#Aggregate confirmed feeding interactions (Interaction = Yes: CI = 5 or Behaviour = looking for food & CI = 4) 
#and reshape to wide format 
data_TR1 <- data_TR %>%
  group_by(Focal_Genus, Hunting_Pressure, Genus, Effort_H_Focal_HP, Tree_ID, Camera_ID, Effort_H_Station, Sampled_Canopy) %>%
  summarise(count = sum(Interaction == "Yes"),.groups = "drop") %>%
  pivot_wider(names_from = Genus, values_from = count, values_fill = 0)

#Upload info about trees correctly sampled but with no primate detections 
#(and thus absent from the main detections dataset), needed to correctly 
#calculate foraging rates for each primate genus per plant genus
CT <- read.csv("Data/Missing_TreeID.csv")

missing_cols <- setdiff(names(data_TR1), names(CT))
CT[missing_cols] <- 0
CT <- CT %>% select(all_of(names(data_TR1)))

data_TR1 <- rbind(data_TR1, CT)

#Correct feeding interaction counts for the proportion of canopy actually sampled (standardized to 100% canopy coverage)
data_TR2 <- data_TR1 %>%
  mutate(across(all_of(c("Callicebus", "Sapajus", "Alouatta", "Aotus", 
                         "Ateles","Pithecia", "Cebus", "Saimiri", "Cacajao",
                         "Lagothrix", "Saguinus")),~ (.x / Sampled_Canopy) * 100))

#Sum corrected interactions per tree and calculate the Relative Abundance Index (RAI)  
#for each primate genus, standardized by the sampling effort
data_TR3 <- data_TR2 %>%
  group_by(Focal_Genus, Hunting_Pressure, Tree_ID, Effort_H_Station) %>%
  summarise(Callicebus = sum(Callicebus),
            Sapajus = sum(Sapajus),
            Alouatta = sum(Alouatta),
            Aotus = sum(Aotus),
            Ateles = sum(Ateles),
            Pithecia = sum(Pithecia),
            Cebus = sum(Cebus),
            Saimiri = sum(Saimiri),
            Cacajao = sum(Cacajao),
            Lagothrix = sum(Lagothrix),
            Saguinus = sum(Saguinus)) %>%
  mutate(RAI_Callicebus = (Callicebus/Effort_H_Station)*100,
         RAI_Sapajus = (Sapajus/Effort_H_Station)*100,
         RAI_Alouatta = (Alouatta/Effort_H_Station)*100,
         RAI_Aotus =(Aotus/Effort_H_Station)*100,
         RAI_Ateles = (Ateles/Effort_H_Station)*100,
         RAI_Pithecia = (Pithecia/Effort_H_Station)*100,
         RAI_Cebus = (Cebus/Effort_H_Station)*100,
         RAI_Saimiri = (Saimiri/Effort_H_Station)*100,
         RAI_Cacajao = (Cacajao/Effort_H_Station)*100,
         RAI_Lagothrix = (Lagothrix/Effort_H_Station)*100,
         RAI_Saguinus = (Saguinus/Effort_H_Station)*100) %>% 
  select(Hunting_Pressure,Focal_Genus, Tree_ID,
         RAI_Callicebus,RAI_Sapajus, RAI_Alouatta,
         RAI_Aotus, RAI_Ateles, RAI_Pithecia, RAI_Cebus,
         RAI_Saimiri,RAI_Cacajao, RAI_Lagothrix,RAI_Saguinus)

#Average RAI across trees to obtain the mean foraging rate of each primate genus 
#per plant genus and hunting-pressure category
data_TR4 <- data_TR3 %>%
  group_by(Focal_Genus, Hunting_Pressure) %>%  
  summarise(RAI_Callicebus = mean(RAI_Callicebus),
            RAI_Sapajus =  mean(RAI_Sapajus),
            RAI_Alouatta =  mean(RAI_Alouatta),
            RAI_Aotus = mean(RAI_Aotus),
            RAI_Ateles =  mean(RAI_Ateles),
            RAI_Pithecia = mean(RAI_Pithecia),
            RAI_Cebus =  mean(RAI_Cebus),
            RAI_Saimiri =  mean(RAI_Saimiri),
            RAI_Cacajao =  mean(RAI_Cacajao),
            RAI_Lagothrix =  mean(RAI_Lagothrix),
            RAI_Saguinus =  mean(RAI_Saguinus))

data_TR <- data_TR4
data_TR <- data_TR %>% ungroup()

#Create a different rai matrix for each treatment:
#Overhunted
TR_OH <- data_TR %>% 
  filter(Hunting_Pressure == "Overhunted") %>%
  select(Focal_Genus, starts_with("RAI_")) %>%
  column_to_rownames(var = "Focal_Genus") %>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Hunted
TR_HU <- data_TR %>% 
  filter(Hunting_Pressure == "Hunted") %>%
  select(Focal_Genus, starts_with("RAI_")) %>%
  column_to_rownames(var = "Focal_Genus")%>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Non_hunted
TR_NH <- data_TR %>% 
  filter(Hunting_Pressure == "Non_hunted") %>%
  select(Focal_Genus, starts_with("RAI_")) %>%
  column_to_rownames(var = "Focal_Genus")%>%
  subset(select = colSums(., na.rm = TRUE) > 0)%>%
  rename_with(~ sub("^RAI_", "", .))

#Calculate the Trophic overlap with Pianka index
TR_OH <- as.data.frame(as.matrix(as.dist(niche.overlap(TR_OH, method = "pianka"))))
diag(TR_OH) <- 1

TR_HU <- as.data.frame(as.matrix(as.dist(niche.overlap(TR_HU, method = "pianka"))))
diag(TR_HU) <- 1

TR_NH <- as.data.frame(as.matrix(as.dist(niche.overlap(TR_NH, method = "pianka"))))
diag(TR_NH) <- 1

#Reorder Trophic matrices
TR_OH <- reorder_matrix(TR_OH, species_order)
TR_HU <- reorder_matrix(TR_HU, species_order)
TR_NH <- reorder_matrix(TR_NH, species_order)

#Clean env.
rm(list = c("CT", "data_TR", "data_TR1", "data_TR2", "data_TR3", "data_TR4", ls(pattern = "^RAI_")))

##### Combined Overlap #####

#Calculate the index of trophic_spatio_temporal overlap
#It is the sum of pairwise Temporal + Spatial + Trophic overlaps
#It gives a value between 0 and 3, where 0 means total segregation and 3 total overlap
C_OH <- as.data.frame(as.matrix(T_OH + S_OH + TR_OH))
C_HU <- as.data.frame(as.matrix(T_HU + S_HU + TR_HU)) 
C_NH <- as.data.frame(as.matrix(T_NH + S_NH + TR_NH)) 

#Calculate the normalized contribution of each variable to the total. 
#(Niche ov. / Total ov.)
#Value between 0 and 1 for each overlap matrix, where 0 means total segregation and 1 total overlap.
#The variable with the highest value contributes the most to the total overlap (between 0 and 3)
C_S_OH <- as.data.frame(S_OH  / (C_OH))
C_T_OH  <- as.data.frame(T_OH  / (C_OH))
C_TR_OH <- as.data.frame(TR_OH / (C_OH))

C_S_HU  <- as.data.frame(S_HU  / (C_HU))
C_T_HU  <- as.data.frame(T_HU  / (C_HU))
C_TR_HU <- as.data.frame(TR_HU / (C_HU))

C_S_NH  <- as.data.frame(S_NH  / (C_NH))
C_T_NH  <- as.data.frame(T_NH  / (C_NH))
C_TR_NH <- as.data.frame(TR_NH / (C_NH))


##### Save Matrices #####
##Temporal
write.csv(T_OH, "Results/Group_Analysis/T_Overhunted.csv", row.names = TRUE)
write.csv(T_HU, "Results/Group_Analysis/T_Hunted.csv", row.names = TRUE)
write.csv(T_NH, "Results/Group_Analysis/T_Non_hunted.csv", row.names = TRUE)


##Spatial
write.csv(S_OH, "Results/Group_Analysis/S_Overhunted.csv", row.names = TRUE)
write.csv(S_HU, "Results/Group_Analysis/S_Hunted.csv", row.names = TRUE)
write.csv(S_NH, "Results/Group_Analysis/S_Non_hunted.csv", row.names = TRUE)


##Trophic
write.csv(TR_OH, "Results/Group_Analysis/TR_Overhunted.csv", row.names = TRUE)
write.csv(TR_HU, "Results/Group_Analysis/TR_Hunted.csv", row.names = TRUE)
write.csv(TR_NH, "Results/Group_Analysis/TR_Non_hunted.csv", row.names = TRUE)


##Combined
#Total overlap
write.csv(C_OH, "Results/Group_Analysis/Tot_Overhunted.csv", row.names = TRUE)
write.csv(C_HU, "Results/Group_Analysis/Tot_Hunted.csv", row.names = TRUE)
write.csv(C_NH, "Results/Group_Analysis/Tot_Non_hunted.csv", row.names = TRUE)

#Dimensions contribution
write.csv(C_S_OH, "Results/Group_Analysis/Contr_S_Overhunted.csv", row.names = TRUE)
write.csv(C_S_HU, "Results/Group_Analysis/Contr_S_Hunted.csv", row.names = TRUE)
write.csv(C_S_NH, "Results/Group_Analysis/Contr_S_Non_hunted.csv", row.names = TRUE)

write.csv(C_T_OH, "Results/Group_Analysis/Contr_T_Overhunted.csv", row.names = TRUE)
write.csv(C_T_HU, "Results/Group_Analysis/Contr_T_Hunted.csv", row.names = TRUE)
write.csv(C_T_NH, "Results/Group_Analysis/Contr_T_Non_hunted.csv", row.names = TRUE)

write.csv(C_TR_OH, "Results/Group_Analysis/Contr_TR_Overhunted.csv", row.names = TRUE)
write.csv(C_TR_HU, "Results/Group_Analysis/Contr_TR_Hunted.csv", row.names = TRUE)
write.csv(C_TR_NH, "Results/Group_Analysis/Contr_TR_Non_hunted.csv", row.names = TRUE)

#Save Combined objects as RDS as well
combined_objects <- list(
  C_OH = C_OH, C_HU = C_HU, C_NH = C_NH,
  C_S_OH = C_S_OH, C_T_OH = C_T_OH, C_TR_OH = C_TR_OH,
  C_S_HU = C_S_HU, C_T_HU = C_T_HU, C_TR_HU = C_TR_HU,
  C_S_NH = C_S_NH, C_T_NH = C_T_NH, C_TR_NH = C_TR_NH
)
saveRDS(combined_objects, "Results/Group_Analysis/combined_objects.rds")



# 02_Network metrics ####
#This script calculates network-level and node-level metrics (strength centrality, modularity,
#and nestedness) for each niche dimension (Temporal, Spatial, Trophic) and the combined overlap
#network, across hunting-pressure treatments (Non-hunted, Hunted, Overhunted). All metrics are
#tested against null models to assess statistical significance.
#Centrality data, module membership, and summary plots are exported for use in the network
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

## Function to store network metrics in a table
metrics_list <- list()

add_metrics <- function(network, hunting, Q_real, Q_null, N_real, N_null) {
  metrics_list[[length(metrics_list) + 1]] <<- data.frame(
    Network = network,
    Hunting_Pressure = hunting,
    Q_obs = Q_real,
    Q_z = (Q_real - mean(Q_null)) / sd(Q_null),
    Q_std = Q_real / max(Q_null),
    N_obs = N_real,
    N_z = (N_real - mean(N_null)) / sd(N_null),
    N_std = N_real / max(N_null)
  )
}

##### Temporal #####
##Upload the matrices
T_OH <- read.csv("Results/Group_Analysis/T_Overhunted.csv", row.names = 1)
T_HU <- read.csv("Results/Group_Analysis/T_Hunted.csv", row.names = 1)
T_NH <- read.csv("Results/Group_Analysis/T_Non_hunted.csv", row.names = 1)

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
    colours = c("#2166AC", "#67A9CF", "#F9FFAF", "#F4A582", "#B2182B"),
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
Z_score_OH #-1.64

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0

##Nestedness
mat <- as.matrix(T_OH)
diag(mat) <- 0
n_sim <- 1000

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

nested_real_OH
Z_score_OH

N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH

add_metrics("Temporal", "OH", Q_real_OH, Q_null_OH, nested_real_OH, nested_null_OH)

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
Z_score_HU #-0.86

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0

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

nested_real_HU
Z_score_HU

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU

add_metrics("Temporal", "HU", Q_real_HU, Q_null_HU, nested_real_HU, nested_null_HU)

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
Z_score_NH #-2.39

Q_ratio_NH <- Q_real_NH / max(Q_null_NH)
Q_ratio_NH #0


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

nested_real_NH
Z_score_NH

N_ratio_NH <- nested_real_NH / max(nested_null_NH)
N_ratio_NH

add_metrics("Temporal", "NH", Q_real_NH, Q_null_NH, nested_real_NH, nested_null_NH)


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

#Combined plot
T_met <- ggarrange(mod, nest, nrow = 1)
T_met


##### Spatial #####
##Calculate Modularity
S_OH <- read.csv("Results/Group_Analysis/S_Overhunted.csv", row.names = 1)
S_HU <- read.csv("Results/Group_Analysis/S_Hunted.csv", row.names = 1)
S_NH <- read.csv("Results/Group_Analysis/S_Non_hunted.csv", row.names = 1)


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
    colours = c("#2166AC", "#67A9CF", "#F9FFAF", "#F4A582", "#B2182B"),
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

Q_real_OH #0.15
Z_score_OH #0.89

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0.63


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

nested_real_OH
Z_score_OH


N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH

add_metrics("Spatial", "OH", Q_real_OH, Q_null_OH, nested_real_OH, nested_null_OH)

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

Q_real_HU #0.05
Z_score_HU #-1.04

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0.17

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

nested_real_HU
Z_score_HU

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU

add_metrics("Spatial", "HU", Q_real_HU, Q_null_HU, nested_real_HU, nested_null_HU)

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
Z_score_NH #3.78
membership(comm_real_over)

#Save the membership for the figure
S_membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(S_membership_NH, "Results/Group_Analysis/S_membership_NH.rds")

Q_ratio_NH <- Q_real_NH / max(Q_null_NH)
Q_ratio_NH #1.06

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

nested_real_NH
Z_score_NH

N_ratio_NH <- nested_real_NH / max(nested_null_NH)
N_ratio_NH

add_metrics("Spatial", "NH", Q_real_NH, Q_null_NH, nested_real_NH, nested_null_NH)

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
TR_OH <- read.csv("Results/Group_Analysis/TR_Overhunted.csv", row.names = 1)
TR_HU <- read.csv("Results/Group_Analysis/TR_Hunted.csv", row.names = 1)
TR_NH <- read.csv("Results/Group_Analysis/TR_Non_hunted.csv", row.names = 1)

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
    colours = c("#2166AC", "#67A9CF", "#F9FFAF", "#F4A582", "#B2182B"),
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

Q_real_OH #0.07
Z_score_OH #0.88
membership(comm_real_over)

Q_ratio_OH <- Q_real_OH / max(Q_null_OH)
Q_ratio_OH #0.51

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

nested_real_OH
Z_score_OH

N_ratio_OH <- nested_real_OH / max(nested_null_OH)
N_ratio_OH

add_metrics("Trophic", "OH", Q_real_OH, Q_null_OH, nested_real_OH, nested_null_OH)

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

Q_real_HU #0.04
Z_score_HU #1.01
membership(comm_real_over)

Q_ratio_HU <- Q_real_HU / max(Q_null_HU)
Q_ratio_HU #0.55

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

nested_real_HU
Z_score_HU

N_ratio_HU <- nested_real_HU / max(nested_null_HU)
N_ratio_HU

add_metrics("Trophic", "HU", Q_real_HU, Q_null_HU, nested_real_HU, nested_null_HU)

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

Q_real_NH #0.14
Z_score_NH #3.4
membership(comm_real_over)

#Save the membership for the figure
membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(membership_NH, "Results/Group_Analysis/TR_membership_NH.rds")

Q_ratio_Non_hunted <- Q_real_NH/ max(Q_null_NH)
Q_ratio_Non_hunted #0.93


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

nested_real_NH
Z_score_NH

nested_ratio_NH <- nested_real_NH/ max(nested_null_NH)
nested_ratio_NH

add_metrics("Trophic", "NH", Q_real_NH, Q_null_NH, nested_real_NH, nested_null_NH)

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
C_OH <- as.matrix(read.csv("Results/Group_Analysis/Tot_Overhunted.csv", row.names = 1))
C_HU <- as.matrix(read.csv("Results/Group_Analysis/Tot_Hunted.csv", row.names = 1))
C_NH <- as.matrix(read.csv("Results/Group_Analysis/Tot_Non_hunted.csv", row.names = 1))

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
    colours = c("#2166AC", "#67A9CF", "#F9FFAF", "#F4A582", "#B2182B"),
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

Q_real_OH # 0.02
Z_score_OH # -0.21

Q_ratio_OH <- Q_real_OH/ max(Q_null_OH)
Q_ratio_OH #0.2

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

nested_real_OH
Z_score_OH

nested_ratio_OH <- nested_real_OH/ max(nested_null_OH)
nested_ratio_OH

add_metrics("Combined", "OH", Q_real_OH, Q_null_OH, nested_real_OH, nested_null_OH)

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

Q_real_HU # 0
Z_score_HU # -0.22

Q_ratio_HU <- Q_real_HU/ max(Q_null_HU)
Q_ratio_HU #0

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

nested_real_HU
Z_score_HU

nested_ratio_HU <- nested_real_HU/ max(nested_null_HU)
nested_ratio_HU

add_metrics("Combined", "HU", Q_real_HU, Q_null_HU, nested_real_HU, nested_null_HU)

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

Q_real_NH  # 0.07
Z_score_NH # 2.79
membership(comm_real_over)

#Save the membership for the figure
C_membership_NH <- data.frame(
  Genus  = names(membership(comm_real_over)),
  module = as.factor(membership(comm_real_over))
)
saveRDS(C_membership_NH, "Results/Group_Analysis/C_membership_NH.rds")

Q_ratio_NH <- Q_real_NH/ max(Q_null_NH)
Q_ratio_NH #0.85

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

nested_real_NH
Z_score_NH

nested_ratio_NH <- nested_real_NH/ max(nested_null_NH)
nested_ratio_NH

add_metrics("Combined", "NH", Q_real_NH, Q_null_NH, nested_real_NH, nested_null_NH)

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


#### Save figures and matrices ####
## Centrality
#Matrices
write.csv(T_strength_null, "Results/Group_Analysis/T_centrality.csv", row.names = TRUE)
write.csv(S_strength_null, "Results/Group_Analysis/S_centrality.csv", row.names = TRUE)
write.csv(TR_strength_null, "Results/Group_Analysis/TR_centrality.csv", row.names = TRUE)

#Plots
centr_plot <- ggarrange(T_centr_plot, S_centr_plot, TR_centr_plot, C_centr_plot, ncol = 1, common.legend = T)

ggsave("Figures/Group_analysis/Centrality.jpg", centr_plot, width = 11, height = 14)
ggsave("Figures/Group_analysis/Centrality.pdf", centr_plot, width = 11, height = 14)

## Modularity + Nestedness
#Temporal
ggsave("Figures/Group_analysis/Metrics_Temporal.pdf", T_met, width = 6, height = 6)

#Spatial
ggsave("Figures/Group_analysis/Metrics_Spatial.pdf", S_met, width = 6, height = 6)

#Trophic
ggsave("Figures/Group_analysis/Metrics_Trophic.pdf", Tr_met, width = 6, height = 6)

#Combined
ggsave("Figures/Group_analysis/Metrics_Combined.pdf", C_met, width = 6, height = 6)


#### Results tables ####
# Formatting helpers: values rounded to 2 decimals, exact zeros shown as "0", z-scores with sign
fmt_val <- function(x) ifelse(round(x, 2) == 0, "0", sprintf("%.2f", x))
fmt_obs_z <- function(obs, z) paste0(fmt_val(obs), " (", sprintf("%+.2f", z), ")")

## Table 1: network metrics (modularity and nestedness)
metrics_raw <- bind_rows(metrics_list) %>%
  mutate(
    Network = factor(Network, levels = c("Temporal", "Spatial", "Trophic", "Combined")),
    Hunting_Pressure = factor(Hunting_Pressure, levels = c("NH", "HU", "OH"))
  ) %>%
  arrange(Network, Hunting_Pressure)

metrics_table <- metrics_raw %>%
  transmute(
    Network = as.character(Network),
    `Hunting Pressure` = as.character(Hunting_Pressure),
    `Observed modularity (z-score)` = fmt_obs_z(Q_obs, Q_z),
    `Standardized modularity` = fmt_val(Q_std),
    `Observed nestedness (z-score)` = fmt_obs_z(N_obs, N_z),
    `Standardized nestedness` = fmt_val(N_std)
  )

# Show the network name only in the first row of each block (as in the manuscript table)
metrics_table$Network[duplicated(metrics_table$Network)] <- ""

metrics_table

## Table 2: strength centrality
C_strength_tab <- C_strength_null
C_strength_tab$Dimension <- "Combined"

centrality_raw <- bind_rows(T_strength_null, S_strength_null, TR_strength_null, C_strength_tab) %>%
  mutate(
    Dimension = factor(Dimension, levels = c("Temporal", "Spatial", "Trophic", "Combined")),
    Hunting_Pressure = factor(as.character(Hunting_Pressure),
                              levels = c("Non_hunted", "Hunted", "Overhunted"),
                              labels = c("NH", "HU", "OH"))
  ) %>%
  arrange(Dimension, Hunting_Pressure, desc(Strength_norm))

centrality_table <- centrality_raw %>%
  transmute(
    Network = as.character(Dimension),
    `Hunting Pressure` = as.character(Hunting_Pressure),
    Genus = Genus,
    `Observed strength` = fmt_val(Strength_obs),
    `Normalized strength` = fmt_val(Strength_norm),
    `z-score` = sprintf("%+.2f", Strength_z)
  )

centrality_table

## Save tables
# Excel: formatted tables + raw numeric values
wb <- createWorkbook()
addWorksheet(wb, "Network_metrics")
writeData(wb, "Network_metrics", metrics_table)
addWorksheet(wb, "Network_metrics_raw")
writeData(wb, "Network_metrics_raw", metrics_raw)
addWorksheet(wb, "Centrality")
writeData(wb, "Centrality", centrality_table)
addWorksheet(wb, "Centrality_raw")
writeData(wb, "Centrality_raw", centrality_raw)
saveWorkbook(wb, "Results/Group_Analysis/Network_tables_group.xlsx", overwrite = TRUE)

# CSV
write.csv(metrics_raw, "Results/Group_Analysis/Network_metrics_group.csv", row.names = FALSE)
write.csv(centrality_raw, "Results/Group_Analysis/Centrality_group.csv", row.names = FALSE)

# Word
doc <- read_docx() %>%
  body_add_par("Network metrics (group-level detections)", style = "heading 1") %>%
  body_add_table(metrics_table, style = "table_template") %>%
  body_add_par("Strength centrality (group-level detections)", style = "heading 1") %>%
  body_add_table(centrality_table, style = "table_template")
print(doc, target = "Results/Group_Analysis/Network_tables_group.docx")


# 03_Networks plot ####
#Network and ternary plots
#This script builds unipartite overlap networks for each niche dimension (Temporal, Spatial, Trophic)
#and hunting-pressure treatment (Non-hunted, Hunted, Overhunted), plus the combined overlap network
#(with module hulls where modularity is significant) and ternary plots showing the relative
#contribution of each dimension to the combined overlap between primate pairs.
#All figures are exported as PDF at the end of the script.

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
TR_centr <- read.csv("Results/Group_Analysis/TR_centrality.csv")
T_centr <- read.csv("Results/Group_Analysis/T_centrality.csv")
S_centr <- read.csv("Results/Group_Analysis/S_centrality.csv")
centr <- rbind(S_centr,T_centr,TR_centr)

global_limits <- range(centr$Strength_norm, na.rm = TRUE)
global_limits

#### Trophic ####
#Upload Trophic overlap matrices
TR_OH <- read.csv("Results/Group_Analysis/TR_Overhunted.csv", row.names = 1)
TR_HU <- read.csv("Results/Group_Analysis/TR_Hunted.csv", row.names = 1)
TR_NH <- read.csv("Results/Group_Analysis/TR_Non_hunted.csv", row.names = 1)

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
membership_NH <- readRDS("Results/Group_Analysis/TR_membership_NH.rds")
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
T_OH <- read.csv("Results/Group_Analysis/T_Overhunted.csv", row.names = 1)
T_HU <- read.csv("Results/Group_Analysis/T_Hunted.csv", row.names = 1)
T_NH <- read.csv("Results/Group_Analysis/T_Non_hunted.csv", row.names = 1)

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
S_OH <- read.csv("Results/Group_Analysis/S_Overhunted.csv", row.names = 1)
S_HU <- read.csv("Results/Group_Analysis/S_Hunted.csv", row.names = 1)
S_NH <- read.csv("Results/Group_Analysis/S_Non_hunted.csv", row.names = 1)

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
  left_join(df_NH %>% select(Genus, Strength_norm), by = c("name" = "Genus"))

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
S_membership_NH <- readRDS("Results/Group_Analysis/S_membership_NH.rds")
nodes_df <- nodes_df %>%
  left_join(df_NH %>% select(Genus, Strength_norm), by = c("name" = "Genus")) %>%
  left_join(S_membership_NH, by = c("name" = "Genus"))

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
combined_objects <- readRDS("Results/Group_Analysis/combined_objects.rds")
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
C_membership_NH <- readRDS("Results/Group_Analysis/C_membership_NH.rds")
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
  ggtitle("Hunted - Trophic") +
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
  ggtitle("Hunted - Trophic") +
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
  ggtitle("Hunted - Trophic") +
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
ggarrange(Temporal1, Spatial1, Trophic1,
          nrow = 4 , legend = "none")

all <- ggarrange(Temporal1, Spatial1, Trophic1,
                 nrow = 3, common.legend = T)
all

ggsave("Figures/Group_analysis/Spatial_networks.pdf", Spatial1, width = 14, height = 6)
ggsave("Figures/Group_analysis/Trophic_networks.pdf", Trophic1, width = 14, height = 6)
ggsave("Figures/Group_analysis/Temporal_networks.pdf", Temporal1, width = 14, height = 6)
ggsave("Figures/Group_analysis/Ternary_plots.pdf", tern, width = 14, height = 6)
ggsave("Figures/Group_analysis/Combined_networks.pdf", net, width = 14, height = 6)
ggsave("Figures/Group_analysis/All_networks.pdf", all, width = 14, height = 14)
ggsave("Figures/Group_analysis/comb_networks.pdf", comb, width = 14, height = 11)

