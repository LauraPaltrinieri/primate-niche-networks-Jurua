#Niche Overlap
#Author: Laura Paltrinieri
#16/04/2026
#Modified: 07/09/2026

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
data <- read.csv("Data/Primates_detections.csv")

#Filter for minimum 10 detections for primate genus in each site
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


library(overlap)
library(ggforce)

ordine_primati <- c("Ateles", "Lagothrix", "Alouatta", "Sapajus", "Cacajao", "Cebus",
                    "Pithecia", "Callicebus", "Aotus", "Saimiri", "Saguinus")

grid_T <- seq(0, 2 * pi, length.out = 97)[-97]
n_risorse_T <- length(grid_T)
min_n <- 10

levins_T <- data_T %>%
  group_by(Primate = Vertebrate_ID, Hunting_Pressure) %>%
  summarise(
    n_det = n(),
    B = if (n() >= min_n) {
      x <- as.numeric(solar_circ)
      d <- densityFit(x, grid_T, getBandWidth(x, kmax = 3))
      p <- d / sum(d)
      1 / sum(p^2)
    } else NA_real_,
    .groups = "drop") %>%
  mutate(BA = (B - 1) / (n_risorse_T - 1),
         Hunting_Pressure = recode(Hunting_Pressure,
                                   "Non_hunted" = "NH", "Hunted" = "HU", "Overhunted" = "OH"),
         Primate = factor(Primate, levels = ordine_primati)) %>%
  arrange(Primate, Hunting_Pressure)

print(levins_T, n = Inf)
cols_hp <- c("NH" = "#85a4c1", "HU" = "#ffd54e", "OH" = "#9e1919")
labs_hp <- c("NH" = "Non-hunted", "HU" = "Hunted", "OH" = "Overhunted")
livelli <- c(ordine_primati, "BA scale")

circ_T <- levins_T %>%
  filter(!is.na(BA)) %>%
  arrange(desc(BA)) %>%
  mutate(id = factor(row_number(), levels = row_number()),
         r = 0.1 + 0.9 * BA,
         Primate = factor(as.character(Primate), levels = livelli)) %>%
  tidyr::crossing(t = seq(0, 2 * pi, length.out = 200)) %>%
  mutate(x = r * cos(t), y = r * sin(t))

ref_BA <- c(0.25, 0.5, 0.75, 1)

ref_circ <- tidyr::crossing(BA = ref_BA, t = seq(0, 2 * pi, length.out = 200)) %>%
  mutate(r = 0.1 + 0.9 * BA, x = r * cos(t), y = r * sin(t),
         Primate = factor("BA scale", levels = livelli))

ref_lab <- data.frame(BA = ref_BA) %>%
  mutate(r = 0.1 + 0.9 * BA, x = r * cos(pi / 4), y = r * sin(pi / 4),
         Primate = factor("BA scale", levels = livelli))

lab_gen <- data.frame(Primate = factor(livelli, levels = livelli),
                      x = 0,
                      y = ifelse(livelli == "BA scale", 1.1, 0.95),
                      face = ifelse(livelli == "BA scale", "plain", "italic"))

fig_levins_T <- ggplot(circ_T, aes(x = x, y = y, group = id,
                                   colour = Hunting_Pressure, fill = Hunting_Pressure)) +
  geom_polygon(linewidth = .6) +
  geom_polygon(data = ref_circ, aes(x = x, y = y, group = BA), inherit.aes = FALSE,
               fill = NA, colour = "grey50", linetype = "dashed", linewidth = .4) +
  geom_label(data = ref_lab, aes(x = x, y = y, label = BA), inherit.aes = FALSE,
             size = 2.8, label.size = 0, fill = "white", label.padding = unit(0.1, "lines")) +
  geom_text(data = lab_gen, aes(x = x, y = y, label = Primate, fontface = face),
            inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~ Primate, ncol = 4, drop = FALSE) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") +
  scale_colour_manual(values = cols_hp, labels = labs_hp, breaks = c("NH", "HU", "OH"), name = "Hunting pressure") +
  scale_fill_manual(values = scales::alpha(cols_hp, 0.3), breaks = c("NH", "HU", "OH"), labels = labs_hp, name = "Hunting pressure") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "bottom")

fig_levins_T

ggsave("Figures/Fig_tempor_niche_breadth.png", fig_levins_T, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("Figures/Fig_tempor_niche_breadth.pdf", fig_levins_T, width = 8, height = 7, dpi = 300, bg = "white")


#Clean env.
rm(hunting_pressure, hp, overlap_matrices, solartime, data_T, hp_data, mat)

##### Spatial Overlap #####
#Spatial data
data_S <- data

#Count the occurrences of each primate genus in each Station_ID
data_S1 <- data_S %>%
  group_by(Tree_ID, Hunting_Pressure, Genus, Effort_H_Station, Focal_Genus, Sampled_Canopy) %>% 
  summarise(count = n(), .groups = "drop") %>%
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

#Average RAI across trees to obtain the mean foraging rate of each primate genus 
#per plant genus and hunting-pressure category
data_S4 <- data_S3 %>%
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


write.csv(data_TR4, "Results/forag_rates.csv", row.names = FALSE)

#Create Table S6 of supporting information
format_rate <- function(x) {
  ifelse(x == 0, "0",
    ifelse(
      abs(x) >= 0.01,
      formatC(x, format = "f", digits = 2),
      formatC(x, format = "g", digits = 2)))
}

Table_S6 <- data_TR4 %>%
  pivot_longer(cols = starts_with("RAI_"),
    names_to = "Primate_Genus",
    names_prefix = "RAI_",
    values_to = "Mean_foraging_rate") %>%
  mutate(Hunting_Pressure = recode(Hunting_Pressure,
                              "Non_hunted" = "NH",
                              "Hunted" = "HU",
                              "Overhunted" = "OH"),
    Mean_foraging_rate = format_rate(Mean_foraging_rate)) %>%
  pivot_wider(
    id_cols = c(Primate_Genus, Hunting_Pressure),
    names_from = Focal_Genus,
    values_from = Mean_foraging_rate) %>%
  select(
    Primate = Primate_Genus,
    `Hunting pressure` = Hunting_Pressure,
    `Pourouma spp.` = Pouroma,
    `Euterpe spp.` = Euterpe,
    `Virola spp.` = Virola,
    `Iryanthera spp.` = Iryanthera,
    `Pouteria spp.` = Pouteria) %>%
  arrange(Primate, factor(`Hunting pressure`, levels = c("NH", "HU", "OH")))

Table_S6


# Levins' niche breadth (standardized, Hurlbert 1978)
tu <- data_TR4 %>%
  pivot_longer(
    cols = starts_with("RAI_"),
    names_to = "Primate",
    names_prefix = "RAI_",
    values_to = "TU") %>%
  mutate(Hunting_Pressure = recode(Hunting_Pressure,
                              "Non_hunted" = "NH",
                              "Hunted" = "HU",
                              "Overhunted" = "OH"),
    Hunting_Pressure = factor(Hunting_Pressure, levels = c("NH", "HU", "OH")))

#Check
#Hunting pressure
table(tu$Hunting_Pressure, useNA = "ifany")

#Focal genus
n_risorse <- n_distinct(tu$Focal_Genus) 

#Calculate niche breadth
levins <- tu %>%
  group_by(Primate, Hunting_Pressure) %>%
  summarise(
    n_piante_usate = sum(TU > 0, na.rm = TRUE),
    B  = ifelse(sum(TU, na.rm = TRUE) > 0,
                1 / sum((TU / sum(TU, na.rm = TRUE))^2, na.rm = TRUE),NA),
    BA = (B - 1) / (n_risorse - 1),.groups = "drop")

# Order the primates
ordine_primati <- c("Ateles", "Lagothrix", "Alouatta", "Sapajus", "Cacajao", "Cebus",
                    "Pithecia", "Callicebus", "Aotus", "Saimiri", "Saguinus")

levins <- levins %>%
  mutate(Primate = factor(Primate, levels = ordine_primati)) %>%
  arrange(Primate, Hunting_Pressure)

print(levins, n = Inf)

#Plot
cols_hp <- c("NH" = "#85a4c1", "HU" = "#ffd54e", "OH" = "#9e1919")
labs_hp <- c("NH" = "Non-hunted", "HU" = "Hunted", "OH" = "Overhunted")

circ_TR <- levins %>%
  filter(!is.na(BA)) %>%
  arrange(desc(BA)) %>%
  mutate(id = factor(row_number(), levels = row_number()),
         r = 0.1 + 0.9 * BA) %>%
  tidyr::crossing(t = seq(0, 2 * pi, length.out = 200)) %>%
  mutate(x = r * cos(t), y = r * sin(t))

lab_gen <- data.frame(Primate = factor(livelli, levels = livelli),
                      x = 0,
                      y = ifelse(livelli == "BA scale", 1.1, 0.95),
                      face = ifelse(livelli == "BA scale", "plain", "italic"))

fig_levins <- ggplot(circ_TR, aes(x = x, y = y, group = id,
                                  colour = Hunting_Pressure, fill = Hunting_Pressure)) +
  geom_polygon(linewidth = .6) +
  geom_polygon(data = ref_circ, aes(x = x, y = y, group = BA), inherit.aes = FALSE,
               fill = NA, colour = "grey50", linetype = "dashed", linewidth = .4) +
  geom_label(data = ref_lab, aes(x = x, y = y, label = BA), inherit.aes = FALSE,
             size = 2.8, label.size = 0, fill = "white", label.padding = unit(0.1, "lines")) +
  geom_text(data = lab_gen, aes(x = x, y = y, label = Primate, fontface = face),
            inherit.aes = FALSE, size = 3.5) +
  facet_wrap(~ Primate, ncol = 4, drop = FALSE) +
  coord_fixed(xlim = c(-1, 1), ylim = c(-1, 1), clip = "off") +
  scale_colour_manual(values = cols_hp, labels = labs_hp, breaks = c("NH", "HU", "OH"), name = "Hunting pressure") +
  scale_fill_manual(values = scales::alpha(cols_hp, 0.3), breaks = c("NH", "HU", "OH"), labels = labs_hp, name = "Hunting pressure") +
  theme_void() +
  theme(strip.text = element_blank(),
        legend.position = "bottom")

fig_levins

ggsave("Figures/Fig_trophic_niche_breadth.png", fig_levins, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("Figures/Fig_trophic_niche_breadth.pdf", fig_levins, width = 8, height = 7, dpi = 300, bg = "white")


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

mean_overlap <- function(m) {
  m <- as.matrix(m)
  v <- m[lower.tri(m)]
  c(mean = round(mean(v), 3), sd = round(sd(v), 3), n_pairs = length(v))
}

rbind(
  Temporal_NH = mean_overlap(T_NH), Temporal_HU = mean_overlap(T_HU), Temporal_OH = mean_overlap(T_OH),
  Spatial_NH  = mean_overlap(S_NH), Spatial_HU  = mean_overlap(S_HU), Spatial_OH  = mean_overlap(S_OH),
  Trophic_NH  = mean_overlap(TR_NH), Trophic_HU = mean_overlap(TR_HU), Trophic_OH = mean_overlap(TR_OH)
)


##### Save Matrices #####
##Temporal
write.csv(T_OH, "Results/T_Overhunted.csv", row.names = TRUE)
write.csv(T_HU, "Results/T_Hunted.csv", row.names = TRUE)
write.csv(T_NH, "Results/T_Non_hunted.csv", row.names = TRUE)


##Spatial
write.csv(S_OH, "Results/S_Overhunted.csv", row.names = TRUE)
write.csv(S_HU, "Results/S_Hunted.csv", row.names = TRUE)
write.csv(S_NH, "Results/S_Non_hunted.csv", row.names = TRUE)


##Trophic
write.csv(TR_OH, "Results/TR_Overhunted.csv", row.names = TRUE)
write.csv(TR_HU, "Results/TR_Hunted.csv", row.names = TRUE)
write.csv(TR_NH, "Results/TR_Non_hunted.csv", row.names = TRUE)


##Combined
#Total overlap
write.csv(C_OH, "Results/Tot_Overhunted.csv", row.names = TRUE)
write.csv(C_HU, "Results/Tot_Hunted.csv", row.names = TRUE)
write.csv(C_NH, "Results/Tot_Non_hunted.csv", row.names = TRUE)

#Dimensions contribution
write.csv(C_S_OH, "Results/Contr_S_Overhunted.csv", row.names = TRUE)
write.csv(C_S_HU, "Results/Contr_S_Hunted.csv", row.names = TRUE)
write.csv(C_S_NH, "Results/Contr_S_Non_hunted.csv", row.names = TRUE)

write.csv(C_T_OH, "Results/Contr_T_Overhunted.csv", row.names = TRUE)
write.csv(C_T_HU, "Results/Contr_T_Hunted.csv", row.names = TRUE)
write.csv(C_T_NH, "Results/Contr_T_Non_hunted.csv", row.names = TRUE)

write.csv(C_TR_OH, "Results/Contr_TR_Overhunted.csv", row.names = TRUE)
write.csv(C_TR_HU, "Results/Contr_TR_Hunted.csv", row.names = TRUE)
write.csv(C_TR_NH, "Results/Contr_TR_Non_hunted.csv", row.names = TRUE)

#Save Combined objects as RDS as well
combined_objects <- list(
  C_OH = C_OH, C_HU = C_HU, C_NH = C_NH,
  C_S_OH = C_S_OH, C_T_OH = C_T_OH, C_TR_OH = C_TR_OH,
  C_S_HU = C_S_HU, C_T_HU = C_T_HU, C_TR_HU = C_TR_HU,
  C_S_NH = C_S_NH, C_T_NH = C_T_NH, C_TR_NH = C_TR_NH
)
saveRDS(combined_objects, "Results/combined_objects.rds")

