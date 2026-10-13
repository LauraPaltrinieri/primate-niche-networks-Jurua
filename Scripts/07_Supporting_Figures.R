####Supporting Information 

### Figure S1 - Behaviours ####
library(dplyr)
library(ggplot2)
#clean_env
rm(list = ls())

# Upload the data
data <- read.csv("Data/Primates_detections.csv")

#Filter for minimum 10 detections for primate genus in each site
data <- data %>%
  add_count(Vertebrate_ID, Hunting_Pressure, name = "n_detections") %>%
  filter(n_detections >= 10)

#Reclassify "Several" as "Foraging" when fruit consumption was confirmed (Confidence_Interaction = 5)
data <- data %>%
  mutate(Behaviour = ifelse(Behaviour == "Several" & Confidence_Interaction %in% 5,
                            "Foraging", Behaviour))

#Order of primate genera (by body size) and behaviours
primate_order <- c("Ateles", "Lagothrix", "Alouatta", "Sapajus", "Cacajao", "Cebus",
                   "Pithecia", "Callicebus", "Aotus", "Saimiri", "Saguinus")

behaviour_order <- c("Foraging", "Looking for food", "Several", "Walking",
                     "Interaction with camera", "Other behaviour")

behaviour_palette <- setNames(hcl.colors(length(behaviour_order), "Heat"), behaviour_order)

#Percentage of each behaviour per genus and hunting pressure
beh_hp <- data %>%
  count(Hunting_Pressure, Genus, Behaviour) %>%
  group_by(Hunting_Pressure, Genus) %>%
  mutate(perc = 100 * n / sum(n), n_tot = sum(n)) %>%
  ungroup() %>%
  mutate(Genus = factor(Genus, levels = primate_order),
         Behaviour = factor(Behaviour, levels = rev(behaviour_order)),
         Hunting_Pressure = factor(Hunting_Pressure,
                                   levels = c("Non_hunted", "Hunted", "Overhunted"),
                                   labels = c("Non-hunted", "Hunted", "Overhunted")))

#Number of detections per genus and hunting pressure
n_lab_hp <- beh_hp %>% distinct(Hunting_Pressure, Genus, n_tot)

#Plot
fig_beh_hp <- ggplot(beh_hp, aes(x = Genus, y = perc, fill = Behaviour)) +
  geom_col(width = 0.75) +
  geom_text(data = n_lab_hp, aes(x = Genus, y = 102, label = paste0("n = ", n_tot)),
            inherit.aes = FALSE, size = 2.8, vjust = 0) +
  facet_wrap(~ Hunting_Pressure, ncol = 1) +
  scale_fill_manual(values = behaviour_palette, breaks = behaviour_order) +
  scale_y_continuous(limits = c(0, 112), breaks = seq(0, 100, 25), expand = c(0, 0)) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = "Detections (%)", fill = "Behaviour") +
  theme_bw() +
  theme(axis.text.x = element_text(face = "italic", angle = 45, hjust = 1))

fig_beh_hp

ggsave("Figures/Barplot_beh.png", fig_beh_hp, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("Figures/Barplot_beh.pdf", fig_beh_hp, width = 8, height = 7, dpi = 300, bg = "white")



### Figure S2 - GLMM ####
library(dplyr)
library(glmmTMB)
library(DHARMa)
library(performance)
library(sjPlot)
library(ggplot2)
library(patchwork)
det <- read.csv("Data/data_glmm.csv", sep = ",")
#Euterpe
m_EUT <- glmmTMB(detected ~ log_BM + log_GS + height + log_effort + canopy +
                   (1 | Site) + (1 | Genus) ,
                 family = binomial, data = filter(det, Focal_Genus == "Euterpe"))
summary(m_EUT)

#Iryanthera
m_IRY <- glmmTMB(detected ~ log_BM + log_GS + height + log_effort + canopy +
                   (1 | Site) + (1 | Genus) ,
                 family = binomial, data = filter(det, Focal_Genus == "Iryanthera"))
summary(m_IRY)

#Pourouma
m_POUR <- glmmTMB(detected ~ log_BM + log_GS + height + log_effort + canopy +
                    (1 | Site) + (1 | Genus) ,
                  family = binomial, data = filter(det, Focal_Genus == "Pouroma"))
summary(m_POUR)

#Pouteria
m_POUT <- glmmTMB(detected ~ log_BM + log_GS + height + log_effort + canopy +
                    (1 | Site) + (1 | Genus) ,
                  family = binomial, data = filter(det, Focal_Genus == "Pouteria"))
summary(m_POUT)

#Virola
m_VIR <- glmmTMB(detected ~ log_BM + log_GS + height + log_effort + canopy +
                   (1 | Site) + (1 | Genus) ,
                 family = binomial, data = filter(det, Focal_Genus == "Virola"))
summary(m_VIR)


performance::r2_nakagawa(m_VIR)
performance::r2_nakagawa(m_EUT)
performance::r2_nakagawa(m_IRY)
performance::r2_nakagawa(m_POUR)
performance::r2_nakagawa(m_POUT)

models <- list(Pourouma   = m_POUR,
               Euterpe    = m_EUT,
               Virola     = m_VIR,
               Iryanthera = m_IRY,
               Pouteria   = m_POUT)

plant_order   <- names(models)
plant_colours <- c(Pourouma = "darkolivegreen", Euterpe = "mediumpurple4", Virola = "sienna2",
                   Iryanthera = "firebrick3", Pouteria = "goldenrod")

pred_labels <- c(log_BM     = "Body mass (log)",
                 log_GS     = "Group size (log)",
                 height     = "Camera height",
                 log_effort = "Sampling effort (log)",
                 canopy     = "Proportion of canopy recorded")

#Odds ratios with 95% CI for each model
coef_df <- bind_rows(lapply(plant_order, function(pg) {
  cf <- as.data.frame(summary(models[[pg]])$coefficients$cond)
  data.frame(Predictor = rownames(cf), Estimate = cf$Estimate, SE = cf$`Std. Error`) %>%
    filter(Predictor != "(Intercept)") %>%
    mutate(Focal_Genus = pg)
})) %>%
  mutate(OR    = exp(Estimate),
         lower = exp(Estimate - 1.96 * SE),
         upper = exp(Estimate + 1.96 * SE),
         Predictor   = factor(pred_labels[Predictor], levels = rev(pred_labels)),
         Focal_Genus = factor(Focal_Genus, levels = rev(plant_order)))

mod <- ggplot(coef_df, aes(x = OR, y = Predictor, colour = Focal_Genus)) +
  geom_hline(yintercept = seq(1.5, 4.5, by = 1), colour = "black", linewidth = 0.3, linetype = "dotted") +
  geom_vline(xintercept = 1, colour = "black") +
  geom_errorbar(aes(xmin = lower, xmax = upper), width = 0, linewidth = 0.7,
                position = position_dodge(width = 0.7), orientation = "y") +
  geom_point(size = 2.5, position = position_dodge(width = 0.7)) +
  scale_x_log10(breaks = c(0.01, 0.1, 1, 10, 100),
                labels = c("0.01", "0.1", "1", "10", "100")) +
  coord_cartesian(xlim = c(0.01, 100)) +
  scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 18)) +
  scale_colour_manual(values = plant_colours, breaks = plant_order) +
  labs(x = "Odds ratio (95% CI)", y = NULL, colour = "Focal plant genus") + 
  guides(colour = guide_legend(nrow = 1, title = NULL)) +
  theme_classic(base_size = 12) +
  theme(axis.text       = element_text(colour = "black"),
        axis.text.y     = element_text(size = 11, lineheight = 0.9),
        legend.position = "top",
        legend.title    = element_text(face = "bold"),
        legend.text     = element_text(face = "italic"))

r2_tab <- bind_rows(lapply(plant_order, function(pg) {
  r2 <- r2_nakagawa(models[[pg]])
  data.frame(Focal_Genus    = pg,
             R2_marginal    = r2$R2_marginal,
             R2_conditional = r2$R2_conditional)
})) %>%
  mutate(Focal_Genus = factor(Focal_Genus, levels = plant_order))

r2_long <- r2_tab %>%
  transmute(Focal_Genus,
            `Fixed effects (marginal R²)` = R2_marginal,
            `Random effects`              = R2_conditional - R2_marginal) %>%
  pivot_longer(-Focal_Genus, names_to = "Component", values_to = "R2") %>%
  mutate(Component = factor(Component, levels = c("Random effects", "Fixed effects (marginal R²)")))

r2 <- ggplot(r2_long, aes(x = Focal_Genus, y = R2, fill = Focal_Genus, alpha = Component)) +
  geom_col(width = 0.6) +
  geom_text(data = r2_tab, aes(x = Focal_Genus, y = R2_conditional,
                               label = sprintf("%.2f", R2_conditional)),
            inherit.aes = FALSE, vjust = -0.4, size = 3.5) +
  geom_text(data = r2_tab, aes(x = Focal_Genus, y = R2_marginal / 2,
                               label = sprintf("%.2f", R2_marginal)),
            inherit.aes = FALSE, colour = "white", size = 3.2) +
  scale_fill_manual(values = plant_colours, guide = "none") +
  scale_alpha_manual(values = c("Fixed effects (marginal R²)" = 1, "Random effects" = 0.35),
                     breaks = c("Fixed effects (marginal R²)", "Random effects"),
                     name = NULL) +
  guides(alpha = guide_legend(override.aes = list(fill = "grey30"))) +
  scale_y_continuous(limits = c(0, 0.7), expand = c(0, 0)) +
  labs(x = NULL, y = expression(R^2)) +
  theme_classic(base_size = 12) +
  theme(axis.text       = element_text(colour = "black"),
        axis.text.x     = element_text(face = "italic"),
        legend.position = "top") + 
  labs(x = "Focal plant genus")


fig_glmm <- (mod | r2) +
  plot_layout(widths = c(1, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(legend.position = "top",
        legend.justification = "left",
        plot.tag = element_text(face = "bold"))

fig_glmm

ggsave("Figures/GLMM_detectability.png", fig_glmm, width = 12, height = 5.5, dpi = 300, bg = "white")
ggsave("Figures/GLMM_detectability.pdf", fig_glmm, width = 12, height = 5.5, bg = "white")



### Figure S3 - %Canopy ####
library(dplyr)
library(ggplot2)
library(ggpubr)

# Load data
crown <- read.csv("Data/Crown_sampled.csv")

# Plant genera (dataset name = label) and colours
genera <- c(Pouroma = "Pourouma", Euterpe = "Euterpe", Virola = "Virola",
            Iryanthera = "Iryanthera", Pouteria = "Pouteria")
colours <- c(Pouroma = "darkolivegreen", Euterpe = "mediumpurple4", Virola = "sienna2",
             Iryanthera = "firebrick3", Pouteria = "goldenrod")

# Histogram for each plant genus
plot_crown <- function(g) {
  ggplot(filter(crown, Focal_Genus == g), aes(x = Sampled_Canopy)) +
    geom_histogram(breaks = seq(0, 100, by = 5), fill = colours[[g]], colour = "white") +
    scale_x_continuous(breaks = seq(0, 100, by = 10)) +
    labs(title = paste0(genera[[g]], " sp."), x = NULL, y = NULL) +
    theme_classic() +
    theme(plot.title = element_text(face = "italic"))
}

# Combined figure
fig <- ggarrange(plotlist = lapply(names(genera), plot_crown))

fig_final <- annotate_figure(fig,
                             left = text_grob("Number of camera deployments", rot = 90, size = 12),
                             bottom = text_grob("Crown sampled (%)", size = 12))

fig_final

# Save
ggsave("Figures/Fig_S_crown_sampled.png", fig_final, width = 10, height = 6, dpi = 300, bg = "white")
ggsave("Figures/Fig_S_crown_sampled.pdf", fig_final, width = 10, height = 6, dpi = 300, bg = "white")


### Figure S4 - SIMPER analysis####

library(vegan)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)

data_TR4 <- read.csv("Results/forag_rates.csv", sep = ",")

#Foraging rates in long format
tu <- data_TR4 %>%
  ungroup() %>%
  pivot_longer(
    cols = starts_with("RAI_"),
    names_to = "Primate",
    names_prefix = "RAI_",
    values_to = "TU"
  ) %>%
  mutate(
    Hunting_Pressure = recode(Hunting_Pressure,
                              "Non_hunted" = "NH",
                              "Hunted" = "HU",
                              "Overhunted" = "OH"),
    Hunting_Pressure = factor(Hunting_Pressure, levels = c("NH", "HU", "OH"))
  )

#Proportion of foraging of each primate on each plant genus, per hunting pressure
prop <- tu %>%
  group_by(Hunting_Pressure, Primate) %>%
  filter(sum(TU, na.rm = TRUE) > 0) %>%
  mutate(p = TU / sum(TU, na.rm = TRUE)) %>%
  ungroup() %>%
  select(Hunting_Pressure, Primate, Focal_Genus, p)

#Trophic overlap modules (non-hunted sites)
primate_modules <- data.frame(
  Primate = c("Ateles", "Lagothrix", "Cacajao", "Pithecia",
              "Alouatta", "Sapajus", "Cebus", "Callicebus", "Aotus", "Saimiri", "Saguinus"),
  Module  = c(1, 1, 1, 1,
              2, 2, 2, 2, 2, 2, 2)
)

#Primate x plant genus matrix (non-hunted sites)
matrix_NH <- prop %>%
  filter(Hunting_Pressure == "NH") %>%
  select(Primate, Focal_Genus, p) %>%
  pivot_wider(names_from = Focal_Genus, values_from = p, values_fill = 0) %>%
  column_to_rownames("Primate")

#Module of each primate, in the same order as the matrix rows
groups <- primate_modules$Module[match(rownames(matrix_NH), primate_modules$Primate)]

# Check: no NA (all primates have a module)
data.frame(Primate = rownames(matrix_NH), Module = groups)

#SIMPER
set.seed(123)
sim <- simper(matrix_NH, group = groups, permutations = 999)
summary(sim)

#Module order used by simper (ava = first module, avb = second module)
comp <- strsplit(names(sim)[1], "_")[[1]]
comp

#Table for the Supplementary Material
Table_SIMPER <- summary(sim)[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("Plant genus") %>%
  transmute(`Plant genus` = paste0(`Plant genus`, " spp."),
            `Average contribution` = round(average, 3),
            `SD` = round(sd, 3),
            !!paste0("Mean proportion module ", comp[1]) := round(ava, 2),
            !!paste0("Mean proportion module ", comp[2]) := round(avb, 2),
            `Cumulative contribution (%)` = round(cumsum * 100, 1),
            `p` = p) %>%
  mutate(`Plant genus` = sub("Pouroma", "Pourouma", `Plant genus`))

Table_SIMPER

#Figure
genus_levels <- rev(Table_SIMPER$`Plant genus`)

fig_simper <- Table_SIMPER %>%
  select(`Plant genus`, `Mean proportion module 1`, `Mean proportion module 2`) %>%
  pivot_longer(
    cols = c(`Mean proportion module 1`, `Mean proportion module 2`),
    names_to = "Module",
    values_to = "Proportion") %>%
  mutate(Module = sub("Mean proportion ", "", Module),
         `Plant genus` = factor(`Plant genus`, levels = genus_levels)) %>%
  ggplot(aes(x = Proportion, y = `Plant genus`, fill = Module)) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(data = Table_SIMPER %>% mutate(`Plant genus` = factor(`Plant genus`, levels = genus_levels)),
            aes(x = pmax(`Mean proportion module 1`, `Mean proportion module 2`) + 0.03,
                y = `Plant genus`,
                label = paste0(round(`Average contribution` * 100, 1), "%")),
            inherit.aes = FALSE,
            hjust = 0,
            size = 3.5) +
  scale_fill_manual(
    values = c("module 1" = "goldenrod",
               "module 2" = "darkolivegreen"),
    labels = c("module 1" = "Module 1",
               "module 2" = "Module 2")) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = "Mean proportion of foraging events",
       y = NULL,
       fill = NULL) +
  theme_classic() +
  theme(axis.text.y = element_text(face = "italic"),
        legend.position = "top")

fig_simper

#Save
ggsave("Figures/Fig_simper.png", fig_simper, width = 8, height = 7, dpi = 300, bg = "white")
ggsave("Figures/Fig_simper.pdf", fig_simper, width = 8, height = 7, dpi = 300, bg = "white")