#Sensitivity analysis independence threshold
library(dplyr)
library(tidyr)

tr15 <- read.csv("Data/tr_15.csv", sep = ",")
tr30 <- read.csv("Data/tr_30.csv", sep = ",")
tr60 <- read.csv("Data/tr_60.csv", sep = ",")

# Number of individuals per tree and primate genus under each threshold
comparison <- bind_rows(
  tr15 = tr15, tr30 = tr30, tr60 = tr60,
  .id = "threshold"
) %>%
  mutate(across(c(Focal_Genus, Tree_ID, Vertebrate_ID), trimws)) %>%
  count(Focal_Genus, Tree_ID, Vertebrate_ID, threshold) %>%
  pivot_wider(names_from = threshold, values_from = n, values_fill = 0) %>%
  select(Focal_Genus, Tree_ID, Vertebrate_ID, tr15, tr30, tr60)

# Only combinations that differ between thresholds
differences <- comparison %>%
  filter(!(tr15 == tr30 & tr30 == tr60)) %>%
  mutate(dif_30 = tr15 - tr30,
         dif_60 = tr15 - tr60) %>%
  arrange(Focal_Genus, Tree_ID, Vertebrate_ID)

differences

# Save
write.csv(differences, "Results/sensitivity_threshold.csv", row.names = FALSE)
