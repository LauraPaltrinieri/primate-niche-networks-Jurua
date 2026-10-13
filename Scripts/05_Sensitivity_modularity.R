#Module stability (sensitivity analysis)
#Author: Laura Paltrinieri

#This script assesses the stability of the module partitions identified by the
#Louvain algorithm in the significantly modular networks (Spatial, Trophic and
#Combined networks in non-hunted forests). The algorithm is run 1000 times with
#different random seeds and each partition is compared with the observed one.
#The observed partitions are also compared with the exact maximum-modularity
#partition (cluster_optimal), which is feasible given the small network size.

# Load libraries
library(igraph)

#clean_env
rm(list = ls())

##### Load data #####
##Observed networks and module partitions (saved in the Network metrics script)
S_mod  <- readRDS("Results/S_modularity_NH.rds")
TR_mod <- readRDS("Results/TR_modularity_NH.rds")
C_mod  <- readRDS("Results/C_modularity_NH.rds")

g_S_NH     <- S_mod$g
comm_S_NH  <- S_mod$comm
g_TR_NH    <- TR_mod$g
comm_TR_NH <- TR_mod$comm
g_C_NH     <- C_mod$g
comm_C_NH  <- C_mod$comm

n_run <- 1000

##### Functions #####
##Relabel modules by order of appearance (so that {1,1,2} and {2,2,1} are equal)
canon <- function(m) match(m, unique(m))

##Stability of the module partition across Louvain runs
check_stability <- function(g, comm, n_run = 1000) {
  ref <- membership(comm)                    # observed partition
  stab <- lapply(1:n_run, function(i) {
    set.seed(i)                              # different seed for each run
    m <- membership(cluster_louvain(g, weights = E(g)$weight))
    data.frame(Q = modularity(g, m, weights = E(g)$weight),
               n_modules = length(unique(m)),
               NMI = compare(ref, m, method = "nmi"),
               part = paste(canon(m), collapse = "-"))
  })
  stab <- do.call(rbind, stab)
  data.frame(Q_mean = round(mean(stab$Q), 3),
             Q_sd = round(sd(stab$Q), 3),
             n_modules = paste(sort(unique(stab$n_modules)), collapse = ", "),
             n_partitions = length(unique(stab$part)),
             NMI_mean = round(mean(stab$NMI), 3),
             perc_identical = round(100 * mean(stab$NMI > 0.999), 1))
}

##Module membership of each partition found and how often it occurs
partition_detail <- function(g, n_run = 1000) {
  runs <- lapply(1:n_run, function(i) {
    set.seed(i)
    canon(membership(cluster_louvain(g, weights = E(g)$weight)))
  })
  keys <- sapply(runs, paste, collapse = "-")
  freq <- sort(table(keys), decreasing = TRUE)
  out <- data.frame(Genus = V(g)$name)
  for (k in names(freq)) {
    perc <- round(100 * freq[[k]] / n_run, 1)
    out[[paste0("Partition (", perc, "%)")]] <- runs[[match(k, keys)]]
  }
  out
}

##Comparison with the exact maximum-modularity partition
check_optimal <- function(g, comm) {
  opt <- cluster_optimal(g, weights = E(g)$weight)
  data.frame(Q_louvain = round(modularity(g, membership(comm), weights = E(g)$weight), 3),
             Q_optimal = round(modularity(opt), 3),
             NMI_optimal = round(compare(membership(comm), membership(opt), method = "nmi"), 3))
}

##### Analyses #####
##Stability across Louvain runs
stability_table <- rbind(
  Spatial_NH  = check_stability(g_S_NH,  comm_S_NH,  n_run),
  Trophic_NH  = check_stability(g_TR_NH, comm_TR_NH, n_run),
  Combined_NH = check_stability(g_C_NH,  comm_C_NH,  n_run)
)
stability_table

##Module membership of each partition
part_S_NH  <- partition_detail(g_S_NH,  n_run)
part_TR_NH <- partition_detail(g_TR_NH, n_run)
part_C_NH  <- partition_detail(g_C_NH,  n_run)

part_S_NH
part_TR_NH
part_C_NH

##Comparison with the maximum-modularity partition
optimal_table <- rbind(
  Spatial_NH  = check_optimal(g_S_NH,  comm_S_NH),
  Trophic_NH  = check_optimal(g_TR_NH, comm_TR_NH),
  Combined_NH = check_optimal(g_C_NH,  comm_C_NH)
)
optimal_table

##### Supplementary table #####
module_stability_table <- data.frame(
  Network = c("Spatial", "Trophic", "Multidimensional"),
  `Observed Q (Louvain)` = optimal_table$Q_louvain,
  `Q across runs (mean ± SD)` = paste0(sprintf("%.3f", stability_table$Q_mean), " ± ",
                                       sprintf("%.3f", stability_table$Q_sd)),
  `Number of partitions` = stability_table$n_partitions,
  `Runs identical to observed partition (%)` = stability_table$perc_identical,
  `Mean NMI` = stability_table$NMI_mean,
  `Maximum Q (exact)` = optimal_table$Q_optimal,
  `NMI with maximum-Q partition` = optimal_table$NMI_optimal,
  check.names = FALSE
)
module_stability_table

##### Export supplementary table to Word #####
library(flextable)
library(officer)

ft <- flextable(module_stability_table)
ft <- theme_booktabs(ft)
ft <- autofit(ft)
ft <- align(ft, j = 2:8, align = "center", part = "all")

doc <- read_docx()
doc <- body_add_par(doc, "Table S7. Stability of module partitions in the significantly modular niche overlap networks of non-hunted forests.", style = "Normal")
doc <- body_add_flextable(doc, ft)
print(doc, target = "Results/Table_S7_Module_stability.docx")

##### Save results #####
write.csv(module_stability_table, "Results/Module_stability_table.csv", row.names = FALSE)
write.csv(part_S_NH,  "Results/Module_partitions_S_NH.csv",  row.names = FALSE)
write.csv(part_TR_NH, "Results/Module_partitions_TR_NH.csv", row.names = FALSE)
write.csv(part_C_NH,  "Results/Module_partitions_C_NH.csv",  row.names = FALSE)
