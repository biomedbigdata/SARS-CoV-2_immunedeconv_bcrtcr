# compute sensitivity & specificity for reconstructed bcr and tcr sequences
library(data.table)
library(dplyr)
library(tidyverse)
library(igraph)


load("data/filtered_bcrtcr_seqs/reconstructed_cdr3_sequences.RData")

full <- cdr3_sequences[[1]]
m50 <- cdr3_sequences[[2]]
m10 <- cdr3_sequences[[3]]

# matrices
distances_full <- as.matrix(fread("data/filtered_bcrtcr_seqs/matrices/all_seqs_dist_mat_10.csv"))
distances_m50 <- as.matrix(fread("data/filtered_bcrtcr_seqs/matrices/all_seqs_50_dist_mat_10.csv"))
distances_m10 <- as.matrix(fread("data/filtered_bcrtcr_seqs/matrices/all_seqs_10_dist_mat_10.csv"))

colnames(distances_full) <- full$CDR3aa
rownames(distances_full) <- full$CDR3aa
colnames(distances_m50) <- m50$CDR3aa
rownames(distances_m50) <- m50$CDR3aa
colnames(distances_m10) <- m10$CDR3aa
rownames(distances_m10) <- m10$CDR3aa

# filter for sequences that were found with both tools
full <- full[Source == "MIXCR_TRUST4"]
m50 <- m50[Source == "MIXCR_TRUST4"]
m10 <- m10[Source == "MIXCR_TRUST4"]

## get the sequences that are unique & not connected to seronegative
# annotate sequences in seronegatives
full[, variant_only := ifelse(!grepl("sero", Groups), T, F)]
m50[, variant_only := ifelse(!grepl("sero", Groups), T, F)]
m10[, variant_only := ifelse(!grepl("sero", Groups), T, F)]

dim(full[variant_only==T])
dim(m50[variant_only==T])
dim(m10[variant_only==T])

# filter matrices
distances_full <- distances_full[full$CDR3aa, full$CDR3aa]
distances_m50 <- distances_m50[m50$CDR3aa, m50$CDR3aa]
distances_m10 <- distances_m10[m10$CDR3aa, m10$CDR3aa]

# create graphs
g_full = graph_from_adjacency_matrix(distances_full, mode = "undirected", weighted = T, diag = F)
g_m50 = graph_from_adjacency_matrix(distances_m50, mode = "undirected", weighted = T, diag = F)
g_m10 = graph_from_adjacency_matrix(distances_m10, mode = "undirected", weighted = T, diag = F)


# get components with seronegative
c_full <- components(g_full)
full[, component := c_full$membership]
c_sero_full <- unique(full[grepl("sero", Groups)]$component)

c_m50 <- components(g_m50)
m50[, component := c_m50$membership]
c_sero_m50 <- unique(m50[grepl("sero", Groups)]$component)

c_m10 <- components(g_m10)
m10[, component := c_m10$membership]
c_sero_m10 <- unique(m10[grepl("sero", Groups)]$component)

full[, unique := ifelse(component %in% c_sero_full, F, T)]
m50[, unique := ifelse(component %in% c_sero_m50, F, T)]
m10[, unique := ifelse(component %in% c_sero_m10, F, T)]

dim(full[unique == T])
dim(m50[unique == T])
dim(m10[unique == T])


# define gold standard seqs 
gold_standard_seqs <- full$CDR3aa

m10_seqs <- m10$CDR3aa

# and compute distances between the two sets of sequences --> calc_distances_gold.ipynb

# write.csv(gold_standard_seqs, file = "data/filtered_bcrtcr_seqs/gold_standard_seqs.csv", row.names = F)
# write.csv(m10_seqs, file = "data/filtered_bcrtcr_seqs/sensitivity_testing_seqs.csv", row.names = F)

gold_vs_m10_distances <- fread("data/filtered_bcrtcr_seqs/matrices/gold_standard_vs_m10.csv")
rownames(gold_vs_m10_distances) <- gold_standard_seqs
colnames(gold_vs_m10_distances) <- m10_seqs

# TODO: define the sets of TP, FP, TN, FN 




