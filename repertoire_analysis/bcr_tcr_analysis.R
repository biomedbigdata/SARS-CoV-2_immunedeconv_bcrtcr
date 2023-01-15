# script for bcr tcr result analysis
# 0. prepare data tables for plotting
# 1. overlap between MIXCR and TRUST4 results
# 2. Graph with all groups, coloring variant_only (for different cutoffs?) and connected components
# 3. Graph with all groups, coloring groups 
# 4. Graph with all groups except Seronegative

library(data.table)
library(dplyr)
library(tidyverse)
library(igraph)

### 0. prepare data tables for plotting

# read the result files
groups <- c("alpha", "alpha_ek", "omicron_ba1", "omicron_ba2", "sero")
groups_short <- c("alpha", "alpha_ek", "ba1", "ba2", "sero")

mixcr_results <- lapply(groups, function(group) 
  fread(paste0("data/filtered_bcrtcr_seqs/mixcr_", group, "_seqs.csv"))[, .(CDR3aa, N)])
trust_results <- lapply(groups, function(group) 
  fread(paste0("data/filtered_bcrtcr_seqs/trust_", group, "_seqs.csv"))[, .(CDR3aa, N)])

names(mixcr_results) <- groups_short
names(trust_results) <- groups_short


# filter the results for sequences that appear in at least half of the samples (rounded, floor)
# and add column for group
num_samples <- c("alpha" = 61, "alpha_ek" = 29, "ba1" = 78, "ba2" = 22, "sero" = 47)

mixcr_results <- lapply(groups_short, function(group){
  mixcr_results[[group]] %>% 
    filter(N >= floor(num_samples[group]/2)) %>%
    add_column(Group = group) %>%
    add_column(Source = "MIXCR")
})

trust_results <- lapply(groups_short, function(group){
  trust_results[[group]] %>% 
    filter(N >= floor(num_samples[group]/2)) %>%
    add_column(Group = group) %>%
    add_column(Source = "TRUST4")
})

# set names again
names(mixcr_results) <- groups_short
names(trust_results) <- groups_short


# combine the results to one data table
mixcr_results_all <- rbindlist(mixcr_results)
trust_results_all <- rbindlist(trust_results)

# remove "out_of_frame" sequence in TRUST4 results
trust_results_all <- trust_results_all[CDR3aa != "out_of_frame"]


# combine the results
results_all <- rbind(mixcr_results_all, trust_results_all)

# create a data table with unique sequences and summarized group and source
seqs_all <- rbindlist(lapply(unique(results_all$CDR3aa), function(seq){
  groups_m = paste(sort(unique(results_all[CDR3aa == seq & Source == "MIXCR"]$Group)), collapse = ", ")
  groups_t = paste(sort(unique(results_all[CDR3aa == seq & Source == "TRUST4"]$Group)), collapse = ", ")
  groups = paste(sort(unique(results_all[CDR3aa == seq]$Group)), collapse = ", ")
  source = paste(sort(unique(results_all[CDR3aa == seq]$Source)), collapse = "_")
  data.table(CDR3aa = seq,
             Group_MIXCR = groups_m,
             Group_TRUST4 = groups_t,
             Groups = groups,
             Source = source)
}))



### Distance matrices

# compute distance matrix (use python script) -> reticulate does not work
# write sequences to file and run with python seperately

# write.csv(seqs_all$CDR3aa, file = "data/filtered_bcrtcr_seqs/all_seqs.csv", row.names = F)

# run calc_dist_mat.ipynb

dist_mat_all_5 <- as.matrix(read.table("data/filtered_bcrtcr_seqs/matrices/all_seqs_dist_mat_5.csv", sep = " "))
dist_mat_all_10 <- as.matrix(read.table("data/filtered_bcrtcr_seqs/matrices/all_seqs_dist_mat_10.csv", sep = " "))
dist_mat_all_15 <- as.matrix(read.table("data/filtered_bcrtcr_seqs/matrices/all_seqs_dist_mat_15.csv", sep = " "))


colnames(dist_mat_all_5) <- seqs_all$CDR3aa
rownames(dist_mat_all_5) <- seqs_all$CDR3aa
colnames(dist_mat_all_10) <- seqs_all$CDR3aa
rownames(dist_mat_all_10) <- seqs_all$CDR3aa
colnames(dist_mat_all_15) <- seqs_all$CDR3aa
rownames(dist_mat_all_15) <- seqs_all$CDR3aa

### 1. overlap between MIXCR and TRUST4 results

# create igraph

# overlap in numbers
num_source <-prop.table(table(seqs_all$Source))

g = graph_from_adjacency_matrix(dist_mat_all_10, mode = "undirected", weighted = T, diag = F)
# run only once for same layout in every graph: 
coords_10 <- layout_nicely(g) 

colors = c("MIXCR"= "#fc8d62", "TRUST4"= "#8da0cb", "MIXCR_TRUST4"= "#66c2a5")
V(g)$color <- sapply(seqs_all$Source, function(group) colors[group])

# with labels
# plot.igraph(g, vertex.size = 6, vertex.label.cex = 0.4, vertex.label.color = "black")

# without labels
plot.igraph(g, vertex.size = 6, vertex.label = NA, layout = coords_10)
legend("bottomright",legend=c("MIXCR only", "TRUST4 only", "MIXCR and TRUST4"), 
       fill=c("#fc8d62","#8da0cb","#66c2a5"), bty = "n")
legend("bottomleft", legend = c(paste0("sequences found by MIXCR only: ", round(num_source["MIXCR"], 2)), 
                                paste0("sequences found by TRUST4 only: ", round(num_source["TRUST4"], 2)), 
                                paste0("sequences found by MIXCR and TRUST4: ", round(num_source["MIXCR_TRUST4"],2))),
       title = "Percentage of nodes per group:", title.adj = 0.19, bty = "n")
title(main = "Comparison of BCR and TCR CDR3 sequences reconstructed by MIXCR and TRUST4")





### 2. Graph with all groups, coloring variant_only (for different cutoffs?) and connected components

# !!!!!!!! IMPORTANT !!!!!!!!!  only keep sequences that are found by mixcr and trust

seqs_both <- seqs_all[Source == "MIXCR_TRUST4"]$CDR3aa
seqs_all <- seqs_all[Source == "MIXCR_TRUST4"]
dist_mat_all_10 <- dist_mat_all_10[seqs_both,seqs_both]
dist_mat_all_5 <- dist_mat_all_5[seqs_both,seqs_both]
dist_mat_all_15 <- dist_mat_all_15[seqs_both,seqs_both]

## cutoff 10
# same graph as above but now different coloring
g = graph_from_adjacency_matrix(dist_mat_all_10, mode = "undirected", weighted = T, diag = F)
coords_10 <- layout_nicely(g)


# add a new column to seqs_all defining whether sero is in groups
# column to combine MIXCR and TRUST groups

seqs_all[, variant_only := ifelse(!grepl("sero", Groups), "true", "false")]
# seqs_all[variant_only == TRUE]    # test

num_variant_only <-table(seqs_all$variant_only)

colors = c("true"= "#8da0cb", "false"= "#f7f7f7")
V(g)$color <- sapply(seqs_all$variant_only, function(group) colors[group])


plot.igraph(g, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_10)
legend("bottomright", legend = c("sequence from infected only", "sequence also in seronegative"),
       fill = c("#8da0cb", "#f7f7f7"), bty = "n")
legend("bottomleft", legend = c(paste0("sequences from infected only: ", num_variant_only["true"]),
                                 paste0("sequences found also in seronegative: ", num_variant_only["false"])),
       title = "Number of nodes per group:",
       title.adj = 0.19, bty = "n")
title("CDR3 sequences found only in samples with COVID-19 infection")


##  analyze clusters/ connected cmponents
g5 <- graph_from_adjacency_matrix(dist_mat_all_5, mode = "undirected", weighted = T, diag = F)
g15 <- graph_from_adjacency_matrix(dist_mat_all_15, mode = "undirected", weighted = T, diag = F)

coords_5 <- layout_nicely(g5)
coords_15 <- layout_nicely(g15)


ccs <- components(g)
ccs5 <- components(g5)
ccs15 <- components(g15)

seqs_all[, component_10 := ccs$membership]
seqs_all[, component_5 := ccs5$membership]
seqs_all[, component_15 := ccs15$membership]

components_sero10 <- unique(seqs_all[grepl("sero", Groups)]$component_10)
components_sero5 <- unique(seqs_all[grepl("sero", Groups)]$component_5)
components_sero15 <- unique(seqs_all[grepl("sero", Groups)]$component_15)


# Graph coloring components without connection to seronegative

# cutoff 5
g5 = graph_from_adjacency_matrix(dist_mat_all_5, mode = "undirected", weighted = T, diag = F)

V(g5)$color <- sapply(seqs_all$component_5, function(component) 
  colors[ifelse(component %in% components_sero5, "false", "true")])


plot.igraph(g5, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_5)
legend("bottomright", legend = c("without connection to seronegative", "with connection to seronegative"),
       fill = c("#8da0cb", "#f7f7f7"), bty = "n")
legend("bottomleft", legend = c(paste0("Number of sequences without connection to a seronegative sequence: ", dim(seqs_all[!component_5 %in% components_sero5])[1])),
       bty = "n")
title("Connected components with and without sequences found also in seronegatives (cutoff = 5)")


# cutoff 10
g = graph_from_adjacency_matrix(dist_mat_all_10, mode = "undirected", weighted = T, diag = F)

colors = c("true"= "#8da0cb", "false"= "#f7f7f7")
V(g)$color <- sapply(seqs_all$component_10, function(component) 
  colors[ifelse(component %in% components_sero10, "false", "true")])


plot.igraph(g, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_10)
legend("bottomright", legend = c("without connection to seronegative", "with connection to seronegative"),
       fill = c("#8da0cb", "#f7f7f7"), bty = "n")
legend("bottomleft", legend = c(paste0("Number of sequences without connection to a seronegative sequence: ", dim(seqs_all[!component_10 %in% components_sero10])[1])),
       bty = "n")
title("Connected components with and without sequences found also in seronegatives (cutoff = 10)")



# cutoff 15
g15 = graph_from_adjacency_matrix(dist_mat_all_15, mode = "undirected", weighted = T, diag = F)

V(g15)$color <- sapply(seqs_all$component_15, function(component) 
  colors[ifelse(component %in% components_sero15, "false", "true")])


plot.igraph(g15, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_15)
legend("bottomright", legend = c("without connection to seronegative", "with connection to seronegative"),
       fill = c("#8da0cb", "#f7f7f7"), bty = "n")
legend("bottomleft", legend = c(paste0("Number of sequences without connection to a seronegative sequence: ", dim(seqs_all[!component_15 %in% components_sero15])[1])),
       bty = "n")
title("Connected components with and without sequences found also in seronegatives (cutoff = 15)")



# TODO: analyze the following sequences

# sequences that are not connected with a sequence that appears in seronegatives
seqs_all[!component_10 %in% components_sero10]
seqs_all[!component_5 %in% components_sero5]
seqs_all[!component_15 %in% components_sero15]








### 3. Graph with all groups, coloring groups 

g = graph_from_adjacency_matrix(dist_mat_all_10, mode = "undirected", weighted = T, diag = F)

colors = c("alpha"= "#d73027", "alpha_ek"= "#fc8d59", "ba1"="#4575b4", "ba2"="#91bfdb", "sero" = "transparent", "any_sero"= "transparent", "any_nonsero" = "#fee090")

# add another column to seqs all
seqs_all[, group_color1 := ifelse(!grepl( ",", Groups), Groups,
                                  ifelse(grepl("sero", Groups), "any_sero", "any_nonsero"))]


num_groups <- table(seqs_all$group_color1)

V(g)$color <- sapply(seqs_all$group_color1, function(group) 
  colors[group])


plot.igraph(g, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_10)
legend("bottomright", legend = c("Alpha only", "Alpha+EK only", "BA.1 only", "BA.2 only",  "more than one variant", "with seronegative" ),
       fill = c("#d73027", "#fc8d59", "#4575b4", "#91bfdb", "#fee090", "transparent"), 
       bty = "n")
legend("bottomleft", legend = c(paste0("Number of sequences in Alpha only: ", num_groups["alpha"]),
                                paste0("Number of sequences in Alpha+EK only: ", num_groups["alpha_ek"]),
                                paste0("Number of sequences in BA.1 only: ", num_groups["ba1"]),
                                paste0("Number of sequences in BA.2 only: ", num_groups["ba2"]),
                                paste0("Number of sequences in more than one variant: ", num_groups["any_nonsero"]),
                                paste0("Number of sequences also in Seronegative: ", num_groups["any_sero"])),
       bty = "n")
title("BCR and TCR sequences in each variant")




# 4. Graph with all groups except Seronegative

seqs_all <- seqs_all[variant_only == "true"]
dist_mat_all_10 <- dist_mat_all_10[seqs_all$CDR3aa,seqs_all$CDR3aa]

g = graph_from_adjacency_matrix(dist_mat_all_10, mode = "undirected", weighted = T, diag = F)
coords_10 <- layout_nicely(g)

colors = c("alpha"= "#d73027", "alpha_ek"= "#fee090", "alpha, alpha_ek"= "#fc8d59",
           "ba1"="#4575b4", "ba2"="#e0f3f8", "ba1, ba2" = "#91bfdb",
           "else" = "#ffffbf")


V(g)$color <- sapply(seqs_all$Groups, function(group) 
  ifelse(grepl("alpha", group) & grepl("ba", group), colors["else"], colors[group]))


num_groups <- table(V(g)$color)

plot.igraph(g, vertex.size = 6, vertex.label = NA, vertex.label.color = "black", layout = coords_10)
legend("bottomright", legend = c("Alpha", "Alpha+EK", "Alpha and Alpha+EK", "BA.1", "BA.2", "BA.1 and BA.2", "any other combination"),
       fill = c("#d73027", "#fee090","#fc8d59","#4575b4", "#e0f3f8", "#91bfdb", "#ffffbf"),
       bty = "n")
legend("bottomleft", legend = c(paste0("Number of sequences in Alpha only: ", num_groups["#d73027"]),
                                paste0("Number of sequences in Alpha+EK only: ", num_groups["#fee090"]),
                                paste0("Number of sequences in Alpha and Alpha+EK: ", num_groups["#fc8d59"]),
                                paste0("Number of sequences in BA.1 only: ", num_groups["#4575b4"]),
                                paste0("Number of sequences in BA.2 only: ", num_groups["#e0f3f8"]),
                                paste0("Number of sequences in BA.1 and BA.2: ", num_groups["#91bfdb"]),
                                paste0("Number of sequences in any oter combination: ", num_groups["#ffffbf"])),
       bty = "n")
title("BCR and TCR sequences in each variant")

