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
m50_seqs <- m50$CDR3aa

# and compute distances between the two sets of sequences --> calc_distances_gold.ipynb

# write.csv(gold_standard_seqs, file = "data/filtered_bcrtcr_seqs/gold_standard_seqs.csv", row.names = F)
# write.csv(m10_seqs, file = "data/filtered_bcrtcr_seqs/sensitivity_testing_seqs.csv", row.names = F)
# write.csv(m50_seqs, file = "data/filtered_bcrtcr_seqs/sensitivity_testing_50m_seqs.csv", row.names = F)

gold_vs_m10_distances <- as.matrix(fread("data/filtered_bcrtcr_seqs/matrices/gold_standard_vs_m10.csv"))
rownames(gold_vs_m10_distances) <- gold_standard_seqs
colnames(gold_vs_m10_distances) <- m10_seqs

gold_vs_m50_distances <- as.matrix(fread("data/filtered_bcrtcr_seqs/matrices/gold_standard_vs_m50.csv"))
rownames(gold_vs_m50_distances) <- gold_standard_seqs
colnames(gold_vs_m50_distances) <- m50_seqs

# define the sets of TP, FP, TN, FN 
gold_true <- full[unique==T]$CDR3aa
m10_true <- m10[unique==T]$CDR3aa
m50_true <- m50[unique==T]$CDR3aa

gold_false <- full[unique==F]$CDR3aa
m10_false <- m10[unique==F]$CDR3aa
m50_false <- m50[unique==F]$CDR3aa

gold_vs_m10_distances[gold_true, m10_true]

gold_vs_m50_distances[gold_true, m50_true]


compute_recall <- function(test_true, test_false, distances){
    tp <- unlist(lapply(test_true, function(seq){
      similar_seqs <- rownames(distances)[distances[,seq] <= 10]
      if (length(similar_seqs[similar_seqs %in% gold_true])>0){
        TRUE
      } else {
        FALSE
      }
    }))
    n_tp <- length(tp[tp])
    
    
    fp <- unlist(lapply(test_true, function(seq){
      similar_seqs <- rownames(distances)[distances[,seq] <= 10]
      if (length(similar_seqs[similar_seqs %in% gold_false])>0){
        TRUE
      } else {
        FALSE
      }
    }))
    n_fp <- length(fp[fp])
    
    tn <- unlist(lapply(test_false, function(seq){
      similar_seqs <- rownames(distances)[distances[,seq] <= 10]
      if (length(similar_seqs[similar_seqs %in% gold_false])>0){
        TRUE
      } else {
        FALSE
      }
    }))
    n_tn <- length(tn[tn])
    
    fn <- unique(unlist(lapply(test_false, function(seq){
      similar_seqs <- rownames(distances)[distances[,seq] <= 10]
      if (length(similar_seqs[similar_seqs %in% gold_true])>0){
        TRUE
      } else {
        FALSE
      }
    })))
    n_fn <- length(fn[fn])
    
    
    # sensitivity: tp/(tp+fn) = recall
    # specificity: tn/(tn+fp)
    print(paste("TP:", n_tp, "FP:", n_fp, "TN:", n_tn, "FN:", n_fn))
    sensitivity <- n_tp/(n_tp+n_fn) 
    specificity <- n_tn/(n_tn+n_fp)
    fpr <- n_fp/(n_fp+n_tn)
    c("sensitivity"=sensitivity, 
      "specificity"=specificity,
      "false positive rate"=fpr)
    
    # confusion matrix
    #TClass <- factor(c("present in Seronegative", "present in Seronegative", "absent in Seronegative", "absent in Seronegative"), levels = c("absent in Seronegative", "present in Seronegative"))
    #PClass <- factor(c("present in Seronegative", "absent in Seronegative", "present in Seronegative", "absent in Seronegative"), levels = c("present in Seronegative", "absent in Seronegative"))
    TClass <- factor(c("present in Seron.", "present in Seron.", "absent in Seron.", "absent in Seron."), levels = c("absent in Seron.", "present in Seron."))
    PClass <- factor(c("present in Seron.", "absent in Seron.", "present in Seron.", "absent in Seron."), levels = c("present in Seron.", "absent in Seron."))
    Y      <- c(n_tn, n_fp, n_fn, n_tp)
    data.frame(TClass, PClass, Y)
    # data.frame(label = c("TP", "FP", "TN", "FN"),
    #            value = c(n_tp, n_fp, n_tn, n_fn))
}
### hmpf
conf_matrix_m10 <- compute_recall(m10_true, m10_false, gold_vs_m10_distances)
conf_matrix_m50 <- compute_recall(m50_true, m50_false, gold_vs_m50_distances)



my_theme2 <- theme(panel.background = element_rect(fill = "white", colour = "white"),
                  panel.grid.major = element_line(colour = "black"),
                  panel.grid.minor = element_line(colour = "black"),
                  text = element_text(size = 12),
                  plot.title = element_text(size = 14),
                  axis.text.y = element_text(angle = 90, vjust = 0, hjust=0.5),
                  plot.margin = margin(t = 10, r = 30, b = 10, l = 30, unit = "pt"))

conf_matrix_m10$color<- c("a", "b", "b", "a")
C <- ggplot(data = conf_matrix_m10, mapping = aes(x = TClass, y = PClass)) +
  geom_tile(aes(fill = color), colour = "black") +
  geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 0.5, size=6) +
  # scale_fill_gradient2(low = "white", high = "dodgerblue") +
  scale_fill_manual(values = c("a" = "palegreen", "b" = "salmon"))+
  my_theme2 + theme(legend.position = "none") +
  labs(x = "full sequencing depth",
       y = "downsampled sequencing depth",
       title = "10 Million")

conf_matrix_m50$color<- c("a", "b", "b", "a")
B <- ggplot(data = conf_matrix_m50, mapping = aes(x = TClass, y = PClass)) +
  geom_tile(aes(fill = color), colour = "black") +
  geom_text(aes(label = sprintf("%1.0f", Y)), vjust = 0.5, size=6) +
  #scale_fill_gradient(low = "white", high = "dodgerblue") +
  scale_fill_manual(values = c("a" = "palegreen", "b" = "salmon"))+
  my_theme2 + theme(legend.position = "none") +
  labs(x = "full sequencing depth",
       y = "downsampled sequencing depth",
       title = "50 Million")


seqs_with_sarscov2 <- c("CYSTDSSGNHRGVF", "CQSYDSSNVVF", "CETWDSNTRVF", "CQQRSNWPPTWTF",
                        "CLQHDNFPLTF", "CQAWDSSVVF", "CSSYTSSSTVF")
seqs_with_sarscov2_10 <- c("CAAWDDSLNGWVF",
                           "CAAWDDSLNGPVF",
                           "CQSADSSGTYVVF",
                           "CQSYDSSLSGSVF",
                           "CGTWDSSLSAGVF",
                           "CLQHNSYPWTF",
                           "CMQATQFPRTF" )


full[, sarscov2 := ifelse(CDR3aa %in% seqs_with_sarscov2, T, F)]

library(tidyr)

# Splitting the Groups column
dt_long <- full[unique==T, .(CDR3aa, sarscov2, Group = unlist(strsplit(Groups, ",\\s*"))), by = .(CDR3aa, sarscov2)]

# Display the transformed data table
print(dt_long)

df <- melt(table(dt_long$Group, dt_long$sarscov2))
df <- df %>%
  mutate(Var1 = case_when(
    Var1 == "alpha"    ~ "Alpha",
    Var1 == "alpha_ek" ~ "Alpha+EK",
    Var1 == "ba1"      ~ "BA.1",
    Var1 == "ba2"      ~ "BA.2",
    TRUE               ~ "0" # Default case if none of the above conditions are met
  ))

A <- ggplot(data = df, mapping = aes(x = Var2, y = Var1)) +
  geom_tile(aes(fill = value), colour = "black") +
  geom_text(aes(label = sprintf("%1.0f", value)), vjust = 0.5, size=6) +
  scale_fill_gradient(low = "white", high = "dodgerblue") +
  my_theme2 + theme(legend.position = "none") +
  labs(y = NULL,
       x = NULL,
       title = "Anti-SARS-CoV-2 hits")


# combine plots into one figure
library(patchwork)
layout <- A | (B/C) 

# Use plot_layout() to specify the relative heights of the rows
layout + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(1, 1))
                                