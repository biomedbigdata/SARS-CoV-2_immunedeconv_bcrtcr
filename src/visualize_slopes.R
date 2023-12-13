# visualize slope p values as heatmaps

library(data.table)
library(pheatmap)

slopes_full <- as.data.frame(fread("results/slopes_p_values.csv"))
rownames(slopes_full) <- paste0(slopes_full$Method, ", ", slopes_full$Cell_type)

mat <- slopes_full %>% select(!Method & !Cell_type)
rownames(mat) <- paste0(slopes_full$Method, ", ", slopes_full$Cell_type)
mat <- apply(mat, 2, as.numeric)

slopes_50m <- as.data.frame(fread("results/slopes_p_values_m50.csv"))
slopes_10m <- as.data.frame(fread("results/slopes_p_values_m10.csv"))
mat_50 <- slopes_50m %>% select(!Method & !Cell_type) %>% apply(., 2, as.numeric)
mat_10 <- slopes_10m %>% select(!Method & !Cell_type) %>% apply(., 2, as.numeric)

mat_total <- cbind(mat, mat_50, mat_10)
rownames(mat_total) <- paste0(slopes_full$Method, ", ", slopes_full$Cell_type)

annotation_df <- slopes_full %>% select(Method)
annotation_c <- data.frame(
  Sequencing_depth = factor(c(rep("full", 4), rep("50 million", 4), rep("10 million", 4)))
)

mat_transformed <- -log10(mat_total)

pheatmap(mat_transformed, cluster_rows = F,
         cluster_cols = F, 
         legend = T,
         annotation_row = annotation_df,
         annotation_names_row = F,
         gaps_row = c(4,8,12),
         gaps_col = c(4,8,12),
         labels_row = c(rep(c("B cell", "Neutrophil", "T cell CD4+", "T cell CD8+"), 4)),
         angle_col = 0,
         main = "P-values (-log10) for slopes")
