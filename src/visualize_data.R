## visualize data, metadata, ...

library(data.table)
library(dplyr)
library(ggplot2)
library(ggpubr)

load(file = "data/combined_metadata_wo_omicron.RData")
#load(file = "data/combined_metadata.RData")


# remove the last excluded samples
meta <- as.data.table(alpha_gamma_om_sero_meta[!(alpha_gamma_om_sero_meta$old_id %in% c("ID43_1st", "ID29_2nd", "ID34_2nd","ID38_1st", "ID38_2nd", "ID38_3rd", "ID53_2nd", "SRR18922948", "SRR18922909", "SRR18922902", "B_425_Seronegative", "B_436_Seronegative", "B_446_Seronegative") |
                                     alpha_gamma_om_sero_meta$sample_id %in% c("Seronegative_ID425_1st", "Seronegative_ID436_1st", "Seronegative_ID446_1st", "Omicron_--_107_1st")),])
meta <- meta[!(is.na(old_id))]
table(meta$group)


meta$num_day <- as.numeric(meta$num_day)
meta$group <- ifelse(meta$group == "Alpha EK", "Alpha+EK", meta$group)

# plot samples per variant over time +++++++++++++++++++++++++++++++++++++++++++
vertical_lines <- c(5.5, 10.5, 15.5, 30.5)
ggplot(meta[group != "Seronegative"], aes(x = num_day)) + geom_bar(stat="count") +
  geom_vline(xintercept = vertical_lines) +
  geom_text(aes(x=1,y=15),label='Day [0-5]') +
  geom_text(x=8,y=12,label='Day [6-10]', angle = 90)+
  geom_text(x=13,y=12,label='Day [11-15]', angle = 90)+
  geom_text(x=22,y=15,label='Day [16-30]')+
  geom_text(x=40,y=15,label='> Day 30')+
  facet_wrap(~group, nrow = 3)+
  labs(x = "Days", y = "Number of samples", title = "Number of samples per group and day of sampling")

meta[group=="Alpha"]


# pca ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
load("data/variants_omicron_ischgl_tpms_prepared.RData")
keep <- colnames(variants_omicron_ischgl_tpms)[colnames(variants_omicron_ischgl_tpms) %in% meta$old_id | colnames(variants_omicron_ischgl_tpms) %in% meta$sample_id]

tpms <- variants_omicron_ischgl_tpms[rowSums(variants_omicron_ischgl_tpms) != 0, keep]
pca <- prcomp(t(tpms), scale = T, center = T)
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
var_explained <- pca_summ$importance[2,] # PC1: 0.87936, PC2: 0.07033

plot_df <- merge(meta, pcs, by="sample_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(x = "PC1 (87.9 %)", y = "PC2 (7.0 %)", title = "PCA of expression values (TPM) for each group")


# color by cohort
load(file = "data/omicron_prepared.RData")
omicron1 <- colnames(omicron_tpms)
meta$cohort <- ifelse(meta$group %in% c("Alpha", "Alpha EK", "Gamma"), "Variants", 
                      ifelse(meta$group == "Seronegative", "Ischgl",
                             ifelse(meta$group == "BA.1" & meta$sample_id %in% omicron1, "Omicron 1 (BA.1 only)", "Omicron 2 (BA.1 & BA.2)")))

ggplot(plot_df[cohort!="Omicron 1 (BA.1 only)"], aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (87.9 %)", y = "PC2 (7.0 %)", title = "PCA of expression values (TPM) for each group")


# check only with one omicron cohort
load(file = "data/variants_ischgl_tpms_prepared.RData") # load file
load(file = "data/omicron_prepared.RData")
test <- cbind(variants_ischgl_tpms, omicron_tpms)

pca <- prcomp(t(test))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
var_explained <- pca_summ$importance[2,] # PC1: 0.88736 PC2 0.06423

plot_df <- merge(meta2, pcs, by.y="sample_id", by.x = "sample_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point() +
  labs(x = "PC1 (88.7 %)", y = "PC2 (6.4 %)", title = "PCA of expression values (TPM) for each group")



## 10 m ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
srr_meta <- fread("data/srr_metadata.csv")
load("data/all_tpms_gene_names_10m.RData")

meta2 <- merge(meta, srr_meta[, .(sample_id, old_id)], by = "sample_id")
meta2$old_id <- meta2$old_id.y

keep <-  meta2$old_id
tpms <- tpms_10m[, keep]
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
var_explained <- pca_summ$importance[2,] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta2, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (75.0 %)", y = "PC2 (18.7 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 10 million)")



# 50m
load("data/all_tpms_gene_names_50m.RData")
keep <- keep[!keep  %in% c("ID68_1st", "SRR13187044")]
tpms <- tpms_50m[, keep]
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
var_explained <- pca_summ$importance[2,] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta2, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (79.8 %)", y = "PC2 (14.1 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 50 million)")


# 7m
load("data/all_tpms_gene_names_7m.RData")
meta3 <- meta
meta3[, old_id := ifelse(group=="Seronegative", as.character(sample_id), old_id)]

tpms <- as.data.frame(tpms_7m)[, keep]
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
pca_summ$importance[2,c(1,2)] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta3, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (74.8 %)", y = "PC2 (19.1 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 7 million)")


# 5m
load("data/all_tpms_gene_names_5m.RData")
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
pca_summ$importance[2,c(1,2)] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta3, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (74.8 %)", y = "PC2 (18.5 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 5 million)")

# 3m
load("data/all_tpms_gene_names_3m.RData")
tpms <- as.data.frame(tpms_3m)[tpms_3m$gene_name!="", keep]
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
pca_summ$importance[2,c(1,2)] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta3, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (74.4 %)", y = "PC2 (18.5 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 3 million)")

# 1m 
load("data/all_tpms_gene_names_1m.RData")
tpms <- as.data.frame(tpms_1m)[, keep]
pca <- prcomp(t(tpms))
pcs <- as.data.frame(pca$x)
pcs$sample_id <- rownames(pcs)

pca_summ <- summary(pca)
pca_summ$importance[2,c(1,2)] # PC1: 0.74991, PC2: 0.18700

plot_df <- merge(meta3, pcs, by.y="sample_id", by.x = "old_id")

ggplot(plot_df, aes(x = PC1, y = PC2, color = cohort)) +
  geom_point() +
  labs(x = "PC1 (74.4 %)", y = "PC2 (18.5 %)", title = "PCA of expression values (TPM) for each group (sequencing depth 1 million)")

 