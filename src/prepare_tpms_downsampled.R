# script to prepare the tpms matrix from the downsampled matrices

library(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)

## bind all tpms from all samples together
all_tpms <- bind_cols(lapply(c("alpha", "alpha_ek", "gamma", "omikron", "sero"), function (variant) {
  files <- list.files(paste0("data/seq_depth_tpms/", variant, "/"), full.names = T)
  variant_df <- bind_cols(lapply(files, function(file) {
    tab <- as.data.frame(fread(file, select = c("Name", "TPM")))
    rownames(tab) <- tab$Name
    tab$Name<- NULL
    n <- unlist(strsplit(basename(file), "_"))
    id <- ifelse(variant %in% c("omikron", "sero"), n[[1]], paste0(n[[1]], "_", unlist(strsplit(basename(file), "_"))[[2]]))
    #print(id)
    colnames(tab) <- id
    tab
  }))
  variant_df
}))

## mapping transcript ids to gene ids
library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## Get all genes defined in Ensembl (version 75):
ids_names <- transcripts(edb, columns=c("tx_id", "gene_id", "gene_name"))
mapping <- unique(cbind(tx_id=ids_names$tx_id, gene_id=ids_names$gene_id, gene_name=ids_names$gene_name))
mapping[,1]
rownames(mapping) <- mapping[, 1]

# get gene names for all transcript ids
# mapped_names <- sapply(rownames(all_tpms), function(x) {
#   y <- unlist(strsplit(x, ".", fixed = T))[[1]]
#   mapping[mapping[,1] == y, 3]
# })
all_tpms$tx_id <- sapply(rownames(all_tpms), function(x) unlist(strsplit(x, ".", fixed=T))[[1]])

indices <- match(all_tpms$tx_id, rownames(mapping))
all_tpms$gene_name <- ifelse(is.na(indices), "", mapping[indices, "gene_name"])
all_tpms$tx_id <- NULL

all_tpms_sum <- all_tpms %>%
  group_by(gene_name) %>%
  summarise_all(sum) 

tpms_10m <- all_tpms_sum

# save tpms as rdata object for easy access
# save(tpms_10m, file = 'data/all_tpms_gene_names_10m.RData')
