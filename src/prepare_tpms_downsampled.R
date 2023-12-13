# script to prepare the tpms matrix from the downsampled matrices

library(data.table)
library(tidyr)
library(tidyverse)
library(dplyr)

## bind all tpms from all samples together
# all_tpms <- bind_cols(lapply(c("alpha", "alpha_ek", "gamma", "omikron", "sero"), function (variant) {
#   files <- list.files(paste0("data/seq_depth_tpms/", variant, "/"), full.names = T)
#   variant_df <- bind_cols(lapply(files, function(file) {
#     tab <- as.data.frame(fread(file, select = c("Name", "TPM")))
#     rownames(tab) <- tab$Name
#     tab$Name<- NULL
#     n <- unlist(strsplit(basename(file), "_"))
#     id <- ifelse(variant %in% c("omikron", "sero"), n[[1]], paste0(n[[1]], "_", unlist(strsplit(basename(file), "_"))[[2]]))
#     #print(id)
#     colnames(tab) <- id
#     tab
#   }))
#   variant_df
# }))


files <- list.files(paste0("data/seq_depth_tpms/M7/"), full.names = T)
df <- bind_cols(lapply(files, function(file) {
    tab <- as.data.frame(fread(file, select = c("Name", "TPM")))
    rownames(tab) <- tab$Name
    tab$Name<- NULL
    n <- unlist(strsplit(basename(file), "_"))
    #print(n)
    id <- ifelse(any(startsWith(n, "SRR")), n[[1]], 
                 ifelse("Seronegative" %in% n, paste0(n[[1]], "_", unlist(strsplit(basename(file), "_"))[[2]], "_", unlist(strsplit(basename(file), "_"))[[3]]), 
                        paste0(n[[1]], "_", unlist(strsplit(basename(file), "_"))[[2]])))
    #print(id)
    colnames(tab) <- id
    tab
  }))
df

all_tpms <- df 

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

all_tpms_sum <- as.data.table(all_tpms_sum)



# map seronegative sample names to reaonable names
load("data/combined_metadata.RData")
meta <- alpha_gamma_om_sero_meta


# Loop through each column name in df
for (col in colnames(all_tpms_sum)) {
  if(any(!is.na(grep("Seronegative", col)))){
    # Find the matching row in meta$A
    match_row = match(col, meta$old_id)
    # print(match_row)
    
    # If a match is found, replace the column name
    if (!is.na(match_row)) {
      print(meta$sample_id[match_row])
      colnames(all_tpms_sum)[colnames(all_tpms_sum) == col] <- as.character(meta$sample_id[match_row])
    }
  }
}

all_tpms_sum <- all_tpms_sum[-1,] # remove "" gene
tpms_7m <- all_tpms_sum


# save tpms as rdata object for easy access
save(tpms_7m, file = 'data/all_tpms_gene_names_7m.RData')



#### remove omicron 1 cohort
load(file = "data/combined_metadata.RData")
# remove the last excluded samples
meta <- as.data.table(alpha_gamma_om_sero_meta[!(alpha_gamma_om_sero_meta$old_id %in% c("ID43_1st", "ID29_2nd", "ID34_2nd","ID38_1st", "ID38_2nd", "ID38_3rd", "ID53_2nd", "SRR18922948", "SRR18922909", "SRR18922902", "B_425_Seronegative", "B_436_Seronegative", "B_446_Seronegative") |
                                                   alpha_gamma_om_sero_meta$sample_id %in% c("Seronegative_ID425_1st", "Seronegative_ID436_1st", "Seronegative_ID446_1st", "Omicron_--_107_1st")),])
meta <- meta[!(is.na(old_id))]
table(meta$group)

meta$num_day <- as.numeric(meta$num_day)
meta$group <- ifelse(meta$group == "Alpha_EK", "Alpha EK", meta$group)

meta[group=="Seronegative"]$old_id <- meta[group=="Seronegative"]$sample_id
load(file = "data/omicron_prepared.RData")
omicron1 <- colnames(omicron_tpms)
meta$cohort <- ifelse(meta$group %in% c("Alpha", "Alpha EK", "Gamma"), "Variants", 
                      ifelse(meta$group == "Seronegative", "Ischgl",
                             ifelse(meta$group == "BA.1" & meta$sample_id %in% omicron1, "Omicron 1 (BA.1 only)", "Omicron 2 (BA.1 & BA.2)")))
meta <- meta[cohort != "Omicron 1 (BA.1 only)"]


# load("data/all_tpms_gene_names_1m.RData")
# tpms_1m_wo_omicron1 <- tpms_1m %>% select(all_of(meta$old_id))
# save(tpms_1m_wo_omicron1, file = "data/all_tpms_gene_names_1m_wo_omicron1.RData")
# 
# load("data/all_tpms_gene_names_3m.RData")
# tpms_3m_wo_omicron1 <- tpms_3m %>% select(all_of(meta$old_id))
# save(tpms_3m_wo_omicron1, file = "data/all_tpms_gene_names_3m_wo_omicron1.RData")
# 
# load("data/all_tpms_gene_names_5m.RData")
# tpms_5m_wo_omicron1 <- tpms_5m %>% select(all_of(meta$old_id))
# save(tpms_5m_wo_omicron1, file = "data/all_tpms_gene_names_5m_wo_omicron1.RData")
# 
# load("data/all_tpms_gene_names_7m.RData")
# tpms_7m_wo_omicron1 <- tpms_7m %>% select(all_of(meta$old_id))
# save(tpms_7m_wo_omicron1, file = "data/all_tpms_gene_names_7m_wo_omicron1.RData")
