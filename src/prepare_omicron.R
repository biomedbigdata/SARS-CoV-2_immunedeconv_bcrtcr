# script to prepare omicron data

library(data.table)
library(tidyr)
library(dplyr)
library(stringr)

## prepare samplesheet
omicron_samples <- fread("~/SARS-CoV-2_immunedeconv_bcrtcr/data/omikron_samplesheet.csv")
head(omicron_samples)

# keep only sample name
omicron_samples <- omicron_samples[, .(sample)]

# split samplename into SRR number and group
omicron_samples <- separate(omicron_samples, 
                         col = "sample",
                         into = c("sample", "group"),
                         sep = "_")

# renaming the groups
omicron_samples[, group := ifelse(group == "COVID-19-Omicron", "Omicron-BA-1", "Healthy")]

# create new IDs (BA1_#)
omicron_samples <- omicron_samples %>% mutate(across("ID", str_replace, "Om", "BA1_"))

# save 
write.table(omicron_samples, file = "data/omicron_samples.tsv", quote = F, 
            sep = "\t", row.names = F)


## prepare tpms and counts
omicron_tpms <- fread("~/SARS-CoV-2_immunedeconv_bcrtcr/data/omikron_tpms.tsv")
omicron_gene_counts <- fread("~/SARS-CoV-2_immunedeconv_bcrtcr/data/omikron_geneCounts.tsv")

omicron_tpms <- as.data.frame(omicron_tpms[!duplicated(gene_name),])
rownames(omicron_tpms) <- omicron_tpms$gene_name
omicron_tpms <- subset(omicron_tpms, select = -c(gene_id, gene_name))

omicron_gene_counts <- as.data.frame(omicron_gene_counts[!duplicated(gene_name),])
rownames(omicron_gene_counts) <- omicron_gene_counts$gene_name
omicron_gene_counts <- subset(omicron_gene_counts, select = -c(gene_id, gene_name))


# renaming
# metadata 
load("data/all_pbmc_metadata.RData") # full_metadata

colnames(omicron_tpms) <- sapply(colnames(omicron_tpms), function(srr){
  raw_srr <- unlist(str_split(srr, "_")[1])
  unique(full_metadata[old_id == raw_srr[1]]$sample_id)[1]
})



# remove Naive samples
omicron_tpms <- omicron_tpms[!is.na(names(omicron_tpms))]


# save prepared data frames
save(omicron_tpms, omicron_gene_counts, file = "data/omicron_prepared.RData")


# add second omicron study
# load("data/variants_omicron_ischgl_tpms_prepared.RData")
omicron2_tpms <- fread("/nfs/data2/covid_hennighausen/omikron_new/output/output/org_exbio_sradownloader_rnaSeq_NfCoreRnaSeq/output/results/star_salmon/salmon.merged.gene_tpm.tsv")


omicron2_tpms <- as.data.frame(omicron2_tpms[!duplicated(gene_name),])
rownames(omicron2_tpms) <- omicron2_tpms$gene_name
omicron2_tpms <- subset(omicron2_tpms, select = -c(gene_id, gene_name))

# metadata 
load("data/all_pbmc_metadata.RData") # full_metadata

# omicron2_meta <- fread("/nfs/data2/covid_hennighausen/omikron_new/_meta/SraRunTable.txt")

colnames(omicron2_tpms) <- sapply(colnames(omicron2_tpms), function(srr){
  unique(full_metadata[old_id == srr]$sample_id)[1]
})




## combine with variants/ischgl df
load(file = "data/variants_ischgl_tpms_prepared.RData") # load file
load(file = "data/omicron_prepared.RData")


variants_ischgl_tpms

variants_omicron_ischgl_tpms <- cbind(variants_ischgl_tpms, omicron_tpms, omicron2_tpms)
save(variants_omicron_ischgl_tpms, file = "data/variants_omicron_ischgl_tpms_prepared.RData")











## prepare batch corrected tpms
bc_tpms <- fread("data/batch_corrected/variants_omicron_ischgl_batch_corrected.csv")
samples <- colnames(bc_tpms)
new_colnames <- c("gene_id",samples[1:178])
bc_tpms <- bc_tpms[, c(2, 4:181)]
colnames(bc_tpms) <- new_colnames

# for gene names
omicron_raw_tpms <- fread("data/omikron_tpms.tsv")
bc_tpms$gene_name <- omicron_raw_tpms$gene_name

# remove duplicated gene names
bc_tpms <- as.data.frame(bc_tpms[!duplicated(gene_name),])
rownames(bc_tpms) <- bc_tpms$gene_name
bc_tpms <- subset(bc_tpms, select = -c(gene_id, gene_name))

# save prepared data
save(bc_tpms, file="data/batch_corrected/variants_omicron_ischgl_bc_prepared.RData")

