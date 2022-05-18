### perform immune deconvolution using the immunedeconv package

# install immunedeconv
# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)


### transform tpms to correct input matrix with HGNC symbols as row names and sample names as column names
# (might skip if data is in correct form already)

## read gene expression matrix for tpms from rds file
## immunedeconv needs only tpms as input
PATH_TO_DIR <- "/nfs/data2/covid_hennighausen/covid_pbmcs_variants/02_nfcore_rnaseq_results/01_inc_cat/star_salmon/"

rds_obj <- readRDS(paste0(PATH_TO_DIR, 'salmon.merged.gene_counts.rds'))
tpms <- assays(rds_obj)$abundance


##  immunedeconv needs HGNC names not ENSG numbers -> mapping using ensembldb package

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("EnsDb.Hsapiens.v75")

library(EnsDb.Hsapiens.v75)
edb <- EnsDb.Hsapiens.v75

## Get all genes defined in Ensembl (version 75):
ids_names <- transcripts(edb, columns=c("gene_id", "gene_name"))
mapping <- unique(cbind(gene_id=ids_names$gene_id, gene_name=ids_names$gene_name))
mapping[,1]

# get gene names for all 
mapped_names <- sapply(rownames(tpms), function(x) mapping[mapping[,1] == x, 2])


# add names to dataframes
# info: some names are duplicated mapping to RNA genes -> remove duplicated TODO: think of a smarter way!!!
tpms <- tpms[!duplicated(mapped_names),]
rownames(tpms) <- mapped_names[!duplicated(mapped_names)]

# save tpms as rdata object for easy access
save(tpms, file = 'data/variants_tpms_gene_names.RData')


################################################################################
### start here if tpm matrix is already a matrix with genes in rows and samples in columns
### with HGNC symbols as row names and samples as column names

load(file = "data/variants_tpms_gene_names.RData") # load file
head(tpms)

## perform deconvolution
# quantiseq
quantiseq_result <- immunedeconv::deconvolute(tpms, "quantiseq")

# timer
# "TIMER uses indication-specific reference profiles. Therefore, you must specify the tumor type when running TIMER" (https://omnideconv.org/immunedeconv/articles/immunedeconv.html)
# timer_result <- immunedeconv::deconvolute(tpms, "timer")

# cibersort
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")
cibersort_result <- immunedeconv::deconvolute(tpms, "cibersort")

# cibersort_abs
cibersort_abs_result <- immunedeconv::deconvolute(tpms, "cibersort_abs")

# mpc_counter
mpc_counter_result <- immunedeconv::deconvolute(tpms, "mpc_counter") # error

# xcell
xcell_result <- immunedeconv::deconvolute(tpms, "xcell")

# epic
epic_result <- immunedeconv::deconvolute(tpms, "epic")

# abis
abis_result <- immunedeconv::deconvolute(tpms, "abis")

# consensus_tme
# needs indications
# consensus_result <- immunedeconv::deconvolute(tpms, "consensus_tme") # error


################################################################################
### VISUALISATION
quantiseq_dt <- as.data.table(quantiseq_result)[, method := 'quantiseq']
cibersort_dt <- as.data.table(cibersort_result)[, method := 'cibersort']
cibersort_abs_dt <- as.data.table(cibersort_abs_result)[, method := 'cibersort_abs']
epic_dt <- as.data.table(epic_result)[, method := 'epic']
xcell_dt <- as.data.table(xcell_result)[, method := 'xcell']
abis_dt <- as.data.table(abis_result)[, method := 'abis']


result_dt <- rbindlist(list(quantiseq_dt, 
                            cibersort_dt, 
                            cibersort_abs_dt, 
                            epic_dt, 
                            xcell_dt,
                            abis_dt
                            ))

result_dt_melt <- melt(result_dt, id.vars = c("cell_type", "method"), variable.name = "sample")


mean_dt <- result_dt_melt[, .(mean_per_type = mean(value)), by = c("cell_type", "method")]

## ggplot approach
# ggplot(mean_dt, aes(x = "", y = mean_per_type, fill = cell_type)) +
#   geom_col(width=1) + coord_polar("y", start=0) + facet_wrap(~method)

mean_dt[, mean_per_type := ifelse(mean_per_type < 0, 0, mean_per_type)]

# par(mfrow=c(3,2))
# par(mfrow=c(1,1))
pie(mean_dt[method == "quantiseq"]$mean_per_type, labels = mean_dt[method == "quantiseq"]$cell_type) 
pie(mean_dt[method == "cibersort"]$mean_per_type, labels = mean_dt[method == "cibersort"]$cell_type)
pie(mean_dt[method == "cibersort_abs"]$mean_per_type, labels = mean_dt[method == "cibersort_abs"]$cell_type)
pie(mean_dt[method == "epic"]$mean_per_type, labels = mean_dt[method == "epic"]$cell_type)
pie(mean_dt[method == "xcell"]$mean_per_type, labels = mean_dt[method == "xcell"]$cell_type)
pie(mean_dt[method == "abis"]$mean_per_type, labels = mean_dt[method == "abis"]$cell_type)


# TODO: 
#   - add info about variants
#   - smarter way of plotting?
#   - benchmarking pipeline






#### notes
methods <- unique(mean_dt$method)

pie_chart <- function(x){
  pie(mean_dt[method == x]$mean_per_type, labels = mean_dt[method == x]$cell_type)
} 

pie_charts <- lapply(methods, function(x){
  plot <- function(c){pie_chart(x)}
  plot
})
plot <- function(x){pie_chart(x)}
plot("quantiseq")