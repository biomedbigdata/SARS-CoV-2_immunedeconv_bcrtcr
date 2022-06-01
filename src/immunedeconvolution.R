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
library(patchwork)


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
save(tpms, file = '../data/variants_tpms_gene_names.RData')


################################################################################
### start here if tpm matrix is already a matrix with genes in rows and samples in columns
### with HGNC symbols as row names and samples as column names

load(file = "../data/variants_tpms_gene_names.RData") # load file
head(tpms)

methods = list("quantiseq",
            "cibersort_abs",
            "xcell",
            "epic" #, "mcp_counter"
            )


## perform deconvolution
# deconv_results <- lapply(methods, function(x) deconvolute(tpms, x))


## perform each method alone
# quantiseq
quantiseq_result <- immunedeconv::deconvolute(tpms, "quantiseq", tumor = FALSE)

# set configurations for cibersort
set_cibersort_binary("cibersort/CIBERSORT.R")
set_cibersort_mat("cibersort/LM22.txt")

# cibersort_abs
cibersort_abs_result <- immunedeconv::deconvolute(tpms, "cibersort_abs")

# mcp_counter
mcp_counter_result <- immunedeconv::deconvolute_mcp_counter(tpms) # somehow works only this way

# xcell
xcell_result <- immunedeconv::deconvolute(tpms, "xcell")

# epic
epic_result <- immunedeconv::deconvolute(tpms, "epic", tumor = FALSE)



################################################################################
### VISUALISATION
# configurations for cibersort
set_cibersort_binary("cibersort/CIBERSORT.R")
set_cibersort_mat("cibersort/LM22.txt")
source('plotting.R') # script and method for visualization of deconvolution results with absolute scores

# load data if not already loaded into work space
if (!exists("tpms")){
  load(file = "data/variants_tpms_gene_names.RData") # load file
}

# load meta data about samples
if (!exists("full_metadata")){
  full_metadata <- fread("data/full_pbmc_metadata_original.csv")
}

create_score_plots("xcell", "plots/xcell_plots/")
create_score_plots("cibersort_abs", "plots/cibersort_abs_plots/")
create_score_plots("mcp_counter", "plots/mcp_counter_plots/")

create_fraction_plot("quantiseq", "plots/")
create_fraction_plot("epic", "plots/")




# TODO: 
#   - add color bar annotation to fraction plots
#   - benchmarking pipeline

