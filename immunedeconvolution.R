## perform immune deconvolution using the immunedeconv package

# install immunedeconv
# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")

library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
library(immunedeconv)
library(tibble)

# read gene expression matrix for counts and tpms from rds file

PATH_TO_DIR <- "/nfs/data2/covid_hennighausen/covid_pbmcs_variants/02_nfcore_rnaseq_results/01_inc_cat/star_salmon/"

rds_obj <- readRDS(paste0(PATH_TO_DIR, 'salmon.merged.gene_counts.rds'))
counts <- assays(rds_obj)$counts
tpms <- assays(rds_obj)$abundance

source("CIBERSORT.R")
immunedeconv::deconvolute(counts, "quantiseq")
