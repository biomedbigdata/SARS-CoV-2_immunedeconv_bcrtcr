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

# read gene expression matrix and transform to a matrix with genes in rows and 
# samples in columns

PATH_TO_EXPRESSION_MAT <- "/nfs/data2/covid_hennighausen/covid_pbmcs_variants/02_nfcore_rnaseq_results/01_inc_cat/star_salmon/salmon.merged.transcript_tpm.tsv"
expr_mat <- fread(PATH_TO_EXPRESSION_MAT)

# remove first column (double the tx)
# TODO: have a look at the files, which ID to use
expr_mat <- expr_mat[, -c("V1")] 
expr_mat <- as.matrix(expr_mat)
expr_mat[1:6, 1:6]
