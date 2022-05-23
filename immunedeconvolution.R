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

methods = list("quantiseq",
            "cibersort_abs",
            "xcell",
            "epic" #, "mpc_counter"
            )

## perform deconvolution
set_cibersort_binary("CIBERSORT.R") # make sure cibersort is known
set_cibersort_mat("LM22.txt")

# deconv_results <- lapply(methods, function(x) deconvolute(tpms, x))


## perform each method alone
# quantiseq
quantiseq_result <- immunedeconv::deconvolute(tpms, "quantiseq")

# set configurations for cibersort
set_cibersort_binary("CIBERSORT.R")
set_cibersort_mat("LM22.txt")

# cibersort_abs
cibersort_abs_result <- immunedeconv::deconvolute(tpms, "cibersort_abs")

# mpc_counter
mpc_counter_result <- immunedeconv::deconvolute_mcp_counter(tpms) # somehow works only this way

# xcell
xcell_result <- immunedeconv::deconvolute(tpms, "xcell")

# epic
epic_result <- immunedeconv::deconvolute(tpms, "epic")



################################################################################
### VISUALISATION

source('plotting.R') # script and method for visualization of deconvolution results with absolute scores

create_score_plots("xcell", "plots/xcell_plots/")
create_score_plots("mpc_counter", "plots/mpc_counter_plots/")
create_score_plots("cibersort_abs", "plots/cibersort_abs_plots/")

create_fraction_plot("quantiseq", "plots/")
create_fraction_plot("epic", "plots/")



# TODO: 
#   - add info about variants
#   - smarter way of plotting?
#   - benchmarking pipeline




## old?

mpc_cell_types <- rownames(mpc_counter_result)
mpc_counter_result_dt <- as.data.table(mpc_counter_result)
mpc_counter_result_dt[, cell_type := mpc_cell_types]

### try visualization as in https://omnideconv.org/immunedeconv/articles/detailed_example.html
quantiseq_result %>%
  gather(sample, fraction, -cell_type) %>%
  # plot as stacked bar chart
  ggplot(aes(x=sample, y=fraction, fill=cell_type)) +
  geom_bar(stat='identity') +
  coord_flip() +
  scale_fill_brewer(palette="Paired") +
  scale_x_discrete(limits = rev(levels(quantiseq_result)))

xcell_result_dt <- as.data.table(xcell_result)
xcell_result_dt[cell_type %in% c("B cell", "B cell memory", "B cell naive", "T cell CD8+", "T cell CD4+ memory", "T cell CD4+ naive")] %>%
  gather(sample, score, -cell_type) %>%
  ggplot(aes(x=sample, y=score, color=cell_type)) +
  geom_point(size=4) +
  facet_wrap(~cell_type, scales="free_x", ncol=3) +
  scale_color_brewer(palette="Paired", guide=FALSE) +
  coord_flip() +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))





