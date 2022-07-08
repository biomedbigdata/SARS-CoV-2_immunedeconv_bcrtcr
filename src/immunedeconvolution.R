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
library(RColorBrewer)


################################################################################
### prepare the data

## transform tpms to correct input matrix with HGNC symbols as row names and sample names as column names
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
# save(tpms, file = '../data/variants_tpms_gene_names.RData')





### prepare data for variants, ischgl and nuns data from csv files
variants_ischgl_tpms <- read.csv("data/variants_ischgl_pbmcs_tpms.csv")
nuns_tpms <- read.csv("data/nuns_pbmcs_tpms.csv")

# remove gene_id, delete duplicated genes and make gene_name to rownames
nuns_tpms <- nuns_tpms[!duplicated(nuns_tpms$gene_name),]
rownames(nuns_tpms) <- nuns_tpms$gene_name
nuns_tpms <- nuns_tpms[, !names(nuns_tpms) %in% c("X", "gene_id", "gene_name")]
save(nuns_tpms, file = "data/nuns_tpms_prepared.RData")

# do the same for variants and ischgl set
variants_ischgl_tpms <- variants_ischgl_tpms[!duplicated(variants_ischgl_tpms$gene_name),]
rownames(variants_ischgl_tpms) <- variants_ischgl_tpms$gene_name
variants_ischgl_tpms <- variants_ischgl_tpms[, !names(variants_ischgl_tpms) %in% c("X", "gene_id", "gene_name", "gene_id.1", "gene_name.1")]

# filter variants_ischgl dataframe for desired samples (only Seronegative or variant)
variants_ischgl_tpms <- cbind(variants_ischgl_tpms[, startsWith(names(variants_ischgl_tpms), "Alpha")],
                                variants_ischgl_tpms[, startsWith(names(variants_ischgl_tpms), "Gamma")],
                                variants_ischgl_tpms[, startsWith(names(variants_ischgl_tpms), "Seronegative")])

save(variants_ischgl_tpms, file = "data/variants_ischgl_tpms_prepared.RData")


# add something to metadata
load(file = "data/all_pbmc_metadata.RData") # loaa metadata
# add day group column
full_metadata[, day_group := ifelse(num_day <= 11, "<= day 11", "> day 11")]
save(full_metadata, file = "data/all_pbmc_metadata.RData")



################################################################################
### start here if tpm matrix is already a matrix with genes in rows and samples in columns
### with HGNC symbols as row names and samples as column names

load(file = "data/variants_ischgl_tpms_prepared.RData") # load file
load(file = "data/nuns_tpms_prepared.RData") # load file
load(file = "data/all_pbmc_metadata.RData") # load metadata
head(nuns_tpms)
head(variants_ischgl_tpms)


## set the output directory and which dataset to use
result_plotting_dir <- "plots/nuns/all_by_group/boxplots/"
dir.create(file.path(result_plotting_dir))
dataset <- "nuns"  # choose one of "nuns", "nuns_naive", "nuns_convalescent", "variants_ischgl"


# set tpms to correct data and whether to use the old id or the new id
if(dataset == "nuns"){
  tpms <- nuns_tpms
} else if (dataset == "nuns_naive") {
  tpms <- nuns_tpms[, startsWith(names(nuns_tpms), "BNT_Naive")]
} else if (dataset == "nuns_convalescent") {
  tpms <- nuns_tpms[, startsWith(names(nuns_tpms), "BNT_Convalescent")]
} else if (dataset == "variants_ischgl") {
  tpms <- variants_ischgl_tpms
}
old_id <- ifelse(startsWith(dataset,"nuns"), TRUE, FALSE)


## perform each of the deconvolution methods 
# quantiseq
quantiseq_result <- immunedeconv::deconvolute(tpms, "quantiseq", tumor = FALSE)

# set configurations for cibersort
set_cibersort_binary("src/cibersort/CIBERSORT.R")
set_cibersort_mat("src/cibersort/LM22.txt")

# cibersort_abs
cibersort_abs_result <- immunedeconv::deconvolute(tpms, "cibersort_abs", tumor = FALSE)

# mcp_counter
mcp_counter_result <- immunedeconv::deconvolute_mcp_counter(tpms) # somehow works only this way

# xcell
xcell_result <- immunedeconv::deconvolute(tpms, "xcell", tumor = FALSE)

# epic
epic_result <- immunedeconv::deconvolute(tpms, "epic", tumor = FALSE)



################################################################################
### VISUALISATION

color_by_sampling = F # parameter whether to color by group (F) or by sampling (T) (for nuns color sampling == by day)
result_plotting_dir <- "plots/nuns/all_by_day/boxplots/"


# configurations for cibersort
set_cibersort_binary("src/cibersort/CIBERSORT.R")
set_cibersort_mat("src/cibersort/LM22.txt")
source('src/plotting.R') # script and method for visualization of deconvolution results with absolute scores


# load meta data about samples
if (!exists("full_metadata")){
  full_metadata <- fread("data/all_pbmc_metadata.csv")
}

# stacked bar plots
create_fraction_plot("quantiseq", result_plotting_dir)
create_fraction_plot("epic", result_plotting_dir)

# scatter plots
create_score_plots("xcell", paste0(result_plotting_dir, "xcell_plots/"), color_by_sampling = color_by_sampling)
create_score_plots("cibersort_abs", paste0(result_plotting_dir,"cibersort_abs_plots/"), color_by_sampling = color_by_sampling)
create_score_plots("mcp_counter", paste0(result_plotting_dir,"mcp_counter_plots/"), color_by_sampling = color_by_sampling)

# heatmaps
create_hm("epic", result_plotting_dir, color_by_sampling = color_by_sampling)
create_hm("xcell", result_plotting_dir, color_by_sampling = color_by_sampling)
create_hm("cibersort_abs", result_plotting_dir, color_by_sampling = color_by_sampling)
create_hm("mcp_counter", result_plotting_dir, color_by_sampling = color_by_sampling)
create_hm("quantiseq", result_plotting_dir, color_by_sampling = color_by_sampling)

# boxplots
color = "group"
reference = ".all." # for the reference group in variants data set
create_boxplot("epic", paste0(result_plotting_dir, "epic/"), color, reference)
create_boxplot("xcell", paste0(result_plotting_dir, "xcell/"), color, reference)
create_boxplot("cibersort_abs", paste0(result_plotting_dir, "cibersort_abs/"), color, reference)
create_boxplot("mcp_counter", paste0(result_plotting_dir, "mcp_counter/"), color, reference)
create_boxplot("quantiseq", paste0(result_plotting_dir, "quantiseq/"), color, reference)

# boxplots for time series
color = "day_group"
variants = dataset == "variants_ischgl" # check which dataset is used, variants needs different metadata
create_boxplot("epic", paste0(result_plotting_dir, "epic/"), color = color, time = T, variants = variants)
create_boxplot("xcell", paste0(result_plotting_dir, "xcell/"), color = color, time = T, variants = variants)
create_boxplot("cibersort_abs", paste0(result_plotting_dir, "cibersort_abs/"), color = color, time = T, variants = variants)
create_boxplot("mcp_counter", paste0(result_plotting_dir, "mcp_counter/"), color = color, time = T, variants = variants)
create_boxplot("quantiseq", paste0(result_plotting_dir, "quantiseq/"), color = color, time = T, variants = variants)




# plot distribution of days for both datasets
load("data/alpha_gamma_time_meta.RData")
ggplot(alpha_gamma_meta, aes(x = num_day, fill = sampling)) +
  geom_histogram(color = "grey", bins = max(alpha_gamma_meta$num_day)) +
  scale_fill_grey()+
  labs(title = "Distribution of days when samples were taken for variants_ischgl data")
ggsave("plots/variants_ischgl_distribution_days.png")

ggplot(full_metadata[group %in% c("Naive", "Convalescent")], aes(x = num_day)) +
  geom_histogram(color = "grey", bins = max(full_metadata[group %in% c("Naive", "Convalescent")]$num_day)) +
  scale_fill_grey()+
  labs(title = "Distribution of days when samples were taken for nuns data")
ggsave("plots/nuns_distribution_days.png")

ggplot(full_metadata[group == "Naive"], aes(x = num_day)) +
  geom_histogram(color = "grey", bins = max(full_metadata[group %in% c("Naive", "Convalescent")]$num_day)) +
  scale_fill_grey()+
  labs(title = "Distribution of days when samples were taken for nuns naive data")
ggsave("plots/nuns_naive_distribution_days.png")

ggplot(full_metadata[group == "Convalescent"], aes(x = num_day)) +
  geom_histogram(color = "grey", bins = max(full_metadata[group %in% c("Naive", "Convalescent")]$num_day)) +
  scale_fill_grey()+
  labs(title = "Distribution of days when samples were taken for nuns naive data")
ggsave("plots/nuns_convalescent_distribution_days.png")


