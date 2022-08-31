###
# script for 
#   1. computing immunedeconv results and combining results from different methods into one dt
#   2. adding metadata to results

load(file = "data/variants_ischgl_tpms_prepared.RData") # load file
load(file = "data/nuns_tpms_prepared.RData") # load file
head(nuns_tpms)
head(variants_ischgl_tpms)



## perform each of the deconvolution methods 
# set configurations for cibersort
set_cibersort_binary("src/cibersort/CIBERSORT.R")
set_cibersort_mat("src/cibersort/LM22.txt")

# quantiseq
quantiseq_result <- as.data.table(immunedeconv::deconvolute(variants_ischgl_tpms, "quantiseq", tumor = FALSE))

# cibersort_abs
cibersort_abs_result <- as.data.table(immunedeconv::deconvolute(variants_ischgl_tpms, "cibersort_abs", tumor = FALSE))

# mcp_counter
mcp_counter_result <- immunedeconv::deconvolute_mcp_counter(variants_ischgl_tpms) # somehow works only this way
cell_types <- rownames(mcp_counter_result)
mcp_counter_result <- as.data.table(mcp_counter_result)
mcp_counter_result[, cell_type := cell_types]

# xcell
xcell_result <- as.data.table(immunedeconv::deconvolute(variants_ischgl_tpms, "xcell", tumor = FALSE))

# epic
epic_result <- as.data.table(immunedeconv::deconvolute(variants_ischgl_tpms, "epic", tumor = FALSE))


quantiseq_result[, method := "quantiseq"]
cibersort_abs_result[, method := "cibersort_abs"]
mcp_counter_result[, method := "mcp_counter"]
xcell_result[, method := "xcell"]
epic_result[, method := "epic"]



# combine and save results
variants_ischgl_results <- rbindlist(list(
  quantiseq_result, 
  cibersort_abs_result,
  mcp_counter_result,
  xcell_result,
  epic_result
), use.names = T)



## same for nuns
# quantiseq
quantiseq_result <- as.data.table(immunedeconv::deconvolute(nuns_tpms, "quantiseq", tumor = FALSE))

# cibersort_abs
cibersort_abs_result <- as.data.table(immunedeconv::deconvolute(nuns_tpms, "cibersort_abs", tumor = FALSE))

# mcp_counter
mcp_counter_result <- immunedeconv::deconvolute_mcp_counter(nuns_tpms) # somehow works only this way
cell_types <- rownames(mcp_counter_result)
mcp_counter_result <- as.data.table(mcp_counter_result)
mcp_counter_result[, cell_type := cell_types]
mcp_counter_result <- as.data.table(mcp_counter_result)

# xcell
xcell_result <- as.data.table(immunedeconv::deconvolute(nuns_tpms, "xcell", tumor = FALSE))

# epic
epic_result <- as.data.table(immunedeconv::deconvolute(nuns_tpms, "epic", tumor = FALSE))

quantiseq_result[, method := "quantiseq"]
cibersort_abs_result[, method := "cibersort_abs"]
mcp_counter_result[, method := "mcp_counter"]
xcell_result[, method := "xcell"]
epic_result[, method := "epic"]

# combine and save results
nuns_results <- rbindlist(list(
  quantiseq_result, 
  cibersort_abs_result,
  mcp_counter_result,
  xcell_result,
  epic_result
), use.names = T)




# melt data tables
variants_ischgl_results <- melt(variants_ischgl_results, 
                                id.vars = c("cell_type", "method"),
                                variable.name = "sample",
                                value_name = "score")

nuns_results <- melt(nuns_results, 
                                id.vars = c("cell_type", "method"),
                                variable.name = "sample",
                                value_name = "score")



################################################################################
### add cibersortx results (result generated at https://cibersortx.stanford.edu))

## variants_ischgl
cibersortx_rel_result <- fread("cibersortx/CIBERSORTx_variants_ischgl_results_relative.csv")
cibersortx_abs_result <- fread("cibersortx/CIBERSORTx_variants_ischgl_results_absolute.csv")

# melt for long data table
cibersortx_rel_result <- melt(cibersortx_rel_result, id.vars = c("Mixture"), variable.name = "cell_type", value.name = "value")
cibersortx_abs_result <- melt(cibersortx_abs_result, id.vars = c("Mixture"), variable.name = "cell_type", value.name = "value")

# rename column "Mixture" to "sample"
names(cibersortx_rel_result)[names(cibersortx_rel_result) == 'Mixture'] <- 'sample'
names(cibersortx_abs_result)[names(cibersortx_abs_result) == 'Mixture'] <- 'sample'

# add method column
cibersortx_rel_result[, method := "cibersortx_rel"]
cibersortx_abs_result[, method := "cibersortx_abs"]

# bind cibersortx results to variants_ischgl results
variants_ischgl_results <- rbindlist(list(variants_ischgl_results,  cibersortx_rel_result, cibersortx_abs_result), use.names = T)
variants_ischgl_results$value <- as.numeric(variants_ischgl_results$value)



## nuns
cibersortx_rel_result <- fread("cibersortx/CIBERSORTx_nuns_results_relative.csv")
cibersortx_abs_result <- fread("cibersortx/CIBERSORTx_nuns_results_absolute.csv")

# melt for long data table
cibersortx_rel_result <- melt(cibersortx_rel_result, id.vars = c("Mixture"), variable.name = "cell_type", value.name = "value")
cibersortx_abs_result <- melt(cibersortx_abs_result, id.vars = c("Mixture"), variable.name = "cell_type", value.name = "value")

# rename column "Mixture" to "sample"
names(cibersortx_rel_result)[names(cibersortx_rel_result) == 'Mixture'] <- 'sample'
names(cibersortx_abs_result)[names(cibersortx_abs_result) == 'Mixture'] <- 'sample'

# add method column
cibersortx_rel_result[, method := "cibersortx_rel"]
cibersortx_abs_result[, method := "cibersortx_abs"]

# bind cibersortx results to variants_ischgl results
nuns_results <- rbindlist(list(nuns_results,  cibersortx_rel_result, cibersortx_abs_result), use.names = T)
nuns_results$value <- as.numeric(nuns_results$value)

################################################################################


# combine with metadata
load(file = "data/all_pbmc_metadata.RData") # load metadata
load(file = "data/alpha_gamma_sero_meta.RData") # load time series metadata for variants

variants_ischgl_results <- merge(variants_ischgl_results, alpha_gamma_sero_meta,
                                 by.x = "sample", by.y = "sample_id", all.x = T, allow.cartesian = T)


nuns_results <- merge(nuns_results, full_metadata, by.x = "sample", by.y = "old_id", 
                      all.x = T, allow.cartesian = T)

# save the results as RData
# save(variants_ischgl_results, file = "data/variants_ischgl_deconv.RData")
# save(nuns_results, file = "data/nuns_deconv.RData")

# test loading
load("data/variants_ischgl_deconv.RData")
load("data/nuns_deconv.RData")



### other additional info and renaming for plotting

# rename day group groups in variants set
variants_ischgl_results[, day_group := ifelse(day_group == "<= day 5", "day [0,5]",
                                              ifelse(day_group == "<= day 10", "day [6,10]",
                                                     ifelse(day_group == "<= day 15", "day [11,15]",
                                                            ifelse(day_group == "<= day 30", "day [16,30]", "> day 30"))))]
# save(variants_ischgl_results, file = "data/variants_ischgl_deconv.RData")



# add controlled vocabulary (cv) for the (important) cell types
variants_ischgl_results$cell_type <- as.character(variants_ischgl_results$cell_type)
variants_ischgl_results[, cv_cell_type := ifelse(startsWith(cell_type, "T cell CD4"), "T cell CD4+",
                                                 ifelse(startsWith(cell_type, "T cells CD4"), "T cell CD4+",
                                                 ifelse(cell_type == "T cells", "T cell CD4+",
                                                 ifelse(startsWith(cell_type, "T cell regulatory"), "T cell CD4+",
                                                        ifelse(startsWith(cell_type, "T cells regulatory"), "T cell CD4+",
                                                 ifelse(cell_type == "CD8 T cells", "T cell CD8+",
                                                 ifelse(startsWith(cell_type, "T cell CD8"), "T cell CD8+",
                                                        ifelse(startsWith(cell_type, "T cells CD8"), "T cell CD8+",
                                                        ifelse(startsWith(cell_type, "B "), "B cell",
                                                               ifelse(startsWith(cell_type, "Neutro"), "Neutrophil",
                                                                                 ifelse(startsWith(cell_type, "NK"), "NK cell", cell_type)))))))))))]
# save(variants_ischgl_results, file = "data/variants_ischgl_deconv.RData")


# add cv for nuns
nuns_results$cell_type <- as.character(nuns_results$cell_type)
nuns_results[, cv_cell_type := ifelse(startsWith(cell_type, "T cell CD4"), "T cell CD4+",
                                      ifelse(startsWith(cell_type, "T cells CD4"), "T cell CD4+",
                                                 ifelse(cell_type == "T cells", "T cell CD4+",
                                                        ifelse(startsWith(cell_type, "T cell regulatory"), "T cell CD4+",
                                                               ifelse(startsWith(cell_type, "T cells regulatory"), "T cell CD4+",
                                                                ifelse(cell_type == "CD8 T cells", "T cell CD8+",
                                                                       ifelse(startsWith(cell_type, "T cell CD8"), "T cell CD8+",
                                                                              ifelse(startsWith(cell_type, "T cells CD8"), "T cell CD8+",
                                                                              ifelse(startsWith(cell_type, "B "), "B cell",
                                                                                     ifelse(startsWith(cell_type, "Neutro"), "Neutrophil",
                                                                                            ifelse(startsWith(cell_type, "NK"), "NK cell", cell_type)))))))))))]
# save(nuns_results, file = "data/nuns_deconv.RData")


# alternative day groups for nuns
nuns_results[, day_group_2 := ifelse(num_day <= 6, "day [0,1]", # only day 0 for Naive
                                     ifelse(num_day <= 11, "day [7,11]", "> day 11"))]
# save(nuns_results, file = "data/nuns_deconv.RData")


## save cibersortx results
# save(variants_ischgl_results, file = "data/variants_ischgl_deconv_with_cibersortx.RData")
# save(nuns_results, file = "data/nuns_deconv_with_cibersortx.RData")


