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

# combine with metadata
load(file = "data/all_pbmc_metadata.RData") # load metadata
load(file = "data/alpha_gamma_time_meta.RData") # load time series metadata for variants

variants_ischgl_results <- merge(variants_ischgl_results, alpha_gamma_sero_meta_test,
                                 by.x = "sample", by.y = "sample_id", all.x = T, allow.cartesian = T)


nuns_results <- merge(nuns_results, full_metadata, by.x = "sample", by.y = "old_id", 
                      all.x = T, allow.cartesian = T)

# save the results as RData
save(variants_ischgl_results, file = "data/variants_ischgl_deconv.RData")
save(nuns_results, file = "data/nuns_deconv.RData")

# test loading
load("data/variants_ischgl_deconv.RData")
load("data/nuns_deconv.RData")

