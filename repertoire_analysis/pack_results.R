# bcr tcr results
cdr3_seqs <- fread("data/filtered_bcrtcr_seqs/CDR3_seqs.csv", header = T)
cdr3_seqs_10m <- fread("data/filtered_bcrtcr_seqs/CDR3_seqs_10.csv", header = T)


cdr3_seqs$V1 <- NULL
cdr3_seqs_10m$V1 <- NULL

cdr3_seqs <- cdr3_seqs[, .(CDR3aa, Groups, Source)]
cdr3_seqs_10m <- cdr3_seqs_10m[, .(CDR3aa, Groups, Source)]

bcrtcr_results <- list(cdr3_seqs, cdr3_seqs_10m)

# save(cdr3_sequences, file = "data/filtered_bcrtcr_seqs/reconstructed_cdr3_sequences.RData")

# deconvolution results
load("data/variants_omicron_ischgl_deconv.RData")
deconvolution_result <- variants_omicron_ischgl_result
deconvolution_result$SRA <- NULL
colnames(deconvolution_result) <- c("sample_id", "cell_type", "method","value", "sampling", "id_num", "sample", "patient", "group", "infected", "day_group", "num_day", "cv_cell_type")


load("data/m10_variants_omicron_ischgl_deconv.RData")
deconvolution_result_10m <- variants_omicron_ischgl_result
deconvolution_result_10m$SRA <- NULL

load("data/m50_variants_omicron_ischgl_deconv.RData")
deconvolution_result_50m <- variants_omicron_ischgl_result
deconvolution_result_50m$SRA <- NULL


deconvolution_results <- list(deconvolution_result, deconvolution_result_50m, deconvolution_result_10m)

save(bcrtcr_results, deconvolution_results, file = "data/SARS_CoV_2_immunedeconv_bcrtcr.RData")

