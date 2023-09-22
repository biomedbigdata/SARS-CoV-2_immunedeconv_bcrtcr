# downsample tpms matrices to 10 M using scuttle

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("scuttle")

library(scuttle)


load("data/variants_omicron_ischgl_tpms_prepared.RData")


downsampled_M10_tpms <- as.matrix(downsampleMatrix(variants_omicron_ischgl_tpms, 0.1, bycol = F))

dim(downsampled_M10_tpms)

save(downsampled_M10_tpms, file = "data/downsampled_M10_tpms.RData")

