# script to simulate different sequncing depths with polyester and save simulated data

if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("polyester")

library(polyester)
?polyester
