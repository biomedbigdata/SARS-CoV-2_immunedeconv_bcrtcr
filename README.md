# SARS-CoV-2_immunedeconv_bcrtcr
Immune deconvolution using immunedeconv (https://github.com/icbi-lab/immunedeconv) and bcr tcr analysis

## Scripts
* mito_qc.R: for looking into fastqc result to find samples with too high number of mitochondrial genes or other abnormalities in distribution of the genes along the chromosomes
* immunedeconvolution.R: 
  - creating input tpm matrix for immune deconvolution methods
  - apply relevant immune deconvolution methods
  - create plots visualizing the results
* plotting.R: contains functions for creating useful data tables from result and visualizing the results
