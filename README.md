# SARS-CoV-2_immunedeconv_bcrtcr
Immune deconvolution using immunedeconv (https://github.com/icbi-lab/immunedeconv) and BCR TCR analysis

## Scripts for QC and Immunedeconvolution (scr folder)
* mito_qc.R: for looking into fastqc result to find samples with too high number of mitochondrial genes or other abnormalities in distribution of the genes along the chromosomes
* qc.R: visualizing and evaluating results (used for nuns PBMC data)
* immunedeconvolution.R: 
  - creating input tpm matrix for immune deconvolution methods
  - apply relevant immune deconvolution methods
  - create simple plots visualizing the results
* plotting.R: contains functions for creating useful data tables from result and visualizing the results
* combine_deconv_results.R: runs the immune deconvolution methods and combines the results with metadata to a large table
* combined_plots.R: visualizing the results cobined in combine_deconv_results.R
* cbc_comparison.R: script to compare deconvolution resuts from combine_deconv_results.R to CBC data
* downsample_tpms_scuttle.R: testing to downsample tpms matrices
* prepare*.R: scripts to reformat and annotate the data
* samplesheet_analysis.R: looking at the metadata and preparing special metadata for Omicron and Gamma
* simulate_seq_depth_polyester.R: testing to downsample sequencing depth with polyester

## Scripts for BCR TCR analysis (repertoire_analysis folder)
* filter_results.R: reading results from MiXCR and TRUST4 and preparing the sequences
* bcr_tcr_analysis.R: using the sequences from filter_results.R, filter to keep sequences that appear in more than half of the samples per each group, read in distance matrices and creating graphs (igraph) and sequence logos, finding unique sequences for pBLAST analysis
* bcrtcr_seq_depth.R: similar to bcr_tcr_analysis.R but with option to include results from different sequencing depths
* pack_results.R: packing all results into an RData object for publication
* seq_depth.R: analyzing the original sequencing depth of our samples
* python/calc_dist_matrix.py: functions to calculate pairwise distace matrix from sequences
* python/calc_dist_mat.ipynb: calculating pairwise distance matrices
* python/cdr3_distances.ipynb: testing to create network graphs and siatnce visualization in python
