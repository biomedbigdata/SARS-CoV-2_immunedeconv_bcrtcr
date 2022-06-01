## script for finding samples with high mitochondrial gene count to exclude in further analysis

library(data.table)

THRESHOLD <- 0.1

mapped_per_contig_variants_dt <- fread('/nfs/data2/covid_hennighausen/covid_pbmcs_variants/02_nfcore_rnaseq_results/01_inc_cat/multiqc/star_salmon/multiqc_data/mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.txt')
mapped_per_contig_ischgl_dt <- fread('/nfs/data3/covid_hennighausen/total_RNA_seq_covid_ischgl_PBMCs/07_nfcore_results/multiqc/star_salmon/multiqc_data/mqc_samtools-idxstats-mapped-reads-plot_Normalised_Counts.txt')

samples_to_remove <- mapped_per_contig_variants_dt[MT>THRESHOLD, Sample]

samples_to_remove <- c(samples_to_remove, mapped_per_contig_ischgl_dt[MT>THRESHOLD, Sample])
samples_to_remove

# find the one sample in variants pbmcs that has only counts at chr 21
problematic <- mapped_per_contig_variants_dt[order(`21`, decreasing = T), Sample][1] 
problematic
