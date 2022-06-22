## script for quality control

library(data.table)
library(ggplot2)
library(tidyr)
library(heatmaply)
library(pheatmap)

# create mapped reads per contig plot and chack if percent of mitochondrial genes is higher than 10%
DIR <- "/nfs/data2/covid_hennighausen/vaccine_studies_nuns/03_nfcore_rnaseq_results/star_salmon/samtools_stats"
filenames <- list.files(DIR, pattern="*.idxstats", full.names=T)
names(filenames) <- basename(filenames)
tables <- lapply(filenames, fread) # read all files at once into a list of data.tables
idxstats_dt <- rbindlist(tables, idcol = 'filename')
colnames(idxstats_dt) <- c("filename", "contig", "contig_length", "n_mapped", "n_unmapped")

idxstats_dt[, percent_mapped := n_mapped/sum(n_mapped), by = filename]

ggplot(idxstats_dt[contig %in% c(c(1:22), "MT", "X", "Y")], aes(x = contig, y = percent_mapped, color = filename)) +
  geom_point() + theme(legend.position="none")

# all samples with MT percent > 10 %
unlist(lapply(idxstats_dt[contig == "MT" & percent_mapped > 0.1]$filename, function(x){
  strsplit(x, ".", fixed=T)[[1]][[1]]
})) 



# deSeq pca plot
DIR <- "/nfs/data2/covid_hennighausen/vaccine_studies_nuns/03_nfcore_rnaseq_results/star_salmon/deseq2_qc/"
pca_dt <- fread(paste0(DIR, "deseq2.pca.vals.txt"))
pc_labels <- colnames(pca_dt)
colnames(pca_dt) <- c("sample", "PC1", "PC2")

# separate samplename for coloring
pca_dt <- unite(separate(pca_dt, col = sample, into = c(NA, "type", "number", "day"), sep = "_"),
                'sample', type:number, sep = "_", remove=F)


ggplot(pca_dt, aes(x = PC1, y = PC2, color = sample)) +
  geom_point() +
  labs(x = pc_labels[[2]], y = pc_labels[[3]]) +
  theme(legend.position = "none")


# get samples that are quite far away
pca_dt <- fread(paste0(DIR, "deseq2.pca.vals.txt"))
colnames(pca_dt) <- c("sample", "PC1", "PC2")
pca_dt[PC2 > 100, sample]



# deSeq distance heatmap
distances <- read.csv(paste0(DIR, "deseq2.sample.dists.txt"), sep = "\t")
rownames(distances) <- distances$sample
distances_mat <- as.matrix(distances[, -1])

pheatmap(distances_mat, cluster_cols = F, cluster_rows = F)
heatmaply(distances_mat, show_dendrogram = c(F,F)) # interactive and clustered


