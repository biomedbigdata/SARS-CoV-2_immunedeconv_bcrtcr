# batch correction for all PBMCs using combat

library(sva)
library(data.table)

version <- "v1"

expr <- fread("data/variants_omicron_ischgl_pbmcs_tpms.csv")
gene_counts <- fread("data/variants_omicron_ischgl_pbmcs_geneCounts.csv")
batches_dt <- fread(paste0('data/batches/batches_', version, '.csv'))


# samples that failed QC or are excluded due to similar reasons to remove from batch sheet
# samples_to_drop_variants <- c('ID43_1st','B.1.351.ID5_3rd','B.1.351.ID1_3rd','ID29_2nd','ID34_2nd','ID38_3rd', 'ID53_2nd','BNT_Aus_20_2','ID38_2nd','ID38_3rd','ID38_1st')
# samples_to_drop_omicron <- c("IDOm18_1st")
# samples_to_drop_ischgl <- c('A_118_Asymptom', 'B_425_Seronegative','B_436_Seronegative','B_446_Seronegative')
# samples_to_drop_nuns <- c("BNT_Convalescent_4_Day7", "BNT_Convalescent_9_Day0", "BNT_Convalescent_9_Day1", "BNT_Convalescent_9_Day7", "BNT_Convalescent_13_Day0", "BNT_Convalescent_2_Day0", "BNT_Convalescent_3_Day0", "BNT_Naive_13_Day34", "BNT_Naive_13_Day41")


# get the correct current batches
# batches_dt <- batches_dt[group %in% c("Seronegative", "Alpha", "Alpha_EK", "Gamma", "Omicron")]
# batches_filtered <- batches_dt[!sample_id %in% c(samples_to_drop_ischgl, samples_to_drop_nuns, samples_to_drop_variants)]


# remove gene id column
gene_ids <- expr$gene_id
expr <- expr[, -1]
gene_counts <- gene_counts[, -1]
# filter expr
expr <- expr %>% select(batches_dt$sample_id)



# create batch ids (match sample to batch)
batchid <- match(batches_dt$batch, unique(batches_dt$batch)) 

# create groups
# groups <- sapply(colnames(expr), function(name) head(unlist(strsplit(name, '_')),1))
groups <- batches_dt$group


## do the correction for tpms
corrected_expr <- ComBat(dat = expr, batch = batchid, par.prior = T, prior.plots = F)
rownames(corrected_expr) <- gene_ids

# write to file
write.table(corrected_expr, 
            file = paste0('data/batch_corrected/corrected_tpms_', version, '.csv'), sep = ",", quote=F, row.names = T)




## do the correction for gene counts
# remove control 4 as only 1 sample for batch rep4
# gene_counts <- gene_counts[, -c('Control_rep4')]
batchid_counts <- batchid
groups_counts <- groups 

corrected_counts <- ComBat_seq(counts = as.matrix(gene_counts), 
                               batch = batchid_counts, 
                               group = groups_counts,
                               full_mod = F)

rownames(corrected_counts) <- gene_ids

# write to file
write.table(corrected_counts, 
            file = paste0('batch_corrected/corrected_gene_counts_', version, '.csv'), quote=F, row.names = T)