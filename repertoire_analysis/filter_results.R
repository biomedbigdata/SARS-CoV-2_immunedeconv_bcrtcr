# script to filter mixcr and trust4 output to perform clustering

library(data.table)
library(dplyr)
library(ggplot2)

# exclude samples
samples_to_exclude_pattern = "ID43_1st|ID29_2nd|ID34_2nd|ID38|ID53_2nd|SRR18922948|SRR18922902|SRR18922909|SRR13187045|SRR13187053|SRR13187055"

filter_samples <- function(filenames, pattern){
  filenames[!grepl(pattern, filenames)]
}

# MIXCR results
alpha_files <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/results_alpha", pattern = ".*clones.txt|.*clones.tsv", full.names = TRUE)
alpha_ek_files <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/results_alpha_ek", pattern = ".*clones.txt|.*clones.tsv", full.names = TRUE)
gamma_files <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/results_gamma", pattern = ".*clones.txt|.*clones.tsv", full.names = TRUE)
omicron_files <- list.files("/nfs/data2/covid_hennighausen/omikron/01_bcrtcr_repertoire/mixcr", pattern = ".*clones.txt|.*clones.tsv", full.names = TRUE)
omicron2_files <- list.files("/nfs/data2/covid_hennighausen/omikron_new/01_bcrtcr_repertoire/results_mixcr/", pattern = ".*clones.tsv", full.names = TRUE)
sero_files <- list.files("/nfs/data3/covid_hennighausen/total_RNA_seq_covid_ischgl_PBMCs/08_bcrtcr_repertoire/mixcr", pattern = ".*clones.txt|.*clones.tsv", full.names = TRUE)

## TRUST4 results
alpha_files_trust <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/trust4/results_alpha", pattern=".*report.tsv", full.names=TRUE)
alpha_ek_files_trust <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/trust4/results_alpha_ek", pattern=".*report.tsv", full.names=TRUE)
gamma_files_trust <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/trust4/results_gamma", pattern=".*report.tsv", full.names=TRUE)
omicron_files_trust <- list.files("/nfs/data2/covid_hennighausen/omikron/01_bcrtcr_repertoire/trust4", pattern=".*report.tsv", full.names=TRUE)
omicron2_files_trust <- list.files("/nfs/data2/covid_hennighausen/omikron_new/01_bcrtcr_repertoire/results_trust", pattern=".*report.tsv", full.names=TRUE)
sero_files_trust <- list.files("/nfs/data3/covid_hennighausen/total_RNA_seq_covid_ischgl_PBMCs/08_bcrtcr_repertoire/trust4", pattern=".*report.tsv", full.names=TRUE)

# filter for excluded samples, MIXCR
alpha_files <- filter_samples(alpha_files, samples_to_exclude_pattern)
alpha_ek_files <- filter_samples(alpha_ek_files, samples_to_exclude_pattern)
omicron_files <- filter_samples(omicron_files, samples_to_exclude_pattern)
omicron2_files <- filter_samples(omicron2_files, samples_to_exclude_pattern)
sero_files <- filter_samples(sero_files, samples_to_exclude_pattern)

# filter for excluded samples, TRUST4
alpha_files_trust <- filter_samples(alpha_files_trust, samples_to_exclude_pattern)
alpha_ek_files_trust <- filter_samples(alpha_ek_files_trust, samples_to_exclude_pattern)
omicron_files_trust <- filter_samples(omicron_files_trust, samples_to_exclude_pattern)
omicron2_files_trust <- filter_samples(omicron2_files_trust, samples_to_exclude_pattern)
sero_files_trust <- filter_samples(sero_files_trust, samples_to_exclude_pattern)


## combine omicron BA.1 and omicron BA.2 together
# extract SRR ids
omicron_meta <- fread("data/omicron_samplings.csv")
omicron_ids <- sapply(omicron_files_trust, function(x) unlist(strsplit(x, "_"))[5])
omicron2_ids <- sapply(omicron2_files_trust, function(x) unlist(strsplit(x, "_"))[7])
omicron_df <- data.table(mixcr_file = omicron_files, trust_file = omicron_files_trust, Run = omicron_ids)
omicron2_df <- data.table(mixcr_file = omicron2_files, trust_file = omicron2_files_trust, Run = omicron2_ids)

# combine and merge with metadata
omicron_total_df <- rbindlist(list(omicron_df, omicron2_df))
omicron_total_merge <- merge(omicron_total_df, omicron_meta, by = "Run", all.x = T)

# get the files
omicron_ba1_files <- omicron_total_merge[grepl("BA.1", Sample)]$mixcr_file
omicron_ba2_files <- omicron_total_merge[grepl("BA.2", Sample)]$mixcr_file
omicron_ba1_files_trust <- omicron_total_merge[grepl("BA.1", Sample)]$trust_file
omicron_ba2_files_trust <- omicron_total_merge[grepl("BA.2", Sample)]$trust_file




# method to combine mixcr results into one result table
# use omicron metadata for sampling info on omicron samples (see code above: omicron_meta)
combine_group_results <- function(files, method = "mixcr", omicron = F){
  id_index <- ifelse(method == "mixcr", 1, 2)
  rbindlist(lapply(files, function(file){
    dt <- fread(file)
    source = basename(file)
    dt[, source := source,]
    run <- strsplit(source, "_")[[1]][id_index]
    dt[, id := run]
    # sampling is either in the filename or from metadata for omicron
    if (omicron) {
      sampling <- omicron_meta[Run == run]$sampling[1]
      dt[, sampling := sampling]
    } 
    else {
      dt[, sampling := strsplit(source, "_")[[1]][id_index + 1]]
    }

    return(dt)
  }), fill=T)
}

# method to filter the sequences nd keep only sequences that appears in more than <cutoff> samples
# exclude first sampling samples if desired
# keep each sequence only once
get_common_aaseqs <- function(dt, cutoff = 5, filter_1st = FALSE){
 
  # rename column with aa cdr3 seq for consistency between mixccr and trust
  names(dt)[names(dt) == 'aaSeqCDR3'] <- 'CDR3aa'
  
  if (filter_1st){
    dt <- dt[sampling != "1st"] # exclude sampling "1st" 
  }
  dt <- dt[, .(id, sampling, CDR3aa)]
  dt <- dt[!duplicated(dt)] # remove duplicated sequences that appear more than once in each sample
  dt <- dt[, .N, by = CDR3aa] # count in how many samples a sequence appears
  dt <- dt[N >= cutoff] # filter to keep only sequences that appear in more than <cutoff> samples
  return(dt[order(N, decreasing=TRUE)]) # sort the dt and return
}


# combine results from all samples in each group
alpha_mixcr <- combine_group_results(alpha_files) 
alpha_ek_mixcr <- combine_group_results(alpha_ek_files) 
gamma_mixcr <- combine_group_results(gamma_files) 
omicron_ba1_mixcr <- combine_group_results(omicron_ba1_files, omicron = T)
omicron_ba2_mixcr <- combine_group_results(omicron_ba2_files, omicron = T)
sero_mixcr <- combine_group_results(sero_files) 


alpha_trust <- combine_group_results(alpha_files_trust, method = "trust")
alpha_ek_trust <- combine_group_results(alpha_ek_files_trust, method = "trust")
gamma_trust <- combine_group_results(gamma_files_trust, method = "trust")
omicron_ba1_trust <- combine_group_results(omicron_ba1_files_trust, method = "trust", omicron = T)
omicron_ba2_trust <- combine_group_results(omicron_ba2_files_trust, method = "trust", omicron = T)
sero_trust <- combine_group_results(sero_files_trust, method = "trust")



# only keep the sequences that appear in more than <cutoff> samples
alpha_mixcr_seqs <- get_common_aaseqs(alpha_mixcr)
alpha_ek_mixcr_seqs <- get_common_aaseqs(alpha_ek_mixcr, cutoff = 3)
gamma_mixcr_seqs <- get_common_aaseqs(gamma_mixcr, cutoff = 2)
omicron_ba1_mixcr_seqs <- get_common_aaseqs(omicron_ba1_mixcr, cutoff = 6)
omicron_ba2_mixcr_seqs <- get_common_aaseqs(omicron_ba2_mixcr, cutoff = 3)
sero_mixcr_seqs <- get_common_aaseqs(sero_mixcr)

alpha_trust_seqs <- get_common_aaseqs(alpha_trust)
alpha_ek_trust_seqs <- get_common_aaseqs(alpha_ek_trust, cutoff = 3)
gamma_trust_seqs <- get_common_aaseqs(gamma_trust, cutoff = 2)
omicron_ba1_trust_seqs <- get_common_aaseqs(omicron_ba1_trust, cutoff = 6)
omicron_ba2_trust_seqs <- get_common_aaseqs(omicron_ba2_trust, cutoff = 3)
sero_trust_seqs <- get_common_aaseqs(sero_trust)



# save all mixcr results
write.csv(alpha_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_alpha_seqs.csv")
write.csv(alpha_ek_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_alpha_ek_seqs.csv")
write.csv(gamma_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_gamma_seqs.csv")
write.csv(omicron_ba1_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_omicron_ba1_seqs.csv")
write.csv(omicron_ba2_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_omicron_ba2_seqs.csv")
write.csv(sero_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/mixcr_sero_seqs.csv")

# save all trust results
write.csv(alpha_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_alpha_seqs.csv")
write.csv(alpha_ek_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_alpha_ek_seqs.csv")
write.csv(gamma_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_gamma_seqs.csv")
write.csv(omicron_ba1_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_omicron_ba1_seqs.csv")
write.csv(omicron_ba2_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_omicron_ba2_seqs.csv")
write.csv(sero_trust_seqs, file = "data/filtered_bcrtcr_seqs/trust_sero_seqs.csv")




## exclude 1st sampling -> no gamma because gamma has only first sampling, no sero because doesnt make sense

alpha_mixcr_seqs <- get_common_aaseqs(alpha_mixcr, filter_1st = T)
alpha_ek_mixcr_seqs <- get_common_aaseqs(alpha_ek_mixcr, cutoff = 3, filter_1st = T)
omicron_ba1_mixcr_seqs <- get_common_aaseqs(omicron_ba1_mixcr, cutoff = 6, filter_1st = T)
omicron_ba2_mixcr_seqs <- get_common_aaseqs(omicron_ba2_mixcr, cutoff = 3, filter_1st = T)

alpha_trust_seqs <- get_common_aaseqs(alpha_trust, filter_1st = T)
alpha_ek_trust_seqs <- get_common_aaseqs(alpha_ek_trust, cutoff = 3, filter_1st = T)
omicron_ba1_trust_seqs <- get_common_aaseqs(omicron_ba1_trust, cutoff = 6, filter_1st = T)
omicron_ba2_trust_seqs <- get_common_aaseqs(omicron_ba2_trust, cutoff = 3, filter_1st = T)

# save all filtered mixcr results
write.csv(alpha_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/mixcr_alpha_seqs.csv")
write.csv(alpha_ek_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/mixcr_alpha_ek_seqs.csv")
write.csv(omicron_ba1_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/mixcr_omicron_ba1_seqs.csv")
write.csv(omicron_ba2_mixcr_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/mixcr_omicron_ba2_seqs.csv")

# save all filtered trust results
write.csv(alpha_trust_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/trust_alpha_seqs.csv")
write.csv(alpha_ek_trust_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/trust_alpha_ek_seqs.csv")
write.csv(omicron_ba1_trust_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/trust_omicron_ba1_seqs.csv")
write.csv(omicron_ba2_trust_seqs, file = "data/filtered_bcrtcr_seqs/excluded_1st/trust_omicron_ba2_seqs.csv")
