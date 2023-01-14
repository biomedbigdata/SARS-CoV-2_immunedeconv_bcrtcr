# script to create metadata summary table variant x time (#)

library(data.table)
library(tidyverse)

alpha_samplesheet <- fread("C://Users/User/LRZ Sync+Share/BA/metadata/geopapers/Supplementary_Table_alpha.csv")

alpha_samplesheet$`discharged sampling` <- as.numeric(alpha_samplesheet$`discharged sampling`)
alpha_samplesheet$`covalescent sampling` <- as.numeric(alpha_samplesheet$`covalescent sampling`)

# all alpha summary
summary(alpha_samplesheet[, .(`hospitalized sampling`, `discharged sampling`, `covalescent sampling`)])
dim(alpha_samplesheet)  
# first sampling:
#    - 45 values from 0-18, mean: 5.422
# second sampling: 
#    - 40 values from 10-39, mean: 20.82
# third sampling:
#    - 12 values from 24-65, mean: 44.75

# include gamma
summary(c(alpha_samplesheet$`hospitalized sampling`, 3, 9, 11))
# 1st sampling: 
#    - 48 values from 0 to 18, mean: 5.562

# split to alpha and alpha_ek
alpha_ek_samplesheet <- alpha_samplesheet[Variant == "Alpha+EK"]
alpha_samplesheet <- alpha_samplesheet[Variant == "Alpha"]


# number of alpha samples per day
summary(alpha_samplesheet[, .(`hospitalized sampling`, `discharged sampling`, `covalescent sampling`)])
dim(alpha_samplesheet)        
# 1st sampling: 
#   - 32 values from 0 to 18, mean: 5.56
# 2nd sampling:
#   - 30 values from 10 to 29, mean: 21.83
# 3rd sampling:
#   - 5 values from 42 to 53, mean: 45.8


summary(alpha_ek_samplesheet[, .(`hospitalized sampling`, `discharged sampling`, `covalescent sampling`)])
dim(alpha_ek_samplesheet)  
# 1st sampling: 
#   - 13 values from 1 to 11, mean: 4.846
# 2nd sampling:
#   - 10 values from 10 to 37, mean: 17.8
# 3rd sampling:
#   - 7 values from 24 to 65, mean: 44



# GAMMA
# 1st sampling: 3 values from 3 to 11, mean: 7.667
# no other sampling


# Omicron?
# 1st sampling: 47 (no info)
# Todo: Time series data

# first experiment
omicron1_samplesheet <- fread("/nfs/data2/covid_hennighausen/omikron/_meta/SraRunTable.csv")
omicron2_samplesheet <- fread("/nfs/data2/covid_hennighausen/omikron_new/_meta/SraRunTable.txt")

table(omicron1_samplesheet$disease_state) # 47 Omicron, 8 Healthy
table(omicron2_samplesheet$omicron_sublineage) # 32 BA.1, 23 BA.2, 1 --

omicron_samplesheet <- rbindlist(list(omicron1_samplesheet, omicron2_samplesheet))
omicron_samplesheet <- omicron_samplesheet[disease_state != "Healthy control"]

omicron_samplesheet <- omicron_samplesheet[, .(Run, `Sample Name`, omicron_sublineage, days_after_positive_pcr_results)]
omicron_samplesheet[, Day := gsub("[^0-9]*", "", days_after_positive_pcr_results)]
omicron_samplesheet$Day <- as.numeric(omicron_samplesheet$Day)

# get number of samples for each group
omicron_samplesheet[grepl("BA.1", omicron_sublineage, fixed = T) & Day %in% c(0:5)]
omicron_samplesheet[grepl("BA.1", omicron_sublineage, fixed = T) & Day %in% c(6:10)]
omicron_samplesheet[grepl("BA.1", omicron_sublineage, fixed = T) & Day %in% c(11:15)]
omicron_samplesheet[grepl("BA.1", omicron_sublineage, fixed = T) & Day %in% c(16:30)]
omicron_samplesheet[grepl("BA.1", omicron_sublineage, fixed = T) & Day > 30]

omicron_samplesheet[grepl("BA.2", omicron_sublineage, fixed = T) & Day %in% c(0:5)]
omicron_samplesheet[grepl("BA.2", omicron_sublineage, fixed = T) & Day %in% c(6:10)]
omicron_samplesheet[grepl("BA.2", omicron_sublineage, fixed = T) & Day %in% c(11:15)]
omicron_samplesheet[grepl("BA.2", omicron_sublineage, fixed = T) & Day %in% c(16:30)]
omicron_samplesheet[grepl("BA.2", omicron_sublineage, fixed = T) & Day > 30]



# read mapping gsm to samplename
omicron_mapping <- fread("data/omicron_GSM_to_name.txt")
omicron_merge <- merge(omicron_samplesheet, omicron_mapping, by.x = "Sample Name", by.y = "GSM")
omicron_merge <- separate(omicron_merge, Sample_name, c("Variant", "ID", "sampling"), sep = "_")
omicron_merge <- unite(omicron_merge, "Sample", c(Variant, omicron_sublineage, ID, sampling), sep = "_", remove = F)

# number of samples for each category
omicron_merge["BA.1" == omicron_sublineage & sampling == "1st"]
omicron_merge["BA.1" == omicron_sublineage & sampling == "2nd"]
omicron_merge["BA.2" == omicron_sublineage & sampling == "1st"]
omicron_merge["BA.2" == omicron_sublineage & sampling == "2nd"]

# for range and mean
summary(omicron_merge["BA.1" == omicron_sublineage & sampling == "1st"]$Day)
summary(omicron_merge["BA.1" == omicron_sublineage & sampling == "2nd"]$Day)
summary(omicron_merge["BA.2" == omicron_sublineage & sampling == "1st"]$Day)
summary(omicron_merge["BA.2" == omicron_sublineage & sampling == "2nd"]$Day)


write.csv(omicron_merge[order(Sample, sampling), .(Run, Sample, ID, sampling, Day)], "data/omicron_samplings.csv",
          row.names = F, quote = F)


# add to full_metadata (all_pbmc_metadata.RData)
omicron_metadata <- omicron_merge[order(Sample, sampling), .(Run, Sample, ID, omicron_sublineage, sampling, Day)] 
colnames(omicron_metadata) <- c("old_id", "sample_id", "patient", "group", "sampling", "num_day")

omicron_metadata[, infected := "infected"]
omicron_metadata[, day_group := ifelse(num_day <= 11, "<= day 11", "> day 11")]

load("data/all_pbmc_metadata.RData")

full_metadata <- rbindlist(list(full_metadata[!group %in% c("BA.1", "BA.2")], omicron_metadata), use.names = T)
save(full_metadata, file = "data/all_pbmc_metadata.RData")

# Seronegative
# 54 samples




# Alpha, Alpha+EK, Gamma, Omicron large table
alpha_samplesheet <- fread("C://Users/User/LRZ Sync+Share/BA/metadata/geopapers/Supplementary_Table_alpha.csv")


alpha_samplesheet <- alpha_samplesheet[, .(Variant, `Sample ID`, `hospitalized sampling`, `discharged sampling`, `covalescent sampling`)]
alpha_samplesheet$`hospitalized sampling` <- as.numeric(alpha_samplesheet$`hospitalized sampling`)
alpha_samplesheet$`discharged sampling` <- as.numeric(alpha_samplesheet$`discharged sampling`)
alpha_samplesheet$`covalescent sampling` <- as.numeric(alpha_samplesheet$`covalescent sampling`)

colnames(alpha_samplesheet) <- c("Variant", "ID", "1st", "2nd", "3rd")


alpha_samplesheet_melt <- melt(alpha_samplesheet, id.vars = c("Variant", "ID"), variable.name = "sampling", value.name = "Day")
alpha_samplesheet_melt <- alpha_samplesheet_melt[!is.na(Day)]
alpha_samplesheet_melt[Variant == "Alpha+EK"]$Variant <- "Alpha_EK"
alpha_samplesheet_unite <- unite(alpha_samplesheet_melt, "Sample", c(Variant, ID), sep = "_")
write.csv(alpha_samplesheet_unite[order(Sample, sampling)], "C://Users/User/LRZ Sync+Share/BA/metadata/alpha_samplings.csv",
          row.names = F, quote = F)

# for summary table
alpha_samplesheet_unite[!grepl("Alpha_EK", Sample, fixed = T) & Day %in% c(16:30)]


gamma_samplesheet <- fread("C://Users/User/LRZ Sync+Share/BA/metadata/geopapers/Supplementary Table 1. patients.csv")

gamma_samplesheet <- gamma_samplesheet[, .(Variant, `Sample ID`, `hospitalized sampling`)]
gamma_samplesheet$`hospitalized sampling` <- as.numeric(gamma_samplesheet$`hospitalized sampling`)

colnames(gamma_samplesheet) <- c("Variant", "ID", "1st")

gamma_samplesheet$Variant <- "Gamma"
gamma_samplesheet_melt <- melt(gamma_samplesheet, id.vars = c("Variant", "ID"), variable.name = "sampling", value.name = "Day")
gamma_samplesheet_melt <- gamma_samplesheet_melt[!is.na(Day)]
gamma_samplesheet_unite <- unite(gamma_samplesheet_melt, "Sample", c(Variant, ID), sep = "_")
write.csv(gamma_samplesheet_unite[order(Sample, sampling)], "C://Users/User/LRZ Sync+Share/BA/metadata/gamma_samplings.csv",
          row.names = F, quote = F)



# TODO: combine days and create new full metadata with correct time data for ALpha, Gamma
