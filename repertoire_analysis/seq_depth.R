# analyze sequencing depth estimates for variants and omicron
library(data.table)
library(ggplot2)

estimates <- fread("bcr_tcr_repertoire/average_seq_depth_estimate_variants.txt")
ggplot(estimates, aes(x=V1)) + geom_histogram()

mean(estimates$V1)
summary(estimates)

estimates_o <- fread("bcr_tcr_repertoire/average_seq_depth_estimate_omikron_new.txt")
ggplot(estimates_o, aes(x=V1)) + geom_histogram()

mean(estimates_o$V1)
summary(estimates_o)


