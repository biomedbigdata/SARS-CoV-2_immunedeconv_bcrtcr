# tests with mixcr output
library(data.table)
library(dplyr)
library(ggplot2)

clones_files <- list.files("/nfs/data2/covid_hennighausen/covid_pbmcs_variants/07_bcrtcr_repertoire/results_alpha", pattern="*clones.txt", full.names=TRUE)

all_clones_alpha <- lapply(clones_files, function(file){
  dt <- fread(file)
  dt[, source := basename(file),]
  dt[, id := strsplit(source, "_")[[1]][1]]
})
names(all_clones_alpha) <- lapply(clones_files, function(file) basename(file))

clones_alpha_bind <- rbindlist(all_clones_alpha)
colnames(clones_alpha_bind)
clones_alpha_bind <- clones_alpha_bind[, .(id, source, cloneId, cloneCount, cloneFraction, targetSequences, nSeqCDR3, aaSeqCDR3)]

# save large dataframe to tsv file
write.csv(clones_alpha_bind, file = "data/alpha_clones.tsv", sep = "\t")


# ID10_1st & ID10_2nd
test_merge12 <-merge(all_clones_alpha[[1]][,.(targetSequences)], all_clones_alpha[[2]][,.(targetSequences)])
overlap12 <- test_merge12[! is.na(cloneId.x) & ! is.na(cloneId.y)]
bind12 <- rbindlist(list(all_clones_alpha[[1]], all_clones_alpha[[2]]))

# ID10_1st & ID11_1st
overlap13 <-merge(all_clones_alpha[[1]], all_clones_alpha[[3]], by = "targetSequences")


# ID11_1st & ID11_2nd
overlap34 <-merge(all_clones_alpha[[3]], all_clones_alpha[[4]], by = "targetSequences")




# overlap in targetSequences/aaCDR3 sequence of all alpha samples
targetSequences <- Reduce(intersect, lapply(all_clones_alpha[c(1:3)], function(dt) dt$aaSeqCDR3))

n_overlaps <- lapply(c(1:length(all_clones_alpha)), function(i) length(Reduce(intersect, lapply(all_clones_alpha[c(1:i)], function(dt) dt$aaSeqCDR3))))
Reduce(intersect, lapply(all_clones_alpha, function(dt) dt$aaSeqCDR3))

n_overlaps_dt <- data.table(
  n = c(1:length(all_clones_alpha)),
  overlap_size = n_overlaps
)
ggplot(n_overlaps_dt[-1,], aes(x = n, y = overlap_size)) + geom_col()


Reduce(intersect, lapply(all_clones_alpha, function(dt) dt$aaSeqCDR3))




# number of clones per sample
n_clones <- lapply(all_clones_alpha, function(dt) nrow(dt))
n_clones

ggplot(data.table(sample = names(n_clones),
                  n_clone = n_clones), 
       aes(x=sample, y = n_clone)) + 
  geom_col() +
  geom_hline(yintercept = 7916.553)


mean(unlist(n_clones))



### try scirpy distance metrics
library(reticulate)
use_condaenv("bcrtcr_python") # environment containing scirpy

sc <- import("scirpy")



### read distance matrix from tcr dist
library(tseries) # for read.matrix()
alpha_distances <- read.matrix("data/alpha_distances.txt")
alpha_count10_distances <- read.matrix("data/alpha_count10_distances.txt")

coords <- cmdscale(alpha_distances)
coords_count10 <- cmdscale(alpha_count10_distances)


coords_dt <- as.data.frame(coords)
names(coords_dt) <- c("x", "y")
coords_dt$sample <- clones_alpha_bind$id[1:10000]
coords_dt$source <- clones_alpha_bind$source[1:10000]

ggplot(coords_count10, aes(x = x, y = y, color = source)) + geom_point()







