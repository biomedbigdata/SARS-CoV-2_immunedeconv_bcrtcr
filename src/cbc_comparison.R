###
# Script to compare deconvolution results with ground truth - CBC data
###

library(data.table)
library(ggplot2)


# read CBC data
cbc_alpha <- fread("data/CBC/cbc_alpha.csv") 
colnames(cbc_alpha) <- unlist(cbc_alpha[2])


cbc_alpha_h <- cbc_alpha[,1:10] # hospitalized samplings
cbc_alpha_d <- cbc_alpha[,c(1:4, 11:16)] # discharged samplings

colnames(cbc_alpha_h) <- c("Variant", "Sample_ID", "Sex", "Age", "CRP",  "PCT", "White_blood_cells", "Neutrophil", "Lymphocyte", "Monocyte")
colnames(cbc_alpha_d) <- c("Variant", "Sample_ID", "Sex", "Age", "CRP",  "PCT", "White_blood_cells", "Neutrophil", "Lymphocyte", "Monocyte")

cbc_alpha_h <- cbc_alpha_h[3:51]
cbc_alpha_d <- cbc_alpha_d[3:51]

# TODO: 
# 1) reformat cbc data
# 2) map to deconv results - check excluded samples!
# 3) compute correlation
# 4) visualize correlation -> as heatmap?

# reformat hospitalized samples
cbc_alpha_h <- melt(cbc_alpha_h, id.vars = c("Variant", "Sample_ID", "Sex", "Age"),
                 variable.name = "cbc_cell_type",
                 value.name = "percentage")
cbc_alpha_h$percentage <- sub(",", ".", cbc_alpha_h$percentage)
cbc_alpha_h$percentage <- as.numeric(cbc_alpha_h$percentage)

colnames(cbc_alpha_h) <- c("Variant", "Sample_ID", "Sex", "Age", "cbc_cell_type", "percentage_h")

# same for discharged samples
cbc_alpha_d <- melt(cbc_alpha_d, id.vars = c("Variant", "Sample_ID", "Sex", "Age"),
                    variable.name = "cbc_cell_type",
                    value.name = "percentage")
cbc_alpha_d$percentage <- sub(",", ".", cbc_alpha_d$percentage)
cbc_alpha_d$percentage <- as.numeric(cbc_alpha_d$percentage)

colnames(cbc_alpha_d) <- c("Variant", "Sample_ID", "Sex", "Age", "cbc_cell_type", "percentage_d")

# merge hospitalized samples and discharged samples into one df
cbc_alpha_full <- merge(cbc_alpha_h, cbc_alpha_d, all = T)

# make id column numeric
cbc_alpha_full$Sample_ID <- as.numeric(cbc_alpha_full$Sample_ID)


# read deconvolution results
load("data/variants_omicron_ischgl_deconv.RData")
deconv_results <- variants_omicron_ischgl_result

cell_types_comp <- c("Neutrophil", "T cell CD4+", "T cell CD8+", "B cell", "Monocyte")

deconv_results <- deconv_results[cv_cell_type %in% cell_types_comp]
deconv_results <- deconv_results[group == "Alpha" | group == "AlphaEK"]

# make wide format
deconv_results <- dcast(deconv_results[,.(id_num, group, method, sampling, cv_cell_type, value)], ...~cv_cell_type, fun.aggregate = sum)

# combine T and B cells into "Lymphocyte"
deconv_results[, Lymphocyte := `B cell`+`T cell CD4+`+`T cell CD8+`]

# make long format again
deconv_results <- melt(deconv_results, id.vars = c("id_num", "group", "method", "sampling"),
                       variable.name = "deconv_cell_type",
                       value.name = "score")

# cast to separate first and second sampling
deconv_results <- dcast(deconv_results, ...~sampling, value.var = "score")


# merge cbc data and deconv results
cbc_deconv <- merge(cbc_alpha_full[, .(Sample_ID, cbc_cell_type, percentage_h, percentage_d)],
      deconv_results,
      by.x = c("Sample_ID", "cbc_cell_type"), by.y = c("id_num", "deconv_cell_type"))


cell_types_comp <- c("Neutrophil", "Lymphocyte", "Monocyte")

# cor_df <- rbindlist(lapply(deconv_methods, function(m){
#   rbindlist(lapply(cell_types_comp, function(c){
#     cor_res <- cor(cbc_deconv[method == m & cbc_cell_type == c, .(percentage_h, `1st`)], method = "spearman")
#     print(cor_res)
#     data.frame(method = m,
#                cell_type = c,
#                correlation = cor_res)
#   }))
# }))

#cor.test(~ percentage_h + `1st`, cbc_deconv[method == "EPIC" & cbc_cell_type == "Neutrophil"],
 #                   method = "spearman")

# cbc_deconv[method == "EPIC" & cbc_cell_type == "Neutrophil"]

# together
cbc_deconv_melt <- melt(cbc_deconv, 
                        measure.vars = c("percentage_h", "percentage_d"), 
                        variable.name = "cbc_time_point", value.name = "percentage")

cbc_deconv_melt <- melt(cbc_deconv_melt, 
                        measure.vars = c("1st", "2nd"), 
                        variable.name = "deconv_time_point", value.name = "estimate")

cbc_deconv_melt <- cbc_deconv_melt[(cbc_time_point == "percentage_h" & deconv_time_point == "1st") | (cbc_time_point == "percentage_d" & deconv_time_point == "2nd")]

cbc_deconv_melt[, method := ifelse(method == "quantiseq", "quanTIseq", ifelse(method == "mcp_counter", "MCP-counter", ifelse(method=="xcell", "xCell", ifelse(method == "epic", "EPIC", method))))]

cbc_deconv_melt$shape <- as.factor(cbc_deconv_melt$group)
ggplot(cbc_deconv_melt, aes(x = percentage, y = estimate, color = cbc_cell_type)) + 
  geom_point() + 
  scale_color_brewer(palette = "Dark2") +
  geom_smooth(method = "lm") +
  facet_wrap(~method, scales = "free") +
  labs(title = "CBC data and deconvolution estimates correlation (sequencing depth 10 million)",
       x = "CBC percentage", y = "Deconvolution estimate",
       color = "Cell type")


cor_df <- rbindlist(lapply(c("epic", "mcp_counter", "quantiseq", "xcell"), function(m){
  rbindlist(lapply(cell_types_comp, function(c){
    print(paste(m, c))
    cor_res <- cor.test(cbc_deconv_melt[method == m & cbc_cell_type == c]$percentage, cbc_deconv_melt[method == m & cbc_cell_type == c]$estimate, method = "spearman")
    print(cor_res)
    data.frame(method = m,
               cell_type = c,
               p_value = cor_res$p.value,
               estimate = cor_res$estimate)
  }))
}))

cor.test(cbc_deconv_melt$percentage, cbc_deconv_melt$estimate, metod = "spearman")



### scatter plot to visualize correlation
# hospitalized / "first" samples
ggplot(cbc_deconv, aes(x = percentage_h, y = `1st`, color = cbc_cell_type, shape = group)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~method, scales = "free") +
  labs(title = "Correlation between CBC data and deconvolution estimates for hospitalized samples",
       x = "CBC percentage", y = "deconvolution estimate")


# discharged / "second" samples
ggplot(cbc_deconv, aes(x = percentage_d, y = `2nd`, color = cbc_cell_type)) + 
  geom_point() + 
  geom_smooth(method = "lm") +
  facet_wrap(~method, scales = "free") +
  labs(title = "Correlation between CBC data and deconvolution estimates for discharged samples",
       x = "CBC percentage", y = "deconvolution estimate")



test <- cor.test(cbc_deconv_melt$percentage, cbc_deconv_melt$estimate, method = "spearman")
test$estimate
