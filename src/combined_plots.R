#########################
# create combined plots # 
#########################
library(data.table)
library(ggplot2)

load("data/variants_ischgl_deconv.RData")
load("data/nuns_deconv.RData")

# conditional plot variants_ischgl
conditional_plot_dt <- variants_ischgl_results[method %in% c("quantiseq", "mcp_counter", "epic", "xcell")]
conditional_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
conditional_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]


conditional_variants_comparisons = list(c("Alpha", "Seronegative"), c("Alpha_EK", "Seronegative"), c("Gamma", "Seronegative"))
ggplot(conditional_plot_dt[cell_type %in% c("B cell", "Neutrophil", "T cells")], 
       aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  stat_compare_means(comparisons = conditional_variants_comparisons) +
  labs(title = "Differences in immune cell abundances between infected and healthy") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# without color because of duplicate information
ggplot(conditional_plot_dt[cell_type %in% c("B cell", "Neutrophil", "T cells")], 
       aes(x = group, y = value)) +
  geom_boxplot() +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  stat_compare_means(comparisons = conditional_variants_comparisons) +
  labs(title = "Differences in immune cell abundances between infected and healthy") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))


# conditional plot nuns
nuns_conditional_plot_dt <- nuns_results[method %in% c("quantiseq", "mcp_counter", "epic", "xcell")]
nuns_conditional_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
nuns_conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
nuns_conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
nuns_conditional_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]

ggplot(nuns_conditional_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
       aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Naive", "Convalescent"))) +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  labs(title = "Differences in immune cell abundances between naive and convalescent") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank())

# same without color
ggplot(nuns_conditional_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
       aes(x = group, y = value)) +
  geom_boxplot() +
  stat_compare_means(comparisons = list(c("Naive", "Convalescent"))) +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  labs(title = "Differences in immune cell abundances between convalescent and naive") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5))

# time series variants_ischgl --> exclude Gamma samples, use epic and mcp_counter as they are best for time series data and the cell types
variants_time_plot_dt <- variants_ischgl_results[method %in% c("mcp_counter", "epic") & group != "Gamma"]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]# variants_time_plot_dt$day_group <- factor(variants_time_plot_dt$day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30" ))

time_variants_comparisons = list(c("Alpha", "Seronegative"), c("Alpha_EK", "Seronegative"))
ggplot(variants_time_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
       aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  facet_grid(rows = vars(cell_type, method), cols = vars(factor(day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30"))), scales = "free") +
  stat_compare_means(comparisons = time_variants_comparisons) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  labs(title = "Difference in immune cell abundance over time after infection")

# same without color (duplicate info)
ggplot(variants_time_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
       aes(x = group, y = value)) +
  geom_boxplot() +
  facet_grid(rows = vars(cell_type, method), cols = vars(factor(day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30"))), scales = "free") +
  stat_compare_means(comparisons = time_variants_comparisons) +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Difference in immune cell abundance over time after infection")

# day_group on x axis, not as facet
ggplot(variants_time_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
       aes(x = factor(day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30")), y = value, fill = group)) +
  geom_boxplot() +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
  labs(title = "Difference in immune cell abundance over time after infection",
       x = "day_group")



# time series nuns
nuns_naive_time_plot_dt <- nuns_results[method %in% c("quantiseq", "epic") & group == "Naive"]
nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]# variants_time_plot_dt$day_group <- factor(variants_time_plot_dt$day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30" ))
# nuns_naive_time_plot_dt[, day_group := ifelse(num_day <= 5, "<= day 5", 
#                                               ifelse(num_day <= 11, "<= day 11",
#                                                      ifelse(num_day <= 40, "<= day 40", "> day 40")))]
# nuns_naive_time_plot_dt$day_group <- factor(nuns_naive_time_plot_dt$day_group, levels = c("<= day 5", "<= day 11", "<= day 40", "> day 40"))

time_nuns_comparisons = list(c("<= day 11", "> day 11"))
ggplot(nuns_naive_time_plot_dt[cell_type %in% c("T cells")], 
       aes(x = day_group, y = value)) +
  geom_boxplot() +
  facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
  stat_compare_means(comparisons = time_nuns_comparisons) +
  labs(title = "Difference in immune cell abundance over time after vaccination")

