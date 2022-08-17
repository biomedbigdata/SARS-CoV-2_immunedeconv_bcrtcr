#########################
# create combined plots # 
#########################
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)

# load("data/variants_ischgl_deconv.RData")
# load("data/nuns_deconv.RData")
load("data/variants_ischgl_deconv_with_cibersortx.RData")
load("data/nuns_deconv_with_cibersortx.RData")

my_theme <- theme(panel.background = element_rect(fill = "white", colour = "grey"),
                  panel.grid.major = element_line(colour = "lightgrey"),
                  panel.grid.minor = element_line(colour = "lightgrey"))



deconv_methods <- c("quantiseq", "mcp_counter", "epic", "xcell", "cibersortx_abs", "cibersortx_rel")

# conditional plot variants_ischgl
conditional_plot_dt <- variants_ischgl_results[method %in% deconv_methods]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]
conditional_plot_dt[, value := ifelse(method == "mcp_counter", log(value), value)]

conditional_variants_comparisons = list(c("Alpha", "Seronegative"), c("Alpha_EK", "Seronegative"), c("Gamma", "Seronegative"))
ggplot(conditional_plot_dt[cv_cell_type %in% c("B cell", "Neutrophil", "T cell CD4+", "T cell CD8+")], 
       aes(x = group, y = value, fill = group)) +
  geom_boxplot() +
  facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
  stat_compare_means(comparisons = conditional_variants_comparisons) +
  labs(title = "Differences in immune cell abundances between infected and healthy") +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  my_theme

# variants conditional plots for each method combined to get free scales
conditional_plots <- lapply(deconv_methods, function(m) {
  ggplot(conditional_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
         aes(x = group, y = value, fill = group)) +
    geom_boxplot() + 
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free_y") +
    stat_compare_means(comparisons = conditional_variants_comparisons, vjust = 1.2) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
    labs(title = m, x = "group") +
    my_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
})

# !!! change number of conditional_plots and heights depending on number of methods used
conditional_figure <- ggarrange(conditional_plots[[1]], 
                                conditional_plots[[2]], 
                                conditional_plots[[3]], 
                                conditional_plots[[4]],
                                conditional_plots[[5]],
                                conditional_plots[[6]],
                    heights = c(1,1,1,1,1,1),
                    labels = names(conditional_plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(conditional_figure, 
                top = text_grob("Differences in immune cell abundances between infected and healthy", size = 16),  
                fig.lab = "A", fig.lab.size = 18, fig.lab.face = "bold")



# without color because of duplicate information
# ggplot(conditional_plot_dt[cell_type %in% c("B cell", "Neutrophil", "T cells")], 
#        aes(x = group, y = value)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
#   stat_compare_means(comparisons = conditional_variants_comparisons) +
#   labs(title = "Differences in immune cell abundances between infected and healthy") +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
#   my_theme


# conditional plot nuns
nuns_conditional_plot_dt <- nuns_results[method %in% deconv_methods]
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
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  my_theme

# same without color
# ggplot(nuns_conditional_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
#        aes(x = group, y = value)) +
#   geom_boxplot() +
#   stat_compare_means(comparisons = list(c("Naive", "Convalescent"))) +
#   facet_grid(rows = vars(method), cols = vars(cell_type), scales = "free") +
#   labs(title = "Differences in immune cell abundances between convalescent and naive") +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
#   my_theme

# time series variants_ischgl --> exclude Gamma samples, use epic and mcp_counter as they are best for time series data and the cell types
variants_time_plot_dt <- variants_ischgl_results[method %in% deconv_methods & group != "Gamma"]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
variants_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]# variants_time_plot_dt$day_group <- factor(variants_time_plot_dt$day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30" ))

# time_variants_comparisons = list(c("Alpha", "Seronegative"), c("Alpha_EK", "Seronegative"))
# ggplot(variants_time_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
#        aes(x = group, y = value, fill = group)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(cell_type, method), cols = vars(factor(day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30"))), scales = "free") +
#   stat_compare_means(comparisons = time_variants_comparisons) +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   labs(title = "Difference in immune cell abundance over time after infection") +
#   my_theme

# same without color (duplicate info)
# ggplot(variants_time_plot_dt[cell_type %in% c("Neutrophil", "B cell", "T cells")], 
#        aes(x = group, y = value)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(cell_type, method), cols = vars(factor(day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30"))), scales = "free") +
#   stat_compare_means(comparisons = time_variants_comparisons) +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
#   labs(title = "Difference in immune cell abundance over time after infection") +
#   my_theme

# day_group on x axis, not as facet
variants_time_plot_dt[, value := ifelse(method == "mcp_counter", log(value), value)]
ggplot(variants_time_plot_dt[cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
       aes(x = factor(day_group, levels = c("day [0,5]", "day [6,10]", "day [11,15]", "day [16,30]", "> day 30")), y = value, color = group)) +
  geom_boxplot() + 
  geom_smooth(method = "lm", aes(group=group), alpha = 0.3) +
  facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
  theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
  labs(title = "Difference in immune cell abundance over time after infection",
       x = "day_group") +
  my_theme


# variants time plots for each method combined to get free scales
plots <- lapply(deconv_methods, function(m) {
  ggplot(variants_time_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
                  aes(x = factor(day_group, levels = c("day [0,5]", "day [6,10]", "day [11,15]", "day [16,30]", "> day 30")), y = value, color = group)) +
    geom_boxplot() + 
    geom_smooth(method = "lm", aes(group=group), alpha = 0.3) +
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free") +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
    labs(title = m, x = "day_group") +
    my_theme
})

# !!! change number of conditional_plots and heights depending on number of methods used
figure <- ggarrange(plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[2]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[3]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[4]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[5]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[6]],
                    heights = c(1,1,1,1,1,1.4),
                    labels = names(plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(figure, 
                top = text_grob("Difference in immune cell abundance over time after infection", size = 16),
                fig.lab = "B", fig.lab.size = 18, fig.lab.face = "bold")





# time series nuns
deconv_methods_nuns <- c("quantiseq", "epic", "mcp_counter", "cibersortx_abs", "cibersortx_rel")
nuns_naive_time_plot_dt <- nuns_results[method %in% deconv_methods_nuns & group == "Naive"]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]


# nuns_naive_time_plot_dt[, day_group := ifelse(num_day <= 5, "<= day 5", 
#                                               ifelse(num_day <= 11, "<= day 11",
#                                                      ifelse(num_day <= 40, "<= day 40", "> day 40")))]
# nuns_naive_time_plot_dt$day_group <- factor(nuns_naive_time_plot_dt$day_group, levels = c("<= day 5", "<= day 11", "<= day 40", "> day 40"))
nuns_naive_time_plot_dt[, value := ifelse(method == "mcp_counter",  log(value), value)]
nuns_naive_time_plot_dt[, value := ifelse(method == "xcell" & cv_cell_type == "Neutrophil", signif(value, digits = 2), value)]
# naive only have day 0 not day 1
nuns_naive_time_plot_dt[, day_group_2 := ifelse(day_group_2 == "day [0,1]", "day 0", day_group_2)]


ggplot(nuns_naive_time_plot_dt[cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
       aes(x = factor(day_group_2, levels = c("day 0", "day [7,11]", "> day 11")), y = value)) +
  geom_boxplot() +
  geom_smooth(method = "lm", aes(group=group), alpha = 0.3) +
  # facet_wrap(~cv_cell_type, ncol = 4, scales = "free") +
  facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
  stat_compare_means(comparisons = time_nuns_comparisons) +
  labs(title = "Difference in immune cell abundance over time after vaccination") +
  my_theme

plots <- lapply(deconv_methods_nuns, function(m) {
  ggplot(nuns_naive_time_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
         aes(x = factor(day_group_2, levels = c("day 0", "day [7,11]", "> day 11")), y = value)) +
    geom_boxplot() +
    geom_smooth(aes(group=group), alpha = 0.3) +
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free") +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
    # facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
    labs(title = m, x = "day_group", y = "value") +
    my_theme
})

# !!! change number of conditional_plots and heights depending on number of methods used
figure <- ggarrange(plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[2]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[3]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                    plots[[4]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                    plots[[5]],
                    heights = c(1,1,1,1,1.4),
                    labels = names(plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(figure, 
                top = text_grob("Difference in immune cell abundance over time after vaccination", size = 16),
                fig.lab = "C", fig.lab.size = 18, fig.lab.face = "bold")

