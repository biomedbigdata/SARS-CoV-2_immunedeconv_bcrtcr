#########################
# create combined plots # 
#########################
library(data.table)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(RColorBrewer)

load("data/variants_omicron_ischgl_10m_wo_omicron1_deconv.RData")
# load("data/nuns_deconv.RData")
# load("data/variants_ischgl_deconv_with_cibersortx.RData")
# load("data/nuns_deconv_with_cibersortx.RData")
variants_omicron_ischgl_result <- variants_omicron_ischgl_10m_wo_omicron1_result

# new deconv names 
variants_omicron_ischgl_result[, method := ifelse(method == "quantiseq", "quanTIseq", ifelse(method == "mcp_counter", "MCP-counter", ifelse(method=="xcell", "xCell", ifelse(method == "epic", "EPIC", ifelse(method == "cibersort_abs", "CIBERSORT", method)))))]
# nuns_results[, method := ifelse(method == "quantiseq", "quanTIseq", ifelse(method == "mcp_counter", "MCP-counter", ifelse(method=="xcell", "xCell", ifelse(method == "epic", "EPIC", method))))]



my_theme <- theme(panel.background = element_rect(fill = "white", colour = "grey"),
                  panel.grid.major = element_line(colour = "lightgrey"),
                  panel.grid.minor = element_line(colour = "lightgrey"),
                  text = element_text(size = 20),
                  plot.title = element_text(size = 20))

variants_omicron_ischgl_colors <- list("alpha" = "#E31A1C", "alpha_ek" = "#FF7F00", "gamma" = "#A6CEE3", "ba1" = "#33A02C", "ba2" = "#B2DF8A", "sero" = "#6A3D9A")

deconv_methods <- c("quanTIseq", "MCP-counter", "EPIC", "xCell")

# conditional plot variants_ischgl
conditional_plot_dt <- variants_omicron_ischgl_result[method %in% deconv_methods & !is.na(group)]
# conditional_plot_dt <- variants_omicron_ischgl_result[num_day < 6]  # filter for first 10 days

# conditional_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
# conditional_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]
conditional_plot_dt[, value := ifelse(method == "MCP-counter", log(value), value)]
conditional_plot_dt[group == "Alpha_EK"]$group <- "Alpha+EK"

conditional_plot_dt$group <- factor(conditional_plot_dt$group, levels = c("Alpha", "Alpha+EK", "Gamma", "BA.1", "BA.2", "Seronegative"))

conditional_variants_comparisons = list(c("Alpha", "Seronegative"), c("Alpha+EK", "Seronegative"), c("Gamma", "Seronegative"), c("BA.1", "Seronegative"), c("BA.2", "Seronegative"))

# ggplot(conditional_plot_dt[cv_cell_type %in% c("B cell", "Neutrophil", "T cell CD4+", "T cell CD8+")], 
#        aes(x = group, y = value, fill = group)) +
#   geom_boxplot() +
#   facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
#   stat_compare_means(comparisons = conditional_variants_comparisons) +
#   labs(title = "Differences in immune cell abundances between infected and healthy") +
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
#   my_theme



colnames(conditional_plot_dt)[colnames(conditional_plot_dt) == "group"] <- "Group"


# variants conditional plots for each method combined to get free scales
conditional_plots <- lapply(deconv_methods, function(m) {
  ggplot(conditional_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
         aes(x = Group, y = value, fill = Group)) +
    geom_boxplot() + 
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free_y") +
    stat_compare_means(comparisons = conditional_variants_comparisons, vjust = 1.2) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
    scale_fill_manual(values = c(variants_omicron_ischgl_colors[["alpha"]], variants_omicron_ischgl_colors[["alpha_ek"]],  variants_omicron_ischgl_colors[["gamma"]],variants_omicron_ischgl_colors[["ba1"]], variants_omicron_ischgl_colors[["ba2"]], variants_omicron_ischgl_colors[["sero"]])) +
    labs(title = m, x = "Group", y = "Value") +
    my_theme + theme(axis.title.x = element_blank(), axis.text.x = element_blank())
})

# !!! change number of conditional_plots and heights depending on number of methods used
conditional_figure <- ggarrange(conditional_plots[[1]], 
                                conditional_plots[[2]], 
                                conditional_plots[[3]], 
                                conditional_plots[[4]],
                               #  conditional_plots[[5]],
                    heights = c(1,1,1,1),
                    labels = names(conditional_plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(conditional_figure, 
                top = text_grob("Differences in immune cell abundances between infected and healthy (sequencing depth 10 million)", size = 22),  
                fig.lab.size = 22, fig.lab.face = "bold")


## get the p-values
unique(ggplot_build(conditional_plots[[1]])$data[[2]]$annotation)

p_val_df <- rbindlist(lapply(c(1:4), function(i){
  m <- deconv_methods[[i]]
  p_values <- (ggplot_build(conditional_plots[[i]])$data[[2]]$annotation)
  print(p_values)
  data.table(
    method = m,
    variant = rep(c(rep("Alpha",3), rep("Alpha_EK",3), rep("Gamma",3), rep("BA.1",3), rep("BA.2",3)), 4),
    cell_type = c(rep("B cell", 15), rep("Neutrophil", 15), rep("T cell CD4+", 15), rep("T cell CD8+", 15)),
    pvalue = p_values
  )
})) %>% unique()

p_val_df[, seq_depth := "full"]

# all_p_val_df <- p_val_df
all_p_val_df <- rbind(all_p_val_df, p_val_df)
all_p_val_df[pvalue == "p < 2.22e-16"]$pvalue <- "2.22e-16"
all_p_val_df$pvalue <- as.numeric(all_p_val_df$pvalue)

# now repeat with seq depth 50 M and 10 M
# write.csv(all_p_val_df, file = "data/conditional_p_values.csv")

all_p_val_df[seq_depth == "50 M"]$seq_depth <- "50 million"
all_p_val_df[seq_depth == "10 M"]$seq_depth <- "10 million"

## plot the pvalues
ggplot(all_p_val_df, aes(x = seq_depth, y = pvalue))+ 
  geom_boxplot() +
  #facet_grid(~variant) +
  labs(title = "P-values for difference between deconvolution results of variants vs seronegative",
       x = "Sequencing depth", y = "P-value") +
  my_theme

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

load("data/variants_omicron_ischgl_50m_wo_omicron1_deconv_time.RData")
#load("data/batch_corrected/variants_omicron_ischgl_10m_corrected_deconv_time.RData")
variants_omicron_ischgl_result <- variants_omicron_ischgl_50m_wo_omicron1_time_result

variants_omicron_ischgl_result[, method := ifelse(method == "quantiseq", "quanTIseq", ifelse(method == "mcp_counter", "MCP-counter", ifelse(method=="xcell", "xCell", ifelse(method == "epic", "EPIC", ifelse(method == "cibersort_abs", "CIBERSORT", method)))))]

# time series variants_ischgl --> exclude Gamma samples, use EPIC and MCPcounter as they are best for time series data and the cell types
variants_time_plot_dt <- variants_omicron_ischgl_result[method %in% deconv_methods & group != "Gamma" & !is.na(group)]
#variants_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
#variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
#variants_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
#variants_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]# variants_time_plot_dt$day_group <- factor(variants_time_plot_dt$day_group, levels = c("<= day 5", "<= day 10", "<= day 15", "<= day 30", "> day 30" ))

variants_time_plot_dt[group == "Alpha_EK"]$group <- "Alpha+EK"

variants_time_plot_dt$group <- factor(variants_time_plot_dt$group, levels = c("Alpha", "Alpha+EK", "Gamma", "BA.1", "BA.2", "Seronegative"))

colnames(variants_time_plot_dt)[colnames(variants_time_plot_dt) == "group"] <- "Group"
variants_time_plot_dt[day_group == "> day 30"]$day_group <- "> Day 30"
variants_time_plot_dt[, day_group := gsub("^day ", "Day ", day_group)]
variants_time_plot_dt <- variants_time_plot_dt[!is.na(day_group)]

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
variants_time_plot_dt[, value := ifelse(method == "MCP-counter", log(value), value)]
# ggplot(variants_time_plot_dt[cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
#        aes(x = factor(day_group, levels = c("day [0,5]", "day [6,10]", "day [11,15]", "day [16,30]", "> day 30")), y = value, color = group)) +
#   geom_boxplot() + 
#   geom_smooth(method = "lm", aes(group=group), alpha = 0.3) +
#   facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
#   theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
#   labs(title = "Difference in immune cell abundance over time after infection",
#        x = "day_group") +
#   my_theme


# variants time plots for each method combined to get free scales
plots <- lapply(deconv_methods, function(m) {
  ggplot(variants_time_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
                  aes(x = factor(day_group, levels = c("Day [0,5]", "Day [6,10]", "Day [11,15]", "Day [16,30]", "> Day 30")), y = value, color = Group)) +
    geom_boxplot() + 
    geom_smooth(method = "lm", aes(group=Group), alpha = 0.3) +
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free") +
    scale_color_manual(values = c(variants_omicron_ischgl_colors[["alpha"]], variants_omicron_ischgl_colors[["alpha_ek"]], variants_omicron_ischgl_colors[["ba1"]], variants_omicron_ischgl_colors[["ba2"]], variants_omicron_ischgl_colors[["sero"]])) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) + # TODO: change color
    labs(title = m, x = "Day Group", y = "Value") +
    my_theme
})

# !!! change number of conditional_plots and heights depending on number of methods used
figure <- ggarrange(plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[2]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[3]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    #plots[[4]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[4]] + theme(axis.title.x = element_blank()),
                    heights = c(1,1,1,1.7),
                    labels = names(plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(figure, 
                top = text_grob("Difference in immune cell abundance over time after infection", size = 22),
                fig.lab.size = 22, fig.lab.face = "bold")


# linear model to get the p values of the slopes
slope_p_values <- rbindlist(lapply(deconv_methods, function(m){
  rbindlist(lapply(c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+"), function(c){
    rbindlist(lapply(c("Alpha", "Alpha+EK", "BA.1", "BA.2"), function(v){
      df <- variants_time_plot_dt[method == m & cv_cell_type == c & Group == v]
      df <-df %>%
        mutate(day_group_num = case_when(
          day_group == "Day [0,5]"    ~ 1,
          day_group == "Day [6,10]" ~ 2,
          day_group == "Day [11,15]"      ~ 3,
          day_group == "Day [16,30]"      ~ 4,
          day_group == "> Day 30"      ~ 5,
          TRUE               ~ 0 # Default case if none of the above conditions are met
        ))
      linearModel <- lm(day_group_num ~ value, data = df)
      slopePValue <- summary(linearModel)$coefficients["value", "Pr(>|t|)"]
      data.table(Method = m, Cell_type = c, Group = v, P.Value = slopePValue)
    }))

  }))
}))

p_value_slopes <- dcast(slope_p_values, Method + Cell_type ~ Group)
# write.csv(p_value_slopes, file = "results/slopes_p_values_m50.csv", row.names = F)




# time series nuns
deconv_methods_nuns <- c("quanTIseq","MCPcounter",  "EPIC", "xCell")
nuns_naive_time_plot_dt <- nuns_results[method %in% deconv_methods_nuns & group == "Naive"]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "Neutrophils", "Neutrophil", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+", "T cells", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "T cell CD4+ (non-regulatory)", "T cells", cell_type)]
# nuns_naive_time_plot_dt[, cell_type := ifelse(cell_type == "B lineage", "B cell", cell_type)]


# nuns_naive_time_plot_dt[, day_group := ifelse(num_day <= 5, "<= day 5", 
#                                               ifelse(num_day <= 11, "<= day 11",
#                                                      ifelse(num_day <= 40, "<= day 40", "> day 40")))]
# nuns_naive_time_plot_dt$day_group <- factor(nuns_naive_time_plot_dt$day_group, levels = c("<= day 5", "<= day 11", "<= day 40", "> day 40"))
nuns_naive_time_plot_dt[, value := ifelse(method == "MCPcounter",  log(value), value)]
nuns_naive_time_plot_dt[, value := ifelse(method == "xCell" & cv_cell_type == "Neutrophil", signif(value, digits = 2), value)]
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

scaleFUN <- function(x) ifelse(x < 8e-18, signif(x,1), x) # function for scaling y axis labels

plots <- lapply(deconv_methods_nuns, function(m) {
  ggplot(nuns_naive_time_plot_dt[method == m & cv_cell_type %in% c("Neutrophil", "B cell", "T cell CD4+", "T cell CD8+")], 
         aes(x = factor(day_group_2, levels = c("day 0", "day [7,11]", "> day 11")), y = value)) +
    geom_boxplot() +
    geom_smooth(aes(group=group), alpha = 0.3) +
    facet_wrap(~cv_cell_type, ncol = 4, scales = "free") +
    scale_y_continuous(labels=scaleFUN) +
    theme(axis.text.x=element_text(angle = 90, vjust = 0.5)) +
    # facet_grid(rows = vars(method), cols = vars(cv_cell_type), scales = "free") +
    labs(title = m, x = "day_group", y = "value") +
    my_theme
})

# !!! change number of conditional_plots and heights depending on number of methods used
figure <- ggarrange(plots[[1]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[2]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()), 
                    plots[[3]] + theme(axis.title.x = element_blank(), axis.text.x = element_blank()),
                    plots[[4]],
                    heights = c(1,1,1,1.7),
                    labels = names(plots),
                    ncol = 1,
                    align = "v",
                    common.legend = T,
                    legend = "bottom")
annotate_figure(figure, 
                top = text_grob("Difference in immune cell abundance over time after vaccination", size = 22),
                fig.lab = "C", fig.lab.size = 22, fig.lab.face = "bold")

