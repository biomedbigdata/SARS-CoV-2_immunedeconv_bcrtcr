## script for creating plots from deconvolution result with scores and fractions
# computes the desired result if not already computed
# creates a plot for every cell type and saves them in the specified dir

library(immunedeconv)
library(ggplot2)
library(data.table)
library(pheatmap)
library(ggpubr)
library(stringr)
library(dplyr)

mytheme <- theme(plot.title = element_text(size = 11))


get_result_table <- function(method){
  # check if results were already computed, otherwise compute
  if (!exists(paste0(method, "_result"))){
    print(paste0(method, "_result is computed"))
    if (method == "mcp_counter"){
      result <- immunedeconv::deconvolute_mcp_counter(tpms)
    } else{
      result <- immunedeconv::deconvolute(tpms, method, tumor = FALSE)
    }
    assign(paste0(method, "_result"), result, env=.GlobalEnv)
  } else {
    print(paste0(method, "_result is loaded from work space"))
    result <- get(paste0(method, "_result"))
  }
  result_dt <- as.data.table(result) # convert to data table for easy access
  
  # add cell type column if result was a matrix with cell types as rownames
  if (!"cell_type" %in% colnames(result_dt)){
    cell_types <- rownames(result)
    result_dt[, cell_type := cell_types]
  }
  
  return(result_dt)
}



# function to create the desired data table format (long/dt or wide/matrix)
create_result_table <- function(method, as_matrix = F, exclude_gamma = F){
  # create the result and get the matrix (cell_type x sample)
  result_dt <- get_result_table(method)
  
  # make long data table
  result_dt_long <- as.data.table(gather(result_dt, sample, score, -c(cell_type)))
  
  # merge result_dt with the meta data containing info about groups
  if (old_id){
    result_dt_long <- merge(result_dt_long, full_metadata,
                            by.x = "sample", by.y = "old_id", all.x = T, allow.cartesian = T)
  } else {
    result_dt_long <- merge(result_dt_long, full_metadata,
                            by.x = "sample", by.y = "sample_id", all.x = T, allow.cartesian = T)
  }
  
  if (exclude_gamma) {
    result_dt_long <- result_dt_long[group != "Gamma"]
  }
  
  result_long_intermediate <<- result_dt_long
  
  # return matrix if data should be in format matrix instead of long data table
  if (as_matrix){
    # convert back to wide format but sample x cell_type
    result_dt_wide <- dcast(result_dt_long, 
                            ... ~ cell_type, value.var = 'score')
    return(result_dt_wide)
  } 
  # else return the long format
  return(result_dt_long[order(group)])
}


# function to create the result table merged with time relevant data for alpha, alpha_ek and gamma plus seronegative for all time groups
create_time_result_table <- function(method, exclude_gamma = F){
  # create the result and get the matrix (cell_type x sample)
  result_dt <- get_result_table(method)
  
  # make long data table
  result_dt_long <- as.data.table(gather(result_dt, sample, score, -c(cell_type)))
  
  # merge result_dt with the meta data containing info about groups
  result_dt_long <- merge(result_dt_long, alpha_gamma_sero_meta_test,
                          by.x = "sample", by.y = "sample_id", all.x = T, allow.cartesian = T)
  result_dt_long <- result_dt_long[group %in% c("Alpha", "Alpha_EK", "Gamma", "Seronegative")]
  if (exclude_gamma) {
    result_dt_long <- result_dt_long[group != "Gamma"]
  }
  
  # return the long format
  return(result_dt_long[order(group)])
}



# funtion to add numeric day column to nuns meta data
add_days_to_meta <- function(meta_dt){
  meta_dt$num_day <- lapply(meta_dt$sampling, function(x) str_extract(x, "[0-9]+"))
  meta_dt
}


# function to create plots for each cell type returned by deconvolution method that returns scores
create_score_plots <- function(method, dir, color_by_sampling = F){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))

  result_dt <- create_result_table(method)
  color_by <- ifelse(color_by_sampling, "sampling", "group")
  
  
  # generate plot for every cell type, save to xcell_plots/<cell type>.png
  for (ct in unique(result_dt$cell_type)){
    # generate plot data table
    plot_dt <- result_dt[cell_type == ct]
    
    # generate plot
    p <- ggplot(plot_dt) +
      geom_point(aes(x=sample, y=score, color=get(color_by)), size=4) +
      scale_color_brewer(palette="Paired") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
      labs(title = paste0(method, " scores for ", ct))
    
    num_samples <- length(unique(plot_dt$sample)) # for plot size
    # create filename
    cell_type <- gsub("\\+", "", ct) # remove "+" from cell type if necessary
    filename <- paste0(gsub(" ", "_", cell_type), ".png") # paste filename from cell type
    # save plot
    # ggsave(filename, plot=p, path=dir, height=num_samples/3.5, width=num_samples/6, units="cm")
    ggsave(filename, plot=p, path=dir, height=num_samples/3.5, units="cm")
  }

}



# function to create and save bar plot from deconvolution method that returns fractions
create_fraction_plot <- function(method, dir){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  
  result_dt <- create_result_table(method)
  colnames(result_dt)[which(names(result_dt) == "score")] <- "fraction"
  
  # generate plot data table
  plot_dt <- result_dt 
  
  # plot as stacked bar chart
  p <- ggplot(plot_dt[order(group)]) +
    geom_bar(aes(x=sample, y=fraction, fill=cell_type), stat='identity') +
    coord_flip() +
    # facet_wrap(~group, scales = "free") +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(result_dt))) +
    labs(title = paste0("cell type fractions for each sample computed with ", method))

  num_samples <- length(unique(plot_dt$sample)) # for plot size
  # create filename
  filename <- paste0(gsub(" ", "_", method), ".png") # paste filename from cell type
  # save plot
  ggsave(filename, plot=p, path=dir, height=num_samples/3.5, width=num_samples/4.5, units="cm")
}


# function to create a sample x cell type matrix with fraction or score values
create_hm <- function(method, dir, exclude_gamma = F, color_by_sampling = F){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  
  result_dt <- create_result_table(method, as_matrix = TRUE, exclude_gamma = exclude_gamma)
  
  # generate data table for the annotation
  if (color_by_sampling){
    # annoColors <<- brewer.pal(length(unique(full_metadata$num_day)), "Paired")
    names(annoColors) <- unique(full_metadata$num_day)

    # create annotation data table
    annotation_dt <- as.data.frame(result_dt[, .(num_day)])
  } else {
    # annoColors <- NULL
    annotation_dt <- as.data.frame(result_dt[, .(group)])
  }
  rownames(annotation_dt) <- result_dt$sample
  # TODO: add annotation color if necessary
  
  # generate plot data table
  plot_mat <- as.matrix(result_dt[, -c("sample", "group", "sampling")])
  rownames(plot_mat) <- result_dt$sample
  class(plot_mat) <- "numeric"

  # filename for saving the heatmap
  filename <- paste0(gsub(" ", "_", method), "_hm.png") # paste filename from method
  
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  # plot as hm
  pheatmap(plot_mat, annotation_row = annotation_dt,
           cluster_cols = FALSE, 
           color = brewer.pal(n = 9, name = "Greys"),
           # annotation_colors = annoColors,
           border_color = NA,
           fontsize_row = 8, 
           main = paste0("cell type ", type, " for each sample computed with ", method),
           filename = paste0(dir, filename), width = ifelse(ncol(plot_mat)<15, ncol(plot_mat), ncol(plot_mat)/3), height = nrow(plot_mat)/8)
}



# function to create a boxplot
create_time_boxplot_nuns <- function(method, dir, color = "group", reference = ".all."){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  
  result_dt <<- create_result_table(method, as_matrix = F, exclude_gamma = F)
  
  p <- ggplot(result_dt, aes(x = cell_type, y = score, color = get(color))) +
    geom_boxplot() +
    # facet_wrap(~cell_type, scales = "free_y") +
    stat_compare_means(label = "p.signif", method = "t.test", hide.ns = T) +
    # stat_compare_means(symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
    labs(title = paste0("Difference in ", type, " from ", method, " for ", dataset, " per cell type by time group")) +
    rotate_x_text() +
    mytheme
  
  filename <- paste0(gsub(" ", "_", method), "_boxplot_", ".png") 
  ggsave(filename, plot=p, path=dir, width = 17, height = 15, units = "cm") # TODO: different sizes!
  # ggsave(filename, plot=p, path=dir, width=ifelse(num_cell_types < 15, num_cell_types*2.4, num_cell_types * 1.6 ), units="cm")
}


# create a boxplot for a method withx = cell types, y = scores coloured by group
create_boxplot <- function(method, dir, reference = ".all."){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  result_dt <<- create_result_table(method, as_matrix = F, exclude_gamma = F)

  title <- paste0("Difference in ", type, " from ", method, " per cell type for", dataset)
  
  p <- ggplot(result_dt, aes(x = cell_type, y = score, color = group)) +
    geom_boxplot() +
    stat_compare_means(label = "p.signif", method = "t.test", ref.group = reference) +
    # facet_grid(. ~ day_group, scales = "free_y") +
    # stat_compare_means(symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
    labs(title = title,
         x = color) +
    rotate_x_text() +
    mytheme
  
  filename <- paste0( gsub(" ", "_", method), "_boxplot", ".png")
  ggsave(filename, plot=p, path=dir, width = 17, height = 15, units = "cm")  # TODO: different sizes!
  # ggsave(filename, plot=p, path=dir, width=ifelse(num_cell_types < 15, num_cell_types*2.4, num_cell_types * 1.6 ), units="cm")
}


# function to create boxplots for a method for each celltype with "color" on the x axis
create_boxplots <- function(method, dir, color = "group", reference = ".all.", time = F, variants = F){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  if (time & variants){
    result_dt <<- create_time_result_table(method, exclude_gamma = F)
  } else {
    result_dt <<- create_result_table(method, as_matrix = F, exclude_gamma = F)
  }
  
  for (ct in unique(result_dt$cell_type)){
    title <- ifelse(time, paste0("Difference in ", ct ," ", type, " from ", method, " over time for ", dataset), 
                    paste0("Difference in ", ct ," ", type, " from ", method, " for ", dataset))
    
    p <- ggplot(result_dt[cell_type == ct], aes(x = get(color), y = score)) +
      geom_boxplot() +
      stat_compare_means(label = "p.signif", method = "t.test", ref.group = reference) +
      # facet_grid(. ~ day_group, scales = "free_y") +
      # stat_compare_means(symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
      labs(title = title,
           x = color) +
      rotate_x_text() +
      mytheme
    
    filename <- paste0(gsub("[ +]", "_", ct), "_boxplot_", gsub(" ", "_", method), ".png") 
    ggsave(filename, plot=p, path=dir, width = 17, height = 15, units = "cm")
    # ggsave(filename, plot=p, path=dir, width=ifelse(num_cell_types < 15, num_cell_types*2.4, num_cell_types * 1.6 ), units="cm")
  }
}

# function to create boxplots annotated with 
create_time_boxplot <- function(method, dir, color = "group", reference = ".all.", variants = F, exclude_gamma = F){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  if (variants){
    result_dt <<- create_time_result_table(method, exclude_gamma = exclude_gamma)
  } else {
    result_dt <<- create_result_table(method, as_matrix = F, exclude_gamma = exclude_gamma)
  }
  for (day in unique(result_dt$day_group)){
    for (ct in unique(result_dt$cell_type)){
      title <- paste0("Difference in ", ct ," ", type, " from ", method, " for day ",day, " from ", dataset)
      
      p <- ggplot(result_dt[cell_type == ct & day_group == day], aes(x = get(color), y = score)) +
        geom_boxplot() +
        stat_compare_means(label = "p.signif", method = "t.test", ref.group = reference) +
        # facet_grid(. ~ day_group, scales = "free_y") +
        # stat_compare_means(symnum.args = list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))) +
        labs(title = title,
             x = color) +
        rotate_x_text() +
        mytheme
      
      filename <- paste0(day, "_", gsub("[ +]", "_", ct), "_boxplot_", gsub(" ", "_", method), ".png") 
      ggsave(filename, plot=p, path=dir, width = 17, height = 15, units = "cm")
      # ggsave(filename, plot=p, path=dir, width=ifelse(num_cell_types < 15, num_cell_types*2.4, num_cell_types * 1.6 ), units="cm")
    }
  }
}



