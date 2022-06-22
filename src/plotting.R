## script for creating plots from deconvolution result with scores and fractions
# computes the desired result if not already computed
# creates a plot for every cell type and saves them in the specified dir

library(immunedeconv)
library(ggplot2)
library(data.table)
library(pheatmap)


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
    result_dt_long <- merge(result_dt_long, full_metadata[, .(old_id, group, sampling)],
                            by.x = "sample", by.y = "old_id", all.x = T, allow.cartesian = T)
  } else {
    result_dt_long <- merge(result_dt_long, full_metadata[, .(sample_id, group)],
                            by.x = "sample", by.y = "sample_id", all.x = T, allow.cartesian = T)
  }
  
  if (exclude_gamma) {
    result_dt_long <- result_dt_long[group != "Gamma"]
  }
  
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
  
  result_dt <<- create_result_table(method, as_matrix = TRUE, exclude_gamma = exclude_gamma)
  
  # generate data table for the annotation
  if (color_by_sampling){
    annotation_dt <- as.data.frame(result_dt[, .(sampling)])
  } else {
    annotation_dt <- as.data.frame(result_dt[, .(group)])
  }
  rownames(annotation_dt) <- result_dt$sample
  # TODO: add annotation color if necessary
  
  # generate plot data table
  plot_mat <- as.matrix(result_dt[, -c("sample", "group", "sampling")])
  rownames(plot_mat) <- result_dt$sample
  class(plot_mat) <- "numeric"
  

  # filename for saving the heatmap
  filename <- paste0(gsub(" ", "_", method), "_hm.png") # paste filename from cell type
  
  type <- ifelse(method %in% c("quantiseq", "epic"), "fractions", "scores")
  # plot as hm
  pheatmap(plot_mat, annotation_row = annotation_dt,
           cluster_cols = FALSE, 
           color = brewer.pal(n = 9, name = "Greys"),
           border_color = NA,
           fontsize_row = 8, 
           main = paste0("cell type ", type, " for each sample computed with ", method),
           filename = paste0(dir, filename), width = ifelse(ncol(plot_mat)<15, ncol(plot_mat), ncol(plot_mat)/3), height = nrow(plot_mat)/8)
}

