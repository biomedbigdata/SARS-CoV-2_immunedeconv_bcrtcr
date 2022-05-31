## script for creating plots from deconvolution result with scores and fractions
# computes the desired result if not already computed
# creates a plot for every cell type and saves them in the specified dir

library(immunedeconv)
library(ggplot2)
library(data.table)


create_result_table <- function(method){
  # check if results were already computed, otherwise compute
  if (!exists(paste0(method, "_result"))){
    print(paste0(method, "_result is computed"))
    if (method == "mpc_counter"){
      result <- immunedeconv::deconvolute_mcp_counter(tpms)
    } else{
      result <- immunedeconv::deconvolute(tpms, method)
    }
    assign(paste0(method, "_result"), result)
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
  
  # make long data table
  result_dt <- as.data.table(gather(result_dt, sample, score, -c(cell_type)))
  
  
  # merge result_dt with the meta data containing info about groups
  result_dt <- merge(result_dt, full_metadata[, .(old_id, group)],
                     by.x = "sample", by.y = "old_id", all.x = T, allow.cartesian = T)
  
  return(result_dt[order(group)])
}


# function to create plots for each cell type returned by deconvolution method that returns scores
create_score_plots <- function(method, dir){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))

  result_dt <- create_result_table(method)
  
  # generate plot for every cell type, save to xcell_plots/<cell type>.png
  for (ct in unique(result_dt$cell_type)){
    # generate plot data table
    plot_dt <- result_dt[cell_type == ct]
    
    # generate plot
    p <- ggplot(plot_dt) +
      geom_point(aes(x=sample, y=score, color=group), size=4) +
      scale_color_brewer(palette="Paired") +
      coord_flip() +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    
    # create filename
    cell_type <- gsub("\\+", "", ct) # remove "+" from cell type if necessary
    filename <- paste0(gsub(" ", "_", cell_type), ".png") # paste filename from cell type
    # save plot
    ggsave(filename, plot=p, path=dir, height=40, width=20, units="cm")
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
    scale_x_discrete(limits = rev(levels(result_dt)))


  # create filename
  filename <- paste0(gsub(" ", "_", method), ".png") # paste filename from cell type
  # save plot
  ggsave(filename, plot=p, path=dir, height=40, width=30, units="cm")
}




