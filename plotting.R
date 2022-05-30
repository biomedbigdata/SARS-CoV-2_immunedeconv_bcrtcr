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
    result <- immunedeconv::deconvolute(tpms, method)
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
  
  return(result_dt)
}


# function to create plots for each cell type returned by deconvolution method that returns scores
create_score_plots <- function(method, dir){
  # create a new dir for result plots if it does not exist yet
  dir.create(file.path(dir))

  result_dt <- create_result_table(method)
  
  # generate plot for every cell type, save to xcell_plots/<cell type>.png
  for (ct in result_dt$cell_type){
    # generate plot data table
    plot_dt <- result_dt[cell_type == ct] %>%
      gather(sample, score, -cell_type)
    plot_dt <- as.data.table(plot_dt)
    plot_dt <- plot_dt[, group := full_metadata[sample==sample_id]$group]
    
    # generate plot
    p <- ggplot(plot_dt) +
      geom_point(aes(x=sample, y=score, color=cell_type), size=4) +
      # facet_wrap(~cell_type, scales="free_x", ncol=3) +
      scale_color_brewer(palette="Paired", guide=FALSE) +
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
  
  # generate plot data table
  plot_dt <- result_dt %>%
    gather(sample, fraction, -cell_type)
  
  # print(head(plot_dt))
  plot_dt <- merge(as.data.table(plot_dt), full_metadata[, .(sample_id, group)], by.x="sample", by.y="sample_id", all.x=T)
  # plot_dt <- plot_dt[, group := full_metadata[sample_id==sample]$group]
  print(plot_dt)
  
  # plot as stacked bar chart
  p <- ggplot(plot_dt) +
    geom_bar(aes(x=sample, y=fraction, fill=cell_type), stat='identity') +
    coord_flip() +
    scale_fill_brewer(palette="Paired") +
    scale_x_discrete(limits = rev(levels(result_dt)))
  
  # create filename
  filename <- paste0(gsub(" ", "_", method), ".png") # paste filename from cell type
  # save plot
  ggsave(filename, plot=p, path=dir, height=40, width=20, units="cm")
}




