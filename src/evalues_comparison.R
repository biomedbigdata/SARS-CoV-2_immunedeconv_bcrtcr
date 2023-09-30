# script to plot Evalues from BCR TCR BLAST analysis
library(ggplot2)
library(data.table)

e_values <- fread("data/e_values.csv")
e_values[, `E-value` := gsub(",", ".", `E-value`)]
e_values$`E-value` <- as.numeric(e_values$`E-value`)

## plot the pvalues
ggplot(e_values, aes(x = `Sequencing depth`, y = `E-value`))+ 
  geom_boxplot() +
  #facet_grid(~variant) +
  labs(title = "E values from BLAST results per sequencing depth",
       x = "Sequencing Depth", y = "E value") +
  my_theme


ggplot(e_values[`E-value` < 1], aes(x = `Sequencing depth`, y = `E-value`))+ 
  geom_boxplot() +
  #facet_grid(~variant) +
  labs(title = "E values from BLAST results per sequencing depth",
       x = "Sequencing Depth", y = "E value") +
  my_theme
