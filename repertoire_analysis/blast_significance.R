# plotting blast significance

library(ggplot2)

M10_signif <- c(4e-05,2e-04,0.002,0.003,0.001,0.002,0.007)

all_signif <- c(6e-05,0.050,0.010,0.27, 2.6, 0.006, 5e-06)


df <- rbindlist(list(
  data.table(type = "10 Mio", values = M10_signif),
  data.table(type="full", values = all_signif)))

ggplot(df[values<1], aes(x=type, y=values)) + geom_boxplot() +
  labs(title = "Significance of BLAST results for total results and downsampled results",
       x = "sequencing depth", y = "BLAST E-value")

