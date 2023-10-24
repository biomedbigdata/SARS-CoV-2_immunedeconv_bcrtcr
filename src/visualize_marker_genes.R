# for plotting the marker genes
library(data.table)
library(ggplot2)
library(readxl)
library(dplyr)
library(reshape2)
library(ComplexUpset)
library(dendextend)
library(tidyverse)

### signature/deconvolution based

quantiseq_signature <- read.csv("data/signatures/quantiseq_TIL10_signature.txt", sep = "\t", row.names = "ID")
colnames(quantiseq_signature) <- c('B cells', 'Macrophages M1', 'Macrophages M2', 'Monocytes', 'Neutrophils', 'NK cells', 'CD4 T cells', 'CD8 T cells', 'Tregs', 'Dendritic cells')

epic_signature <- as.data.frame(read_xlsx("data/signatures/epic_gene_expr_reference.xlsx", range = "A3:G50000"))
epic_signature <- epic_signature %>%
  filter(rowSums(select(., -Gene) != 0) > 0)
rownames(epic_signature) <- epic_signature$Gene
epic_signature$Gene <- NULL
colnames(epic_signature) <- c('B cells','CD4 T cells','CD8 T cells','Monocytes','Neutrophils','NK cells')


# plot quantiseq heatmap
heatmap(as.matrix(quantiseq_signature), Rowv = NA, Colv = NA) + 
  title(main = "quanTIseq signature matrix")


# plot epic heatmap
# epic_signature[apply(epic_signature, 1, function(row) any(row > 1000)), ] %>% dim
# epic_top_315 <- epic_signature[apply(epic_signature, 1, function(row) any(row > 1000)), ]
heatmap(as.matrix(epic_signature), Rowv = F, Colv = F) +
  title(main = "EPIC signature matrix")


# skip the renaming of columns for save for RData 
# merge both data frames
colnames(epic_signature) <- paste(colnames(epic_signature), "(EPIC)")
colnames(quantiseq_signature) <- paste(colnames(quantiseq_signature), "(quanTIseq)")

epic_signature$Gene <- rownames(epic_signature)
quantiseq_signature$Gene <- rownames(quantiseq_signature)

merged_signatures <- merge(quantiseq_signature, epic_signature, by = "Gene")
rownames(merged_signatures) <- merged_signatures$Gene
merged_signatures$Gene <- NULL

heatmap(as.matrix(merged_signatures), na.rm = T, Rowv=F, Colv = F) +
  title(main = "Signatures for EPIC and quanTIseq")


### Marker gene based
mcp_counter_tms <- as.data.frame(read_xlsx("data/signatures/mcpcounter_tm_expression.xlsx", sheet = "Table S2", range = "A2:C550"))
# mcp_counter_tms$Symbol <- sapply(mcp_counter_tms$Symbol, function(x) trimws(unlist(strsplit(x, split = "///", fixed = T))))
mcp_counter_tms$Probeset <- NULL
mcp_counter_tms <- mcp_counter_tms %>% filter(!is.na(Population))

mcp_counter_tms$value <- 1

mcp_matrix <- as.data.frame(dcast.data.table(as.data.table(mcp_counter_tms), Symbol + value ~ Population, fun.aggregate = sum))
mcp_matrix$value <- NULL
# mcp_matrix$Symbol <- sapply(mcp_matrix$Symbol, function(x) trimws(unlist(strsplit(x, split = "///", fixed = T))))

colnames(mcp_matrix) <- c('Symbol', 'B cells', 'CD8 T cells', 'Cytotoxic lymphocytes', 'Endothelial cells', 'Fibroblasts', 'Monocytes', 'Myeloid dendritic cells', 'NK cells', 'Neutrophils', 'T cells')

# mcp_matrix[duplicated(Symbol)]
rownames(mcp_matrix) <- mcp_matrix$Symbol
mcp_matrix <- mcp_matrix[mcp_matrix$Symbol != "NA",]
mcp_matrix$Symbol <- NULL


heatmap(as.matrix(mcp_matrix)) +
  title(main = "Transcriptomic markers for MCP-counter")

## xCell
xcell_tms <- as.data.frame(read_xlsx("data/signatures/xcell_signature.xlsx"))
xcell_tms <- xcell_tms %>% 
  rowwise() %>% 
  mutate(list_col = list(c_across(-c(1,2)))) %>%
  select(1:2, list_col) %>%
  as.data.frame()

# Remove NA from the lists
xcell_tms$list_col <- map(xcell_tms$list_col, ~ .[!is.na(.)])

# split the entries in the first column to keep only the cell type
xcell_tms <- xcell_tms %>%
  mutate(Celltype_Source_ID = str_split(Celltype_Source_ID, pattern = "_", simplify = TRUE)[, 1])

# combine signatures from multiple sources for one cell tye into one entry
xcell_tms <- xcell_tms %>%
  group_by(Celltype_Source_ID) %>%  # Group by the first column
  summarize(
    Genes = list(unlist(list_col))  # Combine the lists
  ) %>%
  ungroup()

colnames(xcell_tms)[colnames(xcell_tms) == "Celltype_Source_ID"] <- "Cell type"

# Expand the tibble so that each gene gets its own row
expanded_data <- xcell_tms %>%
  unnest(Genes)

# Create the contingency table
contingency_table <- with(expanded_data, table(Genes, `Cell type`))

# Convert table to a matrix
xcell_tms <- as.matrix(contingency_table)

## plot
heatmap(xcell_tms) +
  title(main = "Combined transcriptomic markers for xCell")

# make to data frame
df <- data.frame(matrix(unlist(xcell_tms), nrow=nrow(xcell_tms), byrow=TRUE))
colnames(df) <- colnames(xcell_tms)
rownames(df) <- rownames(xcell_tms)
xcell_tms <- df

# create RData for signatures
mcp_counter_markers <- mcp_matrix
xcell_markers <- xcell_tms

signatures <- list(
  "quanTIseq" = quantiseq_signature,
  "EPIC" = epic_signature,
  "MCP-counter" = mcp_counter_markers,
  "xCell" = xcell_markers
)

save(signatures, file = "data/signatures/signatures.RData")
load("data/signatures/signatures.RData")

## merge mcp_counter and xcell
colnames(xcell_tms) <- paste(colnames(xcell_tms), "(xCell)")
xcell_tms$Genes <- row.names(xcell_tms)

colnames(mcp_matrix) <- paste(colnames(mcp_matrix), "(MCP-Counter)")
mcp_matrix$Genes <- row.names(mcp_matrix)

## merging
# Create a function to find matching strings
find_match <- function(gene) {
  matches <- xcell_tms$Genes[str_detect(xcell_tms$Genes, gene)]
  if(length(matches) > 0) return(matches[1]) else return(NA)
}

# Apply the function to table1 and merge
mcp_xcell_matrix <- mcp_matrix %>%
  rowwise() %>%
  mutate(Match = find_match(Genes)) 

# Merge the tables based on the matches
merged_table <- left_join(mcp_xcell_matrix, xcell_tms, by = c("Match" = "Genes"))

mcp_xcell_matrix <- as.data.frame(as.matrix(merged_table))
mcp_xcell_matrix <- mcp_xcell_matrix[!is.na(mcp_xcell_matrix$Match),]
row.names(mcp_xcell_matrix) <- mcp_xcell_matrix$Genes

mcp_xcell_matrix$Match <- NULL
mcp_xcell_matrix$Genes <- NULL

df_numeric <- mcp_xcell_matrix %>%
  mutate_all(as.numeric)

heatmap(as.matrix(df_numeric), na.rm = T, Rowv=F, Colv = F) +
  title(main = "Common transcriptomic markers for MCP-counter and xCell")
