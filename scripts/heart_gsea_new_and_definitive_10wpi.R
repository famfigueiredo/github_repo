library('tidyverse')
library('VennDiagram')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('enrichplot')
library('AnnotationDbi')
library('ReactomePA')

rm(list = setdiff(ls(), grep("res_", ls(), value = TRUE)))

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi'
)

results_files <-
  list.files(pattern = '^res_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')

gsea_simplified_results_dnavaccine_10wpi <- simplify(gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_dnavaccine_10wpi)  # 658 GO terms/pathways

nrow(gsea_simplified_results_dnavaccine_10wpi)  # 220 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_10wpi@result) %>%
  arrange(desc(NES)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_10wpi@result) %>% 
  arrange(NES) %>% 
  top_n(8, wt = desc(NES)) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 200)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_dnavaccine_10wpi$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 10WPI, heart tissue') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(color = 'black', size = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)

y_dnavaccine <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% print(n = 100)

# viewPathway('', readable = T, foldChange = entrez_gene_list)

dnavaccine_pathways <- 
  as_tibble(y_dnavaccine) %>% 
  arrange(NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(dnavaccine_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/dnavaccine_upregulated_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/dnavaccine_upregulated_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")



## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_eomes_vs_conu_10wpi, 'eomes', '10wpi')

gsea_simplified_results_eomes_10wpi <- simplify(gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_eomes_10wpi)  # 3035 GO terms/pathways
nrow(gsea_simplified_results_eomes_10wpi)  # 649 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_eomes_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <- 
#   as_tibble(gsea_simplified_results_eomes_10wpi) %>%
#   filter(NES < 0) %>% 
#   arrange(desc(setSize)) %>% 
#   top_n(10, wt = setSize) %>% # only 5 downregulated terms
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top10_high_nes %>%
  # mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(top10_high_nes$setSize))) +
  scale_x_continuous(limits = c(0, max(top10_high_nes$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
          subtitle = 'EOMES, 10WPI, heart, human orths') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(color = 'black', linewidth = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  )

y_eomes <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(y_eomes) %>% arrange(-NES) %>% print(n = 100)

viewPathway('', readable = T, foldChange = entrez_gene_list)

eomes_pathways <- as_tibble(y_eomes) %>% 
  arrange(-NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(eomes_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/eomes_upregulated_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/eomes_upregulated_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")

### Testing
wrangled_data <-
  improved_data_wrangling(res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')

# Select ortholog_name and log2FC, and remove NA values
ortholog_fc <- wrangled_data %>%
  dplyr::select(ortholog_name, log2FC) %>%
  na.omit()

# Convert gene symbols to Entrez IDs
entrez_ids <-
  bitr(
    ortholog_fc$ortholog_name,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = org.Hs.eg.db
  )


entrez_genes <-
  wrangled_data %>% left_join(entrez_ids,
                              by = c('ortholog_name' = 'SYMBOL'),
                              relationship = 'many-to-many') %>% dplyr::select(ENTREZID, log2FC) %>% na.omit()




## Testing different annotation files

## GATA3 ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_gata3_vs_conu_10wpi, 'gata3', '10wpi')

gsea_simplified_results_gata3_10wpi <- simplify(gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_gata3_10wpi)  # 3122 GO terms/pathways
nrow(gsea_simplified_results_gata3_10wpi)  # 668 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_gata3_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

top10_high_nes %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', breaks = seq(0, 300, by = 50)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', breaks = seq(0, max(top10_high_nes$setSize), by = 50)) +
  scale_x_continuous(limits = c(0, max(top10_high_nes$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
          subtitle = 'GATA3, 10WPI, heart, human orths') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(color = 'black', linewidth = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  )

y_gata3 <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(y_gata3) %>% arrange(-NES) %>% print(n = 100)

viewPathway('', readable = T, foldChange = entrez_gene_list)

gata3_pathways <- as_tibble(y_gata3) %>% 
  arrange(-NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(gata3_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/gata3_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/gata3_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")

### Testing
wrangled_data <-
  improved_data_wrangling(res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')

# Select ortholog_name and log2FC, and remove NA values
ortholog_fc <- wrangled_data %>%
  dplyr::select(ortholog_name, log2FC) %>%
  na.omit()

# Convert gene symbols to Entrez IDs
entrez_ids <-
  bitr(
    ortholog_fc$ortholog_name,
    fromType = 'SYMBOL',
    toType = 'ENTREZID',
    OrgDb = org.Hs.eg.db
  )


entrez_genes <-
  wrangled_data %>% left_join(entrez_ids,
                              by = c('ortholog_name' = 'SYMBOL'),
                              relationship = 'many-to-many') %>% dplyr::select(ENTREZID, log2FC) %>% na.omit()






## IV-HD ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_ivhd_vs_conu_10wpi, 'ivhd', '10wpi')

nrow(gsea_results_ivhd_10wpi)  # 1255 GO terms/pathways

gsea_simplified_results_ivhd_10wpi <- simplify(gsea_results_ivhd_10wpi)  # simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms

nrow(gsea_simplified_results_ivhd_10wpi)  # 345 GO terms/pathways

as_tibble(gsea_simplified_results_ivhd_10wpi) %>% arrange(NES) %>% print(n = 100)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES < 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)
max(low_high_nes_ivhd_10wpi$Count)

low_high_nes_ivhd_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, 200)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_ivhd_10wpi$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivhd_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 10WPI, heart, human orths') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(color = 'black', linewidth = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_ivhd <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                     pvalueCutoff = .2,
                     pAdjustMethod = 'BH',
                     verbose = F)

as_tibble(y_ivhd) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

ivhd_pathways <- as_tibble(y_ivhd) %>% 
  arrange(-NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(ivhd_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/ivhd_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/ivhd_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")


















## IV-LD ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

nrow(gsea_results_ivld_10wpi)  # 1694 GO terms/pathways

gsea_simplified_results_ivld_10wpi <- simplify(gsea_results_ivld_10wpi)  # simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms

nrow(gsea_simplified_results_ivld_10wpi)  # 433 GO terms/pathways

as_tibble(gsea_simplified_results_ivld_10wpi) %>% arrange(NES) %>% print(n = 100)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_ivld_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_ivld_10wpi) %>%
  filter(NES < 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivld_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivld_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, 400)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_ivld_10wpi$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivld_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 10WPI, heart, human orths') +
  theme_bw(base_size = 14) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm'),
    strip.text = element_text(size = 24),
    plot.title = element_text(hjust = .5),
    plot.subtitle = element_text(hjust = .5),
    panel.grid.minor = element_blank(),
    panel.grid = element_line(color = 'black', linewidth = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_ivld <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                     pvalueCutoff = .2,
                     pAdjustMethod = 'BH',
                     verbose = T)

as_tibble(y_ivld) %>% arrange(NES) %>% print(n = 100)

viewPathway('', readable = T, foldChange = entrez_gene_list)

ivld_pathways <- as_tibble(y_ivld) %>% 
  arrange(NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(ivld_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/ivld_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi/pathways/ivld_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <- paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <- paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")


















