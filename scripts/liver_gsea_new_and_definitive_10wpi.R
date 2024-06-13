library('tidyverse')
library('VennDiagram')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('enrichplot')
library('AnnotationDbi')
library('ReactomePA')

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi'
)

load_results('^liver_.*_conu_10wpi')  # regex matching results files

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')

gsea_simplified_results_dnavaccine_10wpi <-
  simplify(gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_dnavaccine_10wpi)  # 653 GO terms/pathways

nrow(gsea_simplified_results_dnavaccine_10wpi)  # 245 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_10wpi@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_10wpi@result) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_10wpi <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 200)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    breaks = seq(0, max(low_high_nes_dnavaccine_10wpi$setSize), by = 50)
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 10WPI, liver, human orths') +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(color = guide_legend(override.aes = list(size = 5)),
             # increase point size in gene count legend
             size = guide_legend(override.aes = list(
               shape = 1,
               fill = NA,
               stroke = .5,
               color = 'red'
             ))) +
  facet_grid(. ~ Regulation)


y_dnavaccine <-
  gsePathway(
    entrez_gene_list,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_dnavaccine) %>% arrange(NES) %>% print(n = 100)

viewPathway(
  'Interleukin-12 family signaling',
  readable = T,
  foldChange = entrez_gene_list
)  # up

liver_dnavaccine_pathways <- as_tibble(y_dnavaccine) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_dnavaccine_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_dnavaccine_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_dnavaccine_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

# Print the Markdown table
cat(markdown_table(data), sep = "\n")



## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_eomes_vs_conu_10wpi, 'eomes', '10wpi')

gsea_simplified_results_eomes_10wpi <-
  simplify(gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_eomes_10wpi)  # 2379 GO terms/pathways
nrow(gsea_simplified_results_eomes_10wpi)  # 551 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_eomes_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_eomes_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_10wpi <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_eomes_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 400)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_eomes_10wpi$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_eomes_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'EOMES, 10WPI, liver, human orths') +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(color = guide_legend(override.aes = list(size = 5)),
             # increase point size in gene count legend
             size = guide_legend(override.aes = list(
               shape = 1,
               fill = NA,
               stroke = .5,
               color = 'red'
             ))) +
  facet_grid(. ~ Regulation)


y_eomes <-
  gsePathway(
    entrez_gene_list,
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_eomes) %>% arrange(NES) %>% print(n = 100)

viewPathway(
  'Interleukin-4 and Interleukin-13 signaling',
  readable = T,
  foldChange = entrez_gene_list
)  # up

liver_eomes_pathways <- as_tibble(y_eomes) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_eomes_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_eomes_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_eomes_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

# Print the Markdown table
cat(markdown_table(data), sep = "\n")




## GATA3 ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_gata3_vs_conu_10wpi, 'gata3', '10wpi')

gsea_simplified_results_gata3_10wpi <-
  simplify(gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_gata3_10wpi)  # 2940 GO terms/pathways
nrow(gsea_simplified_results_gata3_10wpi)  # 627 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_gata3_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_gata3_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%  # only 3 downregulated terms
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_gata3_10wpi <-
  bind_rows(top10_high_nes, bottom10_low_nes)


low_high_nes_gata3_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_gata3_10wpi$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_gata3_10wpi$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_gata3_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'GATA3, 10WPI, liver, human orths') +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(color = guide_legend(override.aes = list(size = 5)),
             # increase point size in gene count legend
             size = guide_legend(override.aes = list(
               shape = 1,
               fill = NA,
               stroke = .5,
               color = 'red'
             ))) +
  facet_grid(. ~ Regulation)


y_gata3 <-
  gsePathway(
    entrez_gene_list,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 1000000,
    verbose = F
  )

as_tibble(y_gata3) %>% arrange(NES) %>% print(n = 100)

liver_gata3_pathways <-
  as_tibble(y_gata3) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_gata3_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_gata3_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_gata3_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

# Print the Markdown table
cat(markdown_table(data), sep = "\n")



## IV-HD ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_ivhd_vs_conu_10wpi, 'ivhd', '10wpi')

gsea_simplified_results_ivhd_10wpi <-
  simplify(gsea_results_ivhd_10wpi)  # simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms

nrow(gsea_results_ivhd_10wpi)  # 514 GO terms/pathways

nrow(gsea_simplified_results_ivhd_10wpi)  # 190 GO terms/pathways

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

low_high_nes_ivhd_10wpi <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(low_high_nes_ivhd_10wpi$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivhd_10wpi$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivhd_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 10WPI, liver, human orths') +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),
    # increase point size in gene count legend
    size = guide_legend(override.aes = list(
      shape = 1,
      fill = NA,
      stroke = .5,
      color = 'red'
    ))) +
  facet_grid(. ~ Regulation)

y_ivhd <-
  gsePathway(
    entrez_gene_list,
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    verbose = F
  )

as_tibble(y_ivhd) %>% arrange(NES) %>% print(n = 100)

liver_ivhd_pathways <-
  as_tibble(y_ivhd) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_ivhd_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_ivhd_pathways.tsv'
)

viewPathway('Interleukin-4 and Interleukin-13 signaling',
            readable = T,
            foldChange = entrez_gene_list)  # down

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_ivhd_pathways.tsv', header = TRUE, sep = "\t")

markdown_table <- function(data) {
  # Get the header
  header <-
    paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <-
    paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <-
    apply(data, 1, function(row) {
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
gsea_formatting(liver_res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

gsea_simplified_results_ivld_10wpi <-
  simplify(gsea_results_ivld_10wpi)

nrow(gsea_simplified_results_ivld_10wpi)  # 21 GO terms/pathways
nrow(gsea_results_ivld_10wpi)  # 11 GO terms/pathways

as_tibble(gsea_simplified_results_ivld_10wpi) %>% arrange(NES) %>% print(n = 100)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_ivld_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_ivld_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivld_10wpi <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivld_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(
    aes(size = setSize),
    shape = 1,
    stroke = 0.2,
    color = 'red'
  ) +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(low_high_nes_ivld_10wpi$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivld_10wpi$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivld_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 10WPI, liver, human orths') +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),
    # increase point size in gene count legend
    size = guide_legend(override.aes = list(
      shape = 1,
      fill = NA,
      stroke = .5,
      color = 'red'
    ))) +
  facet_grid(. ~ Regulation
  )


y_ivld <-
  gsePathway(
    entrez_gene_list,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    verbose = T
  )

as_tibble(y_ivld) %>% arrange(NES) %>% print(n = 100)

viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers',
            readable = T,
            foldChange = entrez_gene_list)  # down

liver_ivld_pathways <-
  as_tibble(y_ivld) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_ivld_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_ivld_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/liver_ivld_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

markdown_table <-
  function(data) {
    # Get the header
    header <-
      paste("|", paste(names(data), collapse = " | "), "|")
    
    # Get the separator line
    separator <-
      paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
    
    # Get the table rows
    rows <-
      apply(data, 1, function(row) {
        paste("|", paste(row, collapse = " | "), "|")
      })
    
    # Combine header, separator, and rows
    c(header, separator, rows)
  }

# Print the Markdown table
cat(markdown_table(data), sep = "\n")
