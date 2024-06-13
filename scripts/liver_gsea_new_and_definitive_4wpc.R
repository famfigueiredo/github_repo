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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc'
)

results_files <-
  list.files(pattern = '^liver_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

gsea_simplified_results_dnavaccine_4wpc <-
  simplify(gsea_results_dnavaccine_4wpc)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_dnavaccine_4wpc)  # 27 GO terms/pathways

nrow(gsea_simplified_results_dnavaccine_4wpc)  # 9 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_4wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_4wpc@result) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_4wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 150)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    breaks = seq(0, max(low_high_nes_dnavaccine_4wpc$setSize), by = 50)
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 4WPC, liver, human orths') +
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
    verbose = F
  )

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% print(n = 100)

liver_dnavaccine_pathways <- as_tibble(y_dnavaccine) %>%
  arrange(-NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_dnavaccine_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_dnavaccine_pathways.tsv'
)

viewPathway('Interleukin-4 and Interleukin-13 signaling',
            readable = T,
            foldChange = entrez_gene_list)


# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_dnavaccine_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

markdown_table <- function(data) {
  # Get the header
  header <-
    paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <-
    paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
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
gsea_formatting(liver_res_eomes_vs_conu_4wpc, 'eomes', '4wpc')

gsea_simplified_results_eomes_4wpc <-
  simplify(gsea_results_eomes_4wpc)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_eomes_4wpc)  # 257 GO terms/pathways
nrow(gsea_simplified_results_eomes_4wpc)  # 118 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_eomes_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_eomes_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only 6 enriched terms

low_high_nes_eomes_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_eomes_4wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_eomes_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_eomes_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_eomes_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'EOMES, 4WPC, liver, human orths') +
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
    verbose = F
  )

as_tibble(y_eomes) %>% arrange(NES) %>% print(n = 100)

# viewPathway(
#   '',
#   readable = T,
#   foldChange = entrez_gene_list
# )

liver_eomes_pathways <- as_tibble(y_eomes) %>%
  arrange(NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>%
  dplyr::select(., Description, NES, setSize, Count)

write_tsv(
  liver_eomes_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_eomes_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_eomes_pathways.tsv',
    header = TRUE,
    sep = "\t"
  )

markdown_table <- function(data) {
  # Get the header
  header <-
    paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <-
    paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Print the Markdown table
cat(markdown_table(data), sep = "\n")

## GATA3 ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(liver_res_gata3_vs_conu_4wpc, 'gata3', '4wpc')

gsea_simplified_results_gata3_4wpc <-
  simplify(gsea_results_gata3_4wpc)  # simplifying GO terms to reduce redundancy

nrow(gsea_results_gata3_4wpc)  # 161 GO terms/pathways
nrow(gsea_simplified_results_gata3_4wpc)  # 72 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_gata3_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_gata3_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%  # only 3 downregulated terms
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_gata3_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)


low_high_nes_gata3_4wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_gata3_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_gata3_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_gata3_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'GATA3, 4WPC, liver, human orths') +
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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_gata3_pathways.tsv'
)

viewPathway(
  'Interleukin-4 and Interleukin-13 signaling',
  readable = T,
  foldChange = entrez_gene_list
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_gata3_pathways.tsv',
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

### Testing
wrangled_data <-
  improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

# Select ortholog_name and log2FC, and remove NA values
ortholog_fc <-
  wrangled_data %>%
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
gsea_formatting(liver_res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')

gsea_simplified_results_ivhd_4wpc <-
  simplify(gsea_results_ivhd_4wpc)  # simplify output from enrichGO and gseGO by removing redundancy of enriched GO terms

nrow(gsea_results_ivhd_4wpc)  # 514 GO terms/pathways

nrow(gsea_simplified_results_ivhd_4wpc)  # 190 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_ivhd_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_ivhd_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_4wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(low_high_nes_ivhd_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivhd_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivhd_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 4WPC, liver, human orths') +
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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_ivhd_pathways.tsv'
)

viewPathway('Interleukin-4 and Interleukin-13 signaling',
            readable = T,
            foldChange = entrez_gene_list)  # down

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_ivhd_pathways.tsv', header = TRUE, sep = "\t")

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
gsea_formatting(liver_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')

gsea_simplified_results_ivld_4wpc <-
  simplify(gsea_results_ivld_4wpc)

nrow(gsea_simplified_results_ivld_4wpc)  # 21 GO terms/pathways
nrow(gsea_results_ivld_4wpc)  # 11 GO terms/pathways

as_tibble(gsea_simplified_results_ivld_4wpc) %>% arrange(NES) %>% print(n = 100)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_ivld_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_ivld_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivld_4wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivld_4wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(low_high_nes_ivld_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivld_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivld_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 4WPC, liver, human orths') +
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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_ivld_pathways.tsv'
)

# Convert to a Markdown table ----
# Read the TSV file
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/pathways/liver_ivld_pathways.tsv',
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
