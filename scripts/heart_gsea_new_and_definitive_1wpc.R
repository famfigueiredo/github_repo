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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_1wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# GSEA ----

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
## added saving feature on 24/06/2024, so I don't have to re-run the formatting function everytime
gsea_formatting(res_dnavaccine_vs_conu_1wpc, 'dnavaccine', '1wpc')
save(gsea_results_dnavaccine_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_dnavaccine_1wpc.RData')

gsea_simplified_results_dnavaccine_1wpc <- simplify(gsea_results_dnavaccine_1wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_dnavaccine_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_dnavaccine_1wpc.RData')
save(entrez_gene_list, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_dnavaccine1wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_1wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_1wpc) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(8, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_1wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_1wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'black') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend') +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend') +
  scale_x_continuous() +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 1WPC, heart tissue') +
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
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)

gseGO_dnavaccine_1wpc <- low_high_nes_dnavaccine_1wpc %>% dplyr::select(Description, setSize, Count, NES) %>% arrange(desc(Count))

write_tsv(gseGO_dnavaccine_1wpc, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gseGO_dnavaccine_1wpc.tsv')

# Convert to a Markdown table ---
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gseGO_dnavaccine_1wpc.tsv', header = TRUE, sep = "\t")
# Print the Markdown table
cat(markdown_table(data), sep = "\n")

# gsePathway
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_dnavaccine1wpc.RData')
y_dnavaccine <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      eps = 1e-300,
                      nPermSimple = 100000,
                      verbose = F)

as_tibble(y_dnavaccine) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

dnavaccine1wpc_pathways <- as_tibble(y_dnavaccine) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine) %>% arrange(NES) %>% filter(., NES < 0) %>% pull(Description)

write_tsv(dnavaccine1wpc_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/dnavaccine1wpc_gsePathways.tsv')

# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/dnavaccine1wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

### significantly differentially regulated genes ###
# gsea formatting starting from a DESeq results table and using only significantly differentially regulated genes
gsea_formatting_significant(res_dnavaccine_vs_conu_1wpc, 'dnavaccine', '1wpc')

gsea_dnavaccine_1wpc_significant <- gseGO(
  geneList = enrichment_gsea_dnavaccine_1wpc,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_dnavaccine_1wpc_significant)  # no enriched terms

## EOMES ----
### all genes ###
gsea_formatting(res_eomes_vs_conu_1wpc, 'eomes', '1wpc')
save(gsea_results_eomes_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_eomes_1wpc.RData')

gsea_simplified_results_eomes_1wpc <- simplify(gsea_results_eomes_1wpc)
save(gsea_simplified_results_eomes_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_eomes_1wpc.RData')
save(entrez_gene_list, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_eomes1wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_eomes_1wpc@result) %>%
  arrange(desc(NES)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_eomes_1wpc@result) %>% 
  arrange(NES) %>% 
  top_n(10, wt = desc(NES)) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_1wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_eomes_1wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend') +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'EOMES, 1WPC, heart tissue') +
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
      size = .05,
      linetype = 2)) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_eomes <- gsePathway(entrez_gene_list,
                pvalueCutoff = .2,
                pAdjustMethod = 'BH',
                verbose = F)

as_tibble(y_eomes) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

eomes_pathways_down <- as_tibble(y_eomes) %>% arrange(NES) %>% filter(., NES < 0) %>% dplyr::select(., Description, NES) 

write_tsv(eomes_pathways_down, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/eomes_downregulated_pathways.tsv')

# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/eomes_downregulated_pathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

### significantly differentially regulated genes ###
gsea_formatting_significant(res_eomes_vs_conu_1wpc, 'eomes', '1wpc')

gsea_eomes_1wpc_significant <- gseGO(
  geneList = enrichment_gsea_eomes_1wpc,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_eomes_1wpc_significant)  # no enriched terms

## GATA3 ----
rm(list = setdiff(ls(), grep("res_", ls(), value = TRUE)))

gsea_formatting(res_gata3_vs_conu_1wpc, 'gata3', '1wpc')

gsea_simplified_results_gata3_1wpc <- simplify(gsea_results_gata3_1wpc)  # simplifying gsea results to remove redundancy of terms

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top20_high_nes <- 
  as_tibble(gsea_results_gata3_1wpc@result) %>%
  arrange(desc(-NES)) %>% 
  top_n(20, wt = -p.adjust) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

top20_high_nes %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend') +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, upregulated genes',
          subtitle = 'GATA3, 1WPC, heart tissue') +
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
      size = .05,
      linetype = 2)) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  )

y_gata3 <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(y_gata3) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

gata3_pathways_down <- as_tibble(y_gata3) %>% arrange(NES) %>% filter(., NES < 0) %>% dplyr::select(., Description, NES) 

write_tsv(gata3_pathways_down, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/gata3_downregulated_pathways.tsv')


# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/gata3_downregulated_pathways.tsv', header = TRUE, sep = "\t")

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

### significantly differentially regulated genes ###
gsea_formatting_significant(res_gata3_vs_conu_1wpc, 'gata3', '1wpc')

gsea_eomes_1wpc_significant <- gseGO(
  geneList = enrichment_gsea_eomes_1wpc,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_eomes_1wpc_significant)  # no enriched terms. very few significantly regulated genes

## IV-HD ----
gsea_formatting(res_ivhd_vs_conu_1wpc, 'ivhd', '1wpc')

gsea_simplified_results_ivhd_1wpc <- simplify(gsea_results_ivhd_1wpc)  # simplifying gsea results to remove redundancy of terms

as_tibble(gsea_results_ivhd_1wpc)
as_tibble(gsea_simplified_results_ivhd_1wpc)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_ivhd_1wpc@result) %>%
  arrange(desc(NES)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_ivhd_1wpc@result) %>% 
  arrange(NES) %>% 
  top_n(8, wt = desc(NES)) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_1wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_1wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend') +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend') +
  scale_x_continuous() +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 1WPC, heart tissue') +
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
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_ivhd <- gsePathway(entrez_gene_list,
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_ivhd) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

ivhd_pathways_down <- as_tibble(y_ivhd) %>% arrange(NES) %>% filter(., NES < 0) %>% dplyr::select(., Description, NES)
ivhd_pathways_up <- as_tibble(y_ivhd) %>% arrange(NES) %>% filter(., NES > 0) %>% dplyr::select(., Description, NES) 

ivhd_pathways <- bind_rows(ivhd_pathways_down, ivhd_pathways_up)

write_tsv(ivhd_pathways_down, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/ivhd_downregulated_pathways.tsv')
write_tsv(ivhd_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/ivhd_pathways.tsv')

# Convert to a Markdown table ---
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/ivhd_pathways.tsv', header = TRUE, sep = "\t")

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

### significantly differentially regulated genes ###
gsea_formatting_significant(res_ivhd_vs_conu_1wpc, 'ivhd', '1wpc')

gsea_ivhd_1wpc_significant <- gseGO(
  geneList = enrichment_gsea_ivhd_1wpc,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_eomes_1wpc_significant)  # no enriched terms.

## IV-LD ----
gsea_formatting(res_ivld_vs_conu_1wpc, 'ivld', '1wpc')

gsea_simplified_results_ivld_1wpc <- simplify(gsea_results_ivld_1wpc)  # simplifying gsea results to remove redundancy of terms

as_tibble(gsea_results_ivld_1wpc)
as_tibble(gsea_simplified_results_ivld_1wpc) %>% arrange(NES) %>% top_n(20, wt = NES)


# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top20_high_nes <- 
  as_tibble(gsea_simplified_results_ivld_1wpc@result) %>%
  arrange(desc(NES)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

top20_high_nes %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(47, 260)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend') +
  scale_x_continuous(limits = c(0, 260), breaks = seq(0, 260, by = 50)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 1WPC, heart tissue') +
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
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  )

y_ivld <- gsePathway(entrez_gene_list,
                     pvalueCutoff = .2,
                     pAdjustMethod = 'BH',
                     verbose = F)

as_tibble(y_ivld) %>% arrange(NES) %>% print(n = 100)
# the only pathway enriched on Reactome was 'Neuronal System'


### significantly differentially regulated genes ###
gsea_formatting_significant(res_ivld_vs_conu_1wpc, 'ivld', '1wpc')

gsea_ivld_1wpc_significant <- gseGO(
  geneList = enrichment_gsea_ivld_1wpc,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_ivld_1wpc_significant)  # no enriched terms. Only one gene enriched in gsea.

