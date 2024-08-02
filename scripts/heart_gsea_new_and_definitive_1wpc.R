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

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_dnavaccine_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_dnavaccine_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_dnavaccine_1wpc.RData')

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
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only 6 terms had NES < 0

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

# Selecting variables for dotplot and sorting by gene count
gseGO_dnavaccine_1wpc <-
  low_high_nes_dnavaccine_1wpc %>% 
  dplyr::select(Description, setSize, Count, NES) %>% 
  arrange(desc(Count))

# Writing out top and bottom terms to .tsv
write_tsv(gseGO_dnavaccine_1wpc, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gseGO_dnavaccine_1wpc.tsv')

# Convert to a Markdown table ---
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gseGO_dnavaccine_1wpc.tsv', header = TRUE, sep = "\t")
# Print the Markdown table
cat(markdown_table(data), sep = "\n")

# gsePathway with 100k permutations in nPermSimple
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_dnavaccine_1wpc.RData')
y_dnavaccine <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      eps = 1e-300,
                      nPermSimple = 100000,
                      verbose = F)

as_tibble(y_dnavaccine) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)
viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = entrez_gene_list)
viewPathway('TNFR1-induced NF-kappa-B signaling pathway', readable = T, foldChange = entrez_gene_list)
viewPathway('Interleukin-6 family signaling', readable = T, foldChange = entrez_gene_list)


dnavaccine_1wpc_pathways <- as_tibble(y_dnavaccine) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine) %>% arrange(NES) %>% filter(., NES < 0) %>% pull(Description)

write_tsv(dnavaccine_1wpc_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/dnavaccine1wpc_gsePathways.tsv')

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
save(eomes_entrez_gene_list_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/heart_entrez_gene_list_eomes_1wpc.RData')

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_eomes_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_eomes_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_eomes_1wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_eomes_1wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_eomes_1wpc) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
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
      linewidth = .05,
      linetype = 2)) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)

genes_69597 <- str_split(low_high_nes_eomes_1wpc %>% filter(ID == 'GO:0006959') %>% pull(core_enrichment),'/')[[1]]
genes_19730 <- str_split(low_high_nes_eomes_1wpc %>% filter(ID == 'GO:0019730') %>% pull(core_enrichment),'/')[[1]]
genes_983695 <- str_split(as_tibble(y_eomes) %>% filter(ID == 'R-HSA-983695') %>% pull(core_enrichment), '/')[[1]]

y_eomes <- setReadable(y_eomes, org.Hs.eg.db)
intersect(genes_69597, genes_19730)
setdiff(genes_69597, genes_19730)

intersect(genes_69597, genes_983695)

# gsePathway with 100k permutations in nPermSimple
y_eomes <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      eps = 1e-300,
                      nPermSimple = 100000,
                      verbose = F)

as_tibble(y_eomes) %>% arrange(NES) %>% print(n = 100)

# viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

eomes_1wpc_pathways <- as_tibble(y_eomes) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

write_tsv(eomes_1wpc_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/eomes_1wpc_gsePathways.tsv')

viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = entrez_gene_list)
viewPathway('ADORA2B mediated anti-inflammatory cytokines production', readable = T, foldChange = entrez_gene_list)


# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/eomes_1wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

# ### significantly differentially regulated genes ###
# gsea_formatting_significant(res_eomes_vs_conu_1wpc, 'eomes', '1wpc')
# 
# gsea_eomes_1wpc_significant <- gseGO(
#   geneList = enrichment_gsea_eomes_1wpc,
#   OrgDb = org.Hs.eg.db,
#   ont = 'BP',
#   pvalueCutoff = 0.05,
#   pAdjustMethod = 'BH',
#   verbose = T)
# 
# as_tibble(gsea_eomes_1wpc_significant)  # no enriched terms

## GATA3 ----

# running gsea and saving files for later use
gsea_formatting(res_gata3_vs_conu_1wpc, 'gata3', '1wpc')
save(gsea_results_gata3_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_gata3_1wpc.RData')

gsea_simplified_results_gata3_1wpc <- simplify(gsea_results_gata3_1wpc)  # simplifying gsea results to remove redundancy of terms
save(gsea_simplified_results_gata3_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_gata3_1wpc.RData')
save(entrez_gene_list, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_gata31wpc.RData')

# loading saved files to re-run analysis
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_gata3_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_gata3_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_gata3_1wpc.RData')

# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_gata3_1wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_gata3_1wpc@result) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only three downregulated pathways

low_high_nes_gata3_1wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_gata3_1wpc %>%
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
      linewidth = .05,
      linetype = 2)) +
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)

# gsePathway with 100k permutations in nPermSimple
y_gata3 <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      eps = 1e-300,
                      nPermSimple = 100000,
                      verbose = F)

as_tibble(y_gata3) %>% arrange(NES) %>% print(n = 100)

gata31wpc_pathways <- as_tibble(y_gata3) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

write_tsv(gata31wpc_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gata31wpc_gsePathways.tsv')

# Convert to a Markdown table ----
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/gata31wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## IV-HD ----
# running gsea and saving files for later use
gsea_formatting(res_ivhd_vs_conu_1wpc, 'ivhd', '1wpc')
save(gsea_results_ivhd_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_ivhd_1wpc.RData')

gsea_simplified_results_ivhd_1wpc <- simplify(gsea_results_ivhd_1wpc)  # simplifying gsea results to remove redundancy of terms
save(gsea_simplified_results_ivhd_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_ivhd_1wpc.RData')
save(entrez_gene_list, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_ivhd_1wpc.RData')

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_ivhd_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_ivhd_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_ivhd_1wpc.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_ivhd_1wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_ivhd_1wpc@result) %>% 
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
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
    panel.grid = element_line(color = 'black', linewidth = .05, linetype = 2)
  ) + 
  guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)

y_ivhd <- gsePathway(
  entrez_gene_list,
  pvalueCutoff = .2,
  pAdjustMethod = 'BH',
  eps = 1e-300,
  nPermSimple = 100000,
  verbose = F
)

as_tibble(y_ivhd) %>% arrange(NES) %>% print(n = 100)

genes_909733 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-909733') %>% pull(core_enrichment), '/')[[1]]
genes_936964 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-936964') %>% pull(core_enrichment), '/')[[1]]
genes_918233 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-918233') %>% pull(core_enrichment), '/')[[1]]
# the genes pulled from core_enrichment do not match the ones in viewPathway. there are more genes in viewPathway.

ivhd_1wpc_pathways <- as_tibble(y_ivhd) %>% arrange(NES) %>% dplyr::select(., Description, NES)
write_tsv(ivhd_1wpc_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/ivhd_1wpc_gsePathways.tsv')

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)
viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = entrez_gene_list)
viewPathway('TRAF3-dependent IRF activation pathway', readable = T, foldChange = entrez_gene_list)


genes_2607 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0032607') %>% pull(core_enrichment),'/')[[1]]
genes_2647 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0032647') %>% pull(core_enrichment), '/')[[1]]
genes_61844 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0061844') %>% pull(core_enrichment), '/')[[1]]

# Convert to a Markdown table ----
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/ivhd_1wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## IV-LD ----
gsea_formatting(res_ivld_vs_conu_1wpc, 'ivld', '1wpc')
save(gsea_results_ivld_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_ivld_1wpc.RData')

gsea_simplified_results_ivld_1wpc <-
  simplify(gsea_results_ivld_1wpc)  # simplifying gsea results to remove redundancy of terms
save(gsea_simplified_results_ivld_1wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_ivld_1wpc.RData')
save(entrez_gene_list, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_ivld_1wpc.RData')

# loading saved files to re-run analysis
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_results_ivld_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/gsea_simplified_results_ivld_1wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/gsea_results_tables/entrez_gene_list_ivld_1wpc.RData')

# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_ivld_1wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_ivld_1wpc@result) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivld_1wpc <-
  bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivld_1wpc$Description

low_high_nes_ivld_1wpc %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 250)) +
  scale_size_continuous('Set size', 
                        range = c(2, 10), 
                        guide = 'legend',
                        limits = c(2, max(low_high_nes_ivld_1wpc$setSize))) +
  scale_x_continuous(limits = c(0, 250)) +
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
    panel.grid = element_line(
      color = 'black',
      linewidth = .05,
      linetype = 2
    )
  ) +
  guides(color = guide_legend(override.aes = list(size = 5)),
         # increase point size in gene count legend
         size = guide_legend(override.aes = list(
           shape = 1,
           fill = NA,
           stroke = .5,
           color = 'red'
         ))) +
  facet_grid(. ~ Regulation)


y_ivld <- gsePathway(
  entrez_gene_list,
  pvalueCutoff = .2,
  pAdjustMethod = 'BH',
  eps = 1e-300,
  nPermSimple = 100000,
  verbose = F
)

as_tibble(y_ivld) %>% arrange(NES) %>% print(n = 100)

ivld_1wpc_pathways <- as_tibble(y_ivld) %>% arrange(NES) %>% dplyr::select(., Description, NES)

write_tsv(ivld_1wpc_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/ivld_1wpc_gsePathways.tsv')

# Convert to a Markdown table ----
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/pathways/ivld_1wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")



