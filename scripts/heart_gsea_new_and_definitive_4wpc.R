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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')
save(gsea_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_dnavaccine_4wpc.RData')

gsea_simplified_results_dnavaccine_4wpc <-
  simplify(gsea_results_dnavaccine_4wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_dnavaccine_4wpc.RData')

dnavaccine_entrez_gene_list_4wpc <- entrez_gene_list
save(dnavaccine_entrez_gene_list_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_dnavaccine_4wpc.RData')

nrow(gsea_results_dnavaccine_4wpc)  # 1222 GO terms/pathways
nrow(gsea_simplified_results_dnavaccine_4wpc)  # 398 GO terms/pathways

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_dnavaccine_4wpc.RData')

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
  scale_color_viridis_c('Gene count', guide = 'legend') +
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
          subtitle = 'DNA vaccine, 4WPC, heart tissue') +
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


y_dnavaccine_4wpc <-
  gsePathway(
    entrez_gene_list,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('Signaling by CSF1 (M-CSF) in myeloid cells', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('GPVI-mediated activation cascade', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('Interleukin-3, Interleukin-5 and GM-CSF signaling', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('Inactivation of CSF3 (G-CSF) signaling', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)
viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = dnavaccine_entrez_gene_list_4wpc)


dnavaccine_4wpc_pathways <- as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% filter(., NES < 0) %>% pull(Description)

write_tsv(dnavaccine_4wpc_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/dnavaccine_1wpc_gsePathways.tsv')

# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/dnavaccine_1wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_eomes_vs_conu_4wpc, 'eomes', '4wpc')
save(gsea_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_eomes_4wpc.RData')

gsea_simplified_results_eomes_4wpc <-
  simplify(gsea_results_eomes_4wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_eomes_4wpc.RData')
save(eomes_entrez_gene_list_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_eomes_4wpc.RData')

nrow(gsea_results_eomes_4wpc)  # 1148 GO terms/pathways
nrow(gsea_simplified_results_eomes_4wpc)  # 367 GO terms/pathways

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_eomes1wpc.RData')

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
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

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
          subtitle = 'EOMES, 4WPC, heart tissue') +
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

eomes_entrez_gene_list_4wpc <- entrez_gene_list

y_eomes_4wpc <-
  gsePathway(
    entrez_gene_list,
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_eomes_4wpc) %>% arrange(NES) %>% print(n = 100)

eomes_4wpc_pathways <- as_tibble(y_eomes_4wpc) %>% arrange(NES) %>% dplyr::select(., Description, NES)

write_tsv(
  eomes_4wpc_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/eomes_4wpc_gsePathways.tsv'
)

# Convert to a Markdown table ----
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/eomes_4wpc_gsePathways.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(data), sep = "\n")

## GATA3 ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_gata3_vs_conu_4wpc, 'gata3', '4wpc')
save(gsea_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_gata3_4wpc.RData')

gsea_simplified_results_gata3_4wpc <-
  simplify(gsea_results_gata3_4wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_gata3_4wpc.RData')

gata3_entrez_gene_list_4wpc <- entrez_gene_list
save(gata3_entrez_gene_list_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_gata3_4wpc.RData')

nrow(gsea_results_gata3_4wpc)  #  796 GO terms/pathways
nrow(gsea_simplified_results_gata3_4wpc)  # 296 GO terms/pathways

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_gata3_4wpc.RData')

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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
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
          subtitle = 'GATA3, 4WPC, heart tissue') +
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


y_gata3_4wpc <-
  gsePathway(
    gata3_entrez_gene_list_4wpc,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_gata3_4wpc) %>% arrange(NES) %>% print(n = 100)
as_tibble(y_gata3_4wpc) %>% filter(NES < 0) %>% print(n = 100)

gata3_4wpc_pathways <- as_tibble(y_gata3_4wpc) %>% 
  arrange(NES) %>% 
  dplyr::select(., Description, NES) %>%
  mutate(NES = sprintf('%.3f', NES))  # format NES to 3 decimal places

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = gata3_entrez_gene_list_4wpc)
viewPathway('Activation of NF-kappaB in B cells', readable = T, foldChange = gata3_entrez_gene_list_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = gata3_entrez_gene_list_4wpc)

write_tsv(
  gata3_4wpc_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/gata3_4wpc_gsePathways.tsv'
)

# Convert to a Markdown table ----
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/gata3_4wpc_gsePathways.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(data), sep = "\n")

## IV-HD ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')
save(gsea_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivhd_4wpc.RData')

gsea_simplified_results_ivhd_4wpc <-
  simplify(gsea_results_ivhd_4wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivhd_4wpc.RData')

ivhd_entrez_gene_list_4wpc <- entrez_gene_list
save(ivhd_entrez_gene_list_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivhd_4wpc.RData')

nrow(gsea_results_ivhd_4wpc)  #   GO terms/pathways
nrow(gsea_simplified_results_ivhd_4wpc)  #  GO terms/pathways

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivhd_4wpc.RData')

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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivhd_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 4WPC, heart tissue') +
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
    
y_ivhd_4wpc <-
  gsePathway(
    ivhd_entrez_gene_list_4wpc,
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_ivhd_4wpc) %>% arrange(NES) %>% print(n = 100)

ivhd_4wpc_pathways <- as_tibble(y_ivhd_4wpc) %>% arrange(NES) %>% dplyr::select(., Description, NES)

ivhd_4wpc_pathways <- as_tibble(y_ivhd_4wpc) %>% 
  arrange(NES) %>% 
  dplyr::select(., Description, NES) %>%
  mutate(NES = sprintf('%.3f', NES))  # format NES to 3 decimal places


write_tsv(
  ivhd_4wpc_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivhd_4wpc_gsePathways.tsv'
)

options(ggrepel.max.overlaps = Inf)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('ROS and RNS production in phagocytes', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Signaling by CSF1 (M-CSF) in myeloid cells', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Cytokine Signaling in Immune system', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('TRAF3-dependent IRF activation pathway', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)


# Convert to a Markdown table ----
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivhd_4wpc_gsePathways.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(data), sep = "\n")


## IV-LD ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_ivld_vs_conu_4wpc, 'ivld', '4wpc')
save(gsea_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivld_4wpc.RData')

gsea_simplified_results_ivld_4wpc <-
  simplify(gsea_results_ivld_4wpc)  # simplifying GO terms to reduce redundancy
save(gsea_simplified_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivld_4wpc.RData')

ivld_entrez_gene_list_4wpc <- entrez_gene_list
save(ivld_entrez_gene_list_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivld_4wpc.RData')

nrow(gsea_results_ivld_4wpc)  # 1543 GO terms/pathways
nrow(gsea_simplified_results_ivld_4wpc)  # 418 GO terms/pathways

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivld_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(gsea_simplified_results_ivld_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_ivld_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivld_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, 300)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 4WPC, heart tissue') +
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

y_ivld_4wpc <-
  gsePathway(
    ivld_entrez_gene_list_4wpc,
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_ivld_4wpc) %>% arrange(NES) %>% print(n = 200)

ivld_4wpc_pathways <- as_tibble(y_ivld_4wpc) %>% arrange(NES) %>% dplyr::select(., Description, NES)

ivld_4wpc_pathways <- as_tibble(y_ivld_4wpc) %>% 
  arrange(NES) %>% 
  dplyr::select(., Description, NES) %>%
  mutate(NES = sprintf('%.3f', NES))  # format NES to 3 decimal places

write_tsv(
  ivld_4wpc_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_4wpc_gsePathways.tsv'
)

options(ggrepel.max.overlaps = Inf)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('ROS and RNS production in phagocytes', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Signaling by CSF1 (M-CSF) in myeloid cells', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('Cytokine Signaling in Immune system', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)
viewPathway('TRAF3-dependent IRF activation pathway', readable = T, foldChange = ivhd_entrez_gene_list_4wpc)


# Convert to a Markdown table ----
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_4wpc_gsePathways.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(data), sep = "\n")
        