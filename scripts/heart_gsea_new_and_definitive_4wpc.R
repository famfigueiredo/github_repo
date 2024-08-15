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
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_dnavaccine_vs_conu_4wpc.RData')

gsea_formatting(heart_res_dnavaccine_vs_conu_4wpc, 'heart', 'dnavaccine', '4wpc')
save(heart_gsea_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_dnavaccine_4wpc.RData')

heart_gsea_simplified_results_dnavaccine_4wpc <-
  simplify(heart_gsea_results_dnavaccine_4wpc)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_dnavaccine_4wpc.RData')

heart_entrez_gene_list_dnavaccine_4wpc <- entrez_gene_list
save(heart_entrez_gene_list_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_dnavaccine_4wpc.RData')

nrow(heart_gsea_results_dnavaccine_4wpc)  # 1423 GO terms/pathways
nrow(heart_gsea_simplified_results_dnavaccine_4wpc)  # 418 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_dnavaccine_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_dnavaccine_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_dnavaccine_4wpc) %>%
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
    heart_entrez_gene_list_dnavaccine_4wpc,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = heart_entrez_gene_list_dnavaccine_4wpc)
viewPathway('Signaling by CSF1 (M-CSF) in myeloid cells', readable = T, foldChange = heart_entrez_gene_list_dnavaccine_4wpc)
viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = heart_entrez_gene_list_dnavaccine_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = heart_entrez_gene_list_dnavaccine_4wpc)
viewPathway('FCERI mediated MAPK activation', readable = T, foldChange = heart_entrez_gene_list_dnavaccine_4wpc)



dnavaccine_4wpc_pathways <- as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine_4wpc) %>% arrange(NES) %>% filter(., NES < 0) %>% pull(Description)

write_tsv(dnavaccine_4wpc_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/dnavaccine_4wpc_gsePathways.tsv')

# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/dnavaccine_4wpc_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_eomes_vs_conu_4wpc.RData')

gsea_formatting(heart_res_eomes_vs_conu_4wpc, 'heart', 'eomes', '4wpc')
save(heart_gsea_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_eomes_4wpc.RData')

heart_gsea_simplified_results_eomes_4wpc <-
  simplify(heart_gsea_results_eomes_4wpc)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_eomes_4wpc.RData')

heart_entrez_gene_list_eomes_4wpc <- entrez_gene_list
save(heart_entrez_gene_list_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_eomes_4wpc.RData')

nrow(heart_gsea_results_eomes_4wpc)  # 1064 GO terms/pathways
nrow(heart_gsea_simplified_results_eomes_4wpc)  # 358 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_eomes1wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_eomes_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_eomes_4wpc) %>%
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
  # scale_color_viridis_c('Gene set')
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

y_eomes_4wpc <-
  gsePathway(
    heart_entrez_gene_list_eomes_4wpc,
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

# Convert to a Markdown table ---
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
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_gata3_vs_conu_4wpc.RData')

gsea_formatting(heart_res_gata3_vs_conu_4wpc, 'heart', 'gata3', '4wpc')
save(heart_gsea_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_gata3_4wpc.RData')

heart_gsea_simplified_results_gata3_4wpc <-
  simplify(heart_gsea_results_gata3_4wpc)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_gata3_4wpc.RData')

heart_entrez_gene_list_gata3_4wpc <- entrez_gene_list
save(heart_entrez_gene_list_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_gata3_4wpc.RData')

nrow(heart_gsea_results_gata3_4wpc)  # 1056 GO terms/pathways
nrow(heart_gsea_simplified_results_gata3_4wpc)  # 351 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_gata31wpc.RData')
# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
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
    heart_entrez_gene_list_gata3_4wpc,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(y_gata3_4wpc) %>% arrange(NES) %>% print(n = 100)

gata3_4wpc_pathways <- as_tibble(y_gata3_4wpc) %>% 
  arrange(NES) %>% 
  dplyr::select(., Description, NES) %>%
  mutate(NES = sprintf('%.3f', NES))  # format NES to 3 decimal places

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = heart_entrez_gene_list_gata3_4wpc)
viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = heart_entrez_gene_list_gata3_4wpc)
viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = heart_entrez_gene_list_gata3_4wpc)

write_tsv(
  gata3_4wpc_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/gata3_4wpc_gsePathways.tsv'
)

# Convert to a Markdown table ---
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
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_ivhd_vs_conu_4wpc.RData')

gsea_formatting(heart_res_ivhd_vs_conu_4wpc, 'heart', 'ivhd', '4wpc')
save(heart_gsea_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivhd_4wpc.RData')

heart_gsea_simplified_results_ivhd_4wpc <-
  simplify(heart_gsea_results_ivhd_4wpc)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivhd_4wpc.RData')

heart_entrez_gene_list_ivhd_4wpc <- entrez_gene_list
save(heart_entrez_gene_list_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivhd_4wpc.RData')

nrow(heart_gsea_results_ivhd_4wpc)  # 1274 GO terms/pathways
nrow(heart_gsea_simplified_results_ivhd_4wpc)  # 380 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivhd_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_ivhd_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_ivhd_4wpc) %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(low_high_nes_ivhd_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivhd_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, 1.1*max(low_high_nes_ivhd_4wpc$Count))) +
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
    
ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_definitive_plots/ivhd_4wpc.png', width = 1000, height = 1023, units = "px", dpi = 72)


y_ivhd_4wpc <-
  gsePathway(
    heart_entrez_gene_list_ivhd_4wpc,
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

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = heart_entrez_gene_list_ivhd_4wpc)
ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivhd_interferon_alpha-beta_signaling.png', width = 1600, height = 1000, units = "px", dpi = 72)

viewPathway('ROS and RNS production in phagocytes', readable = T, foldChange = heart_entrez_gene_list_ivhd_4wpc)
viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = heart_entrez_gene_list_ivhd_4wpc)

# Convert to a Markdown table ---
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
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_ivld_vs_conu_4wpc.RData')

gsea_formatting(heart_res_ivld_vs_conu_4wpc, 'heart', 'ivld', '4wpc')
save(heart_gsea_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivld_4wpc.RData')

heart_gsea_simplified_results_ivld_4wpc <-
  simplify(heart_gsea_results_ivld_4wpc)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivld_4wpc.RData')

heart_entrez_gene_list_ivld_4wpc <- entrez_gene_list
save(heart_entrez_gene_list_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivld_4wpc.RData')

nrow(heart_gsea_results_ivld_4wpc)  # 1521 GO terms/pathways
nrow(heart_gsea_simplified_results_ivld_4wpc)  # 409 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_gsea_simplified_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_results_tables/heart_entrez_gene_list_ivld_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_ivld_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_ivld_4wpc) %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(low_high_nes_ivld_4wpc$Count))) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_ivld_4wpc$setSize))
  ) +
  scale_x_continuous(limits = c(0, 350)) +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/gsea_definitive_plots/ivld_4wpc.png', width = 1000, height = 1023, units = "px", dpi = 72)


y_ivld_4wpc <-
  gsePathway(
    heart_entrez_gene_list_ivld_4wpc,
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

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = heart_entrez_gene_list_ivld_4wpc)
ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_interferon_alpha-beta_signaling.png', width = 1600, height = 1000, units = "px", dpi = 72)

viewPathway('Antigen processing-Cross presentation', readable = T, foldChange = heart_entrez_gene_list_ivld_4wpc)
ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_antigen_processing_cross-presentation.png', width = 1600, height = 1000, units = "px", dpi = 72)

viewPathway('Antigen activates B Cell Receptor (BCR) leading to generation of second messengers', readable = T, foldChange = heart_entrez_gene_list_ivld_4wpc)
ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_antigen_BCR_second_messengers.png', width = 1600, height = 1000, units = "px", dpi = 72)


# Convert to a Markdown table ---
data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/pathways/ivld_4wpc_gsePathways.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(data), sep = "\n")

improved_data_wrangling(heart_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')
#############################
# Experimenting with shrunk data ----























