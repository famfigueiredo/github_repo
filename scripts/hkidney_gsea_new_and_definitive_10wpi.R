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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# GSEA ----
## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
## added saving feature on 24/06/2024, so I don't have to re-run the formatting function everytime
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_dnavaccine_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_dnavaccine_vs_conu_10wpi, 'hkidney', 'dnavaccine', '10wpi')
save(hkidney_gsea_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_dnavaccine_10wpi.RData')

hkidney_gsea_simplified_results_dnavaccine_10wpi <-
  simplify(hkidney_gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_10wpi.RData')

hkidney_entrez_gene_list_dnavaccine_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_10wpi.RData')

nrow(hkidney_gsea_results_dnavaccine_10wpi)  # 1 GO terms/pathways
nrow(hkidney_gsea_simplified_results_dnavaccine_10wpi)  # 1 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <-
  as_tibble(hkidney_gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only 6 terms had NES < 0

low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(low_high_nes_dnavaccine_10wpi$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_dnavaccine_10wpi$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_dnavaccine_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 1WPC, hkidney tissue') +
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
             size = guide_legend(override.aes = list(
               shape = 1,
               fill = NA,
               stroke = .5,
               color = 'red'
             ))) +
  facet_grid(. ~ Regulation)

ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/dnavaccine_10wpi.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)


# Selecting variables for dotplot and sorting by gene count
gseGO_dnavaccine_10wpi <-
  low_high_nes_dnavaccine_10wpi %>% 
  dplyr::select(Description, setSize, Count, NES) %>% 
  arrange(desc(Count))

# Writing out top and bottom terms to .tsv
write_tsv(gseGO_dnavaccine_10wpi, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/gseGO_dnavaccine_10wpi.tsv')

y_dnavaccine_10wpi <- gsePathway(hkidney_entrez_gene_list_dnavaccine_10wpi,
                                pvalueCutoff = .2,
                                pAdjustMethod = 'BH',
                                eps = 1e-300,
                                nPermSimple = 10000,
                                verbose = F)

as_tibble(y_dnavaccine_10wpi) %>% arrange(NES) %>% print(n = 100)

options(ggrepel.max.overlaps = Inf)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = hkidney_entrez_gene_list_dnavaccine_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/dnavaccine_interferon_alpha-beta_signaling.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')

viewPathway('TNFR1-induced NF-kappa-B signaling pathway', readable = T, foldChange = hkidney_entrez_gene_list_dnavaccine_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/dnavaccine_tnfr1-induced_nfkb_signaling.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')

viewPathway('Class I MHC mediated antigen processing & presentation', readable = T, foldChange = hkidney_entrez_gene_list_dnavaccine_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/dnavaccine_classI_mhc_mediated_antigen_p.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')


dnavaccine_10wpi_pathways <- as_tibble(y_dnavaccine_10wpi) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

write_tsv(dnavaccine_10wpi_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/dnavaccine10wpi_gsePathways.tsv')

# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/dnavaccine10wpi_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

### significantly differentially regulated genes ###
# gsea formatting starting from a DESeq results table and using only significantly differentially regulated genes
gsea_formatting_significant(res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')

gsea_dnavaccine_10wpi_significant <- gseGO(
  geneList = enrichment_gsea_dnavaccine_10wpi,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

as_tibble(gsea_dnavaccine_10wpi_significant)  # no enriched terms

## EOMES ----
### all genes ###
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_eomes_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_eomes_vs_conu_10wpi, 'hkidney', 'eomes', '10wpi')
save(hkidney_gsea_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_eomes_10wpi.RData')

hkidney_gsea_simplified_results_eomes_10wpi <-
  simplify(hkidney_gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_eomes_10wpi.RData')

hkidney_entrez_gene_list_eomes_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_eomes_10wpi.RData')

nrow(hkidney_gsea_results_eomes_10wpi)  # 711 GO terms/pathways
nrow(hkidney_gsea_simplified_results_eomes_10wpi)  # 256 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_eomes_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_eomes_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_eomes_10wpi) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_eomes_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(low_high_nes_eomes_10wpi$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_eomes_10wpi$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_eomes_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'EOMES, 1WPC, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/eomes_10wpi.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

# genes_69597 <- str_split(low_high_nes_eomes_10wpi %>% filter(ID == 'GO:0006959') %>% pull(core_enrichment),'/')[[1]]
# genes_19730 <- str_split(low_high_nes_eomes_10wpi %>% filter(ID == 'GO:0019730') %>% pull(core_enrichment),'/')[[1]]
# genes_983695 <- str_split(as_tibble(y_eomes) %>% filter(ID == 'R-HSA-983695') %>% pull(core_enrichment), '/')[[1]]

# y_eomes <- setReadable(y_eomes, org.Hs.eg.db)
# intersect(genes_69597, genes_19730)
# setdiff(genes_69597, genes_19730)
# 
# intersect(genes_69597, genes_983695)

# gsePathway with 100k permutations in nPermSimple
y_eomes_10wpi <- gsePathway(hkidney_entrez_gene_list_eomes_10wpi,
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           eps = 1e-300,
                           nPermSimple = 10000,
                           verbose = F)

as_tibble(y_eomes) %>% arrange(NES) %>% print(n = 100)

eomes_10wpi_pathways <- as_tibble(y_eomes_10wpi) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

write_tsv(eomes_10wpi_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/eomes_10wpi_gsePathways.tsv')

viewPathway('Viral mRNA Translation', readable = T, foldChange = hkidney_entrez_gene_list_eomes_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/eomes_viral_mRNA_translation.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')


# Convert to a Markdown table ---
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/eomes_10wpi_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## GATA3 ----

# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_gata3_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_gata3_vs_conu_10wpi, 'hkidney', 'gata3', '10wpi')
save(hkidney_gsea_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_gata3_10wpi.RData')

hkidney_gsea_simplified_results_gata3_10wpi <-
  simplify(hkidney_gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_gata3_10wpi.RData')

hkidney_entrez_gene_list_gata3_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_gata3_10wpi.RData')

nrow(hkidney_gsea_results_gata3_10wpi)  # 13 GO terms/pathways
nrow(hkidney_gsea_simplified_results_gata3_10wpi)  # 6 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_gata3_10wpi.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_gata3_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_gata3_10wpi) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only three downregulated pathways

low_high_nes_gata3_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_gata3_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(low_high_nes_gata3_10wpi$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_gata3_10wpi$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_gata3_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'GATA3, 1WPC, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/gata3_10wpi.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

# gsePathway with 100k permutations in nPermSimple
y_gata3_10wpi <- gsePathway(hkidney_entrez_gene_list_gata3_10wpi,
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           eps = 1e-300,
                           nPermSimple = 10000,
                           verbose = F)

as_tibble(y_gata3_10wpi) %>% arrange(NES) %>% print(n = 100)

gata3_10wpi_pathways <- as_tibble(y_gata3_10wpi) %>% arrange(NES) %>% dplyr::select(., Description, NES) 

write_tsv(gata3_10wpi_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/gata3_10wpi_gsePathways.tsv')

# Convert to a Markdown table ----
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/gata3_10wpi_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## IV-HD ----
# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_ivhd_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_ivhd_vs_conu_10wpi, 'hkidney', 'ivhd', '10wpi')
save(hkidney_gsea_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivhd_10wpi.RData')

hkidney_gsea_simplified_results_ivhd_10wpi <-
  simplify(hkidney_gsea_results_ivhd_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_10wpi.RData')

hkidney_entrez_gene_list_ivhd_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivhd_10wpi.RData')

nrow(hkidney_gsea_results_ivhd_10wpi)  # 483 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivhd_10wpi)  # 179 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivhd_10wpi.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_10wpi) %>% 
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(low_high_nes_ivhd_10wpi$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_ivhd_10wpi$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_ivhd_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 1WPC, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/ivhd_10wpi.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

y_ivhd_10wpi <- gsePathway(
  hkidney_entrez_gene_list_ivhd_10wpi,
  pvalueCutoff = .2,
  pAdjustMethod = 'BH',
  eps = 1e-300,
  nPermSimple = 10000,
  verbose = F
)

as_tibble(y_ivhd_10wpi) %>% arrange(NES) %>% print(n = 100)

# genes_909733 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-909733') %>% pull(core_enrichment), '/')[[1]]
# genes_936964 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-936964') %>% pull(core_enrichment), '/')[[1]]
# genes_918233 <- str_split(as_tibble(y_ivhd) %>% filter(ID == 'R-HSA-918233') %>% pull(core_enrichment), '/')[[1]]
# the genes pulled from core_enrichment do not match the ones in viewPathway. there are more genes in viewPathway.

ivhd_10wpi_pathways <- as_tibble(y_ivhd_10wpi) %>% arrange(NES) %>% dplyr::select(., Description, NES)
write_tsv(ivhd_10wpi_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/ivhd_10wpi_gsePathways.tsv')

viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = hkidney_entrez_gene_list_ivhd_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/ivhd_activation_irf3-irf7_mediated_tbk1_ikke.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')


viewPathway('TRAF3-dependent IRF activation pathway', readable = T, foldChange = hkidney_entrez_gene_list_ivhd_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/ivhd_traf3-dependent_irf_activation.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')


viewPathway('TRAF6 mediated IRF7 activation', readable = T, foldChange = hkidney_entrez_gene_list_ivhd_10wpi)
ggsave(filename = 
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/ivhd_traf6_mediated_irf7_activation.png', 
       width = 1600, 
       height = 1000, 
       units = "px", 
       dpi = 72, 
       bg = 'white')


# genes_2607 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0032607') %>% pull(core_enrichment),'/')[[1]]
# genes_2647 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0032647') %>% pull(core_enrichment), '/')[[1]]
# genes_61844 <- str_split(bottom10_low_nes %>% filter(ID == 'GO:0061844') %>% pull(core_enrichment), '/')[[1]]

# Convert to a Markdown table ----
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/ivhd_10wpi_gsePathways.tsv', header = TRUE, sep = "\t")
cat(markdown_table(data), sep = "\n")

## IV-LD ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_ivld_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_ivld_vs_conu_10wpi, 'hkidney', 'ivld', '10wpi')
save(hkidney_gsea_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivld_10wpi.RData')

hkidney_gsea_simplified_results_ivld_10wpi <-
  simplify(hkidney_gsea_results_ivld_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivld_10wpi.RData')

hkidney_entrez_gene_list_ivld_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivld_10wpi.RData')

nrow(hkidney_gsea_results_ivld_10wpi)  # 10 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivld_10wpi)  # 7 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivld_10wpi.RData')

# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(hkidney_gsea_simplified_results_ivld_10wpi@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_ivld_10wpi@result) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(low_high_nes_ivld_10wpi$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_ivld_10wpi$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_ivld_10wpi$Count))) +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 1WPC, hkidney tissue') +
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
         size = guide_legend(override.aes = list(
           shape = 1,
           fill = NA,
           stroke = .5,
           color = 'red'
         ))) +
  facet_grid(. ~ Regulation)

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/ivld_10wpi.png', width = 1000, height = 1023, units = "px", dpi = 72)


y_ivld_10wpi <- gsePathway(
  hkidney_entrez_gene_list_ivld_10wpi,
  pvalueCutoff = .2,
  pAdjustMethod = 'BH',
  eps = 1e-300,
  nPermSimple = 100000,
  verbose = F
)  # no enriched terms


