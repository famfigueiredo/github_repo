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
  clusterProfiler::simplify(hkidney_gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_10wpi.RData')

hkidney_entrez_gene_list_dnavaccine_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_10wpi.RData')

nrow(hkidney_gsea_results_dnavaccine_10wpi)  # 2 GO terms/pathways
nrow(hkidney_gsea_simplified_results_dnavaccine_10wpi)  # 2 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_10wpi.RData')

## EOMES ----
### all genes ###
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_eomes_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_eomes_vs_conu_10wpi, 'hkidney', 'eomes', '10wpi')
save(hkidney_gsea_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_eomes_10wpi.RData')

hkidney_gsea_simplified_results_eomes_10wpi <-
  clusterProfiler::simplify(hkidney_gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_eomes_10wpi.RData')

hkidney_entrez_gene_list_eomes_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_eomes_10wpi.RData')

nrow(hkidney_gsea_results_eomes_10wpi)  # 0 GO terms/pathways
nrow(hkidney_gsea_simplified_results_eomes_10wpi)  # 0 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_eomes_10wpi.RData')

## GATA3 ----

# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_gata3_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_gata3_vs_conu_10wpi, 'hkidney', 'gata3', '10wpi')
save(hkidney_gsea_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_gata3_10wpi.RData')

hkidney_gsea_simplified_results_gata3_10wpi <-
  clusterProfiler::simplify(hkidney_gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_gata3_10wpi.RData')

hkidney_entrez_gene_list_gata3_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_gata3_10wpi.RData')

nrow(hkidney_gsea_results_gata3_10wpi)  # 12 GO terms/pathways
nrow(hkidney_gsea_simplified_results_gata3_10wpi)  # 12 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_gata3_10wpi.RData')

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
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top10_high_nes %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(top10_high_nes$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(top10_high_nes$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(top10_high_nes$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, upregulated genes',
          subtitle = 'GATA3, 10WPI, hkidney tissue') +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_definitive_plots/gata3_10wpi_upregulated.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)


## IV-HD ----
# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_ivhd_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_ivhd_vs_conu_10wpi, 'hkidney', 'ivhd', '10wpi')
save(hkidney_gsea_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivhd_10wpi.RData')

hkidney_gsea_simplified_results_ivhd_10wpi <-
  clusterProfiler::simplify(hkidney_gsea_results_ivhd_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_10wpi.RData')

hkidney_entrez_gene_list_ivhd_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivhd_10wpi.RData')

nrow(hkidney_gsea_results_ivhd_10wpi)  # 1940 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivhd_10wpi)  # 524 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivhd_10wpi.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top20_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_10wpi) %>% 
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top20_high_nes %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(top20_high_nes$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(top20_high_nes$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(top20_high_nes$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, upregulated genes',
          subtitle = 'IV-HD, 10WPI, hkidney tissue') +
  theme_bw(base_size = 20) +
  theme(
    text = element_text(family = 'Times New Roman'),
    legend.title = element_text(size = 20),
    legend.text = element_text(size = 18),
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

## IV-LD ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/hkidney_res_ivld_vs_conu_10wpi.RData')

gsea_formatting(hkidney_res_ivld_vs_conu_10wpi, 'hkidney', 'ivld', '10wpi')
save(hkidney_gsea_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivld_10wpi.RData')

hkidney_gsea_simplified_results_ivld_10wpi <-
  clusterProfiler::simplify(hkidney_gsea_results_ivld_10wpi)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivld_10wpi.RData')

hkidney_entrez_gene_list_ivld_10wpi <- entrez_gene_list
save(hkidney_entrez_gene_list_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivld_10wpi.RData')

nrow(hkidney_gsea_results_ivld_10wpi)  # 0 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivld_10wpi)  # 0 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_gsea_simplified_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivld_10wpi.RData')

# No enriched pathways

as_tibble(hkidney_gsea_results_dnavaccine_10wpi)
