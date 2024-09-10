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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# GSEA ----
## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
## added saving feature on 24/06/2024, so I don't have to re-run the formatting function everytime
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/hkidney_res_dnavaccine_vs_conu_4wpc.RData')

gsea_formatting(hkidney_res_dnavaccine_vs_conu_4wpc, 'hkidney', 'dnavaccine', '4wpc')
save(hkidney_gsea_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_dnavaccine_4wpc.RData')

hkidney_gsea_simplified_results_dnavaccine_4wpc <-
  simplify(hkidney_gsea_results_dnavaccine_4wpc)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_4wpc.RData')

hkidney_entrez_gene_list_dnavaccine_4wpc <- entrez_gene_list
save(hkidney_entrez_gene_list_dnavaccine_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_4wpc.RData')

nrow(hkidney_gsea_results_dnavaccine_4wpc)  # 100 GO terms/pathways
nrow(hkidney_gsea_simplified_results_dnavaccine_4wpc)  # 51 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_dnavaccine_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <-
  as_tibble(hkidney_gsea_simplified_results_dnavaccine_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_dnavaccine_4wpc) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # only 6 terms had NES < 0

low_high_nes_dnavaccine_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(low_high_nes_dnavaccine_4wpc$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_dnavaccine_4wpc$setSize), by = 25)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_dnavaccine_4wpc$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 4WPC, head-kidney tissue') +
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
         '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_definitive_plots/dnavaccine_4wpc.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

## EOMES ----
### all genes ###
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/hkidney_res_eomes_vs_conu_4wpc.RData')

gsea_formatting(hkidney_res_eomes_vs_conu_4wpc, 'hkidney', 'eomes', '4wpc')
save(hkidney_gsea_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_eomes_4wpc.RData')

hkidney_gsea_simplified_results_eomes_4wpc <-
  simplify(hkidney_gsea_results_eomes_4wpc)  # Couldn't simplify. Error in as.list.default(X) : no method for coercing this S4 class to a vector
save(hkidney_gsea_simplified_results_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_eomes_4wpc.RData')

hkidney_entrez_gene_list_eomes_4wpc <- entrez_gene_list
save(hkidney_entrez_gene_list_eomes_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_eomes_4wpc.RData')

nrow(hkidney_gsea_results_eomes_4wpc)  # 711 GO terms/pathways
nrow(hkidney_gsea_simplified_results_eomes_4wpc)  # 256 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_eomes_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_eomes_4wpc.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_results_eomes_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_results_eomes_4wpc) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

top10_high_nes %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.1*max(top10_high_nes$Count))) +
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
          subtitle = 'EOMES, 4WPC, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_definitive_plots/eomes_4wpc_upregulated.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

## GATA3 ----

# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/hkidney_res_gata3_vs_conu_4wpc.RData')

gsea_formatting(hkidney_res_gata3_vs_conu_4wpc, 'hkidney', 'gata3', '4wpc')
save(hkidney_gsea_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_gata3_4wpc.RData')

hkidney_gsea_simplified_results_gata3_4wpc <-
  simplify(hkidney_gsea_results_gata3_4wpc)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_gata3_4wpc.RData')

hkidney_entrez_gene_list_gata3_4wpc <- entrez_gene_list
save(hkidney_entrez_gene_list_gata3_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_gata3_4wpc.RData')

nrow(hkidney_gsea_results_gata3_4wpc)  # 13 GO terms/pathways
nrow(hkidney_gsea_simplified_results_gata3_4wpc)  # 6 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_gata3_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_gata3_4wpc.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_gata3_4wpc) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated pathways

low_high_nes_gata3_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

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
          subtitle = 'GATA3, 4wpc, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_definitive_plots/gata3_4wpc_upregulated.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

## IV-HD ----
# running gsea and saving files for later use
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/hkidney_res_ivhd_vs_conu_4wpc.RData')

gsea_formatting(hkidney_res_ivhd_vs_conu_4wpc, 'hkidney', 'ivhd', '4wpc')
save(hkidney_gsea_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_ivhd_4wpc.RData')

hkidney_gsea_simplified_results_ivhd_4wpc <-
  simplify(hkidney_gsea_results_ivhd_4wpc)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_4wpc.RData')

hkidney_entrez_gene_list_ivhd_4wpc <- entrez_gene_list
save(hkidney_entrez_gene_list_ivhd_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_ivhd_4wpc.RData')

nrow(hkidney_gsea_results_ivhd_4wpc)  # 217 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivhd_4wpc)  # 100 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_ivhd_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_ivhd_4wpc.RData')


# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(hkidney_gsea_simplified_results_ivhd_4wpc) %>% 
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(low_high_nes_ivhd_4wpc$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_ivhd_4wpc$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_ivhd_4wpc$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-HD, 4wpc, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_definitive_plots/ivhd_4wpc.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

## IV-LD ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/hkidney_res_ivld_vs_conu_4wpc.RData')

gsea_formatting(hkidney_res_ivld_vs_conu_4wpc, 'hkidney', 'ivld', '4wpc')
save(hkidney_gsea_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_ivld_4wpc.RData')

hkidney_gsea_simplified_results_ivld_4wpc <-
  simplify(hkidney_gsea_results_ivld_4wpc)  # simplifying GO terms to reduce redundancy
save(hkidney_gsea_simplified_results_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_ivld_4wpc.RData')

hkidney_entrez_gene_list_ivld_4wpc <- entrez_gene_list
save(hkidney_entrez_gene_list_ivld_4wpc, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_ivld_4wpc.RData')

nrow(hkidney_gsea_results_ivld_4wpc)  # 0 GO terms/pathways
nrow(hkidney_gsea_simplified_results_ivld_4wpc)  # 0 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_gsea_simplified_results_ivld_4wpc.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_results_tables/hkidney_entrez_gene_list_ivld_4wpc.RData')

# convert the gsea results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <-
  as_tibble(hkidney_gsea_simplified_results_ivld_4wpc@result) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(hkidney_gsea_simplified_results_ivld_4wpc@result) %>%
  filter(NES < 0) %>%
  arrange(setSize) %>%
  top_n(10, wt = NES) %>%
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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 1.2*max(low_high_nes_ivld_4wpc$Count))) +
  scale_size_continuous(
    'Set size', 
    range = c(2, 10), 
    guide = 'legend', 
    breaks = seq(0, max(low_high_nes_ivld_4wpc$setSize), by = 100)) +
  scale_x_continuous(limits = c(0, 1.1 * max(low_high_nes_ivld_4wpc$Count))) +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 4wpc, hkidney tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/gsea_definitive_plots/ivld_4wpc.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 72)

