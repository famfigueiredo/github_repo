library('tidyverse')
library('VennDiagram')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('enrichplot')
library('AnnotationDbi')
library('ReactomePA')

rm(list = setdiff(ls(), grep("^spleen_", ls(), value = TRUE)))

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi'
)

results_files <-
  list.files(pattern = '^spleen_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/spleen_res_dnavaccine_vs_conu_10wpi.RData')

gsea_formatting(spleen_res_dnavaccine_vs_conu_10wpi, 'spleen', 'dnavaccine', '10wpi')
save(spleen_gsea_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_dnavaccine_10wpi.RData')

spleen_gsea_simplified_results_dnavaccine_10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy
save(spleen_gsea_simplified_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_dnavaccine_10wpi.RData')

spleen_entrez_gene_list_dnavaccine_10wpi <- entrez_gene_list
save(spleen_entrez_gene_list_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_dnavaccine_10wpi.RData')

nrow(spleen_gsea_results_dnavaccine_10wpi)  # 1 GO terms/pathways
nrow(spleen_gsea_simplified_results_dnavaccine_10wpi)  # 1 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_dnavaccine_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_dnavaccine_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% # only 5 downregulated terms
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 250)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', breaks = seq(0, max(low_high_nes_dnavaccine_10wpi$setSize), by = 50)) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 10WPI, spleen, human orths') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_definitive_plots/dnavaccine_10wpi.png',
       width = 1100,
       height = 1000,
       units = 'px',
       dpi = 100)

## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/spleen_res_eomes_vs_conu_10wpi.RData')

gsea_formatting(spleen_res_eomes_vs_conu_10wpi, 'spleen', 'eomes', '10wpi')
save(spleen_gsea_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_eomes_10wpi.RData')

spleen_gsea_simplified_results_eomes_10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy
save(spleen_gsea_simplified_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_eomes_10wpi.RData')

spleen_entrez_gene_list_eomes_10wpi <- entrez_gene_list
save(spleen_entrez_gene_list_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_eomes_10wpi.RData')

nrow(spleen_gsea_results_eomes_10wpi)  # 2616 GO terms/pathways
nrow(spleen_gsea_simplified_results_eomes_10wpi)  # 602 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_eomes_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_eomes_10wpi.RData')
# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(spleen_gsea_simplified_results_eomes_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_results_eomes_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% # only 5 downregulated terms
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_eomes_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

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
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 300)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    breaks = seq(0, max(low_high_nes_eomes_10wpi$setSize), by = 50)
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_eomes_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated and upregulated genes',
          subtitle = 'EOMES, 10WPI, spleen, human orths') +
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
             ))) +  # show only borders in set size legend) 
             facet_grid(. ~ Regulation)

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_definitive_plots/eomes_10wpi.png',
       width = 1100,
       height = 1000,
       units = 'px',
       dpi = 100)

## GATA3 ----

### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/spleen_res_gata3_vs_conu_10wpi.RData')

gsea_formatting(spleen_res_gata3_vs_conu_10wpi, 'spleen', 'gata3', '10wpi')
save(spleen_gsea_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_gata3_10wpi.RData')

spleen_gsea_simplified_results_gata3_10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy
save(spleen_gsea_simplified_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_gata3_10wpi.RData')

spleen_entrez_gene_list_gata3_10wpi <- entrez_gene_list
save(spleen_entrez_gene_list_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_gata3_10wpi.RData')

nrow(spleen_gsea_results_gata3_10wpi)  # 1751 GO terms/pathways
nrow(spleen_gsea_simplified_results_gata3_10wpi)  # 481 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_gata3_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_gata3_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(spleen_gsea_simplified_results_gata3_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(spleen_gsea_simplified_results_gata3_10wpi) %>%
  filter(NES < 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_gata3_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)


low_high_nes_gata3_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', breaks = seq(0, 300, by = 50)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', breaks = seq(0, max(low_high_nes_gata3_10wpi$setSize), by = 50)) +
  scale_x_continuous(limits = c(0, max(low_high_nes_gata3_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated and upregulated genes',
          subtitle = 'GATA3, 10WPI, spleen, human orths') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_definitive_plots/gata3_10wpi.png', 
       width = 1000, 
       height = 1023, 
       units = "px", 
       dpi = 100)

## IV-HD ----

### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/spleen_res_ivhd_vs_conu_10wpi.RData')

gsea_formatting(spleen_res_ivhd_vs_conu_10wpi, 'spleen', 'ivhd', '10wpi')
save(spleen_gsea_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_ivhd_10wpi.RData')

spleen_gsea_simplified_results_ivhd_10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_ivhd_10wpi)  # simplifying GO terms to reduce redundancy
save(spleen_gsea_simplified_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_ivhd_10wpi.RData')

spleen_entrez_gene_list_ivhd_10wpi <- entrez_gene_list
save(spleen_entrez_gene_list_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivhd_10wpi.RData')

nrow(spleen_gsea_results_ivhd_10wpi)  # 1949 GO terms/pathways
nrow(spleen_gsea_simplified_results_ivhd_10wpi)  # 518 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_ivhd_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivhd_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(spleen_gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(spleen_gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES < 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_ivhd_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_ivhd_10wpi %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_ivhd_10wpi$Count))) +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_definitive_plots/ivhd_10wpi.png',
       width = 1100,
       height = 1000,
       units = 'px',
       dpi = 100)

## IV-LD ----

### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/spleen_res_ivld_vs_conu_10wpi.RData')

gsea_formatting(spleen_res_ivld_vs_conu_10wpi, 'spleen', 'ivld', '10wpi')
save(spleen_gsea_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_ivld_10wpi.RData')

spleen_gsea_simplified_results_ivld_10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_ivld_10wpi)  # simplifying GO terms to reduce redundancy
save(spleen_gsea_simplified_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_ivld_10wpi.RData')

spleen_entrez_gene_list_ivld_10wpi <- entrez_gene_list
save(spleen_entrez_gene_list_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivld_10wpi.RData')

nrow(spleen_gsea_results_ivld_10wpi)  # 575 GO terms/pathways
nrow(spleen_gsea_simplified_results_ivld_10wpi)  # 232 GO terms/pathways

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_gsea_simplified_results_ivld_10wpi.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivld_10wpi.RData')

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(spleen_gsea_simplified_results_ivld_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(spleen_gsea_simplified_results_ivld_10wpi) %>%
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
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(low_high_nes_ivld_10wpi$setSize))) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_ivld_10wpi$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_ivld_10wpi$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'IV-LD, 10WPI, spleen, human orths') +
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


ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_definitive_plots/ivld_10wpi.png',
       width = 1100,
       height = 1000,
       units = 'px',
       dpi = 100)


















