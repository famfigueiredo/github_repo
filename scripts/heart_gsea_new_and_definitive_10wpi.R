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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi'
)

load_results('^heart_.*_conu_10wpi')

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/heart_res_dnavaccine_vs_conu_10wpi.RData')

gsea_formatting(res_dnavaccine_vs_conu_10wpi, 'heart', 'dnavaccine', '10wpi')
save(heart_gsea_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_results_dnavaccine_10wpi.RData')

heart_gsea_simplified_results_dnavaccine_10wpi <-
  clusterProfiler::simplify(heart_gsea_results_dnavaccine_10wpi)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_simplified_results_dnavaccine_10wpi.RData')

heart_entrez_gene_list_dnavaccine_10wpi <- entrez_gene_list
save(heart_entrez_gene_list_dnavaccine_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_dnavaccine_10wpi.RData')

nrow(heart_gsea_results_dnavaccine_10wpi)  # 681 GO terms/pathways
nrow(heart_gsea_simplified_results_dnavaccine_10wpi)  # 227 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <-
  as_tibble(heart_gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_dnavaccine_10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

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
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 200)) +
  scale_size_continuous(
    'Set size',
    range = c(2, 10),
    guide = 'legend',
    limits = c(2, max(low_high_nes_dnavaccine_10wpi$setSize))
  ) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_10wpi$Count))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 10WPI, heart tissue') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_definitive_plots/dnavaccine_10wpi.png', width = 1000, height = 1023, units = "px", dpi = 72)

## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/heart_res_eomes_vs_conu_10wpi.RData')

gsea_formatting(res_eomes_vs_conu_10wpi, 'heart', 'eomes', '10wpi')
save(heart_gsea_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_results_eomes_10wpi.RData')

heart_gsea_simplified_results_eomes_10wpi <-
  clusterProfiler::simplify(heart_gsea_results_eomes_10wpi)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_simplified_results_eomes_10wpi.RData')

heart_entrez_gene_list_eomes_10wpi <- entrez_gene_list
save(heart_entrez_gene_list_eomes_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_eomes_10wpi.RData')

nrow(heart_gsea_results_eomes_10wpi)  # 2680 GO terms/pathways
nrow(heart_gsea_simplified_results_eomes_10wpi)  # 615 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top20_high_nes <- 
  as_tibble(heart_gsea_simplified_results_eomes_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_results_eomes_10wpi@result) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top20_high_nes %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(top10_high_nes$Count))) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(top10_high_nes$setSize))) +
  scale_x_continuous(limits = c(0, max(top10_high_nes$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
          subtitle = 'EOMES, 10WPI, heart, human orths') +
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
    color = guide_legend(override.aes = list(size = 5)),
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))
  )
  facet_grid(.~Regulation)

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_definitive_plots/eomes_10wpi_upregulated.png', width = 1000, height = 1023, units = "px", dpi = 72)
  
## GATA3 ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE
  
### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/heart_res_gata3_vs_conu_10wpi.RData')

gsea_formatting(res_gata3_vs_conu_10wpi, 'heart', 'gata3', '10wpi')
save(heart_gsea_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_results_gata3_10wpi.RData')

heart_gsea_simplified_results_gata3_10wpi <-
  clusterProfiler::simplify(heart_gsea_results_gata3_10wpi)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_simplified_results_gata3_10wpi.RData')

heart_entrez_gene_list_gata3_10wpi <- entrez_gene_list
save(heart_entrez_gene_list_gata3_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_gata3_10wpi.RData')

nrow(heart_gsea_results_gata3_10wpi)  # 3075 GO terms/pathways
nrow(heart_gsea_simplified_results_gata3_10wpi)  # 658 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top20_high_nes_gata3 <- 
  as_tibble(heart_gsea_simplified_results_gata3_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_results_gata3_10wpi) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = NES) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top20_high_nes_gata3 %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'legend', breaks = seq(0, 400, by = 50)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', breaks = seq(0, max(top20_high_nes_gata3$setSize), by = 50)) +
  scale_x_continuous(limits = c(0, max(top20_high_nes_gata3$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
          subtitle = 'GATA3, 10WPI, heart, human orths') +
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_definitive_plots/gata3_10wpi_upregulated.png', width = 1000, height = 1023, units = "px", dpi = 72)

## IV-HD ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/heart_res_ivhd_vs_conu_10wpi.RData')

gsea_formatting(res_ivhd_vs_conu_10wpi, 'heart', 'ivhd', '10wpi')
save(heart_gsea_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_results_ivhd_10wpi.RData')

heart_gsea_simplified_results_ivhd_10wpi <-
  clusterProfiler::simplify(heart_gsea_results_ivhd_10wpi)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_simplified_results_ivhd_10wpi.RData')

heart_entrez_gene_list_ivhd_10wpi <- entrez_gene_list
save(heart_entrez_gene_list_ivhd_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_ivhd_10wpi.RData')

nrow(heart_gsea_results_ivhd_10wpi)  # 1779 GO terms/pathways
nrow(heart_gsea_simplified_results_ivhd_10wpi)  # 498 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top20_high_nes <- 
  as_tibble(heart_gsea_simplified_results_ivhd_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_results_ivhd_10wpi) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = NES) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))
# 
# low_high_nes_dnavaccine_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top20_high_nes %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, 300)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(top20_high_nes$setSize))) +
  scale_x_continuous(limits = c(0, max(top20_high_nes$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
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

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_definitive_plots/ivhd_10wpi_upregulated.png', width = 1000, height = 1023, units = "px", dpi = 72)
  

## IV-LD ----

### all genes ###
# gsea formatting starting from a DESeq results table
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/heart_res_ivld_vs_conu_10wpi.RData')

gsea_formatting(res_ivld_vs_conu_10wpi, 'heart', 'ivld', '10wpi')
save(heart_gsea_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_results_ivld_10wpi.RData')

heart_gsea_simplified_results_ivld_10wpi <-
  clusterProfiler::simplify(heart_gsea_results_ivld_10wpi)  # simplifying GO terms to reduce redundancy
save(heart_gsea_simplified_results_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_gsea_simplified_results_ivld_10wpi.RData')

heart_entrez_gene_list_ivld_10wpi <- entrez_gene_list
save(heart_entrez_gene_list_ivld_10wpi, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_ivld_10wpi.RData')

nrow(heart_gsea_results_ivld_10wpi)  # 2714 GO terms/pathways
nrow(heart_gsea_simplified_results_ivld_10wpi)  # 621 GO terms/pathways

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top20_high_nes_ivld <- 
  as_tibble(heart_gsea_simplified_results_ivld_10wpi) %>%
  filter(NES > 0) %>% 
  arrange(desc(setSize)) %>% 
  top_n(20, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <- 
#   as_tibble(heart_gsea_simplified_results_ivld_10wpi) %>%
#   filter(NES < 0) %>% 
#   arrange(desc(setSize)) %>% 
#   top_n(10, wt = setSize) %>% 
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))
# 
# low_high_nes_ivld_10wpi <- bind_rows(top10_high_nes, bottom10_low_nes)

top20_high_nes_ivld %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  scale_color_viridis_c('Gene count', guide = 'colourbar', limits = c(2, max(top20_high_nes_ivld$Count))) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(top20_high_nes_ivld$setSize))) +
  scale_x_continuous(limits = c(0, max(top20_high_nes_ivld$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, upregulated genes',
          subtitle = 'IV-LD, 10WPI, heart, human orths') +
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
    color = guide_legend(override.aes = list(size = 5)),
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))
  ) +
  facet_grid(. ~ Regulation)


ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_definitive_plots/heart_ivld_10wpi_upregulated.png', width = 1000, height = 1023, units = "px", dpi = 72)
















