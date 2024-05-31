library('tidyverse')
library('VennDiagram')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('enrichplot')
library('AnnotationDbi')
library(ReactomePA)

rm(list = setdiff(ls(), grep("res_", ls(), value = TRUE)))

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

gsea_simplified_results_dnavaccine_4wpc <- simplify(gsea_results_dnavaccine_4wpc)  # simplifying GO terms to reduce redundancy
as_tibble(gsea_simplified_results_dnavaccine_4wpc)
# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES 
top10_high_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_4wpc@result) %>%
  arrange(desc(NES)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_4wpc@result) %>% 
  arrange(NES) %>% 
  top_n(8, wt = desc(NES)) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, 250)) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, 387)) +
  scale_x_continuous(limits = c(0, 250)) +
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
    panel.grid = element_line(color = 'black', size = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_dnavaccine <- gsePathway(entrez_gene_list,
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_dnavaccine) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

dnavaccine_pathways <- as_tibble(y_dnavaccine) %>% arrange(NES) %>% filter(., NES < 0) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine) %>% arrange(NES) %>% filter(., NES < 0) %>% pull(Description)

write_tsv(dnavaccine_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/dnavaccine_downregulated_pathways.tsv')


# Convert to a Markdown table ---
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/dnavaccine_downregulated_pathways.tsv', header = TRUE, sep = "\t")

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
# gsea formatting starting from a DESeq results table and using only significantly differentially regulated genes
gsea_formatting_significant(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

gsea_dnavaccine_4wpc_significant <- gseGO(
  geneList = enrichment_gsea_dnavaccine_4wpc,
  keyType = 'SYMBOL',
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T)

simplified_dnavaccine_4wpc_significant <- simplify(gsea_dnavaccine_4wpc_significant)

as_tibble(gsea_dnavaccine_4wpc_significant)
as_tibble(simplified_dnavaccine_4wpc_significant)

bottom20_low_nes <-
  as_tibble(simplified_dnavaccine_4wpc_significant) %>% 
  arrange(NES) %>% 
  top_n(20, wt = -NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

color.params = list(foldChange = enrichment_gsea_dnavaccine_4wpc, edge = TRUE)
cex.params = list(category_label = 0.6, gene_label = 0.6)

cnetplot(
  simplified_dnavaccine_4wpc_significant,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = T,
  showCategory = c('immune effector process', 'leukocyte activation', 'T cell activation'),
  max.overlaps = 500
)


bottom20_low_nes %>% slice(1:10) %>% pull(Description)

x_dnavaccine <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(x_dnavaccine) %>% arrange(NES) %>% print(n = 100)

viewPathway('Interleukin-4 and Interleukin-13 signaling', readable = T, foldChange = entrez_gene_list)

# Interleukin-4 and Interleukin-13 signaling
# Interferon gamma signaling
# FCGR3A-mediated phagocytosis
# Signaling by the B Cell Receptor (BCR)  



### Testing
wrangled_data <-
  improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

# Select ortholog_name and log2FC, and remove NA values
ortholog_fc <- wrangled_data %>%
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


