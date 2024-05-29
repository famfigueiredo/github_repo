library('tidyverse')
library('VennDiagram')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('enrichplot')
library('AnnotationDbi')
library(ReactomePA)

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc'
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
gsea_formatting <- function(results_table, treatment, sampling_point) {
  # Install and load required packages
  required_packages <- c('dplyr', 'gprofiler2', 'clusterProfiler', 'org.Hs.eg.db')
  installed_packages <- rownames(installed.packages())
  
  for (pkg in required_packages) {
    if (!(pkg %in% installed_packages)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
  
  # Convert rownames to column 'ensembl'
  results_df <- tibble::rownames_to_column(as.data.frame(results_table), var = 'ensembl')
  
  # Convert salmon genes to human orthologs
  orthologs <- gorth(
    query = rownames(results_table),
    source_organism = 'ssalar',
    target_organism = 'hsapiens',
    mthreshold = 1,
    filter_na = TRUE
  )
  
  # Select relevant variables and join with ortholog data
  merged_df <- results_df %>%
    left_join(orthologs, by = c('ensembl' = 'input')) %>%
    dplyr::select(
      ensembl,
      ortholog_name,
      ortholog_ensg,
      log2FoldChange,
      padj,
      description
    ) %>%
    na.omit()
  
  # Order genes by fold change
  ordered_df <- merged_df[order(-merged_df$log2FoldChange), ]
  
  # Prepare matrix for GSEA
  gene_list <- ordered_df$log2FoldChange
  names(gene_list) <- ordered_df$ortholog_name
  
  # Run GSEA
  gsea_results <- gseGO(
    gene_list,
    keyType = 'SYMBOL',
    OrgDb = org.Hs.eg.db,
    ont = 'BP',
    pvalueCutoff = 0.05,
    pAdjustMethod = 'BH',
    verbose = T,
    eps = 1e-300
  )
  
  # Assign the results to a variable including treatment and sampling_point in the name
  results_name <- paste0('gsea_results_', treatment, '_', sampling_point)
  assign(results_name, gsea_results, envir = .GlobalEnv)
  
  # Print GSEA results
  print(gsea_results)
  
  return(gsea_results)
}

gsea_formatting(res_dnavaccine_vs_conu_1wpc, 'dnavaccine', '1wpc')

gsea_simplified_results_dnavaccine_1wpc <- simplify(gsea_results_dnavaccine_1wpc)

# Convert the GSEA results to a tibble and retrieve top 10 highest and lowest NES
top10_high_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_1wpc@result) %>%
  arrange(desc(NES)) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_dnavaccine_1wpc@result) %>% 
  arrange(NES) %>% 
  top_n(8, wt = desc(NES)) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_1wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_dnavaccine_1wpc %>%
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
    panel.grid = element_line(color = 'black', size = .05, linetype = 2)
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5))  # show only borders in set size legend
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

write_tsv(dnavaccine_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/dnavaccine_downregulated_pathways.tsv')


# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/dnavaccine_downregulated_pathways.tsv', header = TRUE, sep = "\t")

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
gsea_formatting_significant <- function(results_table, treatment, sampling_point) {
  # Wrangle the data using the improved_data_wrangling function
  wrangled_data <- improved_data_wrangling(results_table, treatment, sampling_point)
  
  # Select ortholog_name and log2FC, and remove NA values
  ortholog_fc <- wrangled_data %>%
    dplyr::select(ortholog_name, log2FC) %>%
    na.omit()
  
  # Convert gene symbols to Entrez IDs
  entrez_ids <- bitr(ortholog_fc$ortholog_name, fromType = 'SYMBOL', toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
  
  # Merge ortholog_fc with entrez_ids and select relevant columns
  enrichment <- ortholog_fc %>%
    left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>%
    dplyr::select(ENTREZID, log2FC)
  
  # Prepare the enrichment_gsea matrix
  enrichment_gsea <- enrichment$log2FC
  names(enrichment_gsea) <- enrichment$ENTREZID
  
  # Create a dynamic name for the enrichment_gsea object
  results_name <- paste0("enrichment_gsea_", treatment, "_", sampling_point)
  assign(results_name, enrichment_gsea, envir = .GlobalEnv)
  
  return(enrichment_gsea)
}

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
gsea_formatting(res_eomes_vs_conu_1wpc, 'eomes', '1wpc')

gsea_simplified_results_eomes_1wpc <- simplify(gsea_results_eomes_1wpc)

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

as_tibble(y_eomes) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Activation of IRF3, IRF7 mediated by TBK1, IKKε (IKBKE)', readable = T, foldChange = entrez_gene_list)

eomes_pathways <- as_tibble(y_eomes) %>% arrange(NES) %>% filter(., NES < 0) %>% dplyr::select(., Description, NES) 

write_tsv(eomes_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/eomes_downregulated_pathways.tsv')


# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc/pathways/eomes_downregulated_pathways.tsv', header = TRUE, sep = "\t")

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

as_tibble(gsea_simplified_results_gata3_1wpc)

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
  ) + guides(
    color = guide_legend(override.aes = list(size = 5)),  # increase point size in gene count legend
    size = guide_legend(override.aes = list(shape =1, fill = NA, stroke = .5))  # show only borders in set size legend
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

# Convert to a Markdown table ----
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






