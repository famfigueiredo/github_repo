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


y_dnavaccine <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

dnavaccine_pathways <- as_tibble(y_dnavaccine) %>% arrange(-NES) %>% filter(., NES > 0) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% filter(., NES > 0) %>% pull(Description)

write_tsv(dnavaccine_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/dnavaccine_upregulated_pathways.tsv')


# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/dnavaccine_upregulated_pathways.tsv', header = TRUE, sep = "\t")

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
  geneList = significant_enrichment_gsea_dnavaccine_4wpc,
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

# rm(list = ls(pattern = '^enrichment*'))

bottom20_low_nes %>% slice(1:10) %>% pull(Description)

x_dnavaccine <- gsePathway(entrez_gene_list,
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(x_dnavaccine) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Interleukin-4 and Interleukin-13 signaling', readable = T, foldChange = entrez_gene_list)

# Some pathways of interest:
# Interleukin-4 and Interleukin-13 signaling
# Interferon gamma signaling
# FCGR3A-mediated phagocytosis
# Signaling by the B Cell Receptor (BCR)  

significant_dnavaccine_4wpc_pathways <- as_tibble(simplified_dnavaccine_4wpc_significant) %>%  arrange(-NES) %>% dplyr::select(., Description, NES)

write_tsv(significant_dnavaccine_4wpc_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/significant_dnavaccine_pathways.tsv')


# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/significant_dnavaccine_pathways.tsv', header = TRUE, sep = "\t")

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




## EOMES ----
### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_eomes_vs_conu_4wpc, 'eomes', '4wpc')

gsea_simplified_results_eomes_4wpc <- simplify(gsea_results_eomes_4wpc)  # simplifying GO terms to reduce redundancy
as_tibble(gsea_simplified_results_eomes_4wpc) %>% arrange(NES) %>% print(n = 100)

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

low_high_nes_eomes_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_eomes_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_eomes_4wpc$Count))) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_eomes_4wpc$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_eomes_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'EOMES, 4WPC, heart, human orths') +
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


y_eomes <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_eomes) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

eomes_pathways <- as_tibble(y_eomes) %>% 
  arrange(-NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(eomes_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/eomes_upregulated_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/eomes_upregulated_pathways.tsv', header = TRUE, sep = "\t")

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




## Testing different annotation files
### Zebrafish ----
BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)

gsea_formatting(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

gsea_simplified_results_dnavaccine_4wpc <- simplify(gsea_results_dnavaccine_4wpc)

results_df <-
  tibble::rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_4wpc), var = 'ensembl')

# Convert salmon genes to human orthologs
orthologs <- gorth(
  query = results_df$ensembl,
  source_organism = 'ssalar',
  target_organism = 'drerio',
  mthreshold = 1,
  filter_na = TRUE
)

# Select relevant variables and join with ortholog data
merged_df <- results_df %>%
  left_join(orthologs, by = c('ensembl' = 'input')) %>%
  dplyr::select(ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                description) %>%
  na.omit()

# Order genes by fold change
ordered_df <- merged_df[order(-merged_df$log2FoldChange),]

# Prepare matrix for GSEA
gene_list <- ordered_df$log2FoldChange
names(gene_list) <- ordered_df$ortholog_name

# Prepare matrix for gsePathway
ordered_entrez <-
  bitr(ordered_df$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Dr.eg.db)  # 8.91% of input gene IDs are fail to map...
entrez_genes <-
  ordered_df %>% left_join(
    ordered_entrez,
    by = c('ortholog_name' = 'SYMBOL'),
    relationship = 'many-to-many'
  ) %>% dplyr::select(ENTREZID, log2FoldChange)
distinct_genes <-
  entrez_genes %>% distinct(ENTREZID, .keep_all = T)
entrez_gene_list <- distinct_genes$log2FoldChange
names(entrez_gene_list) <- distinct_genes$ENTREZID

# Run GSEA
gsea_results <- gseGO(
  gene_list,
  keyType = 'SYMBOL',
  OrgDb = org.Dr.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T,
  eps = 1e-300
)

as_tibble(gsea_results)
gsea_simplified_results_zebrafish <- simplify(gsea_results)



top10_high_nes <- 
  as_tibble(gsea_simplified_results_zebrafish) %>%
  filter(NES > 0) %>% 
  arrange(-setSize) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_zebrafish) %>%
  filter(NES < 0) %>% 
  arrange(-setSize) %>% 
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_4wpc <- bind_rows(bottom10_low_nes, top10_high_nes)


low_high_nes_dnavaccine_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(1, 250)) +
  scale_size_continuous('Set size', range = c(2, 15), guide = 'legend', limits = c(2, max(low_high_nes_dnavaccine_4wpc$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_4wpc$Count) * 1.1)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 4WPC, zebrafish orthologs') +
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
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)


y_dnavaccine <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                           organism = 'zebrafish',
                           pvalueCutoff = .2,
                           pAdjustMethod = 'BH',
                           verbose = F)

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Signaling by Rho GTPases, Miro GTPases and RHOBTB3', readable = T, organism = 'zebrafish', foldChange = entrez_gene_list)

dnavaccine_pathways <- as_tibble(y_dnavaccine) %>% arrange(-NES) %>% filter(., NES > 0) %>% dplyr::select(., Description, NES) 

as_tibble(y_dnavaccine) %>% arrange(-NES) %>% filter(., NES > 0) %>% pull(Description)


### Atlantic salmon ----
library(AnnotationHub)

ah <- AnnotationHub()

# query(ah, c('OrgDb', 'Salmo salar'))
# AnnotationHub with 1 record
# snapshotDate(): 2024-04-29
# names(): AH114250
# $dataprovider: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/
# $species: Salmo salar
# $rdataclass: OrgDb
# $rdatadateadded: 2023-10-20
# $title: org.Salmo_salar.eg.sqlite
# $description: NCBI gene ID based annotations about Salmo salar
# $taxonomyid: 8030
# $genome: NCBI genomes
# $sourcetype: NCBI/UniProt
# $sourceurl: ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/, ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmappi...
# $sourcesize: NA
# $tags: c("NCBI", "Gene", "Annotation") 
# retrieve record with 'object[["AH114250"]]' 

sasa <- query(ah, c('OrgDb', 'Salmo salar'))[[1]]  # querying AnnotationHub for an OrgDb object for Atlantic salmon


# Select relevant variables
merged_df <- results_df %>%
  dplyr::select(ensembl,
                log2FoldChange,
                padj) %>%
  na.omit()

# Order genes by fold change
ordered_df <- merged_df[order(-merged_df$log2FoldChange),]

head(ordered_df)

# Prepare matrix for GSEA
gene_list <- ordered_df$log2FoldChange
names(gene_list) <- ordered_df$ensembl

# Prepare matrix for gsePathway
ordered_entrez <-
  bitr(ordered_df$ensembl, 'ENSEMBL', 'ENTREZID', OrgDb = sasa)  # 8.91% of input gene IDs are fail to map...

head(ordered_df)
head(ordered_entrez)


entrez_genes <-
  ordered_df %>% left_join(
    ordered_entrez,
    by = c('ensembl' = 'ENSEMBL'),
    relationship = 'many-to-many'
  ) %>% dplyr::select(ENTREZID, log2FoldChange)

head(entrez_genes)

distinct_genes <-
  entrez_genes %>% distinct(ENTREZID, .keep_all = T)

entrez_gene_list <- distinct_genes$log2FoldChange
names(entrez_gene_list) <- distinct_genes$ENTREZID

# Run GSEA
gsea_results <- gseGO(
  gene_list,
  keyType = 'ENSEMBL',
  OrgDb = sasa,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T,
  eps = 1e-300
)

as_tibble(gsea_results)
gsea_simplified_results_salmon <- simplify(gsea_results)


top10_high_nes <- 
  as_tibble(gsea_simplified_results_salmon) %>%
  filter(NES > 0) %>% 
  arrange(-setSize) %>% 
  top_n(10, wt = NES) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <- 
  as_tibble(gsea_simplified_results_salmon) %>%
  filter(NES < 0) %>% 
  arrange(-setSize) %>% 
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_dnavaccine_4wpc <- bind_rows(bottom10_low_nes, top10_high_nes)


low_high_nes_dnavaccine_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(1, max(low_high_nes_dnavaccine_4wpc$Count))) +
  scale_size_continuous('Set size', range = c(2, 15), guide = 'legend', limits = c(2, max(low_high_nes_dnavaccine_4wpc$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_dnavaccine_4wpc$Count) * 1.1)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'DNA vaccine, 4WPC, Atlantic salmon') +
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
    size = guide_legend(override.aes = list(shape = 1, fill = NA, stroke = .5, color = 'red'))  # show only borders in set size legend
  ) +
  facet_grid(. ~ Regulation)



## GATA3 ----

### all genes ###
# gsea formatting starting from a DESeq results table
gsea_formatting(res_gata3_vs_conu_4wpc, 'gata3', '4wpc')

gsea_simplified_results_gata3_4wpc <- simplify(gsea_results_gata3_4wpc)  # simplifying GO terms to reduce redundancy

as_tibble(gsea_simplified_results_gata3_4wpc) %>% arrange(NES) %>% print(n = 100)

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
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_gata3_4wpc <- bind_rows(top10_high_nes, bottom10_low_nes)

low_high_nes_gata3_4wpc %>%
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(size = setSize), shape = 1, stroke = 0.2, color = 'red') +
  geom_point(aes(color = Count, size = Count), shape = 16) +
  # scale_color_viridis_c('Gene set') +
  scale_color_viridis_c('Gene count', guide = 'legend', limits = c(2, max(low_high_nes_gata3_4wpc$Count))) +
  scale_size_continuous('Set size', range = c(2, 10), guide = 'legend', limits = c(2, max(low_high_nes_gata3_4wpc$setSize))) +
  scale_x_continuous(limits = c(0, max(low_high_nes_gata3_4wpc$Count * 1.1))) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +  
  ggtitle('GSEA, downregulated vs upregulated genes',
          subtitle = 'GATA3, 4WPC, heart, human orths') +
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


y_gata3 <- gsePathway(entrez_gene_list,  # the gsea_formatting function removes the duplicates from this object
                      pvalueCutoff = .2,
                      pAdjustMethod = 'BH',
                      verbose = F)

as_tibble(y_gata3) %>% arrange(-NES) %>% print(n = 100)

viewPathway('Interferon alpha/beta signaling', readable = T, foldChange = entrez_gene_list)

gata3_pathways <- as_tibble(y_gata3) %>% 
  arrange(-NES) %>%   
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  dplyr::select(., Description, NES, setSize, Count) 

write_tsv(gata3_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/gata3_pathways.tsv')

# Convert to a Markdown table ----
# Read the TSV file
data <- read.delim('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc/pathways/gata3_pathways.tsv', header = TRUE, sep = "\t")

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




## Testing different annotation files
