suppressPackageStartupMessages({
  library('tidyverse')
  library('VennDiagram')
  library('clusterProfiler')
  library('gprofiler2')
  library('org.Hs.eg.db')
  library('enrichplot')
  library('AnnotationDbi')
})

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# # Get names of all objects in global environment with a regex pattern matching 'res_.*_conu_4wpc'
# res_objects <- ls(pattern = "res_.*_4wpc")
# 
# # Apply significant_genes function to each res object
# for (res_object in res_objects) {
#   # Get the object from the global environment
#   result <- get(res_object)
#   
#   # Apply significant_genes function
#   significant_result <- significant_genes(result)
#   
#   # Assign the result back to the global environment
#   assign(paste0("significant_", res_object),
#          significant_result,
#          envir = .GlobalEnv)
# }
# 
# rm(list.data,
#    result,
#    i,
#    res_object,
#    res_objects,
#    results_files,
#    significant_result)
# 
# rm(list = ls(pattern = '^res_.'))  # removing results files from Global

# GSEA with significantly differentially regulated genes ----
## DNA vaccine
improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

string_test <-
  results_dnavaccine_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()  # selecting ortholog names and fold change

entrez_ids <-
  bitr(string_test$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)  # converting symbols to entrez for gseaGO

enrichment <-
  string_test %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)  # merging entrez ids with symbols and keeping just the entrez id

enrichment_gsea <- enrichment$log2FC  # preparing matrix for gseGO

names(enrichment_gsea) <- enrichment$ENTREZID  # preparing matrix for gseGO

gseGO(
  geneList = enrichment_gsea,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
) %>% as_tibble() -> gsea_dnavaccine_4wpc


gsea_dnavaccine_4wpc %>% 
  top_n(20, wt = abs(NES)) %>% 
  mutate(Regulation = ifelse(NES > 0, 'Upregulated', 'Downregulated')) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(color = Count, size = Count)) +
  scale_color_viridis_c('Gene count', guide = 'legend') +
  scale_size_continuous('Gene count', range = c(2, 15), guide = 'legend') +
  xlim(40, 200) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
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
    plot.subtitle = element_text(hjust = .5)
  ) +
  facet_grid(. ~ Regulation)

## EOMES
improved_data_wrangling(res_eomes_vs_conu_4wpc, 'eomes', '4wpc')

string_test_eomes <-
  results_eomes_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <-
  bitr(string_test_eomes$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_eomes <-
  string_test_eomes %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

enrichment_gsea_eomes <- enrichment_eomes$log2FC

names(enrichment_gsea_eomes) <- enrichment_eomes$ENTREZID

gsea_eomes <- gseGO(
  geneList = enrichment_gsea_eomes,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

gsea_eomes %>%
  fortify(., showCategory = 20) %>%
  as_tibble() %>%
  mutate_at('.sign', str_replace, 'suppressed', 'Downregulated') %>% mutate_at('.sign', str_replace, 'activated', 'Upregulated') %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_viridis_c('Adjusted p-value', guide = guide_colorbar(reverse = TRUE)) +
  guides(size = guide_bins(show.limits = TRUE)) +
  scale_size_continuous('Gene count', range = c(2, 10)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  theme(text = element_text(size = 10,
                            family = 'Arial Narrow')) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
          subtitle = 'EOMES, 4WPC, heart tissue') +
  theme_bw(base_size = 10) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm')
  ) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(size = 24))

## GATA 3
improved_data_wrangling(res_gata3_vs_conu_4wpc, 'gata3', '4wpc')

string_test_gata3 <-
  results_gata3_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <-
  bitr(string_test_gata3$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_gata3 <-
  string_test_gata3 %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

enrichment_gsea_gata3 <- enrichment_gata3$log2FC

names(enrichment_gsea_gata3) <- enrichment_gata3$ENTREZID

gsea_gata3 <- gseGO(
  geneList = enrichment_gsea_gata3,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

gsea_gata3 %>%
  fortify(., showCategory = 20) %>%
  as_tibble() %>%
  mutate_at('.sign', str_replace, 'suppressed', 'Downregulated') %>% mutate_at('.sign', str_replace, 'activated', 'Upregulated') %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_viridis_c('Adjusted p-value', guide = guide_colorbar(reverse = TRUE)) +
  guides(size = guide_bins(show.limits = TRUE)) +
  scale_size_continuous('Gene count', range = c(2, 10)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  theme(text = element_text(size = 10,
                            family = 'Arial Narrow')) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
          subtitle = 'GATA 3, 4WPC, heart tissue') +
  theme_bw(base_size = 10) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm')
  ) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(size = 24))


## IV-HD
improved_data_wrangling(res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')

string_test_ivhd <-
  results_ivhd_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <-
  bitr(string_test_eomes$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_ivhd <-
  string_test_ivhd %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

enrichment_gsea_ivhd <- enrichment_ivhd$log2FC

names(enrichment_gsea_ivhd) <- enrichment_ivhd$ENTREZID

gsea_ivhd <- gseGO(
  geneList = enrichment_gsea_ivhd,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

gsea_ivhd %>%
  fortify(., showCategory = 20) %>%
  as_tibble() %>%
  mutate_at('.sign', str_replace, 'suppressed', 'Downregulated') %>% mutate_at('.sign', str_replace, 'activated', 'Upregulated') %>%
  mutate(Description = fct_reorder(Description, Count)) %>%
  ggplot(aes(Count, Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_viridis_c('Adjusted p-value', guide = guide_colorbar(reverse = TRUE)) +
  guides(size = guide_bins(show.limits = TRUE)) +
  scale_size_continuous('Gene count', range = c(2, 10)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  theme(text = element_text(size = 10,
                            family = 'Arial Narrow')) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
          subtitle = 'IV-HD, 4WPC, heart tissue') +
  theme_bw(base_size = 10) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm')
  ) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(size = 24))


## cnetplots with GSEA output data ----

# List of treatments
treatments <- c("dnavaccine", "eomes", "gata3", "ivhd", "ivld")

# Initialize an empty list to store data frames
result_list <- list()

# Loop through each treatment and apply the function, also adding a column with treatment information
for (treatment in treatments) {
  obj <- paste0("res_", treatment, "_vs_conu_4wpc")
  result_list[[treatment]] <-
    improved_data_wrangling(get(obj), treatment = treatment, sampling_point = "4wpc")
}

rm(i,
   obj,
   objects,
   ora_down,
   ora_up,
   results_files,
   treatment,
   treatments)

# Vertically merge data frames from all treatments and remove NA
merged_df <- do.call(rbind, result_list) %>% na.omit()

# Reduce and intersect the ID column in result_list to find common differentially regulated genes among all treatments
reduced_intersected_ids <-
  Reduce(intersect, lapply(result_list, function(df)
    df$ID))

# Print the reduced and intersected IDs
print(reduced_intersected_ids)

# Find the common IDs between reduced_intersected_ids and merged_df. This way we are keeping only the common regulated salmon genes
common_ids <- intersect(merged_df$ID, reduced_intersected_ids)

# Subset merged_df to keep only rows with common IDs
# merged_subset <- merged_df[merged_df$ID %in% common_ids,]
merged_subset <-
  subset(merged_df, ID %in% common_ids)  # tidy version

merged_df[duplicated(merged_df$ortholog_name),] %>% arrange(ortholog_name) %>% print(n = 100)  # identifying duplicates in merged_df

# Splitting between up and downregulated genes, and sorting by gene name and adjusted p-value
upregulated_gsea <-
  merged_subset %>% dplyr::filter(log2FC > 0) %>% arrange(ortholog_name, adjusted_p.val)

downregulated_gsea <-
  merged_subset %>% dplyr::filter(log2FC < 0) %>% arrange(ortholog_name, adjusted_p.val)


# Keeping just one of the ortholog copies, the one with the lowest adjusted p-value
unique_upregulated <-
  upregulated_gsea[!duplicated(upregulated_gsea$ortholog_ensg),] %>% arrange(ortholog_name)
unique_downregulated <-
  downregulated_gsea[!duplicated(downregulated_gsea$ortholog_ensg),] %>% arrange(ortholog_name)

# Merging down and upregulated dataframes after removing duplicates
enrichment_unique_common <-
  rbind(unique_downregulated, unique_upregulated) %>% dplyr::arrange(desc(log2FC))

# Arranging matrix for GSEA
enrichment_gsea_common <- enrichment_unique_common$log2FC

names(enrichment_gsea_common) <-
  enrichment_unique_common$ortholog_ensg

length(enrichment_gsea_common)
# Running GSEA on common regulated genes
gsea_common <- gseGO(
  geneList = enrichment_gsea_common,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENSEMBL',
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  eps = 0,
  verbose = T
)

length(enrichment_gsea_common)

as_tibble(gsea_common)

gsea_common %>%
  fortify(., showCategory = 20) %>%
  as_tibble() %>%
  mutate_at('.sign', str_replace, 'suppressed', 'Downregulated') %>% mutate_at('.sign', str_replace, 'activated', 'Upregulated') %>% 
  mutate(Description = fct_reorder(Description, p.adjust)) %>% 
  ggplot(aes(Count, Description)) +
  geom_point(aes(color = p.adjust, size = Count)) +
  scale_color_viridis_c('Adjusted p-value', guide = guide_colorbar(reverse = TRUE)) +
  guides(size = guide_bins(show.limits = TRUE)) +
  scale_size_continuous('Gene count', range = c(2, 10)) +
  scale_y_discrete() +
  xlab('Gene count') +
  ylab(NULL) +
  theme(text = element_text(size = 10,
                            family = 'Arial Narrow')) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
          subtitle = 'Common genes, 4WPC, heart tissue') +
  theme_bw(base_size = 10) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm')
  ) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(size = 24))

edox <- setReadable(gsea_common, 'org.Hs.eg.db', 'ENSEMBL')


cnetplot(
  edox,
  color.params = list(foldChange = enrichment_gsea_common),
  circular = T,
  colorEdge = TRUE,
  cex_category = 1,
  cex_gene = 1,
  cex_label_category = 0.8,
  cex_label_gene = 0.8
)

## cnetplots at 4WPC ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

## DNA vaccine
orth_hs <- gorth(
  query = rownames(res_dnavaccine_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_4wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

head(res)

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_dnavaccine_4wpc <- as_tibble(gse)
gse_tibble_dnavaccine_4wpc %>% arrange(., NES) %>% slice_head(n = 4) %>% pull(Description)

# edox <- setReadable(gse, org.Hs.eg.db, 'auto')

## use new way of specifying visualization options
color.params = list(foldChange = gene_list, edge = TRUE)
cex.params = list(category_label = 0.6, gene_label = 0.6)

top_downregulated_dnavaccine_4wpc <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('neutrophil activation', 'B cell receptor signaling pathway', 'interleukin-6 production'),
  max.overlaps =  500
)

top_downregulated_dnavaccine_4wpc + ggtitle("Top 3 downregulated immune pathways\nDNA vaccine, 4WPC") + theme(plot.title = element_text(hjust= 0.5)) + theme(text = element_text(family = "Times New Roman"))


downregulated_dnavaccine_4wpc <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('immune system process', 'regulation of immune system process', 'immune response'),
  max.overlaps =  500
)

library(AnnotationDbi)
# Get the GO terms related to "immune system process" - NOT WORKING
immune_go_terms <- AnnotationDbi::select(org.Hs.eg.db, 
                                         keys = 'GO:0002376', 
                                         keytype = "GO", 
                                         columns = c("GO", 'ONTOLOGYALL'))
?AnnotationDbi::select


BiocManager::install("GO.db")
library(GO.db)
# Get all child terms of GO:0002376 (immune system process)
# Load the GOBPOFFSPRING dataset
as.list(GOBPOFFSPRING)

# Get the child terms for GO:0002376 (immune system process)
child_terms <- GOBPOFFSPRING[["GO:0002376"]]

# Get details of the child terms
child_terms_details <- AnnotationDbi::select(GO.db, 
                                             keys = child_terms, 
                                             keytype = "GOID", 
                                             columns = c("GOID", "TERM", "ONTOLOGY", 'DEFINITION')) %>% as_tibble()

# View the child terms and their details
as_tibble(child_terms_details)


columns(GO.db)
# View the child terms and their details
print(child_terms_details)

keytypes(org.Hs.eg.db)
# View the GO terms
print(immune_go_terms)


## GATA3
orth_hs <- gorth(
  query = rownames(res_gata3_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_gata3_vs_conu_4wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_gata3_4wpc <- as_tibble(gse)

gse_tibble_gata3_4wpc %>% arrange(., NES) %>% slice_head(n = 4) %>% pull(Description)


# edox <- setReadable(gse, org.Hs.eg.db, 'auto')

top_downregulated_gata3_4wpc <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('humoral immune response', 'immunoglobulin mediated immune response', 'B cell mediated immunity'),
  max.overlaps =  1000
)

top_downregulated_gata3_4wpc + ggtitle("Top 3 downregulated immune pathways\nGATA3, 4WPC") + theme(plot.title = element_text(hjust= 0.5)) + theme(text = element_text(family = "Times New Roman"))


## EOMES
orth_hs <- gorth(
  query = rownames(res_eomes_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_eomes_vs_conu_4wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_eomes_4wpc <- as_tibble(gse)
gse_tibble_eomes_4wpc %>% arrange(., NES) %>% slice_head(n = 4) %>% pull(Description)


top_downregulated_eomes_4wpc <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('humoral immune response mediated by circulating immunoglobulin', 'complement activation, classical pathway', 'immunoglobulin mediated immune response'),
  max.overlaps =  1000
)

top_downregulated_eomes_4wpc + ggtitle("Top downregulated immune pathways\nEOMES 4WPC") + theme(plot.title = element_text(hjust= 0.5)) + theme(text = element_text(family = "Times New Roman"))

## IV-HD
orth_hs <- gorth(
  query = rownames(res_ivhd_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_ivhd_vs_conu_4wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_ivhd <- as_tibble(gse)

# edox <- setReadable(gse, org.Hs.eg.db, 'auto')

top_downregulated_ivhd <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('adaptive immune response', 'B cell mediated immunity', 'lymphocyte mediated immunity'),
  max.overlaps =  500
)

top_downregulated_ivhd + ggtitle("Top downregulated IV-HD, 4WPC") + theme(plot.title = element_text(hjust= 0.5))


## IV-LD
orth_hs <- gorth(
  query = rownames(res_ivld_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_ivld_vs_conu_4wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_ivld <- as_tibble(gse)

# edox <- setReadable(gse, org.Hs.eg.db, 'auto')

top_downregulated_ivld <- cnetplot(
  gse,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('regulation of antigen processing and presentation', 'positive regulation of interferon-alpha production', 'positive regulation of CD4-positive, alpha-beta T cell activation'),
  max.overlaps =  500
)

top_downregulated_ivhd + ggtitle("Top downregulated IV-LD, 4WPC") + theme(plot.title = element_text(hjust= 0.5))

## cnetplots at 1WPC ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_1wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

### DNA vaccine ----
orth_hs <- gorth(
  query = rownames(res_dnavaccine_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_1wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

head(res)

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse_dnavaccine <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

gse_tibble_dnavaccine_1wpc <- as_tibble(gse_dnavaccine)

## use new way of specifying visualization options
color.params = list(foldChange = gene_list, edge = TRUE)
cex.params = list(category_label = 0.6, gene_label = 0.6)

top_downregulated_dnavaccine_1wpc <- cnetplot(
  gse_dnavaccine,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('detection of external stimulus', 'cellular response to external stimulus','activation of immune response'),
  max.overlaps =  500
)

top_downregulated_dnavaccine_1wpc + ggtitle("Top downregulated DNA vaccine, 1WPC") +theme(plot.title = element_text(hjust= 0.5))

### GATA3 ----
orth_hs <- gorth(
  query = rownames(res_gata3_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)  # converting ensembl gene IDs to human orthologs

res <- tibble::rownames_to_column(as.data.frame(res_gata3_vs_conu_1wpc), var = 'ensembl')  # setting rownames (gene names) as a variable column

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()  # joining gorth output with results table

res <- res[order(-res$log2FoldChange),]  # ordering by fold-change

gene_list <- res$log2FoldChange  # creating matrix with fold-change values
names(gene_list) <- res$ortholog_name  # adding gene names to fold-change

gse_gata3 <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)  # running gsea

gse_tibble_gata3_1wpc <- as_tibble(gse_gata3)  # converting enrichResult table to tibble


top_downregulated_gata3_1wpc <- cnetplot(
  gse_gata3,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('positive regulation of response to external stimulus', 'negative regulation of response to external stimulus', 'viral genome replication'),
  max.overlaps =  500
)

top_downregulated_gata3_1wpc + ggtitle("Top downregulated GATA3, 1WPC") + theme(plot.title = element_text(hjust= 0.5))


### EOMES ----
orth_hs <- gorth(
  query = rownames(res_eomes_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_eomes_vs_conu_1wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse_eomes <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

# gseaplot(gse, geneSetID = 4)
gse_tibble_eomes <- as_tibble(gse_eomes)

# edox <- setReadable(gse, org.Hs.eg.db, 'auto')

top_downregulated_eomes_1wpc <- cnetplot(
  gse_eomes,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('positive regulation of response to external stimulus', 'negative regulation of response to external stimulus', 'regulation of response to external stimulus'),
  max.overlaps =  500
)

top_downregulated_eomes_1wpc + ggtitle("Top downregulated EOMES, 1WPC") + theme(plot.title = element_text(hjust= 0.5))

## IV-HD ----
orth_hs <- gorth(
  query = rownames(res_ivhd_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_ivhd_vs_conu_1wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse_ivhd <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

gse_tibble_ivhd <- as_tibble(gse_ivhd)

top_downregulated_ivhd_1wpc <- cnetplot(
  gse_ivhd,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('response to extracellular stimulus', 'detection of stimulus', 'cellular response to abiotic stimulus'),
  max.overlaps =  500
)

top_downregulated_ivhd_1wpc + ggtitle("Top downregulated IV-HD, 1WPC") + theme(plot.title = element_text(hjust= 0.5))

## IV-LD ----
orth_hs <- gorth(
  query = rownames(res_ivld_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

res <- tibble::rownames_to_column(as.data.frame(res_ivld_vs_conu_1wpc), var = 'ensembl')

res <- res %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>%
  dplyr::select(.,
                ensembl,
                ortholog_name,
                ortholog_ensg,
                log2FoldChange,
                padj,
                ortholog_ensg,
                description) %>% 
  na.omit()

res <- res[order(-res$log2FoldChange),]

gene_list <- res$log2FoldChange
names(gene_list) <- res$ortholog_name

gse_ivld <- gseGO(gene_list,
             keyType = 'SYMBOL',
             OrgDb = org.Hs.eg.db,
             eps = 1e-300)

gse_tibble_ivld <- as_tibble(gse_ivld)

top_downregulated_ivld <- cnetplot(
  gse_ivld,
  color.params = color.params,
  cex.params = cex.params,
  circular = T,
  colorEdge = TRUE,
  showCategory = c('detection of external stimulus', 'negative regulation of response to external stimulus', 'positive regulation of response to external stimulus'),
  max.overlaps =  500
)

top_downregulated_ivld + ggtitle("Top downregulated IV-LD, 1WPC") + theme(plot.title = element_text(hjust= 0.5))



### Histoscore ----
histoscore_heart_4wpc <- read_delim('histoscore_hear_4wpc.tsv', col_names = c('treatment', 'histoscore'), delim = '\t', show_col_types = F)


histoscore_heart_4wpc %>% 
  group_by(treatment) %>% 
  summarise(
    mean = mean(histoscore, na.rm = T),
    sd = sd(histoscore, na.rm = T),
    n = n(),
    se = sd/sqrt(n)
  ) %>% 
  ggplot(., aes(x = treatment, y = mean)) +
  geom_boxplot()

install.packages('viridis')
install.packages('hrbrthemes')
library(hrbrthemes)
library(viridis)

histoscore_heart_4wpc %>% na.omit() %>%
  ggplot(., aes(x = factor(treatment, level = c('CONU', 'IVLD', 'IVHD', 'DNAvaccine', 'EOMES', 'GATA3')), y = histoscore, fill = treatment)) +
  geom_boxplot() +
  theme_ipsum() +
  theme(
    legend.position = 'none',
    plot.title = element_text(size = 11)
  ) +
  scale_fill_viridis(discrete = T, alpha = 0.6, option = 'B') +
  ggtitle('Histoscore at heart 4WPC') +
  ylab('') + xlab('')
  



### Retrieving GO annotations per gene SYMBOL ----

# GO terms associated with "GATA3"
gene_symbol_gata3 <- "GATA3"

# Get the GO annotations for the gene
go_terms_gata3 <- select(org.Hs.eg.db, keys = gene_symbol_gata3, keytype = "SYMBOL", columns = c("GO"),
                   keytypes = "SYMBOL")

# Print the retrieved GO annotations
print(go_terms)

gse_tibble_gata3_1wpc %>% left_join(., go_terms_gata3, by = c('ID' = 'GO')) %>% na.omit() %>% arrange(NES)

head(go_terms_gata3)
head(gse_tibble_gata3_1wpc)

gene_symbol_il4 <- 'IL4'
go_terms_il4 <- select(org.Hs.eg.db, keys = gene_symbol_il4, keytype = "SYMBOL", columns = c("GO"),
                       keytypes = "SYMBOL")

gse_tibble_gata3_1wpc %>% left_join(., go_terms_il4, by = c('ID' = 'GO')) %>% na.omit()

head(go_terms_il4)

gene_symbol_il13 <- 'IL13'
go_terms_il13 <- select(org.Hs.eg.db, keys = gene_symbol_il13, keytype = "SYMBOL", columns = c("GO"),
                       keytypes = "SYMBOL")

gse_tibble_gata3_1wpc %>% left_join(., go_terms_il13, by = c('ID' = 'GO')) %>% na.omit()


## Looking for STAT6 in the results table for GATA3
rownames_to_column(as.data.frame(res_gata3_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000080589'))

## Looking for RSAD2 in the results table for all treatments. It is present in all of them.
rownames_to_column(as.data.frame(res_eomes_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000048046'))
rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000048046'))
rownames_to_column(as.data.frame(res_gata3_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000048046'))
rownames_to_column(as.data.frame(res_ivhd_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000048046'))
rownames_to_column(as.data.frame(res_ivld_vs_conu_1wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000048046'))

## go terms associated with RSAD2
gene_symbol_rsad2 <- 'RSAD2'
go_terms_rsad2 <- select(org.Hs.eg.db, keys = gene_symbol_rsad2, keytype = "SYMBOL", columns = c("GO"),
                        keytypes = "SYMBOL")

gse_tibble_gata3_1wpc %>% left_join(., go_terms_rsad2, by = c('ID' = 'GO')) %>% na.omit()
gse_tibble_eomes %>% left_join(., go_terms_rsad2, by = c('ID' = 'GO')) %>% na.omit()
gse_tibble_ivhd %>% left_join(., go_terms_rsad2, by = c('ID' = 'GO')) %>% na.omit()
gse_tibble_ivld %>% left_join(., go_terms_rsad2, by = c('ID' = 'GO')) %>% na.omit()
gse_tibble_dnavaccine_1wpc %>% left_join(., go_terms_rsad2, by = c('ID' = 'GO')) %>% na.omit()

## goterms associated with MHC I (HLA in humans)
gene_symbol_hlaa <- 'HLA-A'
go_terms_hlaa <- select(org.Hs.eg.db, keys = gene_symbol_hlaa, keytype = "SYMBOL", columns = c("GO"),
                       keytypes = "SYMBOL")

gse_tibble_dnavaccine_1wpc %>% left_join(., go_terms_hlaa, by = c('ID' = 'GO')) %>% na.omit()

gene_symbol_ifng <- 'IFNG'

go_terms_ifng <- select(org.Hs.eg.db, keys = gene_symbol_ifng, keytype = "SYMBOL", columns = c("GO"),
                        keytypes = "SYMBOL")


gse_tibble_dnavaccine_1wpc %>% left_join(., go_terms_ifng, by = c('ID' = 'GO')) %>% na.omit() %>% pull(Description)

## At 4wpc, all the immune related pathways have negative NES
### A negative NES indicates that the gene set is over represented among the genes 
### that are negatively correlated with the phenotype 
### (e.g., highly expressed in the opposite condition or less expressed in the condition of interest).
gse_tibble_dnavaccine_4wpc %>% left_join(., go_terms_ifng, by = c('ID' = 'GO')) %>% na.omit() %>% pull(Description)

gene_symbol_ifna <- 'IFNA1'

go_terms_ifna <- select(org.Hs.eg.db, keys = gene_symbol_ifna, keytype = "SYMBOL", columns = c("GO"),
                        keytypes = "SYMBOL")


gse_tibble_dnavaccine_4wpc %>% left_join(., go_terms_ifna, by = c('ID' = 'GO')) %>% na.omit() %>% pull(Description)
gse_tibble_dnavaccine_1wpc %>% left_join(., go_terms_ifna, by = c('ID' = 'GO')) %>% na.omit() %>% pull(Description)

gse_tibble_dnavaccine_4wpc 


head(results_dnavaccine_4wpc)

as_tibble(res_dnavaccine_vs_conu_4wpc) %>% filter(str_detect())

### Looking for MHC II
rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_4wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000004635'))

### Looking for ACTB
rownames_to_column(as.data.frame(res_dnavaccine_vs_conu_4wpc), var = 'ensembl') %>% as_tibble() %>% filter(str_detect(ensembl, 'ENSSSAG00000116649'))




## Testing cnetplot on exclusive genes ----

# List of objects
objects <- ls(pattern = "^res_.*_conu_4wpc")

# List of treatments
treatments <- c("dnavaccine", "eomes", "eomes", "ivhd", "ivld")

# Loop through each treatment and apply the function
for (treatment in treatments) {
  obj <- paste0("res_", treatment, "_vs_conu_4wpc")
  improved_data_wrangling(get(obj), treatment = treatment, sampling_point = "4wpc")
}

### GATA 3
eomes_exclusives <- setdiff(
  results_eomes_4wpc$ortholog_ensg,
  c(
    results_dnavaccine_4wpc$ortholog_ensg,
    results_eomes_4wpc$ortholog_ensg,
    results_ivhd_4wpc$ortholog_ensg,
    results_ivld_4wpc$ortholog_ensg
  )
)

eomes_exclusives_all <-
  results_eomes_4wpc %>% filter(ortholog_ensg %in% eomes_exclusives) %>% na.omit() %>% arrange()

# Arranging matrix for GSEA
eomes_gsea_exclusive <- eomes_exclusives_all$log2FC

names(eomes_gsea_exclusive) <-
  eomes_exclusives_all$ortholog_ensg

# Running GSEA on common regulated genes
gsea_eomes <- gseGO(
  geneList = eomes_gsea_exclusive,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENSEMBL',
  ont = 'BP',
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

edox <- setReadable(gsea_common, 'org.Hs.eg.db', 'ENSEMBL')


gsea_common_cnet <-
  cnetplot(edox,
           foldChange = enrichment_gsea_common,
           circular = T,
           colorEdge = TRUE)


### EOMES
eomes_exclusives <- setdiff(
  results_eomes_4wpc$ortholog_ensg,
  c(
    results_dnavaccine_4wpc$ortholog_ensg,
    results_gata3_4wpc$ortholog_ensg,
    results_ivhd_4wpc$ortholog_ensg,
    results_ivld_4wpc$ortholog_ensg
  )
)

eomes_exclusives_all <-
  results_eomes_4wpc %>% filter(ortholog_ensg %in% eomes_exclusives) %>% na.omit() %>% arrange()

# Arranging matrix for GSEA
eomes_gsea_exclusive <- eomes_exclusives_all$log2FC

names(eomes_gsea_exclusive) <-
  eomes_exclusives_all$ortholog_ensg

length(eomes_gsea_exclusive)

# Running GSEA on common regulated genes
gsea_eomes <- gseGO(
  geneList = eomes_gsea_exclusive,
  OrgDb = org.Hs.eg.db,
  keyType = 'ENSEMBL',
  ont = 'BP',
  minGSSize = 10,
  maxGSSize = 10000,
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

edox <- setReadable(gsea_common, 'org.Hs.eg.db', 'ENSEMBL')

gsea_common_cnet <-
  cnetplot(edox,
           foldChange = enrichment_gsea_common,
           circular = T,
           colorEdge = TRUE)
