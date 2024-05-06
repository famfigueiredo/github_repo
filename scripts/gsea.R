suppressPackageStartupMessages({
  library('tidyverse')
  library('VennDiagram')
  library('clusterProfiler')
  library('gprofiler2')
  library('org.Hs.eg.db')
  library('enrichplot')
})


## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

immune_related_GOterms <-
  read.table(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/goTerm_lists/immune_related_GOterms.tsv',
    header = T,
    sep = '\t'
  )  # loading dataframe containing immune related GO terms

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_4wpc'
)

results_files <-
  list.files(pattern = '^res_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}


# Get names of all objects in global environment with a regex pattern matching 'res_.*_conu_4wpc'
res_objects <- ls(pattern = "res_.*_4wpc")

# Apply significant_genes function to each res object
for (res_object in res_objects) {
  # Get the object from the global environment
  result <- get(res_object)
  
  # Apply significant_genes function
  significant_result <- significant_genes(result)
  
  # Assign the result back to the global environment
  assign(paste0("significant_", res_object),
         significant_result,
         envir = .GlobalEnv)
}

rm(list.data,
   result,
   i,
   res_object,
   res_objects,
   results_files,
   significant_result)
rm(list = ls(pattern = '^res_.'))  # removing results files from Global

# DNA vaccine ----

improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

string_test <-
  results_dnavaccine_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <-
  bitr(string_test$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment <-
  string_test %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

enrichment_gsea <- enrichment$log2FC

names(enrichment_gsea) <- enrichment$ENTREZID

gsea_DNAvacc <- gseGO(
  geneList = enrichment_gsea,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

gsea_DNAvacc %>%
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
  theme(text = element_text(size = 30,
                            family = 'Arial Narrow')) +
  ggtitle('GSEA, downregulated vs upregulated DEGs',
          subtitle = 'DNA vaccine, 4WPC, heart tissue') +
  theme_bw(base_size = 14) +
  theme(
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(1, 'cm')
  ) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(size = 24))

# EOMES ----
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

# GATA 3 ----
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


# IV-HD ----
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




## Testing cnetplots with GSEA output data ----

# List of objects containing all results files matching the regex pattern
objects <- ls(pattern = "^res_.*_conu_4wpc")

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

## Testing cnet plot on exclusive genes ----

# List of objects
objects <- ls(pattern = "^res_.*_conu_4wpc")

# List of treatments
treatments <- c("dnavaccine", "eomes", "eomes", "ivhd", "ivld")

# Loop through each treatment and apply the function
for (treatment in treatments) {
  obj <- paste0("res_", treatment, "_vs_conu_4wpc")
  improved_data_wrangling(get(obj), treatment = treatment, sampling_point = "4wpc")
}

### GATA 3 ----
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


### EOMES ----
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
