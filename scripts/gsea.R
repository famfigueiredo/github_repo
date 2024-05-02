suppressPackageStartupMessages({
  library('tidyverse')
  library('VennDiagram')
  library('clusterProfiler')
  library('gprofiler2')
  library('org.Hs.eg.db')
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
  list.files(pattern = 'res_.*')  # regex matching results files
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
  assign(paste0("significant_", res_object), significant_result, envir = .GlobalEnv)
}

rm(list.data, result, i, res_object, res_objects, results_files, significant_result)
rm(list = ls(pattern = '^res_.'))  # removing results files from Global


improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

string_test <- results_dnavaccine_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

write.table(string_test, 'string_test.tsv', sep = '\t', row.names = F, quote = F)

# select significant genes from DESeq2 results table
a <-
  significant_genes(res_dnavaccine_vs_conu_4wpc)

# convert ssalar gene IDs to human orthologs
orth_hs <- gorth(
  query = a$ID,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)




sig_genes <-
  b %>% left_join(gene_id, by = c('id' = 'entrezgene_accession')) %>%
  mutate(complete_id = if_else(is.na(entrezgene_id), id, entrezgene_id)) %>%
  mutate(complete_id = gsub('LOC', '', complete_id)) %>%
  dplyr::select(.,
                complete_id,
                log2FoldChange,
                padj) %>%
  dplyr::arrange(., desc(log2FoldChange))

enrichment_gsea <<-
  sig_genes$log2FoldChange  # converting to matrix for GSEA

names(enrichment_gsea) <<-
  sig_genes$complete_id  # adding gene ids for GSEA


string_test <- results_dnavaccine_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

enrichment_gsea <- string_test$log2FC
names(enrichment_gsea) <- string_test$ortholog_name

gsea_DNAvacc <- gseGO(
  geneList = enrichment_gsea,
  OrgDb = org.Hs.eg.db,
  ont = 'BP',
  pvalueCutoff = 0.05,
  pAdjustMethod = 'BH',
  verbose = T
)

entrez_ids <- bitr(string_test$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)




enrichment <- string_test %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

enrichment_gsea <- enrichment$log2FC
names(enrichment_gsea) <- enrichment$ENTREZID

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
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, 'cm'),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm')) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(
    size = 24
  ))



