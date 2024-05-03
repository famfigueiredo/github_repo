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

# DNA vaccine ----

improved_data_wrangling(res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')

string_test <- results_dnavaccine_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <- bitr(string_test$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment <- string_test %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

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
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, 'cm'),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm')) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(
    size = 24
  ))

# EOMES ----  
improved_data_wrangling(res_eomes_vs_conu_4wpc, 'eomes', '4wpc')

string_test_eomes <- results_eomes_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <- bitr(string_test_eomes$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_eomes <- string_test_eomes %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

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
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, 'cm'),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm')) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(
    size = 24
  ))

# GATA 3 ----  
improved_data_wrangling(res_gata3_vs_conu_4wpc, 'gata3', '4wpc')

string_test_gata3 <- results_gata3_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <- bitr(string_test_eomes$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_gata3 <- string_test_gata3 %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

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
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, 'cm'),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm')) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(
    size = 24
  ))


# IV-HD ----  
improved_data_wrangling(res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')

string_test_ivhd <- results_ivhd_4wpc %>% dplyr::select(ortholog_name, log2FC) %>% na.omit()

entrez_ids <- bitr(string_test_eomes$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)

enrichment_ivhd <- string_test_ivhd %>% left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>% dplyr::select(ENTREZID, log2FC)

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
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        legend.key.size = unit(1, 'cm'),
        legend.position = 'right',
        legend.key.height = unit(1, 'cm')) +
  facet_grid(. ~ .sign) +
  theme(strip.text = element_text(
    size = 24
  ))

read.gmt('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/cytoscape/Supplementary_Table5_hsapiens.pathways.NAME.gmt')

library(DOSE)
data(geneList)
de <- names(geneList)[abs(geneList) > 2]

edo <- enrichDGN(de)

mutate(edo, qscore = -log(p.adjust, base=10)) %>% 
  barplot(x="qscore")

edo2 <- gseDO(geneList)
dotplot(edo, showCategory=30) + ggtitle("dotplot for ORA")
dotplot(edo2, showCategory=30) + ggtitle("dotplot for GSEA")



## convert gene ID to Symbol
edox <- setReadable(edo, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(edox, foldChange=geneList)
## categorySize can be scaled by 'pvalue' or 'geneNum'
p2 <- cnetplot(edox, categorySize="pvalue", foldChange=geneList)
gsea_ivhd_cnet <- cnetplot(edox, foldChange=enrichment_gsea_ivhd, circular = F, colorEdge = TRUE, layout = 'star') 
cowplot::plot_grid(p1, p2, p3, ncol=3, labels=LETTERS[1:3], rel_widths=c(.8, .8, 1.2))
?cnetplot


edox <- setReadable(gsea_ivhd, 'org.Hs.eg.db', 'ENTREZID')
cnetplot()


gsea_ivhd@result$Description


