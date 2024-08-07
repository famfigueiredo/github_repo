## NOT UPDATED WITH NEW LOG2FC AND P-VALUE ##

# Over Representation Analysis, 6WPC vs 10WPI ----
## Loading results ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions. It uses setdiff to find the subset of objects in the global environment
# (as returned by ls()) that don't have mode function (as returned by lsf.str())

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi'
)

results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

immune_related_GOterms <- read.table('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/immune_and_viral_related_GOterms.tsv', 
                                     header = F, sep = '\t')  # loading dataframe containing immune related GO terms

colnames(immune_related_GOterms) <- c('goterms', 'ontology', 'description')  # renaming columns

immune_related_GOterms <- immune_related_GOterms[!duplicated(immune_related_GOterms$goterms),]  # removing duplicated GO terms

nrow(immune_related_GOterms)
## Loading functions ----

# Significant genes grabs a DESeq2 result table, subsets the genes with padj < 0, and selects the ID, log2FoldChange, and padj columns. Also arranges log2FC in a descending manner.
significant_genes <- function(results_files) {
  b <- as.data.frame(subset(results_files, padj < 0.05)) %>%
    rownames_to_column(var = 'ID') %>%
    as_tibble()
  
  sig_genes <- b %>%
    dplyr::select(ID, log2FC = log2FoldChange, adjusted_p.val = padj) %>%
    dplyr::arrange(desc(log2FC))
  
  return(sig_genes)
}

# improved data wrangling function selects significantly regulated genes, converts ssalar IDs to hsapiens orthologs, and joins both
improved_data_wrangling <- function(results_table, treatment, sampling_point) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  
  if (!requireNamespace("gprofiler2", quietly = TRUE)) {
    install.packages("gprofiler2")
  }
  
  # Load the required packages
  library(dplyr)
  library(gprofiler2)
  
  # select significant genes from DESeq2 results table
  a <-
    significant_genes(results_table)
  
  # convert ssalar gene IDs to human orthologs
  orth_hs <- gorth(
    query = a$ID,
    source_organism = 'ssalar',
    target_organism = 'hsapiens',
    mthreshold = 1,
    filter_na = T
  )
  
  # join significant genes table with human ortholog names
  results <- a %>% left_join(orth_hs, by = c('ID' = 'input')) %>%
    dplyr::select(.,
                  ID,
                  ortholog_name,
                  log2FC,
                  adjusted_p.val,
                  ortholog_ensg,
                  description)
  
  # getting gene lists for ORA
  ora_up <<- 
    results %>% drop_na() %>% dplyr::filter(log2FC > 0) %>%  pull(ortholog_ensg)
  
  ora_down <<-
    results %>% drop_na() %>% dplyr::filter(log2FC < 0) %>%  pull(ortholog_ensg)
  
  # create results name
  results_name <-
    paste('results', treatment, sampling_point, sep = '_')
  
  # assign the results data frame to the dynamic name
  assign(results_name, results, envir = .GlobalEnv)
  
  print(get(results_name))
}

# combine and label pathways merges both downregulated and upregulated dataframes into one
combine_and_label_pathways <- function(downreg_tibble, upreg_tibble) {
  # Add a "regulation" column to downregulated tibble
  df_downregulated <- as_tibble(downreg_tibble) %>%
    mutate(regulation = "downregulated")
  
  # Add a "regulation" column to upregulated tibble
  df_upregulated <- as_tibble(upreg_tibble) %>%
    mutate(regulation = "upregulated")
  
  # Concatenate the two tibbles
  final_tibble <- bind_rows(df_downregulated, df_upregulated)
  
  return(final_tibble)
}

# filter rows by GO term intersects the above merged dataframe with the immune related GO terms and returns a dataframe containing those
filter_rows_by_GO_term <- function(df1, df2, id_column_name) {
  # Perform the intersection
  common_GO_terms <- intersect(df1$ID, df2[[id_column_name]])
  
  # Filter rows based on common GO terms
  filtered_rows <- df1 %>% 
    filter(ID %in% common_GO_terms)
  
  return(filtered_rows)
}

## CONU ----
improved_data_wrangling(res_conu_6wpc_vs_10wpi, 'conu', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
conu_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(conu_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/conu_pathways.xlsx', colNames = T, asTable = F)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 3))  %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 4) and upregulated (level 3) DEGs (human orthologs) \
    CONU @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))


## pTagRFP ----
improved_data_wrangling(res_ptagrfp_6wpc_vs_10wpi, 'ptagrfp', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )


egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ptagrfp_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(ptagrfp_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/ptagrfp_pathways.xlsx', colNames = T, asTable = F)


downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))


upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 3) DEGs (human orthologs), first 15 terms\
    pTagRFP @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

## IV-LD ----
improved_data_wrangling(res_ivld_6wpc_vs_10wpi, 'ivld', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ivld_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(ivld_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/ivld_pathways.xlsx', colNames = T, asTable = F)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 4) DEGs (human orthologs)\
    IV-LD @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))



## IV-HD ----
improved_data_wrangling(res_ivhd_6wpc_vs_10wpi, 'ivhd', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ivhd_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(ivhd_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/ivhd_pathways.xlsx', colNames = T, asTable = F)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 4) DEGs (human orthologs)\
    IV-HD @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

## DNA vaccine ----
improved_data_wrangling(res_dnavaccine_6wpc_vs_10wpi, 'dnavaccine', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
dnavaccine_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(dnavaccine_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/dnavaccine_pathways.xlsx', colNames = T, asTable = F)


downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 4) DEGs (human orthologs)\
    DNA vaccine @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))


## EOMES ----
improved_data_wrangling(res_eomes_6wpc_vs_10wpi, 'eomes', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
eomes_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(eomes_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/eomes_pathways.xlsx', colNames = T, asTable = F) 

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 4) DEGs (human orthologs)\
    EOMES @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))


## GATA3 ----
improved_data_wrangling(res_gata3_6wpc_vs_10wpi, 'gata3', '6wpc_vs_10wpi')

## ORA ##
egoORA_downreg <-
  enrichGO(
    gene = ora_down,
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    minGSSize = 10,
    maxGSSize = 1000,
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

egoORA_upreg <-
  enrichGO(
    gene = ora_up,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <- combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
gata3_regulated_pathways <- filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(eomes_regulated_pathways, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/gata3_pathways.xlsx', colNames = T, asTable = F)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           size = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))


upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    size = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (level 3) and upregulated (level 5) DEGs (human orthologs)\
    GATA3 @ 6 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 8)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 6),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

# Summary ----
conu_regulated_pathways <- as_tibble(conu_regulated_pathways) %>% 
  mutate(treatment = 'conu')

ptagrfp_regulated_pathways <- as_tibble(ptagrfp_regulated_pathways) %>% 
  mutate(treatment = 'ptagrfp')

ivld_regulated_pathways <- as_tibble(ivld_regulated_pathways) %>% 
  mutate(treatment = 'ivld')

ivhd_regulated_pathways <- as_tibble(ivhd_regulated_pathways) %>% 
  mutate(treatment = 'ivhd')

eomes_regulated_pathways <- as_tibble(eomes_regulated_pathways) %>% 
  mutate(treatment = 'eomes')

gata3_regulated_pathways <- as_tibble(gata3_regulated_pathways) %>% 
  mutate(treatment = 'gata3')

dnavaccine_regulated_pathways <- as_tibble(dnavaccine_regulated_pathways) %>% 
  mutate(treatment = 'dnavaccine')


regulated_pathways_6wpc_10wpi <- bind_rows(mget(ls(pattern = '_regulated_pathways')))


write.xlsx(regulated_pathways_6wpc_10wpi, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi/regulated_pathways_6wpc_10wpi.xlsx', colNames = T, asTable = F)






