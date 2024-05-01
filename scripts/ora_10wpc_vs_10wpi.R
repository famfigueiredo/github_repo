# Over Representation Analysis, 10WPC vs 10WPI ----
## Loading results ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi'
)

results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

immune_related_GOterms <- read.table('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/immune_related_GOterms.tsv', 
                                     header = F, sep = '\t')  # loading dataframe containing immune related GO terms

colnames(immune_related_GOterms) <- c('goterms', 'ontology', 'description')  # renaming columns

immune_related_GOterms <- immune_related_GOterms[!duplicated(immune_related_GOterms$goterms),]  # removing duplicated GO terms

nrow(immune_related_GOterms)
## Loading functions ----
source('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/functions_data-wrangling_march24.R')

## CONU ----
improved_data_wrangling(res_conu_10wpc_vs_10wpi, 'conu', '10wpc_vs_10wpi')

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

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
conu_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  conu_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_conu_pathways.xlsx',
  colNames = T,
  asTable = F
)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (lvl. 5) and upregulated (lvl. 3) DEGs (human orthologs) \
    CONU @ 10 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/CONU.png', plot = last_plot())

## pTagRFP ----
improved_data_wrangling(res_ptagrfp_10wpc_vs_10wpi, 'ptagrfp', '10wpc_vs_10wpi')

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
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ptagrfp_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  ptagrfp_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_ptagrfp_pathways.xlsx',
  colNames = T,
  asTable = F
)

# no enriched terms in downregulated
downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(egoORA_upreg) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = T
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  ggtitle('ORA upregulated (lvl. 5) DEGs (human orthologs) \
    pTagRFP @ 10 WPC vs 10 WPI, heart tissue') +
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/pTagRFP.png', plot = last_plot())

## IV-LD ----
improved_data_wrangling(res_ivld_10wpc_vs_10wpi, 'ivld', '10wpc_vs_10wpi')

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
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ivld_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  ivld_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_ivld_pathways.xlsx',
  colNames = T,
  asTable = F
)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (lvl. 4) and upregulated (lvl. 4) DEGs (human orthologs) \
    IV-LD @ 10 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/IV-LD.png', plot = last_plot())

 ## IV-HD ----
improved_data_wrangling(res_ivhd_10wpc_vs_10wpi, 'ivhd', '10wpc_vs_10wpi')

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
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
ivhd_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  ivhd_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_ivhd_pathways.xlsx',
  colNames = T,
  asTable = F
)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (lvl. 5) and upregulated (lvl. 4) DEGs (human orthologs)\
    IV-HD @ 10 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/IV-HD.png', plot = last_plot())

## DNA vaccine ----
improved_data_wrangling(res_dnavaccine_10wpc_vs_10wpi, 'dnavaccine', '10wpc_vs_10wpi')

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

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
dnavaccine_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  dnavaccine_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_dnavaccine_pathways.xlsx',
  colNames = T,
  asTable = F
)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()


regulation_ora <- (downreg / upreg)
regulation_ora +
  plot_layout(guides = 'collect') +
  plot_annotation(
    title = 'ORA downregulated (lvl. 3) and upregulated (4) DEGs (human orthologs)\
    DNA vaccine @ 10 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/DNA vaccine.png', plot = last_plot())

## EOMES ----
improved_data_wrangling(res_eomes_10wpc_vs_10wpi, 'eomes', '10wpc_vs_10wpi')

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

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
eomes_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  eomes_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_eomes_pathways.xlsx',
  colNames = T,
  asTable = F
)

downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()

upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
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
    title = 'ORA downregulated (lvl. 4) and upregulated (lvl. 5) DEGs (human orthologs)\
    EOMES @ 10 WPC vs 10 WPI, heart tissue',
    theme = theme(plot.title = element_text(hjust = .5, vjust = 2)),
    tag_levels = 'A'
  ) &
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'bottom',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/EOMES.png', plot = last_plot())

## GATA3 ----
improved_data_wrangling(res_gata3_10wpc_vs_10wpi, 'gata3', '10wpc_vs_10wpi')

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

regulated_pathways <-
  combine_and_label_pathways(egoORA_downreg, egoORA_upreg)
gata3_regulated_pathways <-
  filter_rows_by_GO_term(regulated_pathways, immune_related_GOterms, 'goterms')
write.xlsx(
  eomes_regulated_pathways,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_gata3_pathways.xlsx',
  colNames = T,
  asTable = F
)

# no enriched terms
downreg <-
  as_tibble(gofilter(egoORA_downreg, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous(name = NULL) +
  # geom_text(aes(label = Count), hjust = -.5, colour = 'black') +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic()
# ggtitle('downregulated') +
# theme(plot.title = element_text(hjust = 0.5, size = 8))


upreg <-
  as_tibble(gofilter(egoORA_upreg, level = 3)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = FALSE
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Oranges',
                       name = 'Adjusted \n p-value',
                       guide = guide_colorbar(reverse = TRUE)) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  ggtitle('ORA upregulated (lvl. 3) DEGs (human orthologs) \
    GATA3 @ 10 WPC vs 10 WPI, heart tissue') +
  theme(text = element_text(size = 10,
                            family = 'serif'),
        axis.text.y = element_text(size = 10)) &
  theme(
    plot.tag = element_blank(),
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 10),
    legend.key.size = unit(1, 'cm'),
    legend.position = 'right',
    legend.key.height = unit(0.5, 'cm')
  ) +
  theme(plot.margin = margin(0.5, 0.5, 0, 0, 'cm'))

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/human_orthologs/10wpc_vs_10wpi/GATA3.png', plot = last_plot())

# Summary ----
conu_regulated_pathways <- as_tibble(conu_regulated_pathways) %>%
  mutate(treatment = 'conu')

ptagrfp_regulated_pathways <-
  as_tibble(ptagrfp_regulated_pathways) %>%
  mutate(treatment = 'ptagrfp')

ivld_regulated_pathways <- as_tibble(ivld_regulated_pathways) %>%
  mutate(treatment = 'ivld')

ivhd_regulated_pathways <- as_tibble(ivhd_regulated_pathways) %>%
  mutate(treatment = 'ivhd')

eomes_regulated_pathways <- as_tibble(eomes_regulated_pathways) %>%
  mutate(treatment = 'eomes')

gata3_regulated_pathways <- as_tibble(gata3_regulated_pathways) %>%
  mutate(treatment = 'gata3')

dnavaccine_regulated_pathways <-
  as_tibble(dnavaccine_regulated_pathways) %>%
  mutate(treatment = 'dnavaccine')

regulated_pathways_10wpc_10wpi <-
  bind_rows(mget(ls(pattern = '_regulated_pathways')))

write.xlsx(
  regulated_pathways_10wpc_10wpi,
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi/2404_regulated_pathways_10wpc_10wpi.xlsx',
  colNames = T,
  asTable = F
)