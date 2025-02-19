load('~/Documents/PhD/Papers/Paper III/data/RData/sampleTables/sampleTable.RData')

sampleTable <- droplevels(sampleTable)
table(sampleTable$treatment, sampleTable$samplingPoint, sampleTable$tissue)
sampleTable_heart_10wpi <- sampleTable %>% filter(tissue == 'h' & treatment %in% c('conu', 'ivld', 'eomes', 'gata3') & samplingPoint == '10wpi')

table(sampleTable_heart_10wpi$treatment, sampleTable_heart_10wpi$samplingPoint, sampleTable_heart_10wpi$tissue)

sampleTable_heart_10wpi <- droplevels(sampleTable_heart_10wpi)

levels(sampleTable_heart_10wpi$treatment) <- c('conu', 'ivld', 'eomes', 'gata3')

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_heart_10wpi,
  directory = directory,
  design = ~ treatment
)

collapsed_dds_10wpi <- collapseReplicates(dds,
                                          groupby = dds$n,
                                          run = dds$lane)

keep <-
  rowSums(counts(collapsed_dds_10wpi)) >= 10  # removing low count genes (<10)
collapsed_dds_10wpi <-
  collapsed_dds_10wpi[keep, ]

as.data.frame(colData(collapsed_dds_10wpi))

levels(collapsed_dds_10wpi$samplingPoint)
levels(collapsed_dds_10wpi$treatment)

ddsDGE_heart10wpi <-
  DESeq(collapsed_dds_10wpi, parallel = T)

resultsNames(ddsDGE_heart10wpi)

res_shrunk_ivld_vs_conu <- lfcShrink(ddsDGE_heart10wpi, coef = "treatment_ivld_vs_conu", type = "ashr")

res_shrunk_eomes_vs_conu <- lfcShrink(ddsDGE_heart10wpi, coef = "treatment_eomes_vs_conu", type = "ashr")

res_shrunk_gata3_vs_conu <- lfcShrink(ddsDGE_heart10wpi, coef = "treatment_gata3_vs_conu", type = "ashr")

gsea_formatting(res_shrunk_ivld_vs_conu, 'shrunk', 'ivld', '10wpi')
gsea_formatting(res_shrunk_eomes_vs_conu, 'shrunk', 'eomes', '10wpi')
gsea_formatting(res_shrunk_gata3_vs_conu, 'shrunk', 'gata3', '10wpi')

shrunk_simplified_ivld_10wpi <- clusterProfiler::simplify(shrunk_gsea_results_ivld_10wpi)

shrunk_simplified_eomes_10wpi <- clusterProfiler::simplify(shrunk_gsea_results_eomes_10wpi)

shrunk_simplified_gata3_10wpi <- clusterProfiler::simplify(shrunk_gsea_results_gata3_10wpi)





# IV-LD ####
nrow(shrunk_gsea_results_ivld_10wpi)  # 2822 GO terms/pathways
nrow(shrunk_simplified_ivld_10wpi)  # 587 GO terms/pathways


top10_nes_ivld <- 
  as_tibble(shrunk_simplified_ivld_10wpi) %>%
  filter(NES > 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(desc(NES))


bottom10_nes_ivld <- 
  as_tibble(shrunk_simplified_ivld_10wpi) %>%
  filter(NES < 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(NES)

nes_ivld_10wpi <- bind_rows(top10_nes_ivld, bottom10_nes_ivld)

# EOMES ####
nrow(shrunk_gsea_results_eomes_10wpi)  # 2536 GO terms/pathways
nrow(shrunk_simplified_eomes_10wpi)  # 564 GO terms/pathways


top10_nes_eomes <- 
  as_tibble(shrunk_simplified_eomes_10wpi) %>%
  filter(NES > 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(desc(NES))


bottom10_nes_eomes <- 
  as_tibble(shrunk_simplified_eomes_10wpi) %>%
  filter(NES < 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(NES)

nes_eomes_10wpi <- bind_rows(top10_nes_eomes, bottom10_nes_eomes)

# EOMES ####
nrow(shrunk_gsea_results_gata3_10wpi)  # 2798 GO terms/pathways
nrow(shrunk_simplified_gata3_10wpi)  # 607 GO terms/pathways


top10_nes_gata3 <- 
  as_tibble(shrunk_simplified_gata3_10wpi) %>%
  filter(NES > 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(desc(NES))


bottom10_nes_gata3 <- 
  as_tibble(shrunk_simplified_gata3_10wpi) %>%
  filter(NES < 0) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length)) %>% 
  mutate(geneRatio = Count/setSize) %>% 
  top_n(10, wt = NES) %>% 
  arrange(NES)

nes_gata3_10wpi <- bind_rows(top10_nes_gata3, bottom10_nes_gata3)

intersect(nes_gata3_10wpi$Description, nes_ivld_10wpi$Description)






















