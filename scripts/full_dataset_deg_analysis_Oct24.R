#### Loading packages ####
suppressPackageStartupMessages({
  library('tidyverse')
  library('apeglm')
  library('DESeq2')
  library('LSD')
  library('BiocParallel')
  library('ashr')
  library('ggrepel')
  library('ggpmisc')
  library('DataCombine')
  library('org.Hs.eg.db')
  library('org.Dr.eg.db')
  library('ggrepel')
  register(MulticoreParam(10))
})

# Creating count matrix ####
# directory containing HTSeq count files
directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts'

sampleFiles <- grep('readcount', list.files(directory), value = T)

sampleTable <-
  data.frame(
    filename = sampleFiles,
    n = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[1]),
    treatment = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[2]),
    samplingPoint = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[3]),
    tissue = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[4]),
    lane = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[5]),
    readcount = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[6])
  )  # splitting sample names into different columns


head(sampleTable, n = 100)
summary(sampleTable)

sampleTable$tissue <-
  gsub('[0-9]+', '', sampleTable$tissue)  # removing numerical items from tissue column

sampleTable$treatment <-
  gsub('-', '', sampleTable$treatment) # removing dashes from treatment column

# defining factors
sampleTable$treatment <- factor(sampleTable$treatment)
sampleTable$samplingPoint <- factor(sampleTable$samplingPoint)
sampleTable$tissue <- factor(sampleTable$tissue)
sampleTable$lane <- factor(sampleTable$lane)

sampleTable$sample <-
  paste(
    rownames(sampleTable),
    as.character(sampleTable$treatment),
    as.character(sampleTable$samplingPoint),
    as.character(sampleTable$tissue),
    sep = '_'
  )

sampleTable <-
  as.data.frame(sampleTable) %>% dplyr::select(sample, filename, n, treatment, samplingPoint, tissue, lane)  # selecting columns of interest

sampleTable <- sampleTable %>%
  filter(!(tissue == 'hkr' |
             tissue == 'hr' |
             tissue == 'sr' |
             treatment == 'tbet'))
head(sampleTable, n = 1000) %>% filter(., lane == 'L4' | lane == 'L6')
head(sampleTable, n = 100)

save(sampleTable, file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable.RData')
# load(file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable.RData')

# Creating DESeqDataSet object and modelling with ~treatment * samplingPoint ####
load('~/Documents/PhD/Papers/Paper III/data/RData/sampleTables/sampleTable.RData')

sampleTable <- droplevels(sampleTable)
table(sampleTable$treatment, sampleTable$samplingPoint, sampleTable$tissue)
sampleTable_heart <- sampleTable %>% filter(tissue == 'h' & treatment %in% c('conu', 'ivld', 'eomes', 'gata3'))
sampleTable_heart <- droplevels(sampleTable_heart)
levels(sampleTable_heart$treatment)
sampleTable$samplingPoint <- relevel(sampleTable$samplingPoint, ref = '10wpi')
levels(sampleTable$samplingPoint)

dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_heart,
  directory = directory,
  design = ~ treatment * samplingPoint
)

keep <-
  rowSums(counts(dds)) >= 10  # removing low count genes (<10)
dds <-
  dds[keep, ]

as.data.frame(colData(dds))

collapsed_dds <- collapseReplicates(dds,
                   groupby = dds$n,
                   run = dds$lane)

collapsed_dds$samplingPoint <- relevel(collapsed_dds$samplingPoint, ref = '10wpi')
levels(collapsed_dds$samplingPoint)
levels(collapsed_dds$treatment)


ddsDGE_ensembl_heart <-
  DESeq(collapsed_dds, parallel = T)

resultsNames(ddsDGE_ensembl_heart)

save(ddsDGE_ensembl_heart, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/dds_DGE_ensembl_heart.RData')

save(ddsDGE_ensembl_fulldataset, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_ensembl_fulldataset.RData')

# Creating DESeqDataSet object and modelling with ~treatment * tissue at 10wpi ####
sampleTable_heartSpleen <- sampleTable %>% filter(samplingPoint == '10wpi' & tissue %in% c('h', 's') & treatment %in% c('conu', 'ivld', 'eomes', 'gata3'))

dds_10wpi <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_heartSpleen,
  directory = directory,
  design = ~ treatment * tissue
)

keep <-
  rowSums(counts(dds_10wpi)) >= 10  # removing low count genes (<10)
dds_10wpi <-
  dds_10wpi[keep, ]

as.data.frame(colData(dds_10wpi))

collapsed_dds_10wpi <- collapseReplicates(dds_10wpi,
                                    groupby = dds_10wpi$n,
                                    run = dds_10wpi$lane)

ddsDGE_heartSpleen_10wpi <-
  DESeq(collapsed_dds_10wpi, parallel = T)

resultsNames(ddsDGE_heartSpleen_10wpi)
## Exploratory analysis ####
# Plotting dispersion estimation to ensure that the assumption that most genes are not differentially expressed holds
DESeq2::plotDispEsts(ddsDGE_heartSpleen_10wpi)

# Transformation to stabilize variance across the mean through *variance stabilizing transformation*
vst_counts <- vst(ddsDGE_heartSpleen_10wpi, blind = T)

PCA <-
  plotPCA(
    vst_counts,
    intgroup = c('treatment', 'tissue'),
    returnData = T
  )

PCA <- PCA %>% 
  filter(treatment != 'ptagrfp')  # removing all ptagRFP samples for the PCA


percentVar <- round(100 * attr(PCA, "percentVar"))

library(RSkittleBrewer)
plotSkittles()
wildberry <- RSkittleBrewer('wildberry')


PCA$treatment <- factor(PCA$treatment,
                        levels = c('conu', 'ivld', 'eomes', 'gata3'))


PCA <- PCA %>% 
  filter(name != '639_dnavaccine_1wpc_h')  # removing heart outlier

## PCA

heartSpleen_10wpi_dataset_pca <- PCA %>%
  ggplot(.,
         aes(
           x = PC2,
           y = PC1,
           color = treatment,
           shape = tissue,
         )) +
  geom_point(size = 2) +
  # stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = .7) +
  ggtitle("Heart and Spleen PCA, treatment and tissue as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'eomes', 'gata3')
  ) +
  # scale_shape_manual(
  #   name = "tissue type",
  #   # Custom legend title
  #   values = c(16, 17, 18, 15),
  #   labels = c("heart", "head-kidney", "liver", "spleen")  # Replace with actual tissue types
  # ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))   

ggsave(filename = '~/Desktop/PCAs_full-dataset/full_dataset_pca.png', plot = full_dataset_pca)

spleen_pca <- PCA %>% subset(tissue == 's' & treatment %in% c('conu', 'ivld', 'eomes', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.2) +
  ggtitle("Spleen PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'eomes', 'gata3')
  ) +
  # scale_shape_manual(
  #   name = "sampling point",
  #   # Custom legend title
  #   values = c(16, 17, 18),
  #   labels = c("10 wpi", "4 wpc", "6 wpc")  # Replace with actual tissue types
  # ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )


heart_pca <- PCA %>% subset(tissue == 'h' & treatment %in% c('conu', 'ivld', 'eomes', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.2) +
  ggtitle("Heart PCA, treatment as factor") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'eomes', 'gata3')
  ) +
  # scale_shape_manual(
  #   name = "sampling point",
  #   # Custom legend title
  #   values = c(16, 17, 18),
  #   labels = c("10 wpi", "4 wpc", "6 wpc")  # Replace with actual tissue types
  # ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )



# Creating DESeqDataSet object and modelling with ~treatment * tissue at 4wpc ####
sampleTable_heartSpleen_4wpc <- sampleTable %>% filter(samplingPoint == '4wpc' & tissue %in% c('h', 's') & treatment %in% c('conu', 'ivld', 'eomes', 'gata3'))

dds_4wpc <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_heartSpleen_4wpc,
  directory = directory,
  design = ~ treatment * tissue
)

keep <-
  rowSums(counts(dds_4wpc)) >= 10  # removing low count genes (<10)
dds_4wpc <-
  dds_4wpc[keep, ]

as.data.frame(colData(dds_4wpc))

collapsed_dds_4wpc <- collapseReplicates(dds_4wpc,
                                          groupby = dds_4wpc$n,
                                          run = dds_4wpc$lane)

ddsDGE_heartSpleen_4wpc <-
  DESeq(collapsed_dds_4wpc, parallel = T)

resultsNames(ddsDGE_heartSpleen_4wpc)
## Exploratory analysis ####
# Plotting dispersion estimation to ensure that the assumption that most genes are not differentially expressed holds
DESeq2::plotDispEsts(ddsDGE_heartSpleen_4wpc)

# Transformation to stabilize variance across the mean through *variance stabilizing transformation*
vst_counts <- vst(ddsDGE_heartSpleen_4wpc, blind = T)

PCA <-
  plotPCA(
    vst_counts,
    intgroup = c('treatment', 'tissue'),
    returnData = T
  )

PCA <- PCA %>% 
  filter(treatment != 'ptagrfp')  # removing all ptagRFP samples for the PCA


percentVar <- round(100 * attr(PCA, "percentVar"))

library(RSkittleBrewer)
plotSkittles()
wildberry <- RSkittleBrewer('wildberry')


PCA$treatment <- factor(PCA$treatment,
                        levels = c('conu', 'ivld', 'eomes', 'gata3'))


PCA <- PCA %>% 
  filter(name != '639_dnavaccine_1wpc_h')  # removing heart outlier

## PCA

heartSpleen_10wpi_dataset_pca <- PCA %>%
  ggplot(.,
         aes(
           x = PC2,
           y = PC1,
           color = treatment,
           shape = tissue,
         )) +
  geom_point(size = 2) +
  # stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = .7) +
  ggtitle("Heart and Spleen PCA, treatment and tissue as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'eomes', 'gata3')
  ) +
  # scale_shape_manual(
  #   name = "tissue type",
  #   # Custom legend title
  #   values = c(16, 17, 18, 15),
  #   labels = c("heart", "head-kidney", "liver", "spleen")  # Replace with actual tissue types
  # ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))   

ggsave(filename = '~/Desktop/PCAs_full-dataset/full_dataset_pca.png', plot = full_dataset_pca)

spleen_pca <- PCA %>% subset(tissue == 's' & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment,
                shape = samplingPoint,
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.2) +
  ggtitle("Spleen PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
  ) +
  scale_shape_manual(
    name = "sampling point",
    # Custom legend title
    values = c(16, 17, 18),
    labels = c("10 wpi", "4 wpc", "6 wpc")  # Replace with actual tissue types
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )



########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
## Treatment contrasts within sampling points - HEART - 10wpi. Modelled with ~treatment * samplingPoint ----

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/dds_DGE_ensembl_heart.RData')
resultsNames(ddsDGE_ensembl_heart)
levels(ddsDGE_ensembl_heart$samplingPoint)
levels(ddsDGE_ensembl_heart$treatment)
ddsDGE_ensembl_heart <- relevel(ddsDGE_ensembl_heart$samplingPoint, ref = '10wpi')



res_shrunk_eomes_10wpi <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatment_eomes_vs_conu', type = 'ashr', parallel = T)
res_shrunk_gata3_10wpi <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatment_gata3_vs_conu', type = 'ashr', parallel = T)
res_shrunk_ivld_10wpi <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatment_ivld_vs_conu', type = 'ashr', parallel = T)

res_shrunk_eomes_4wpc <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatmenteomes.samplingPoint4wpc', type = 'ashr', parallel = T)
res_shrunk_gata3_4wpc <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatmentgata3.samplingPoint4wpc', type = 'ashr', parallel = T)
res_shrunk_ivld_4wpc <- lfcShrink(ddsDGE_ensembl_heart, coef = 'treatmentivld.samplingPoint4wpc', type = 'ashr', parallel = T)


########################################################################################################
########################################################################################################
########################################################################################################
########################################################################################################
## Exploratory analysis ####

# Plotting dispersion estimation to ensure that the assumption that most genes are not differentially expressed holds
DESeq2::plotDispEsts(ddsDGE_ensembl_fulldataset)

# Transformation to stabilize variance across the mean through *variance stabilizing transformation*
vst_counts <- vst(ddsDGE_ensembl_fulldataset, blind = T)

PCA <-
  plotPCA(
    vst_counts,
    intgroup = c('samplingPoint', 'treatment', 'tissue'),
    returnData = T
  )

PCA <- PCA %>% 
  filter(treatment != 'ptagrfp')  # removing all ptagRFP samples for the PCA


percentVar <- round(100 * attr(PCA, "percentVar"))

library(RSkittleBrewer)
plotSkittles()
original <- RSkittleBrewer('original') %>% .[. != 'yellow1']
wildberry <- RSkittleBrewer('wildberry')
mm <- RSkittleBrewer('M&M')

PCA$samplingPoint <- factor(PCA$samplingPoint,
                            levels = c('10wpi', '1wpc', '4wpc', '6wpc', '10wpc'))

PCA$treatment <- factor(PCA$treatment,
                        levels = c('conu', 'ptagrfp', 'ivld', 'ivhd', 'eomes', 'gata3', 'dnavaccine'))


PCA %>% 
  filter(tissue == 'h', PC1 > 50, PC2 > 50)

PCA <- PCA %>% 
  filter(name != '639_dnavaccine_1wpc_h')  # removing heart outlier

## PCA

full_dataset_pca <- PCA %>% subset(treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>%
  ggplot(.,
         aes(
           x = PC2,
           y = PC1,
           color = samplingPoint,
           shape = tissue,
         )) +
  geom_point(size = 2) +
  # stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = .7) +
  ggtitle("Full dataset PCA, tissue and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'sampling point',
    values = wildberry,
    labels = c('10 wpi', '1 wpc', '4 wpc', '6 wpc', '10 wpc')
  ) +
  scale_shape_manual(
    name = "tissue type",
    # Custom legend title
    values = c(16, 17, 18, 15),
    labels = c("heart", "head-kidney", "liver", "spleen")  # Replace with actual tissue types
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(color = guide_legend(order = 1), 
         shape = guide_legend(order = 2))   

ggsave(filename = '~/Desktop/PCAs_full-dataset/full_dataset_pca.png', plot = full_dataset_pca)

spleen_pca <- PCA %>% subset(tissue == 's' & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment,
                shape = samplingPoint,
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.2) +
  ggtitle("Spleen PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
  ) +
  scale_shape_manual(
    name = "sampling point",
    # Custom legend title
    values = c(16, 17, 18),
    labels = c("10 wpi", "4 wpc", "6 wpc")  # Replace with actual tissue types
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )

ggsave(filename = '~/Desktop/PCAs_full-dataset/spleen_pca.png', plot = spleen_pca)


hkidney_pca <- PCA %>% subset(tissue == 'hk' & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment,
                shape = samplingPoint,
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.2) +
  ggtitle("Head-kidney PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
  ) +
  scale_shape_manual(
    name = "sampling point",
    # Custom legend title
    values = c(16, 17, 18),
    labels = c("10 wpi", "4 wpc", "6 wpc")  # Replace with actual tissue types
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )

ggsave(filename = '~/Desktop/PCAs_full-dataset/hkidney_pca.png', plot = hkidney_pca)


liver_pca <- PCA %>% subset(tissue == 'l' & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>% 
  ggplot(., aes(x = PC2,
                y = PC1,
                color = treatment,
                shape = samplingPoint,
  )) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.4) +
  ggtitle("Liver PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'treatment',
    values = wildberry,
    labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
  ) +
  scale_shape_manual(
    name = "sampling point",
    # Custom legend title
    values = c(16, 17),
    labels = c("10 wpi", "4 wpc")  # Replace with actual tissue types
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 2),   # Ensure treatment legend is first
    shape = guide_legend(order = 1)    # Ensure sampling point legend is second
  )


ggsave(filename = '~/Desktop/PCAs_full-dataset/liver_pca.png', plot = liver_pca)

modified_wildberry <- c(wildberry, "red3", "purple4", "darkorange1")
# spleen_hk_pca <- PCA %>% subset(tissue %in% c('s', 'hk') & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>% 
#   ggplot(., aes(x = PC2,
#                 y = PC1,
#                 color = samplingPoint,
#                 shape = treatment,
#   )) +
#   geom_point(size = 2) +
#   xlab(paste0("PC2: ", percentVar[1], "% variance")) +
#   ylab(paste0("PC1: ", percentVar[2], "% variance")) +
#   coord_fixed(ratio = 1.4) +
#   stat_ellipse(aes(group = samplingPoint, color = samplingPoint), linetype = 'dashed') +
#   ggtitle("Spleen and head-kidney PCA, treatment and sampling point as factors") +
#   theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
#   scale_color_manual(
#     name = 'sampling point',
#     values = modified_wildberry,
#     labels = c('10 wpi', '4 wpc', '6 wpc')
#   ) +
#   scale_shape_manual(
#     name = 'treatment',
#     # Custom legend title
#     values = c(16, 17, 18, 0, 4),
#     labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
#   ) +
#   # scale_color_manual(
#   #   name = 'sampling point (ellipse)',
#   #   values = c('10 wpi' = 'red3', '4 wpc' = 'purple4', '6 wpc' = 'darkorange1')
#   # ) +
#   theme(
#     plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
#     panel.grid.minor = element_blank(),
#     plot.title = element_text(hjust = .5)
#   ) +
#   guides(
#     color = guide_legend(order = 3),   
#     shape = guide_legend(order = 1),
#     # color_guide = guide_legend(order = 2)
#   )
#   
# ggsave(filename = '~/Desktop/PCAs_full-dataset/spleen_hk_pca.png', plot = spleen_hk_pca)

PCA$samplingPoint <- as.factor(PCA$samplingPoint)
PCA$treatment <- as.factor(PCA$treatment)

spleen_hk_pca <- PCA %>% 
  subset(tissue %in% c('s', 'hk') & treatment %in% c('conu', 'eomes', 'ivld', 'ivhd', 'gata3')) %>%
  ggplot(aes(x = PC2, y = PC1, color = samplingPoint, shape = treatment)) +
  geom_point(size = 2) +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1.4) +
  stat_ellipse(aes(group = samplingPoint, color = samplingPoint), linetype = 'dashed') +  # Ensure grouping is clear
  ggtitle("Spleen and head-kidney PCA, treatment and sampling point as factors") +
  theme_linedraw(base_size = 14, base_family = 'Times New Roman') +
  scale_color_manual(
    name = 'sampling point',
    values = wildberry,
    labels = c('10 wpi', '4 wpc', '6 wpc')
  ) +
  scale_shape_manual(
    name = 'treatment',
    values = c(16, 17, 18, 0, 4),
    labels = c('conu', 'ivld', 'ivhd', 'eomes', 'gata3')
  ) +
  theme(
    plot.margin = grid::unit(c(2, 3, 2, 3), 'mm'),
    panel.grid.minor = element_blank(),
    plot.title = element_text(hjust = .5)
  ) +
  guides(
    color = guide_legend(order = 3),
    shape = guide_legend(order = 1)
  )

spleen_hk_pca
ggsave(filename = '~/Desktop/PCAs_full-dataset/spleen_hk_pca.png', plot = spleen_hk_pca)



# Summarize the data by tissue, treatment, and samplingPoint
summary_table <- sampleTable %>%
  group_by(treatment, tissue, samplingPoint) %>%
  summarise(count = n()) %>%
  arrange(treatment, tissue, samplingPoint)

# View the result
print(summary_table, n = Inf)

sampleTable %>% subset(., tissue == 's') %>%  # spleen ptagrfp was not sequenced by Alex at 4wpc
  dplyr::select(sample, treatment, samplingPoint, tissue) %>% 
  group_by(treatment, tissue, samplingPoint) %>% 
  summarise(count = n()) %>% 
  arrange(samplingPoint, treatment, tissue)

sampleTable %>% filter(tissue == 's' & treatment == 'ptagrfp')



# Creating DESeqDataSet object and modelling with ~group ####
## heart
sampleTable <- sampleTable %>% 
  filter(sample != '639_dnavaccine_1wpc_h')  # removing heart outlier identified in exploratory analysis


sampleTable_heart <- sampleTable %>%
  filter(tissue == 'h' &
           treatment %in% c('conu', 'ptagrfp', 'ivld', 'ivhd', 'eomes', 'gata3'))

sampleTable_heart %>%
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint)) -> sampleTable_heart_group 

# formatting variables
sampleTable_heart_group$group <-
  as.factor(sampleTable_heart_group$group)
sampleTable_heart_group <-
  droplevels(sampleTable_heart_group)  # dropping the reference

# creating dds object
ddsGroup_heart <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_heart_group,
                             directory = directory,
                             design = ~group)

# removing low count genes (<10)
keep <-
  rowSums(counts(ddsGroup_heart)) >= 10  

ddsGroup_heart <-
  ddsGroup_heart[keep, ]

as.data.frame(colData(ddsGroup_heart))

# collapsing replicates (if applicable)
collapsed_heart <- collapseReplicates(ddsGroup_heart,
                                       groupby = ddsGroup_heart$n,
                                       run = ddsGroup_heart$lane)

## spleen
sampleTable_spleen <- sampleTable %>%
  filter(tissue == 's' &
           treatment %in% c('conu', 'ptagrfp', 'ivld', 'ivhd', 'eomes', 'gata3'))

sampleTable_spleen %>%
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint)) -> sampleTable_spleen_group 


# formatting variables
sampleTable_spleen_group$group <-
  as.factor(sampleTable_spleen_group$group)
sampleTable_spleen_group <-
  droplevels(sampleTable_spleen_group)  # dropping the reference

# creating dds object
ddsGroup_spleen <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_spleen_group,
                             directory = directory,
                             design = ~group)

# removing low count genes (<10)
keep <-
  rowSums(counts(ddsGroup_spleen)) >= 10  
ddsGroup_spleen <-
  ddsGroup_spleen[keep, ]

as.data.frame(colData(ddsGroup_spleen))

# collapsing replicates (if applicable)
collapsed_spleen <- collapseReplicates(ddsGroup_spleen,
                                       groupby = ddsGroup_spleen$n,
                                       run = ddsGroup_spleen$lane)


ddsDGE_group_heart <-
  DESeq(collapsed_heart, parallel = T)

save(ddsDGE_group_heart, file = '~/Documents/PhD/Papers/Paper III/data/RData/ddsDGE_group_heart.RData')
save(sampleTable_heart_group, file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_heart_group.RData')

ddsDGE_group_spleen <-
  DESeq(collapsed_spleen, parallel = T)

save(ddsDGE_group_spleen, file = '~/Documents/PhD/Papers/Paper III/data/RData/ddsDGE_group_spleen.RData')
save(sampleTable_spleen_group, file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_spleen_group.RData')


# Loading data ----
load(file = '~/Documents/PhD/Papers/Paper III/data/RData/sampleTable_heart_group.RData')  # heart sampleTable
load(file = '~/Documents/PhD/Papers/Paper III/data/RData/ddsDGE_group_spleen.RData')  # spleen sampleTable
load(file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_heart.RData')  # heart DESeq model
load(file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_spleen.RData')  # spleen DESeq model

## Treatment contrasts within sampling points - HEART - 10wpi ----
### IVLD
heart_res_ivld_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### IVHD
heart_res_ivhd_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### GATA3 
heart_res_gata3_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### EOMES
heart_res_eomes_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)


## Treatment contrasts within sampling points - HEART - 4wpc ----

### IVLD
heart_res_ivld_vs_conu_4wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### IVHD
heart_res_ivhd_vs_conu_4wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### GATA3 
heart_res_gata3_vs_conu_4wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### EOMES
heart_res_eomes_vs_conu_4wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### GATA3 vs EOMES at 4wpc
heart_res_gata3_vs_eomes_4wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.4wpc', 'eomes.4wpc'), type = 'ashr', parallel = T)

## Treatment contrasts within sampling points - HEART - 6wpc ----

### IVLD
heart_res_ivld_vs_conu_6wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### IVHD
heart_res_ivhd_vs_conu_6wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### GATA3 
heart_res_gata3_vs_conu_6wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### EOMES
heart_res_eomes_vs_conu_6wpc <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)



## Treatment contrasts within sampling points - SPLEEN - 10wpi ----

### IVLD
spleen_res_ivld_vs_conu_10wpi <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### IVHD
spleen_res_ivhd_vs_conu_10wpi <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### GATA3 
spleen_res_gata3_vs_conu_10wpi <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### EOMES
spleen_res_eomes_vs_conu_10wpi <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)



## Treatment contrasts within sampling points - SPLEEN - 4wpc ----

### IVLD
spleen_res_ivld_vs_conu_4wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### IVHD
spleen_res_ivhd_vs_conu_4wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### GATA3 
spleen_res_gata3_vs_conu_4wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### EOMES
spleen_res_eomes_vs_conu_4wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

## Treatment contrasts within sampling points - SPLEEN - 6wpc ----

### IVLD
spleen_res_ivld_vs_conu_6wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### IVHD
spleen_res_ivhd_vs_conu_6wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### GATA3 
spleen_res_gata3_vs_conu_6wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)

### EOMES
spleen_res_eomes_vs_conu_6wpc <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.6wpc', 'conu.6wpc'), type = 'ashr', parallel = T)


## Treatment contrasts within treatments ACROSS sampling points - HEART ####

### 4 wpc vs 10 wpi
heart_res_conu_vs_conu_4v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'conu.4wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

heart_res_ivld_vs_ivld_4v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.4wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

heart_res_ivhd_vs_ivhd_4v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.4wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

heart_res_eomes_vs_eomes_4v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.4wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

heart_res_gata3_vs_gata3_4v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.4wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

summary(heart_res_gata3_vs_gata3_4v10)

### 6 wpc vs 10 wpi
heart_res_conu_vs_conu_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'conu.6wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

heart_res_ivld_vs_ivld_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.6wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

heart_res_ivhd_vs_ivhd_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.6wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

heart_res_eomes_vs_eomes_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.6wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

heart_res_gata3_vs_gata3_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.6wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

summary(heart_res_gata3_vs_gata3_6v10)


## Treatment contrasts within treatments ACROSS sampling points - SPLEEN ####
### 4 wpc vs 10 wpi
spleen_res_conu_vs_conu_4v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'conu.4wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivld_vs_ivld_4v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.4wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivhd_vs_ivhd_4v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.4wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

spleen_res_eomes_vs_eomes_4v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.4wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

spleen_res_gata3_vs_gata3_4v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.4wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

### 6 wpc vs 10 wpi
spleen_res_conu_vs_conu_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'conu.6wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivld_vs_ivld_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.6wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivhd_vs_ivhd_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.6wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

spleen_res_eomes_vs_eomes_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.6wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

spleen_res_gata3_vs_gata3_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.6wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Papers/Paper III/data/results_tables'
)

obj <- ls(pattern = '_res_.*_4v10')  # regex pattern matching files containing res
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


# Differential Gene Expression Analysis ----
## Install and load the required packages ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("clusterProfiler", "tidyverse", "pathview", "ReactomePA"))

suppressPackageStartupMessages({
library('clusterProfiler')
library('tidyverse')
library('pathview')
library('ReactomePA')
})

renv::snapshot()


## CONU vs CONU @ 4 wpc vs 10 wpi - HEART ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/results_tables'
)

results_files <-
  list.files(pattern = '^heart_.*_4v10')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

gsea_formatting(heart_res_conu_vs_conu_4v10, 'heart', 'conu', '4v10')

heart_gsea_simplified_conu4v10 <-
  clusterProfiler::simplify(heart_gsea_results_conu_4v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_conu4v10 <- entrez_gene_list

nrow(heart_gsea_results_conu_4v10)  # 2422 GO terms/pathways
nrow(heart_gsea_simplified_conu4v10)  # 547 GO terms/pathways

as_tibble(heart_gsea_simplified_conu4v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_conu4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_conu4v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes_conu4v10 <-
  bind_rows(top10_high_nes, bottom10_low_nes)

gsea_4v10_heart_conu <- as_tibble(low_high_nes_conu4v10) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

# Install the writexl package if you don't have it
if (!requireNamespace("writexl", quietly = TRUE)) {
  install.packages("writexl")
}

# Load the writexl library
library(writexl)

write_xlsx(gsea_4v10_heart_conu, '~/Documents/PhD/Papers/Paper III/data/gsea_4v10_heart_conu.xlsx')

## EOMES vs EOMES @ 4 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_eomes_vs_eomes_4v10, 'heart', 'eomes', '4v10')

heart_gsea_simplified_eomes4v10 <-
  clusterProfiler::simplify(heart_gsea_results_eomes_4v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_eomes4v10 <- entrez_gene_list

nrow(heart_gsea_results_eomes_4v10)  # 484 GO terms/pathways
nrow(heart_gsea_simplified_eomes4v10)  # 255 GO terms/pathways

as_tibble(heart_gsea_simplified_eomes4v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_eomes4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))


gsea_4v10_heart_eomes <- as_tibble(top10_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)


# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_eomes4v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))
# 
# low_high_nes_eomes4v10 <-
#   bind_rows(top10_high_nes, bottom10_low_nes)


## GATA3 vs GATA3 @ 4 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_gata3_vs_gata3_4v10, 'heart', 'gata3', '4v10')

heart_gsea_simplified_gata34v10 <-
  clusterProfiler::simplify(heart_gsea_results_gata3_4v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_gata34v10 <- entrez_gene_list

nrow(heart_gsea_results_gata3_4v10)  # 153 GO terms/pathways
nrow(heart_gsea_simplified_gata34v10)  # 96 GO terms/pathways

as_tibble(heart_gsea_simplified_gata34v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_gata34v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_gata34v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))


gsea_4v10_heart_gata3 <- as_tibble(top10_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

## Very few differentially regulated genes in GATA3

## IV-LD vs IV-LD @ 4 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_ivld_vs_ivld_4v10, 'heart', 'ivld', '4v10')

heart_gsea_simplified_ivld4v10 <-
  clusterProfiler::simplify(heart_gsea_results_ivld_4v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_ivld4v10 <- entrez_gene_list

nrow(heart_gsea_results_ivld_4v10)  # 2003 GO terms/pathways
nrow(heart_gsea_simplified_ivld4v10)  # 481 GO terms/pathways

as_tibble(heart_gsea_simplified_ivld4v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_ivld4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_ivld4v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

gsea_4v10_heart_ivld <- as_tibble(top10_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '^heart*.*_gsea_.*4v10.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '*_gene_list_*.*4v10')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
############################
setwd(
  '~/Documents/PhD/Papers/Paper III/data/results_tables'
)

results_files <-
  list.files(pattern = '^spleen_.*_4v10')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}
## CONU vs CONU @ 4 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_conu_vs_conu_4v10, 'spleen', 'conu', '4v10')

spleen_gsea_simplified_conu4v10 <-
  clusterProfiler::simplify(spleen_gsea_results_conu_4v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_conu4v10 <- entrez_gene_list

nrow(spleen_gsea_results_conu_4v10)  # 1615 GO terms/pathways
nrow(spleen_gsea_simplified_conu4v10)  # 408 GO terms/pathways

as_tibble(spleen_gsea_simplified_conu4v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_conu4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_conu4v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

gsea_4v10_spleen_conu <- as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

## EOMES vs EOMES @ 4 wpc vs 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_eomes_vs_eomes_4v10, 'spleen', 'eomes', '4v10')

spleen_gsea_simplified_eomes4v10 <-
  clusterProfiler::simplify(spleen_gsea_results_eomes_4v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_eomes4v10 <- entrez_gene_list

nrow(spleen_gsea_results_eomes_4v10)  # 82 GO terms/pathways
nrow(spleen_gsea_simplified_eomes4v10)  # 39 GO terms/pathways

as_tibble(spleen_gsea_simplified_eomes4v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_eomes4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_eomes4v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

gsea_4v10_spleen_eomes <- as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)
## GATA3 vs GATA3 @ 4 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_gata3_vs_gata3_4v10, 'spleen', 'gata3', '4v10')

spleen_gsea_simplified_gata34v10 <-
  clusterProfiler::simplify(spleen_gsea_results_gata3_4v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_gata34v10 <- entrez_gene_list

nrow(spleen_gsea_results_gata3_4v10)  # 482 GO terms/pathways
nrow(spleen_gsea_simplified_gata34v10)  # 190 GO terms/pathways

as_tibble(spleen_gsea_simplified_gata34v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_gata34v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_gata34v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))


low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

gsea_4v10_spleen_gata3 <- as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

## IV-LD vs IV-LD @ 4 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_ivld_vs_ivld_4v10, 'spleen', 'ivld', '4v10')

spleen_gsea_simplified_ivld4v10 <-
  clusterProfiler::simplify(spleen_gsea_results_ivld_4v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_ivld4v10 <- entrez_gene_list

nrow(spleen_gsea_results_ivld_4v10)  # 1646 GO terms/pathways
nrow(spleen_gsea_simplified_ivld4v10)  # 441 GO terms/pathways

as_tibble(spleen_gsea_simplified_ivld4v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_ivld4v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_ivld4v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

gsea_4v10_spleen_ivld <- as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)


# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^spleen*.*_gsea_.*4v10.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '*_gene_list_*.*4v10')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
#############################



## CONU vs CONU @ 6 wpc vs 10 wpi - HEART ----
gsea_formatting(heart_res_conu_vs_conu_6v10, 'heart', 'conu', '6v10')

heart_gsea_simplified_conu6v10 <-
  clusterProfiler::simplify(heart_gsea_results_conu_6v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_conu6v10 <- entrez_gene_list

nrow(heart_gsea_results_conu_6v10)  # 2426 GO terms/pathways
nrow(heart_gsea_simplified_conu6v10)  # 719 GO terms/pathways

as_tibble(heart_gsea_simplified_conu6v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_conu6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_conu6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)
## EOMES vs EOMES @ 6 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_eomes_vs_eomes_6v10, 'heart', 'eomes', '6v10')

heart_gsea_simplified_eomes6v10 <-
  clusterProfiler::simplify(heart_gsea_results_eomes_6v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_eomes6v10 <- entrez_gene_list

nrow(heart_gsea_results_eomes_6v10)  # 2426 GO terms/pathways
nrow(heart_gsea_simplified_eomes6v10)  # 719 GO terms/pathways

as_tibble(heart_gsea_simplified_eomes6v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_eomes6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_eomes6v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>% 
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

## GATA3 vs GATA3 @ 6 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_gata3_vs_gata3_6v10, 'heart', 'gata3', '6v10')

heart_gsea_simplified_gata36v10 <-
  clusterProfiler::simplify(heart_gsea_results_gata3_6v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_gata36v10 <- entrez_gene_list

nrow(heart_gsea_results_gata3_6v10)  # 58 GO terms/pathways
nrow(heart_gsea_simplified_gata36v10)  # 39 GO terms/pathways

as_tibble(heart_gsea_simplified_gata36v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_gata36v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_gata36v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

summary(heart_res_gata3_vs_conu_10wpi)
summary(heart_res_gata3_vs_conu_6wpc)
summary(heart_res_gata3_vs_gata3_6v10)

sampleTable %>% filter(., treatment == 'gata3' & samplingPoint == '6wpc' & tissue == 's')

## Very few differentially regulated genes in GATA3

## IV-LD vs IV-LD @ 6 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_ivld_vs_ivld_6v10, 'heart', 'ivld', '6v10')

heart_gsea_simplified_ivld6v10 <-
  clusterProfiler::simplify(heart_gsea_results_ivld_6v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_ivld6v10 <- entrez_gene_list

nrow(heart_gsea_results_ivld_6v10)  # 1476 GO terms/pathways
nrow(heart_gsea_simplified_ivld6v10)  # 425 GO terms/pathways

as_tibble(heart_gsea_simplified_ivld6v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_ivld6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_ivld6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-HD vs IV-HD @ 6 wpc vs 10 wpi - HEART ----

gsea_formatting(heart_res_ivhd_vs_ivhd_6v10, 'heart', 'ivhd', '6v10')

heart_gsea_simplified_ivhd6v10 <-
  clusterProfiler::simplify(heart_gsea_results_ivhd_6v10)  # simplifying GO terms to reduce redundancy

heart_gene_list_ivhd6v10 <- entrez_gene_list

nrow(heart_gsea_results_ivhd_6v10)  # 1676 GO terms/pathways
nrow(heart_gsea_simplified_ivhd6v10)  # 468 GO terms/pathways

as_tibble(heart_gsea_simplified_ivhd6v10)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_ivhd6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# bottom10_low_nes <-
#   as_tibble(heart_gsea_simplified_ivhd6v10) %>%
#   filter(NES < 0) %>%
#   arrange(desc(setSize)) %>%
#   top_n(10, wt = setSize) %>%
#   mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

# low_high_nes <-
#   bind_rows(top10_high_nes, bottom10_low_nes)


# Saving results files
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '.*_gsea_.*6v10.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

############################
## CONU vs CONU @ 6 wpc vs 10 wpi - SPLEEN ----
  
gsea_formatting(spleen_res_conu_vs_conu_6v10, 'spleen', 'conu', '6v10')

spleen_gsea_simplified_conu6v10 <-
  clusterProfiler::simplify(spleen_gsea_results_conu_6v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_conu6v10 <- entrez_gene_list

nrow(spleen_gsea_results_conu_6v10)  # 2284 GO terms/pathways
nrow(spleen_gsea_simplified_conu6v10)  # 517 GO terms/pathways

as_tibble(spleen_gsea_simplified_conu6v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_conu6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_conu6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)


## EOMES vs EOMES @ 6 wpc vs 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_eomes_vs_eomes_6v10, 'spleen', 'eomes', '6v10')

spleen_gsea_simplified_eomes6v10 <-
  clusterProfiler::simplify(spleen_gsea_results_eomes_6v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_eomes6v10 <- entrez_gene_list

nrow(spleen_gsea_results_eomes_6v10)  # 1786 GO terms/pathways
nrow(spleen_gsea_simplified_eomes6v10)  # 449 GO terms/pathways

as_tibble(spleen_gsea_simplified_eomes6v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_eomes6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_eomes6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## GATA3 vs GATA3 @ 6 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_gata3_vs_gata3_6v10, 'spleen', 'gata3', '6v10')

spleen_gsea_simplified_gata36v10 <-
  clusterProfiler::simplify(spleen_gsea_results_gata3_6v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_gata36v10 <- entrez_gene_list

nrow(spleen_gsea_results_gata3_6v10)  # 1705 GO terms/pathways
nrow(spleen_gsea_simplified_gata36v10)  # 443 GO terms/pathways

as_tibble(spleen_gsea_simplified_gata36v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_gata36v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_gata36v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))


low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-LD vs IV-LD @ 6 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_ivld_vs_ivld_6v10, 'spleen', 'ivld', '6v10')

spleen_gsea_simplified_ivld6v10 <-
  clusterProfiler::simplify(spleen_gsea_results_ivld_6v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_ivld6v10 <- entrez_gene_list

nrow(spleen_gsea_results_ivld_6v10)  # 1646 GO terms/pathways
nrow(spleen_gsea_simplified_ivld6v10)  # 441 GO terms/pathways

as_tibble(spleen_gsea_simplified_ivld6v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_ivld6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_ivld6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-HD vs IV-HD @ 6 wpc vs 10 wpi - SPLEEN ----

gsea_formatting(spleen_res_ivhd_vs_ivhd_6v10, 'spleen', 'ivhd', '6v10')

spleen_gsea_simplified_ivhd6v10 <-
  clusterProfiler::simplify(spleen_gsea_results_ivhd_6v10)  # simplifying GO terms to reduce redundancy

spleen_gene_list_ivhd6v10 <- entrez_gene_list

nrow(spleen_gsea_results_ivhd_6v10)  # 906 GO terms/pathways
nrow(spleen_gsea_simplified_ivhd6v10)  # 320 GO terms/pathways

as_tibble(spleen_gsea_simplified_ivhd6v10)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_ivhd6v10) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_ivhd6v10) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)


# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^spleen*.*_gsea_.*6v10.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
  
# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '*_gene_list_*.*6v10')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
#############################
## EOMES vs CONU @ 10 wpi - HEART ----
gsea_formatting(heart_res_eomes_vs_conu_10wpi, 'heart', 'EOMESconu', '10wpi')

heart_gsea_simplified_EOMESconu10wpi <-
  clusterProfiler::simplify(heart_gsea_results_EOMESconu_10wpi)  # simplifying GO terms to reduce redundancy

heart_gene_list_EOMESconu_10wpi <- entrez_gene_list

nrow(heart_gsea_results_EOMESconu_10wpi)  # 1 GO terms/pathways
nrow(heart_gsea_simplified_EOMESconu10wpi)  # 1 GO terms/pathways

as_tibble(heart_gsea_simplified_EOMESconu10wpi)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_EOMESconu10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_EOMESconu10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## GATA3 vs CONU @ 10 wpi - HEART ----
gsea_formatting(heart_res_gata3_vs_conu_10wpi, 'heart', 'GATA3conu', '10wpi')

heart_gsea_simplified_GATA3conu10wpi <-
  clusterProfiler::simplify(heart_gsea_results_GATA3conu_10wpi)  # simplifying GO terms to reduce redundancy

heart_gene_list_GATA3conu_10wpi <- entrez_gene_list

nrow(heart_gsea_results_GATA3conu_10wpi)  # 1 GO terms/pathways
nrow(heart_gsea_simplified_GATA3conu10wpi)  # 1 GO terms/pathways

as_tibble(heart_gsea_simplified_GATA3conu10wpi)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_gata3CONU10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_gata3CONU10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)
## IV-LD vs CONU @ 10 wpi - HEART ----
gsea_formatting(heart_res_ivld_vs_conu_10wpi, 'heart', 'IVLDconu', '10wpi')

heart_gsea_simplified_IVLDconu10wpi <-
  clusterProfiler::simplify(heart_gsea_results_IVLDconu_10wpi)  # simplifying GO terms to reduce redundancy

heart_gene_list_IVLDconu_10wpi <- entrez_gene_list

nrow(heart_gsea_results_IVLDconu_10wpi)  # no enriched terms
nrow(heart_gsea_simplified_IVLDconu10wpi)  # no enriched terms

as_tibble(heart_gsea_simplified_IVLDconu10wpi)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_ivldCONU10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_ivldCONU10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)
## IV-HD vs CONU @ 10 wpi - HEART ----
gsea_formatting(heart_res_ivhd_vs_conu_10wpi, 'heart', 'IVHDconu', '10wpi')

heart_gsea_simplified_IVHDconu10wpi <-
  clusterProfiler::simplify(heart_gsea_results_IVHDconu_10wpi)  # simplifying GO terms to reduce redundancy

heart_gene_list_IVHDconu_10wpi <- entrez_gene_list

nrow(heart_gsea_results_IVHDconu_10wpi)  # no enriched terms
nrow(heart_gsea_simplified_IVHDconu10wpi)  # no enriched terms

as_tibble(heart_gsea_simplified_IVHDconu10wpi)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_ivhdCONU10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_ivhdCONU10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^heart*.*_gsea_.*10wpi.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^heart_gene_list_*.*10wpi')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
##############################
## EOMES vs CONU @ 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_eomes_vs_conu_10wpi, 'spleen', 'EOMESconu', '10wpi')

spleen_gsea_simplified_EOMESconu10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_EOMESconu_10wpi)  # simplifying GO terms to reduce redundancy

spleen_gene_list_EOMESconu_10wpi <- entrez_gene_list

nrow(spleen_gsea_results_EOMESconu_10wpi)  # 785 GO terms/pathways
nrow(spleen_gsea_simplified_EOMESconu10wpi)  # 262 GO terms/pathways

as_tibble(spleen_gsea_simplified_EOMESconu10wpi)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_EOMESconu10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-  
  as_tibble(spleen_gsea_simplified_EOMESconu10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)


## GATA3 vs CONU @ 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_gata3_vs_conu_10wpi, 'spleen', 'GATA3conu', '10wpi')

spleen_gsea_simplified_GATA3conu10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_GATA3conu_10wpi)  # simplifying GO terms to reduce redundancy

spleen_gene_list_GATA3conu_10wpi <- entrez_gene_list

nrow(spleen_gsea_results_GATA3conu_10wpi)  # 46 GO terms/pathways
nrow(spleen_gsea_simplified_GATA3conu10wpi)  # 34 GO terms/pathways

as_tibble(spleen_gsea_simplified_GATA3conu10wpi)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_GATA3conu10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_GATA3conu10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

## IV-LD vs CONU @ 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_ivld_vs_conu_10wpi, 'spleen', 'IVLDconu', '10wpi')

spleen_gsea_simplified_IVLDconu10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_IVLDconu_10wpi)  # simplifying GO terms to reduce redundancy

spleen_gene_list_IVLDconu_10wpi <- entrez_gene_list

nrow(spleen_gsea_results_IVLDconu_10wpi)  # 2554 enriched terms/pathways
nrow(spleen_gsea_simplified_IVLDconu10wpi)  # 1599 enriched terms/pathways

as_tibble(spleen_gsea_simplified_IVLDconu10wpi)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_IVLDconu10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_IVLDconu10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-HD vs CONU @ 10 wpi - SPLEEN ----
gsea_formatting(spleen_res_ivhd_vs_conu_10wpi, 'spleen', 'IVHDconu', '10wpi')

spleen_gsea_simplified_IVHDconu10wpi <-
  clusterProfiler::simplify(spleen_gsea_results_IVHDconu_10wpi)  # simplifying GO terms to reduce redundancy

spleen_gene_list_IVHDconu_10wpi <- entrez_gene_list

nrow(spleen_gsea_results_IVHDconu_10wpi)  # no enriched terms
nrow(spleen_gsea_simplified_IVHDconu10wpi)  # no enriched terms

as_tibble(spleen_gsea_simplified_IVHDconu10wpi)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_IVHDconu10wpi) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_IVHDconu10wpi) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))


summary(spleen_res_ivhd_vs_conu_10wpi)
# out of 36905 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 34, 0.092%
# LFC < 0 (down)     : 8, 0.022%
# outliers [1]       : 139, 0.38%
# low counts [2]     : 8567, 23%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

summary(spleen_res_ivld_vs_conu_10wpi)
# out of 36905 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4, 0.011%
# LFC < 0 (down)     : 3, 0.0081%
# outliers [1]       : 139, 0.38%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

obj_to_remove <- ls(pattern = 'gsea_results|gsea_simplified')
gene_list_to_remove <- ls(pattern = 'gene_list*.*10wpi')
rm(list = obj_to_remove)
rm(list = gene_list_to_remove)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^spleen*.*_gsea_.*10wpi.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '*_gene_list_*.*10wpi')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  



##############################
## EOMES vs CONU @ 6 wpc - HEART ----
gsea_formatting(heart_res_eomes_vs_conu_6wpc, 'heart', 'EOMESconu', '6wpc')

heart_gsea_simplified_EOMESconu6wpc <-
  clusterProfiler::simplify(heart_gsea_results_EOMESconu_6wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_EOMESconu_6wpc <- entrez_gene_list

nrow(heart_gsea_results_EOMESconu_6wpc)  # 917 GO terms/pathways
nrow(heart_gsea_simplified_EOMESconu6wpc)  # 316 GO terms/pathways

as_tibble(heart_gsea_simplified_EOMESconu6wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_EOMESconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_EOMESconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## GATA3 vs CONU @ 6 wpc - HEART ----
gsea_formatting(heart_res_gata3_vs_conu_6wpc, 'heart', 'GATA3conu', '6wpc')

heart_gsea_simplified_GATA3conu6wpc <-
  clusterProfiler::simplify(heart_gsea_results_GATA3conu_6wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_GATA3conu_6wpc <- entrez_gene_list

nrow(heart_gsea_results_GATA3conu_6wpc)  # 236 GO terms/pathways
nrow(heart_gsea_simplified_GATA3conu6wpc)  # 134 GO terms/pathways

as_tibble(heart_gsea_simplified_GATA3conu6wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_GATA3conu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_GATA3conu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-LD vs CONU @ 6 wpc - HEART ----
gsea_formatting(heart_res_ivld_vs_conu_6wpc, 'heart', 'IVLDconu', '6wpc')

heart_gsea_simplified_IVLDconu6wpc <-
  clusterProfiler::simplify(heart_gsea_results_IVLDconu_6wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_IVLDconu_6wpc <- entrez_gene_list

nrow(heart_gsea_results_IVLDconu_6wpc)  # no enriched terms
nrow(heart_gsea_simplified_IVLDconu6wpc)  # no enriched terms

as_tibble(heart_gsea_simplified_IVLDconu6wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_IVLDconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_IVLDconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)
## IV-HD vs CONU @ 6 wpc - HEART ----
gsea_formatting(heart_res_ivhd_vs_conu_6wpc, 'heart', 'IVHDconu', '6wpc')

heart_gsea_simplified_IVHDconu6wpc <-
  clusterProfiler::simplify(heart_gsea_results_IVHDconu_6wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_IVHDconu_6wpc <- entrez_gene_list

nrow(heart_gsea_results_IVHDconu_6wpc)  # no enriched terms
nrow(heart_gsea_simplified_IVHDconu6wpc)  # no enriched terms

as_tibble(heart_gsea_simplified_IVHDconu6wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_IVHDconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_IVHDconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^heart*.*_gsea_.*6wpc.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^heart_gene_list_*.*6wpc')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  
##############################
## EOMES vs CONU @ 6 wpc - SPLEEN ----
gsea_formatting(spleen_res_eomes_vs_conu_6wpc, 'spleen', 'EOMESconu', '6wpc')

spleen_gsea_simplified_EOMESconu6wpc <-
  clusterProfiler::simplify(spleen_gsea_results_EOMESconu_6wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_EOMESconu_6wpc <- entrez_gene_list

nrow(spleen_gsea_results_EOMESconu_6wpc)  # 0 GO terms/pathways
nrow(spleen_gsea_simplified_EOMESconu6wpc)  # 0 GO terms/pathways

as_tibble(spleen_gsea_simplified_EOMESconu6wpc)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_EOMESconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-  
  as_tibble(spleen_gsea_simplified_EOMESconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)


## GATA3 vs CONU @ 6 wpc - SPLEEN ----
gsea_formatting(spleen_res_gata3_vs_conu_6wpc, 'spleen', 'GATA3conu', '6wpc')

spleen_gsea_simplified_GATA3conu6wpc <-
  clusterProfiler::simplify(spleen_gsea_results_GATA3conu_6wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_GATA3conu_6wpc <- entrez_gene_list

nrow(spleen_gsea_results_GATA3conu_6wpc)  # 46 GO terms/pathways
nrow(spleen_gsea_simplified_GATA3conu6wpc)  # 34 GO terms/pathways

as_tibble(spleen_gsea_simplified_GATA3conu6wpc)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_GATA3conu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_GATA3conu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-LD vs CONU @ 6 wpc - SPLEEN ----
gsea_formatting(spleen_res_ivld_vs_conu_6wpc, 'spleen', 'IVLDconu', '6wpc')

spleen_gsea_simplified_IVLDconu6wpc <-
  clusterProfiler::simplify(spleen_gsea_results_IVLDconu_6wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_IVLDconu_6wpc <- entrez_gene_list

nrow(spleen_gsea_results_IVLDconu_6wpc)  # 28 enriched terms/pathways
nrow(spleen_gsea_simplified_IVLDconu6wpc)  # 23 enriched terms/pathways

as_tibble(spleen_gsea_simplified_IVLDconu6wpc)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_IVLDconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_IVLDconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

## IV-HD vs CONU @ 6 wpc - SPLEEN ----
gsea_formatting(spleen_res_ivhd_vs_conu_6wpc, 'spleen', 'IVHDconu', '6wpc')

spleen_gsea_simplified_IVHDconu6wpc <-
  clusterProfiler::simplify(spleen_gsea_results_IVHDconu_6wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_IVHDconu_6wpc <- entrez_gene_list

nrow(spleen_gsea_results_IVHDconu_6wpc)  # no enriched terms
nrow(spleen_gsea_simplified_IVHDconu6wpc)  # no enriched terms

as_tibble(spleen_gsea_simplified_IVHDconu6wpc)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_IVHDconu6wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(spleen_gsea_simplified_IVHDconu6wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)


summary(spleen_res_ivhd_vs_conu_6wpc)
# out of 36905 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 34, 0.092%
# LFC < 0 (down)     : 8, 0.022%
# outliers [1]       : 139, 0.38%
# low counts [2]     : 8567, 23%
# (mean count < 1)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

summary(spleen_res_ivld_vs_conu_6wpc)
# out of 36905 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 4, 0.011%
# LFC < 0 (down)     : 3, 0.0081%
# outliers [1]       : 139, 0.38%
# low counts [2]     : 0, 0%
# (mean count < 0)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

obj_to_remove <- ls(pattern = 'gsea_results|gsea_simplified')
gene_list_to_remove <- ls(pattern = 'gene_list*.*6wpc')
rm(list = obj_to_remove)
rm(list = gene_list_to_remove)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '^spleen*.*_gsea_.*6wpc.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results'
)

obj <- ls(pattern = '*_gene_list_*.*6wpc')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  




load('~/Documents/PhD/Papers/Paper III/data/gsea_results/6wpc vs 10wpi/spleen_gsea_simplified_eomes6v10.RData')
load('~/Documents/PhD/Papers/Paper III/data/gsea_results/6wpc vs 10wpi/spleen_gsea_simplified_gata36v10.RData')
load('~/Documents/PhD/Papers/Paper III/data/gsea_results/6wpc vs 10wpi/spleen_gsea_simplified_conu6v10.RData')

as_tibble(spleen_gsea_simplified_eomes6v10)
inner_join(as_tibble(spleen_gsea_simplified_eomes6v10), as_tibble(spleen_gsea_simplified_gata36v10), by = 'Description')

conu_6v10 <- as_tibble(spleen_gsea_simplified_conu6v10) %>% dplyr::select(Description, setSize, NES, core_enrichment)
eomes_6v10 <- as_tibble(spleen_gsea_simplified_eomes6v10) %>% dplyr::select(Description, setSize, NES, core_enrichment)
gata_6v10 <- as_tibble(spleen_gsea_simplified_gata36v10) %>% dplyr::select(Description, setSize, NES, core_enrichment)


inner_join(eomes_6v10, conu_6v10, by = 'Description') %>% dplyr::select(Description, NES.x, NES.y) %>% head(n = 20)


load('~/Documents/PhD/Papers/Paper III/data/gsea_results/6 wpc/spleen_gsea_results_EOMESconu_6wpc.RData')
load('~/Documents/PhD/Papers/Paper III/data/gsea_results/6 wpc/spleen_gsea_results_GATA3conu_6wpc.RData')

EOMESconu_6wpc <- as_tibble(spleen_gsea_results_EOMESconu_6wpc)


EOMESvGATA3 <- intersect(eomes_6v10$Description, gata_6v10$Description)

unique_EOMESvGATA3 <- setdiff(EOMESvGATA3,
                              c(conu_6v10$Description))


eomes_6v10 %>% filter(Description == 'regulation of immune response')
gata_6v10 %>% filter(Description == 'regulation of immune response')

# pull genes from pathway, remove duplicates, sort, and print 1 per line
eomes_6v10 %>% filter(Description == 'regulation of immune response') %>% 
  pull(core_enrichment) %>% 
  strsplit(split = '/') %>% 
  unlist() %>% 
  unique() %>% 
  sort() %>% 
  cat(., sep = '\n')

eomes_6v10_genes <- eomes_6v10 %>% filter(Description == 'regulation of immune response') %>% 
  pull(core_enrichment) %>% 
  strsplit(split = '/') %>% 
  unlist() %>% 
  unique() %>% 
  sort()

gata_6v10_genes <- gata_6v10 %>% filter(Description == 'regulation of immune response') %>% 
  pull(core_enrichment) %>% 
  strsplit(split = '/') %>% 
  unlist() %>% 
  unique() %>% 
  sort()

setdiff.Vector(eomes_6v10_genes, gata_6v10_genes) %>% unlist() %>% cat(., sep = '\n')







### Exporting data for Hetron ----
as_tibble(heart_gsea_simplified_EOMESconu6wpc)
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables/heart_res_eomes_vs_conu_10wpi.RData')
improved_data_wrangling(heart_res_eomes_vs_conu_10wpi, 'EOMESconu', '10wpi')
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables/heart_res_eomes_vs_conu_6wpc.RData')
improved_data_wrangling(heart_res_eomes_vs_conu_6wpc, 'EOMESconu', '6wpc')
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables/heart_res_eomes_vs_eomes_6v10.RData')
improved_data_wrangling(heart_res_eomes_vs_eomes_6v10, 'EOMESeomes', '6wpcVS10wpi')
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables/heart_res_conu_vs_conu_6v10.RData')
improved_data_wrangling(heart_res_conu_vs_conu_6v10, 'CONUconu', '6wpcVS10wpi')
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/gsea_results/6wpc vs 10wpi/heart_gsea_simplified_eomes6v10.RData')
as_tibble(heart_gsea_simplified_eomes6v10)
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/gsea_results/6wpc vs 10wpi/heart_gsea_simplified_conu6v10.RData')

library(readr)
write_csv(results_EOMESconu_10wpi, '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/results_EOMESconu_10wpi.csv')
write_csv(results_EOMESconu_6wpc, '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/results_EOMESconu_6wpc.csv')
write_csv(results_EOMESeomes_6wpcVS10wpi, '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/results_EOMESeomes_6wpcVS10wpi.csv')
write_csv(results_CONUconu_6wpcVS10wpi, '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/results_CONUconu_6wpcVS10wpi.csv')
write_csv(as_tibble(heart_gsea_simplified_conu6v10), '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/heart_gsea_conu_6wpcVS10wpi.csv')
write_csv(as_tibble(heart_gsea_simplified_eomes6v10), '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/heart_gsea_eomes_6wpcVS10wpi.csv')
write_csv(as_tibble(heart_gsea_simplified_EOMESconu6wpc), '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/data_for_Hetron/heart_gsea_EOMESconu_6wpc.csv')






load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/gsea_results/10 wpi/simplified/heart_gsea_simplified_EOMESconu10wpi.RData')
load('/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/gsea_results/6 wpc/simplified/heart_gsea_simplified_EOMESconu6wpc.RData')

as_tibble(heart_gsea_simplified_EOMESconu10wpi)
as_tibble(heart_gsea_simplified_EOMESconu6wpc)



# Find objects that are of class DESeqResults
obj_to_remove <- ls()[sapply(ls(), function(x) inherits(get(x), "DESeqResults"))]

# Remove those objects from the environment
rm(list = obj_to_remove)




### Testing zebrafish annotation Db ----
zebrafish_data_wrangling <-
  function(results_table) {
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
      target_organism = 'drerio',
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
                    pvalue,
                    ortholog_ensg,
                    description)
    
    # create results name based on input table name
    results_name <- paste0(deparse(substitute(results_table)), '_wrangled')
    
    # assign the results data frame to the dynamic name
    assign(results_name, results, envir = .GlobalEnv)
  }


is.na(results_EOMESconu_drerio_6wpc$ortholog_name) %>% sum()
is.na(results_EOMESconu_6wpc$ortholog_name) %>% sum()



if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("org.Dr.eg.db")
library(org.Dr.eg.db)

gsea_formatting_drerio <-
  function(results_table, tissue, treatment, sampling_point) {
    # Install and load required packages
    required_packages <-
      c('dplyr', 'gprofiler2', 'clusterProfiler', 'org.Hs.eg.db')
    installed_packages <- rownames(installed.packages())
    
    for (pkg in required_packages) {
      if (!(pkg %in% installed_packages)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
    
    # Convert rownames to column 'ensembl'
    results_df <-
      tibble::rownames_to_column(as.data.frame(results_table), var = 'ensembl')
    
    # Convert salmon genes to human orthologs
    orthologs <- gorth(
      query = rownames(results_table),
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
      bitr(ordered_df$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Dr.eg.db)
    entrez_genes <-
      ordered_df %>% left_join(
        ordered_entrez,
        by = c('ortholog_name' = 'SYMBOL'),
        relationship = 'many-to-many'
      ) %>% dplyr::select(ENTREZID, log2FoldChange)
    distinct_genes <-
      entrez_genes %>% distinct(ENTREZID, .keep_all = T)
    entrez_gene_list <<- distinct_genes$log2FoldChange
    names(entrez_gene_list) <<- distinct_genes$ENTREZID
    
    # Run GSEA
    gsea_results <- gseGO(
      gene_list,
      keyType = 'SYMBOL',
      OrgDb = org.Dr.eg.db,
      ont = 'BP',
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      verbose = T,
      eps = 1e-300,
      nPermSimple = 10000
    )
    
    # Assign the results to a variable including treatment and sampling_point in the name
    results_name <-
      paste0(tissue, '_', 'gsea_results_', treatment, '_', sampling_point)
    assign(results_name, gsea_results, envir = .GlobalEnv)
    
    return(gsea_results)
  }


# gsea_formatting(heart_res_eomes_vs_conu_6wpc, 'heart', 'EOMESconu', '6wpc')
# 
# heart_gsea_simplified_EOMESconu6wpc <-
#   clusterProfiler::simplify(heart_gsea_results_EOMESconu_6wpc)  # simplifying GO terms to reduce redundancy
# 
# heart_gene_list_EOMESconu_6wpc <- entrez_gene_list
# 
# nrow(heart_gsea_results_EOMESconu_6wpc)  # 917 GO terms/pathways
# nrow(heart_gsea_simplified_EOMESconu6wpc)  # 316 GO terms/pathways

gsea_formatting_drerio(heart_res_eomes_vs_conu_6wpc, 'heart', 'EOMESconu_drerio', '6wpc')

heart_gsea_simplified_EOMESconu_drerio_6wpc <-
  clusterProfiler::simplify(heart_gsea_results_EOMESconu_drerio_6wpc)  # simplifying GO terms to reduce redundancy


as_tibble(heart_gsea_results_EOMESconu_drerio_6wpc) %>% arrange(NES) %>% print(n = Inf)
as_tibble(heart_gsea_simplified_EOMESconu_drerio_6wpc) %>% arrange(NES) %>% print(n = Inf)

nrow(heart_gsea_simplified_EOMESconu6wpc)  # 917 GO terms/pathways
nrow(heart_gsea_simplified_EOMESconu_drerio_6wpc)  # 316 GO terms/pathways



rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE



install.packages('VennDiagram')
library(VennDiagram)

# Helper function to display Venn diagram
display_venn <- function(x, ...) {
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}




# Load results files
setwd(
  '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables'
)
results_files <-
  list.files(pattern = '^spleen_.*_conu_6wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# applying data_wrangling function to all DESeqResults tables in the Global Environment
list2env(setNames(
  lapply(ls()[sapply(ls(), 
                     function(x) inherits(get(x), 'DESeqResults'))], 
         function(x) zebrafish_data_wrangling(get(x))), 
  paste0(ls()[sapply(ls(), 
                     function(x) inherits(get(x), 'DESeqResults'))], '_wrangled')), 
         envir = .GlobalEnv)

rm(list=ls(pattern=".*6wpc$"))  # removing objects ending in 6wpc

## downregulated
a <- list(
  A = spleen_res_eomes_vs_conu_6wpc_wrangled[spleen_res_eomes_vs_conu_6wpc_wrangled$log2FC < 0, ]$ID,
  B = spleen_res_gata3_vs_conu_6wpc_wrangled[spleen_res_gata3_vs_conu_6wpc_wrangled$log2FC < 0, ]$ID,
  C = spleen_res_ivhd_vs_conu_6wpc_wrangled[spleen_res_ivhd_vs_conu_6wpc_wrangled$log2FC < 0, ]$ID,
  D = spleen_res_ivld_vs_conu_6wpc_wrangled[spleen_res_ivld_vs_conu_6wpc_wrangled$log2FC < 0, ]$ID
)

names(a) <-
  c('EOMES', 'GATA3', 'IV-HD', 'IV-LD')

display_venn(
  a,
  fill = c('#cdb4db', '#bde0fe', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.10, 0.10),
  cat.pos = c(360, 10, 340, 400)
)

## upregulated
b <- list(
  A = spleen_res_eomes_vs_conu_6wpc_wrangled[spleen_res_eomes_vs_conu_6wpc_wrangled$log2FC > 0, ]$ID,
  B = spleen_res_gata3_vs_conu_6wpc_wrangled[spleen_res_gata3_vs_conu_6wpc_wrangled$log2FC > 0, ]$ID,
  C = spleen_res_ivhd_vs_conu_6wpc_wrangled[spleen_res_ivhd_vs_conu_6wpc_wrangled$log2FC > 0, ]$ID,
  D = spleen_res_ivld_vs_conu_6wpc_wrangled[spleen_res_ivld_vs_conu_6wpc_wrangled$log2FC > 0, ]$ID
)

names(b) <-
  c('EOMES', 'GATA3', 'IV-HD', 'IV-LD')

display_venn(
  b,
  fill = c('#cdb4db', '#bde0fe', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.10, 0.10),
  cat.pos = c(360, 10, 340, 400)
)



setwd(
  '/Users/ffi007/Library/CloudStorage/OneDrive-UiTOffice365/Documents/PhD/Papers/Paper III/data/results_tables'
)
results_files <-
  list.files(pattern = '^spleen_.*_6v10')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# applying data_wrangling function to all DESeqResults tables in the Global Environment
list2env(setNames(
  lapply(ls()[sapply(ls(), 
                     function(x) inherits(get(x), 'DESeqResults'))], 
         function(x) zebrafish_data_wrangling(get(x))), 
  paste0(ls()[sapply(ls(), 
                     function(x) inherits(get(x), 'DESeqResults'))], '_wrangled')), 
  envir = .GlobalEnv)

rm(list=ls(pattern=".*6v10$"))  # removing objects ending in 6wpc

## downregulated
a <- list(
  A = spleen_res_conu_vs_conu_6v10_wrangled[spleen_res_conu_vs_conu_6v10_wrangled$log2FC < 0, ]$ID,
  B = spleen_res_eomes_vs_eomes_6v10_wrangled[spleen_res_eomes_vs_eomes_6v10_wrangled$log2FC < 0, ]$ID,
  C = spleen_res_gata3_vs_gata3_6v10_wrangled[spleen_res_gata3_vs_gata3_6v10_wrangled$log2FC < 0, ]$ID,
  D = spleen_res_ivhd_vs_ivhd_6v10_wrangled[spleen_res_ivhd_vs_ivhd_6v10_wrangled$log2FC < 0, ]$ID,
  E = spleen_res_ivld_vs_ivld_6v10_wrangled[spleen_res_ivld_vs_ivld_6v10_wrangled$log2FC < 0, ]$ID
)

names(a) <-
  c('CONU', 'EOMES', 'GATA3', 'IV-HD', 'IV-LD')

display_venn(
  a,
  fill = c('#cdb4db', '#bde0fe', '#ccd5ae', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.22, 0.20,.20),
  cat.pos = c(360, 360, 250, 90, 360)
)

## upregulated
b <- list(
  A = spleen_res_conu_vs_conu_6v10_wrangled[spleen_res_conu_vs_conu_6v10_wrangled$log2FC > 0, ]$ID,
  B = spleen_res_eomes_vs_eomes_6v10_wrangled[spleen_res_eomes_vs_eomes_6v10_wrangled$log2FC > 0, ]$ID,
  C = spleen_res_gata3_vs_gata3_6v10_wrangled[spleen_res_gata3_vs_gata3_6v10_wrangled$log2FC > 0, ]$ID,
  D = spleen_res_ivhd_vs_ivhd_6v10_wrangled[spleen_res_ivhd_vs_ivhd_6v10_wrangled$log2FC > 0, ]$ID,
  E = spleen_res_ivld_vs_ivld_6v10_wrangled[spleen_res_ivld_vs_ivld_6v10_wrangled$log2FC > 0, ]$ID
)

names(b) <-
  c('CONU', 'EOMES', 'GATA3', 'IV-HD', 'IV-LD')

display_venn(
  b,
  fill = c('#cdb4db', '#bde0fe', '#ccd5ae', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.22, 0.20,.20),
  cat.pos = c(360, 360, 250, 90, 360)
)

renv::snapshot()









## EOMES vs CONU @ 4 wpc - SPLEEN ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/results_tables'
)

results_files <-
  list.files(pattern = '^spleen_.*_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

gsea_formatting(spleen_res_eomes_vs_conu_4wpc, 'spleen', 'EOMESconu', '4wpc')

spleen_gsea_simplified_EOMESconu4wpc <-
  clusterProfiler::simplify(spleen_gsea_results_EOMESconu_4wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_EOMESconu_4wpc <- entrez_gene_list

nrow(spleen_gsea_results_EOMESconu_4wpc)  # 51 GO terms/pathways
nrow(spleen_gsea_simplified_EOMESconu4wpc)  # 29 GO terms/pathways

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_EOMESconu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

as_tibble(top10_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, core_enrichment, Count)

## GATA3 vs CONU @ 4 wpc - SPLEEN ----
gsea_formatting(spleen_res_gata3_vs_conu_4wpc, 'spleen', 'GATA3conu', '4wpc')

spleen_gsea_simplified_GATA3conu4wpc <-
  clusterProfiler::simplify(spleen_gsea_results_GATA3conu_4wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_GATA3conu_4wpc <- entrez_gene_list

nrow(spleen_gsea_results_GATA3conu_4wpc)  # 80 GO terms/pathways
nrow(spleen_gsea_simplified_GATA3conu4wpc)  # 47 GO terms/pathways

as_tibble(spleen_gsea_simplified_GATA3conu4wpc)

top10_high_nes <-
  as_tibble(spleen_gsea_simplified_GATA3conu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

as_tibble(top10_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, core_enrichment, Count)

## IV-LD vs CONU @ 4 wpc - SPLEEN ----
gsea_formatting(spleen_res_ivld_vs_conu_4wpc, 'spleen', 'IVLDconu', '4wpc')

spleen_gsea_simplified_IVLDconu4wpc <-
  clusterProfiler::simplify(spleen_gsea_results_IVLDconu_4wpc)  # simplifying GO terms to reduce redundancy

spleen_gene_list_IVLDconu_4wpc <- entrez_gene_list

nrow(spleen_gsea_results_IVLDconu_4wpc)  # 19 enriched terms/pathways
nrow(spleen_gsea_simplified_IVLDconu4wpc)  # 14 enriched terms/pathways

top20_high_nes <-
  as_tibble(spleen_gsea_simplified_IVLDconu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(20, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

as_tibble(top20_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '^spleen*.*_gsea_.*4wpc.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '^spleen_gene_list_*.*4wpc')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

## EOMES vs CONU @ 4 wpc - HEART ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/results_tables'
)

results_files <-
  list.files(pattern = '^heart_.*_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

gsea_formatting(heart_res_eomes_vs_conu_4wpc, 'heart', 'EOMESconu', '4wpc')

heart_gsea_simplified_EOMESconu4wpc <-
  clusterProfiler::simplify(heart_gsea_results_EOMESconu_4wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_EOMESconu_4wpc <- entrez_gene_list

nrow(heart_gsea_results_EOMESconu_4wpc)  # 322 GO terms/pathways
nrow(heart_gsea_simplified_EOMESconu4wpc)  # 162 GO terms/pathways

as_tibble(heart_gsea_simplified_EOMESconu4wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_EOMESconu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-  
  as_tibble(heart_gsea_simplified_EOMESconu4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # 1 downregulated pathway

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)

## GATA3 vs CONU @ 4 wpc - HEART ----
gsea_formatting(heart_res_gata3_vs_conu_4wpc, 'heart', 'GATA3conu', '4wpc')

heart_gsea_simplified_GATA3conu4wpc <-
  clusterProfiler::simplify(heart_gsea_results_GATA3conu_4wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_GATA3conu_4wpc <- entrez_gene_list

nrow(heart_gsea_results_GATA3conu_4wpc)  # 949 GO terms/pathways
nrow(heart_gsea_simplified_GATA3conu4wpc)  # 331 GO terms/pathways

as_tibble(heart_gsea_simplified_GATA3conu4wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_GATA3conu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_GATA3conu4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)


## IV-LD vs CONU @ 4 wpc - HEART ----
gsea_formatting(heart_res_ivld_vs_conu_4wpc, 'heart', 'IVLDconu', '4wpc')

heart_gsea_simplified_IVLDconu4wpc <-
  clusterProfiler::simplify(heart_gsea_results_IVLDconu_4wpc)  # simplifying GO terms to reduce redundancy

heart_gene_list_IVLDconu_4wpc <- entrez_gene_list

nrow(heart_gsea_results_IVLDconu_4wpc)  # 1469 enriched terms/pathways
nrow(heart_gsea_simplified_IVLDconu4wpc)  # 428 enriched terms/pathways

as_tibble(heart_gsea_simplified_IVLDconu4wpc)

top10_high_nes <-
  as_tibble(heart_gsea_simplified_IVLDconu4wpc) %>%
  filter(NES > 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = NES) %>%
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))

bottom10_low_nes <-
  as_tibble(heart_gsea_simplified_IVLDconu4wpc) %>%
  filter(NES < 0) %>%
  arrange(desc(setSize)) %>%
  top_n(10, wt = setSize) %>% 
  mutate(Count = sapply(strsplit(as.character(core_enrichment), '/'), length))  # no downregulated terms

low_high_nes <-
  bind_rows(top10_high_nes, bottom10_low_nes)

as_tibble(low_high_nes) %>% dplyr::select(ID, Description, NES, p.adjust, Count)


# Saving results files ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '^heart*.*_gsea_.*4wpc.*')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  

# Saving gene lists ----
setwd(
  '~/Documents/PhD/Papers/Paper III/data/gsea_results/4 wpc'
)

obj <- ls(pattern = '^heart_gene_list_*.*4wpc')  # regex pattern matching files containing gsea
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}  


### Removing objects  ----
obj_to_remove <- ls(pattern = 'gsea_results|gsea_simplified')
gene_list_to_remove <- ls(pattern = 'gene_list*.*4wpc')
rm(list = obj_to_remove)
rm(list = gene_list_to_remove)





### WGCNA ----
install.packages('WGCNA')
library(WGCNA)
library(DESeq2)


normalized_counts <- counts(ddsDGE_group_heart, normalized=TRUE)

str(sampleTable_heart_group)























