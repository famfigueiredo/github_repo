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
  library('ggrepel')
  register(MulticoreParam(10))
})


# directory containing HTSeq count files
directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/full_dataset/readcounts'

# Creating count matrix ####
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
# sampleTable$treatment <-
#   relevel(sampleTable$treatment, ref = 'conu')

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


# Creating DESeqDataSet object and modelling with ~treatment + samplingPoint ####
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~ treatment + samplingPoint + tissue
)

keep <-
  rowSums(counts(dds)) >= 10  # removing low count genes (<10)
dds <-
  dds[keep, ]

as.data.frame(colData(dds))

collapsed_dds <- collapseReplicates(dds,
                   groupby = dds$n,
                   run = dds$lane)

colData(collapsed_dds) 

print(as.data.frame(colData(collapsed_dds)), max = nrow(colData(collapsed_dds)))

options(max.print = 1000)  # Set this to a higher number or Inf to print everything
print(as.data.frame(colData(collapsed_dds)))


ddsDGE_ensembl_fulldataset <-
  DESeq(collapsed_dds, parallel = T)

save(ddsDGE_ensembl_fulldataset, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_ensembl_fulldataset.RData')

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

save(ddsDGE_group_heart, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_heart.RData')
load(file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_heart.RData')
load(file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_spleen.RData')

ddsDGE_group_spleen <-
  DESeq(collapsed_spleen, parallel = T)

save(ddsDGE_group_spleen, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_group_spleen.RData')


## Treatment contrasts within sampling points - HEART - 10wpi ----
### IVLD
heart_res_ivld_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### IVHD
heart_res_ivhd_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### GATA3 
heart_res_gata3_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)

### EOMES
heart_res_eomes_vs_conu_10wpi <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.10wpi', 'conu.10wpi'), type = 'ashr', parallel = T)


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
### 6 wpc vs 10 wpi
heart_res_conu_vs_conu_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'conu.6wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

heart_res_ivld_vs_ivld_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivld.6wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

heart_res_ivhd_vs_ivhd_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'ivhd.6wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

heart_res_eomes_vs_eomes_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'eomes.6wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

heart_res_gata3_vs_gata3_6v10 <- lfcShrink(ddsDGE_group_heart, contrast = c('group', 'gata3.6wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

summary(heart_res_gata3_vs_gata3_6v10)


## Treatment contrasts within treatments ACROSS sampling points - SPLEEN
### 6 wpc vs 10 wpi
spleen_res_conu_vs_conu_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'conu.6wpc', 'conu.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivld_vs_ivld_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivld.6wpc', 'ivld.10wpi'), type = 'ashr', parallel = T)

spleen_res_ivhd_vs_ivhd_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'ivhd.6wpc', 'ivhd.10wpi'), type = 'ashr', parallel = T)

spleen_res_eomes_vs_eomes_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'eomes.6wpc', 'eomes.10wpi'), type = 'ashr', parallel = T)

spleen_res_gata3_vs_gata3_6v10 <- lfcShrink(ddsDGE_group_spleen, contrast = c('group', 'gata3.6wpc', 'gata3.10wpi'), type = 'ashr', parallel = T)

summary(spleen_res_conu_vs_conu_6v10)







