# Loading packages ####
suppressPackageStartupMessages({
  library('tidyverse')
  library('apeglm')
  library('DESeq2')
  library('LSD')
  library('BiocParallel')
  library('ExploreModelMatrix')
  library('limma')
  library('ashr')
  library('ggrepel')
  library('ggpmisc')
  library('DataCombine')
  library('org.Hs.eg.db')
  library('ggrepel')
  register(MulticoreParam(8))
})


# directory containing HTSeq count files
directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/readcounts/ensembl_htseq-count/non-heart-ensembl-reads'


# Creating count matrix ####
sampleFiles <- grep('readcount', list.files(directory), value = T)

sampleTable <-
  data.frame(
    n = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[1]),
    filename = sampleFiles,
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

sampleTable$tissue <-
  gsub('[0-9]+', '', sampleTable$tissue)  # removing numerical items from tissue column

sampleTable$treatment <-
  gsub('-', '', sampleTable$treatment) # removing dashes from treatment column

# defining factors
sampleTable$treatment <- factor(sampleTable$treatment)
sampleTable$samplingPoint <- factor(sampleTable$samplingPoint)
sampleTable$tissue <- factor(sampleTable$tissue)
sampleTable$lane <- factor(sampleTable$lane)
sampleTable$treatment <-
  relevel(sampleTable$treatment, ref = 'conu')

sampleTable$sample <-
  paste(
    rownames(sampleTable),
    as.character(sampleTable$treatment),
    as.character(sampleTable$samplingPoint),
    as.character(sampleTable$tissue),
    sep = '_'
  )

sampleTable <-
  as.data.frame(sampleTable) %>% dplyr::select(sample, filename, treatment, samplingPoint, tissue, lane)  # selecting columns of interest

sampleTable <- sampleTable %>%
  filter(
    !(
      tissue == 'hkr' |
        tissue == 'hr' |
        tissue == 'sr' |
        tissue == 'h' |
        samplingPoint == '10wpc' |
        treatment == 'ptagrfp' | treatment == 'tbet'
    )
  )
head(sampleTable)
nrow(sampleTable)
summary(sampleTable)

# Removing PCA outliers from sampleTable
sampleTable <- sampleTable %>%
  filter(
    sample != '25_ivhd_6wpc_s' &
      sample != '105_ivhd_4wpc_l' &
      sample != '208_ivld_10wpi_l' &
      sample != '275_ivld_10wpi_hk' &
      sample != '285_eomes_4wpc_hk' &
      sample != '303_gata3_10wpi_hk'
  )

# Creating DESeqDataSet object and modelling with ~treatment + samplingPoint ####
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~ treatment + samplingPoint + tissue
)

keep <-
  rowSums(counts(dds)) >= 10  # removing low count genes (<10)
dds <-
  dds[keep,]

ddsDGE_ensembl_lymphoid <-
  DESeq(dds, parallel = T)

# colData_df <- as.data.frame(colData(ddsDGE_ensembl_lymphoid))
# 
# options(max.print = 99999999)
# colData_df %>% print(n = Inf, na.print = "NA")
# 
# colData(ddsDGE_ensembl_lymphoid) <- as.data.frame(colData(ddsDGE_ensembl_lymphoid)) %>%
#   filter(
#     sample != '25_ivhd_6wpc_s' &
#       sample != '105_ivhd_4wpc_l' &
#       sample != '208_ivld_10wpi_l' &
#       sample != '275_ivld_10wpi_hk' &
#       sample != '285_eomes_4wpc_hk' & 
#       sample != '303_gata3_10wpi_hk'
#   )
# 
# nrow(colData_df)
# nrow(colData_filtered)
# 
# colData(ddsDGE_ensembl_lymphoid) <- colData_filtered

save(ddsDGE_ensembl_lymphoid, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_ensembl_lymphoid.RData')

## Exploratory analysis ####

# Plotting dispersion estimation to ensure that the assumption that most genes are not differential expressed holds
DESeq2::plotDispEsts(ddsDGE_ensembl_lymphoid)

vst_counts <- vst(ddsDGE_ensembl_lymphoid, blind = T)

vst_counts_spleen <- vst_counts[, vst_counts$tissue %in% c('s')]
vst_counts_liver <- vst_counts[, vst_counts$tissue %in% c('l')]
vst_counts_hkidney <- vst_counts[, vst_counts$tissue %in% c('hk')]

PCA <-
  plotPCA(vst_counts_spleen,
          intgroup = c('treatment', 'samplingPoint'),
          returnData = T)

percentVar <- round(100 * attr(PCA, "percentVar"))

## PCA
# options(ggrepel.max.overlaps = Inf)
pca_spleen <- ggplot(PCA, aes(
  x = PC2,
  y = PC1,
  color = treatment,
  shape = samplingPoint
)) +
  geom_point(size = 2) +
  # stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed(ratio = 1) +
  ggtitle("spleen  and treatment/samplingPoint PCA") +
  theme_linedraw(base_size = 10, base_family = 'Times New Roman') +
  theme(plot.margin = margin(10, 10, 50, 10))  # Adjust the bottom margin to increase height
# scale_color_viridis_d(name = 'samplingPoint',
#                       labels = c('Heart',
#                                  'Head-kidney',
#                                  'Liver',
#                                  'Spleen')) +


pca_spleen
pca_liver
pca_hkidney

# Spleen group ----
sampleTable_grouped_spleen <- sampleTable %>%
  filter(tissue == 's') %>% 
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint))  # creating the treatment.samplingPoint grouping
# OR
# sampleTable %>% dplyr::mutate(group = interaction(treatment, samplingPoint)) -> sampleTable

# formatting variables
sampleTable_grouped_spleen$group <-
  as.factor(sampleTable_grouped_spleen$group)
sampleTable_grouped_spleen <-
  droplevels(sampleTable_grouped_spleen)  # dropping the reference

# creating dds object from htseq count output files
ddsGroup_spleen <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_grouped_spleen,
                             directory = directory,
                             design = ~ group)

keep <-
  rowSums(counts(ddsGroup_spleen)) >= 10  # removing low count genes (<10)
ddsGroup_spleen <-
  ddsGroup_spleen[keep, ]


ddsDGE_grouped_spleen <-
  DESeq(ddsGroup_spleen, parallel = T)

save(
  ddsDGE_grouped_spleen,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_spleen.RData'
)  # Saving DESeqDataSet object

save(
  sampleTable_grouped_spleen,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_grouped_spleen.Rda'
)



# Liver group ----
sampleTable_grouped_liver <- sampleTable %>%
  filter(tissue == 'l') %>% 
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint))  # creating the treatment.samplingPoint grouping
# OR
# sampleTable %>% dplyr::mutate(group = interaction(treatment, samplingPoint)) -> sampleTable

summary(sampleTable_grouped_liver)
# formatting variables
sampleTable_grouped_liver$group <-
  as.factor(sampleTable_grouped_liver$group)
sampleTable_grouped_liver <-
  droplevels(sampleTable_grouped_liver)  # dropping the reference

# creating dds object from htseq count output files
ddsGroup_liver <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_grouped_liver,
                             directory = directory,
                             design = ~ group)

keep <-
  rowSums(counts(ddsGroup_liver)) >= 10  # removing low count genes (<10)
ddsGroup_liver <-
  ddsGroup_liver[keep, ]

ddsDGE_grouped_liver <-
  DESeq(ddsGroup_liver, parallel = T)

save(
  ddsDGE_grouped_liver,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_liver.RData'
)  # Saving DESeqDataSet object

save(
  sampleTable_grouped_liver,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_grouped_liver.Rda'
)


# Head-kidney group ----
sampleTable_grouped_hkidney <- sampleTable %>%
  filter(tissue == 'hk') %>% 
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint))  # creating the treatment.samplingPoint grouping
# OR
# sampleTable %>% dplyr::mutate(group = interaction(treatment, samplingPoint)) -> sampleTable

summary(sampleTable_grouped_hkidney)

# formatting variables
sampleTable_grouped_hkidney$group <-
  as.factor(sampleTable_grouped_hkidney$group)
sampleTable_grouped_hkidney <-
  droplevels(sampleTable_grouped_hkidney)  # dropping the reference

# creating dds object from htseq count output files
ddsGroup_hkidney <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_grouped_hkidney,
                             directory = directory,
                             design = ~ group)

keep <-
  rowSums(counts(ddsGroup_hkidney)) >= 10  # removing low count genes (<10)
ddsGroup_hkidney <-
  ddsGroup_hkidney[keep, ]

ddsDGE_grouped_hkidney <-
  DESeq(ddsGroup_hkidney, parallel = T)

save(
  ddsDGE_grouped_hkidney,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_hkidney.RData'
)  # Saving DESeqDataSet object

save(
  sampleTable_grouped_hkidney,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_grouped_hkidney.Rda'
)
