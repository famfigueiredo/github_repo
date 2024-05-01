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
  library('patchwork')
  library('ggrepel')
  library('ggVennDiagram')
  library('ggpmisc')
  library('hrbrthemes')
  library('DataCombine')
  library('ReactomePA')
  library('clusterProfiler')
  library('gprofiler2')
  library('org.Hs.eg.db')
  register(MulticoreParam(10))
})

# Load functions and databases ####
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/functions_data-wrangling_march24.R'
)

# directory containing HTSeq count files
directory <-
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/ensembl_htseq-count'

# Creating count matrix ####
sampleFiles <- grep('readcount', list.files(directory), value = T)

sampleTable <-
  data.frame(
    treatment = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[1]),
    fileName = sampleFiles,
    samplingPoint = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[2]),
    tissue = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[3]),
    lane = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[4]),
    readcount.txt = sapply(strsplit(sampleFiles, '_', fixed = T), function(x)
      x[5])
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
  relevel(sampleTable$treatment, ref = 'ptagrfp')

# adding sample name to rowname
sampleTable$sample <-
  paste(
    rownames(sampleTable),
    as.character(sampleTable$treatment),
    as.character(sampleTable$samplingPoint),
    as.character(sampleTable$tissue),
    sep = '_'
  )

# changing variables to uppercase to facilitate reading
# sampleTable[,3] <- as.factor(toupper(sampleTable[,3]))
# sampleTable[,4] <- as.factor(toupper(sampleTable[,4]))
# sampleTable[,5] <- as.factor(toupper(sampleTable[,5]))
# sampleTable[,6] <- as.factor(toupper(sampleTable[,6]))

sampleTable <-
  as.data.frame(sampleTable) %>% dplyr::select(sample, fileName, treatment, samplingPoint, tissue, lane)  # selecting columns of interest

head(sampleTable)
nrow(sampleTable)
summary(sampleTable)

# finding and removing outlier in PCA #
sampleTable %>% filter(treatment == 'dnavaccine' &
                         samplingPoint == '1wpc')
sampleTable %>% filter(sample == '53_dnavaccine_1wpc_h')

sampleTable <-
  sampleTable %>% filter(.,
                         sample != '53_dnavaccine_1wpc_h')

save(sampleTable, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/sampleTable.Rda')

# Creating DESeqDataSet object and modelling with ~treatment + samplingPoint ####
dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable,
  directory = directory,
  design = ~ treatment + samplingPoint
)

keep <-
  rowSums(counts(dds)) >= 10  # removing low count genes (<10)
dds <-
  dds[keep, ]

ddsDGE_heart_ensembl <-
  DESeq(dds, parallel = T)

save(ddsDGE_heart_ensembl, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/ddsDGE_heart_ensembl.RData')

## Exploratory analysis ####
load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/ddsDGE_heart_ensembl.RData'
)
# Plotting dispersion estimation to ensure that the assumption that most genes are not differential expressed holds
DESeq2::plotDispEsts(ddsDGE_heart_ensembl)

# Transformation to stabilize variance across the mean through *variance stabilizing transformation*
vst_counts <- vst(ddsDGE_heart_ensembl, blind = T)

PCA <-
  plotPCA(
    vst_counts,
    intgroup = c('samplingPoint', 'treatment'),
    returnData = T
  )

percentVar <- round(100 * attr(PCA, "percentVar"))

## PCA
pca <- ggplot(PCA, aes(
  x = PC2,
  y = PC1,
  color = samplingPoint,
  label = treatment
)) +
  geom_point(size = 1) +
  geom_label_repel(
    data = subset(PCA, samplingPoint == '1wpc'),
    aes(PC2, PC1, label = treatment),
    family = 'mono',
    size = 3
  ) +
  stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("Heart tissue, treatment and sampling point PCA") +
  theme_linedraw(base_size = 10, base_family = 'Arial Narrow') +
  # scale_color_viridis_d(name = 'samplingPoint',
  #                       labels = c('Heart',
  #                                  'Head-kidney',
  #                                  'Liver',
  #                                  'Spleen')) +
  theme(plot.margin = grid::unit(c(0, 0, 0, 2), 'mm'),
        panel.grid.minor = element_blank())

options(ggrepel.max.overlaps = Inf)

pca

# Creating DESeqDataSet object and modelling with ~group ####
load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/sampleTable.Rda'
)

sampleTable %>%
  dplyr::mutate(group = paste(treatment, sep = '.', samplingPoint)) -> sampleTable_group_ensembl  # creating the treatment.samplingPoint grouping
# OR
# sampleTable %>% dplyr::mutate(group = interaction(treatment, samplingPoint)) -> sampleTable

# formatting variables
sampleTable_group_ensembl$group <-
  as.factor(sampleTable_group_ensembl$group)
sampleTable_group_ensembl <-
  droplevels(sampleTable_group_ensembl)  # dropping the reference

# creating dds object from htseq count output files
ddsGroup <-
  DESeqDataSetFromHTSeqCount(sampleTable = sampleTable_group_ensembl,
                             directory = directory,
                             design = ~ group)

keep <-
  rowSums(counts(ddsGroup)) >= 10  # removing low count genes (<10)
ddsGroup <-
  ddsGroup[keep, ]


ddsGroup_ensembl <-
  DESeq(ddsGroup, parallel = T)

save(
  ddsGroup_ensembl,
  file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/ddsGroup_ensembl.RData'
)  # Saving DESeqDataSet object

save(sampleTable_group_ensembl,
     file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/sampleTable_group_ensembl.Rda')

## Exploratory analysis, part II ####
vst_counts <- vst(ddsGroupDGE, blind = T)

PCA <-
  plotPCA(vst_counts,
          intgroup = c('group'),
          returnData = T)

percentVar <- round(100 * attr(PCA, "percentVar"))

# Plotting
pca <- ggplot(PCA, aes(x = PC2,
                       y = PC1,
                       color = group)) +
  geom_point(size = 1) +
  stat_ellipse() +
  xlab(paste0("PC2: ", percentVar[1], "% variance")) +
  ylab(paste0("PC1: ", percentVar[2], "% variance")) +
  coord_fixed() +
  ggtitle("Heart tissue, treatment and sampling point PCA") +
  theme_linedraw(base_size = 10, base_family = 'Arial Narrow') +
  theme(plot.margin = grid::unit(c(0, 0, 0, 2), 'mm'),
        panel.grid.minor = element_blank())

options(ggrepel.max.overlaps = Inf)

pca


# Loading model data and extracting results ####
load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/ddsGroup_ensembl.RData'
)

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/sampleTable_group_ensembl.Rda'
)

summary(sampleTable_group_ensembl)
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions


# 1 WPC vs 10 WPI data analysis ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_1wpc_vs_10wpi'
)

## Loading results files ----
results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

# Create an empty list to store results
all_significant_genes <- list()

# Assuming results_files is a list of object names
for (name in results_files) {
  # Remove the file extension ".RData" if present
  name <- sub(".RData$", "", name)
  
  # Retrieve the object from the global environment
  results <- get(name, envir = .GlobalEnv)
  
  # Call the significant_genes function and store the result
  result_genes <- significant_genes(results)
  
  # Store the result in the list
  all_significant_genes[[name]] <- result_genes
}

### Renaming dataframes ----
names(all_significant_genes) <-
  c('CONU',
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD',
    'pTagRFP')

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'CONU',
      'DNA vaccine',
      'EOMES',
      'GATA3',
      'IV-HD',
      'IV-LD',
      'pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

### Reordering factors to plot ----
deg_regulation_summary$Treatment <-
  factor(
    deg_regulation_summary$Treatment,
    levels = c(
      'CONU',
      'pTagRFP',
      'IV-LD',
      'IV-HD',
      'DNA vaccine',
      'EOMES',
      'GATA3'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Set3') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 1 WPC vs 10 WPI') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    linetype = 'dotted',
    size = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2)

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within treatments_across sampling points/DGE_1wpc_vs_10wpi.png', plot = last_plot())

# 4 WPC vs 10 WPI data analysis ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc_vs_10wpi'
)
results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

# Create an empty list to store results
all_significant_genes <- list()

# Assuming results_files is a list of object names
for (name in results_files) {
  # Remove the file extension ".RData" if present
  name <- sub(".RData$", "", name)
  
  # Retrieve the object from the global environment
  results <- get(name, envir = .GlobalEnv)
  
  # Call the significant_genes function and store the result
  result_genes <- significant_genes(results)
  
  # Store the result in the list
  all_significant_genes[[name]] <- result_genes
}

### Renaming dataframes ----
names(all_significant_genes) <-
  c('CONU',
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD',
    'pTagRFP')

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents
### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'CONU',
      'DNA vaccine',
      'EOMES',
      'GATA3',
      'IV-HD',
      'IV-LD',
      'pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

### Reordering factors to plot ----
deg_regulation_summary$Treatment <-
  factor(
    deg_regulation_summary$Treatment,
    levels = c(
      'CONU',
      'pTagRFP',
      'IV-LD',
      'IV-HD',
      'DNA vaccine',
      'EOMES',
      'GATA3'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Set1') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 4 WPC vs 10 WPI') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    linetype = 'dotted',
    size = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2)

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within treatments_across sampling points/DGE_4wpc_vs_10wpi.png', plot = last_plot())

# 6 WPC vs 10 WPI data analysis ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi'
)

## Loading results files ----
results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

# Create an empty list to store results
all_significant_genes <- list()

# Assuming results_files is a list of object names
for (name in results_files) {
  # Remove the file extension ".RData" if present
  name <- sub(".RData$", "", name)
  
  # Retrieve the object from the global environment
  results <- get(name, envir = .GlobalEnv)
  
  # Call the significant_genes function and store the result
  result_genes <- significant_genes(results)
  
  # Store the result in the list
  all_significant_genes[[name]] <- result_genes
}

### Renaming dataframes ----
names(all_significant_genes) <-
  c('CONU',
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD',
    'pTagRFP')

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'CONU',
      'DNA vaccine',
      'EOMES',
      'GATA3',
      'IV-HD',
      'IV-LD',
      'pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))
  
### Reordering factors to plot ----
deg_regulation_summary$Treatment <-
  factor(
    deg_regulation_summary$Treatment,
    levels = c(
      'CONU',
      'pTagRFP',
      'IV-LD',
      'IV-HD',
      'DNA vaccine',
      'EOMES',
      'GATA3'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Set2') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 6 WPC vs 10 WPI') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    linetype = 'dotted',
    size = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2)

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within treatments_across sampling points/DGE_6wpc_vs_10wpi.png', plot = last_plot())


# 10 WPC vs 10 WPI data analysis ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions. It uses setdiff to find the subset of objects in the global environment
# (as returned by ls()) that don't have mode function (as returned by lsf.str())

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi'
)

## Loading results files ----
results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

# Create an empty list to store results
all_significant_genes <- list()

# Assuming results_files is a list of object names
for (name in results_files) {
  # Remove the file extension ".RData" if present
  name <- sub(".RData$", "", name)
  
  # Retrieve the object from the global environment
  results <- get(name, envir = .GlobalEnv)
  
  # Call the significant_genes function and store the result
  result_genes <- significant_genes(results)
  
  # Store the result in the list
  all_significant_genes[[name]] <- result_genes
}

### Renaming dataframes ----
names(all_significant_genes) <-
  c('CONU',
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD',
    'pTagRFP')

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'CONU',
      'DNA vaccine',
      'EOMES',
      'GATA3',
      'IV-HD',
      'IV-LD',
      'pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

### Reordering factors to plot ----
deg_regulation_summary$Treatment <-
  factor(
    deg_regulation_summary$Treatment,
    levels = c(
      'CONU',
      'pTagRFP',
      'IV-LD',
      'IV-HD',
      'DNA vaccine',
      'EOMES',
      'GATA3'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Accent') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 3
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 10 WPC vs 10 WPI') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 5.5, 6.5),
    linetype = 'dotted',
    size = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2)

ggsave('~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within treatments_across sampling points/DGE_10wpc_vs_10wpi.png', plot = last_plot())
