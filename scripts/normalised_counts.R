BiocManager::install("tximport")
library(tximport)
BiocManager::install("DESeq2")
library(DESeq2)
install.packages('tidyverse')
BiocManager::install("apeglm")
BiocManager::install("org.Hs.eg.db")


# Loading packages ####
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

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_heart.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_group_heart.Rda')


# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = countData,
                              colData = colData,
                              design = ~ 1)  # No condition variable


dds <- DESeqDataSetFromHTSeqCount(
  sampleTable = sampleTable_group_ensembl,
  directory = directory,
  design = ~ 1
)

dds <- estimateSizeFactors(dds)
normalized_counts <- counts(dds, normalized = TRUE)

filtered_counts <- normalized_counts[rowSums(normalized_counts) > 0, ]

filtered_counts['ENSSSAG00000062937', ]



filtered_counts_df <- rownames_to_column(as.data.frame(filtered_counts), var = "gene")
head(filtered_counts_df)


BiocManager::install("ReactomePA")
library(ReactomePA)

orth_hs <- gorth(
  query = filtered_counts_df$gene,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

BiocManager::install("biomaRt")
library('biomaRt')

# Connect to the Ensembl BioMart
mart <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

# Query BioMart to get gene names corresponding to Ensembl IDs
gene_info <- getBM(attributes = c("ensembl_gene_id", "external_gene_name"),  # Attributes you want
                   filters = "ensembl_gene_id",                              # Filtering by Ensembl gene IDs
                   values = filtered_counts_df$gene,                                     # List of Ensembl gene IDs
                   mart = mart)                                              # BioMart dataset

gene_info %>% filter(grepl('ENSSSAG00000062937', ensembl_gene_id))  # ef1a
filtered_counts_df %>% filter(grepl('ENSSSAG00000062937', gene))  # ef1a
filtered_counts_df %>% filter(grepl('ENSSSAG00000015324', gene))  # tbx21
filtered_counts_df %>% filter(grepl('ENSSSAG00000065097', gene))  # gata3
filtered_counts_df %>% filter(grepl('ENSSSAG00000058121', gene))  # irf1
filtered_counts_df %>% filter(grepl('ENSSSAG00000040910', gene))  # tlr3
filtered_counts_df %>% filter(grepl('ENSSSAG00000076873', gene))  # tlr9
filtered_counts_df %>% filter(grepl('ENSSSAG00000119673', gene))  # ddx58


colnames(filtered_counts_df)

conu_10wpc_h_COUNTS <- filtered_counts_df %>% 
  dplyr::select(., matches('gene|conu_10wpc_h')) %>% 
  mutate(average = rowMeans(across(contains("conu_10wpc_h")))) %>% 
  dplyr::select(gene, average) %>% 
  format(.$average, digits = 1, nsmall = 0, scientific = F)

conu_10wpi_h_COUNTS <- filtered_counts_df %>% 
  dplyr::select(., matches('gene|conu_10wpi_h')) %>% 
  mutate(average = rowMeans(across(contains("conu_10wpi_h")))) %>% 
  dplyr::select(gene, average) %>% 
  format(.$average, digits = 1, nsmall = 0, scientific = F)

conu_1wpc_h_COUNTS <- filtered_counts_df %>% 
  dplyr::select(., matches('gene|conu_1wpc_h')) %>% 
  mutate(average = rowMeans(across(contains("conu_1wpc_h")))) %>% 
  dplyr::select(gene, average) %>% 
  format(.$average, digits = 1, nsmall = 0, scientific = F)

conu_1wpc_h_COUNTS <- filtered_counts_df %>% 
  dplyr::select(., matches('gene|conu_1wpc_h')) %>% 
  mutate(average = rowMeans(across(contains("conu_1wpc_h")))) %>% 
  dplyr::select(gene, average) %>% 
  format(.$average, digits = 1, nsmall = 0, scientific = F)


groups <- gsub('[0-9]+_', '', colnames(filtered_counts_df)) %>% unique(.)  # removing numerical items from tissue column
print(groups)



# Define the function
process_group <- function(group_name) {
  # Dynamically select columns that contain the current group_name
  group_counts <- filtered_counts_df %>%
    dplyr::select(gene, matches(group_name)) %>%
    mutate(!!group_name := rowMeans(across(contains(group_name)))) %>%
    dplyr::select(gene, !!group_name) %>%
    mutate(!!group_name := format(!!sym(group_name), digits = 1, nsmall = 0, scientific = FALSE))
  
  return(group_counts)
}

result_list <- lapply(groups, process_group)
final_result <- Reduce(function(x, y) merge(x, y, by = "gene"), result_list)









