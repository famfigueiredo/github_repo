library(knitr)
library(kableExtra)
library(tidyverse)
library(gprofiler2)

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

### spleen ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_4wpc'
)

results_files <-
  list.files(pattern = '^spleen_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(spleen_res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')
improved_data_wrangling(spleen_res_eomes_vs_conu_4wpc, 'eomes', '4wpc')
improved_data_wrangling(spleen_res_gata3_vs_conu_4wpc, 'gata3', '4wpc')
improved_data_wrangling(spleen_res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')
improved_data_wrangling(spleen_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')

rm(list  = ls(pattern = '^spleen'))

# downregulated
b <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC < 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC < 0, ]$ID,
  C = results_gata3_4wpc[results_gata3_4wpc$log2FC < 0, ]$ID,
  D = results_ivhd_4wpc[results_ivhd_4wpc$log2FC < 0, ]$ID,
  E = results_ivld_4wpc[results_ivld_4wpc$log2FC < 0, ]$ID
)

# add treatment names
names(b) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
b_clean <- lapply(b, function(x) x[!is.na(x)])
b <- b_clean

png('venn_downregulated_spleen.png', width = 800, height = 800, res = 100)
display_venn(
  b_clean,
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
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(b, length)), col.names = c('count'))

# common genes
EOMESvGATA3 <- intersect(b$EOMES, b$GATA3)

unique_EOMESvGATA3 <- setdiff(EOMESvGATA3,
                              c(b$`DNA vaccine`, b$IVHD, b$IVLD))


downregulated_EOMESvGATA3 <- gorth(
  query = unique_EOMESvGATA3,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs GATA3' )


EOMESvDNAvacc <- intersect(b$EOMES, b$`DNA vaccine`)

unique_EOMESvDNAvacc <- setdiff(EOMESvDNAvacc,
                                c(b$GATA3, b$IVHD, b$IVLD))

downregulated_EOMESvDNAvacc <- gorth(
  query = unique_EOMESvDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs DNA vaccine' )


DNAvaccvIVLD <- intersect(b$`DNA vaccine`, b$IVLD)

unique_DNAvaccvIVLD <- setdiff(DNAvaccvIVLD,
                               c(b$GATA3, b$IVHD, b$EOMES))

downregulated_DNAvaccvIVLD <- gorth(
  query = unique_DNAvaccvIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs IV-LD' )

downregulated_common_genes <- bind_rows(downregulated_EOMESvGATA3, downregulated_EOMESvDNAvacc, downregulated_DNAvaccvIVLD)

table_md <- downregulated_common_genes %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_4wpc/venn_diagrams/downregulated_common_genes.md')

na_orthologs_vector <- as.vector(downregulated_common_genes$ssalar_ensembl[downregulated_common_genes$hsapiens_ortholog == "N/A"])

# Query for InterPro domains associated with the gene with not hsapiens ortholog
ensembl <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro_description', 'description'),
  filters = 'ensembl_gene_id',
  values = na_orthologs_vector,
  mart = ensembl
)

print(interpro_results)

paralog_results <- getBM(
  attributes = c('ensembl_gene_id', 'ssalar_paralog_canonical_transcript_protein', 'ssalar_paralog_ensembl_gene', 'description'),
  filters = 'ensembl_gene_id',
  values = na_orthologs_vector,
  mart = ensembl
)

# List available attributes for this dataset
attributes <- listAttributes(ensembl)

# Search for UniProt-related attributes
attributes[grep("interpro", attributes$name), ]


# upregulated
c <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC > 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC > 0, ]$ID,
  C = results_gata3_4wpc[results_gata3_4wpc$log2FC > 0, ]$ID,
  D = results_ivhd_4wpc[results_ivhd_4wpc$log2FC > 0, ]$ID,
  E = results_ivld_4wpc[results_ivld_4wpc$log2FC > 0, ]$ID
)

# add treatment names
names(c) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
c_clean <- lapply(c, function(x) x[!is.na(x)])
c <- c_clean

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

png('venn_upregulated_spleen.png', width = 800, height = 800, res = 100)
display_venn(
  c,
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
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

# common genes
GATA3vDNAvacc <- intersect(c$`DNA vaccine`, c$GATA3)

unique_GATA3vDNAvacc <- setdiff(GATA3vDNAvacc,
                                c(c$EOMES, c$IVHD, c$IVLD))


upregulated_GATA3vDNAvacc <- gorth(
  query = unique_GATA3vDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs DNA vaccine' )


GATA3vIVLD <- intersect(c$GATA3, c$IVLD)

unique_GATA3vIVLD <- setdiff(GATA3vIVLD,
                             c(c$`DNA vaccine`, c$EOMES, c$IVHD))

upregulated_GATA3vIVLD <- gorth(
  query = unique_GATA3vIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs IV-LD' )


EOMESvDNAvacc <- intersect(c$`DNA vaccine`, c$EOMES)

unique_EOMESvDNAvacc <- setdiff(EOMESvDNAvacc,
                               c(c$GATA3, c$IVHD, c$IVLD))

upregulated_EOMESvDNAvacc <- gorth(
  query = unique_EOMESvDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs EOMES')



EOMESvGATA3 <- intersect(c$GATA3, c$EOMES)

unique_EOMESvGATA3 <- setdiff(EOMESvGATA3,
                                c(c$`DNA vaccine`, c$IVHD, c$IVLD))

upregulated_EOMESvGATA3 <- gorth(
  query = unique_EOMESvGATA3,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs EOMES')




upregulated_common_genes <- bind_rows(upregulated_GATA3vDNAvacc, upregulated_GATA3vIVLD, upregulated_EOMESvDNAvacc, upregulated_EOMESvGATA3)


na_orthologs_vector_upregulated <- as.vector(upregulated_common_genes$ssalar_ensembl[upregulated_common_genes$hsapiens_ortholog == "N/A"])

# Query for InterPro domains associated with the gene with not hsapiens ortholog
ensembl <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro_description', 'description'),
  filters = 'ensembl_gene_id',
  values = na_orthologs_vector_upregulated,
  mart = ensembl
)

print(interpro_results)

unique <- interpro_results[!duplicated(interpro_results$ensembl_gene_id), ]

# Rename the column in 'unique' to match 'upregulated_common_genes'
unique_renamed <- unique %>%
  rename(ssalar_ensembl = ensembl_gene_id)

# Merge the data frames
merged_df <- upregulated_common_genes %>%
  left_join(unique_renamed %>% select(ssalar_ensembl, interpro_description, description), 
            by = "ssalar_ensembl")

# Update the 'description' column
merged_df <- merged_df %>%
  mutate(description = case_when(
    !is.na(description.x) & description.x != "N/A" ~ description.x,
    is.na(interpro_description) ~ "novel gene: lncRNA",
    TRUE ~ paste("InterPro:", interpro_description)
  )) %>%
  select(ssalar_ensembl, hsapiens_ortholog, hsapiens_ensembl, description, intersection)

# View the updated merged data frame
print(merged_df)


table_md <- merged_df %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_4wpc/venn_diagrams/upregulated_common_genes.md')


### heart ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(heart_res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')
improved_data_wrangling(heart_res_eomes_vs_conu_4wpc, 'eomes', '4wpc')
improved_data_wrangling(heart_res_gata3_vs_conu_4wpc, 'gata3', '4wpc')
improved_data_wrangling(heart_res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')
improved_data_wrangling(heart_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')

rm(list  = ls(pattern = '^heart_'))

# downregulated ----
b <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC < 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC < 0, ]$ID,
  C = results_gata3_4wpc[results_gata3_4wpc$log2FC < 0, ]$ID,
  D = results_ivhd_4wpc[results_ivhd_4wpc$log2FC < 0, ]$ID,
  E = results_ivld_4wpc[results_ivld_4wpc$log2FC < 0, ]$ID
)

# add treatment names
names(b) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
b_clean <- lapply(b, function(x) x[!is.na(x)])
b <- b_clean

png('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/venn_diagrams/venn_downregulated_heart.png', width = 800, height = 800, res = 100)
display_venn(
  b_clean,
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
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(b, length)), col.names = c('count'))

# common genes
DNAvaccvIVLD <- intersect(b$`DNA vaccine`, b$IVLD)

unique_DNAvaccvIVLD <- setdiff(DNAvaccvIVLD,
                               c(b$GATA3, b$IVHD, b$EOMES))

downregulated_DNAvaccvIVLD <- gorth(
  query = unique_DNAvaccvIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs IV-LD' )


DNAvaccvEOMES <- intersect(b$`DNA vaccine`, b$EOMES)

unique_DNAvaccvEOMES <- setdiff(DNAvaccvEOMES,
                               c(b$GATA3, b$IVHD, b$IVLD))

downregulated_DNAvaccvEOMES <- gorth(
  query = unique_DNAvaccvEOMES,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs EOMES')


IVLDvEOMES <- intersect(b$IVLD, b$EOMES)

unique_IVLDvEOMES <- setdiff(IVLDvEOMES,
                                c(b$GATA3, b$IVHD, b$`DNA vaccine`))

downregulated_IVLDvEOMES <- gorth(
  query = unique_IVLDvEOMES,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'IV-LD vs EOMES')




DNAvaccvGATA3 <- intersect(b$`DNA vaccine`, b$GATA3)

unique_DNAvaccvGATA3 <- setdiff(DNAvaccvGATA3,
                                c(b$EOMES, b$IVHD, b$IVLD))

downregulated_DNAvaccvGATA3 <- gorth(
  query = unique_DNAvaccvGATA3,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs GATA3')



IVLDvGATA3 <- intersect(b$IVLD, b$GATA3)

unique_IVLDvGATA3 <- setdiff(IVLDvGATA3,
                             c(b$EOMES, b$IVHD, b$`DNA vaccine`))

downregulated_IVLDvGATA3 <- gorth(
  query = unique_IVLDvGATA3,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'IV-LD vs GATA3')


GATA3vEOMES <- intersect(b$EOMES, b$GATA3)

unique_GATA3vEOMES <- setdiff(GATA3vEOMES,
                             c(b$IVLD, b$IVHD, b$`DNA vaccine`))

downregulated_GATA3vEOMES <- gorth(
  query = unique_GATA3vEOMES,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs GATA3')


downregulated_common_genes_heart <- bind_rows(downregulated_DNAvaccvEOMES, 
                                              downregulated_DNAvaccvGATA3, 
                                              downregulated_DNAvaccvIVLD, 
                                              downregulated_GATA3vEOMES, 
                                              downregulated_IVLDvEOMES, 
                                              downregulated_IVLDvGATA3)

c

# View the updated merged data frame
print(merged_df)

table_md <- merged_df %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/venn_diagrams/downregulated_common_genes.md')

unique_GATA3vEOMES

# upregulated ----
c <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC > 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC > 0, ]$ID,
  C = results_gata3_4wpc[results_gata3_4wpc$log2FC > 0, ]$ID,
  D = results_ivhd_4wpc[results_ivhd_4wpc$log2FC > 0, ]$ID,
  E = results_ivld_4wpc[results_ivld_4wpc$log2FC > 0, ]$ID
)

# add treatment names
names(c) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
c_clean <- lapply(c, function(x) x[!is.na(x)])
c <- c_clean

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

png('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/venn_diagrams/venn_upregulated_heart.png', width = 800, height = 800, res = 100)
display_venn(
  c,
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
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

# common genes
EOMES_GATA3 <- intersect(c$EOMES, c$GATA3)

unique_EOMES_GATA3 <- setdiff(EOMES_GATA3, 
                              c(c$IVHD, c$`DNA vaccine`, c$IVLD))

upregulated_EOMES_GATA3 <- gorth(
  query = unique_EOMES_GATA3,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs GATA3')


EOMES_IVLD <- intersect(c$EOMES, c$IVLD)

unique_EOMES_IVLD <- setdiff(EOMES_IVLD,
                             c(c$`DNA vaccine`, c$GATA3, c$IVHD))


upregulated_EOMES_IVLD <- gorth(
  query = unique_EOMES_IVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs IV-LD')



GATA3_IVLD <- intersect(c$GATA3, c$IVLD)

unique_GATA3_IVLD <- setdiff(GATA3_IVLD,
                             c(c$`DNA vaccine`, c$EOMES, c$IVHD))

upregulated_GATA3_IVLD <- gorth(
  query = unique_GATA3_IVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs IV-LD')


EOMES_DNAvacc <- intersect(c$EOMES, c$`DNA vaccine`)

unique_EOMES_DNAvacc <- setdiff(EOMES_DNAvacc,
                             c(c$GATA3, c$IVLD, c$IVHD))

upregulated_EOMES_DNAvacc <- gorth(
  query = unique_EOMES_DNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs DNA vaccine')


GATA3_DNAvacc <- intersect(c$GATA3, c$`DNA vaccine`)

unique_GATA3_DNAvacc <- setdiff(GATA3_DNAvacc,
                                c(c$EOMES, c$IVLD, c$IVHD))

upregulated_GATA3_DNAvacc <- gorth(
  query = unique_GATA3_DNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs DNA vaccine')


IVLD_DNAvacc <- intersect(c$IVLD, c$`DNA vaccine`)

unique_IVLD_DNAvacc <- setdiff(IVLD_DNAvacc,
                                c(c$EOMES, c$GATA3, c$IVHD))

upregulated_IVLD_DNAvacc <- gorth(
  query = unique_IVLD_DNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'IV-LD vs DNA vaccine')


upregulated_common_genes <- bind_rows(upregulated_EOMES_DNAvacc, upregulated_EOMES_GATA3, upregulated_EOMES_IVLD, upregulated_GATA3_IVLD, upregulated_GATA3_DNAvacc, upregulated_IVLD_DNAvacc)

na_orthologs_vector_upregulated <- as.vector(upregulated_common_genes$ssalar_ensembl[upregulated_common_genes$hsapiens_ortholog == "N/A"])

# Query for InterPro domains associated with the gene with not hsapiens ortholog
ensembl <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro_description', 'description'),
  filters = 'ensembl_gene_id',
  values = na_orthologs_vector_upregulated,
  mart = ensembl
)

print(interpro_results)

unique <- interpro_results[!duplicated(interpro_results$ensembl_gene_id), ]

# Rename the column in 'unique' to match 'upregulated_common_genes'
unique_renamed <- unique %>%
  rename(ssalar_ensembl = ensembl_gene_id)

# Merge the data frames
merged_df <- upregulated_common_genes %>%
  left_join(unique_renamed %>% select(ssalar_ensembl, interpro_description, description), 
            by = "ssalar_ensembl")

# Update the 'description' column
merged_df <- merged_df %>%
  mutate(description = case_when(
    !is.na(description.x) & description.x != "N/A" ~ description.x,
    is.na(interpro_description) ~ "novel gene: lncRNA",
    TRUE ~ paste("InterPro:", interpro_description)
  )) %>%
  select(ssalar_ensembl, hsapiens_ortholog, hsapiens_ensembl, description, intersection)


table_md <- merged_df %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/venn_diagrams/upregulated_common_genes.md')











### liver ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc'
)

results_files <-
  list.files(pattern = '^liver_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(liver_res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')
improved_data_wrangling(liver_res_eomes_vs_conu_4wpc, 'eomes', '4wpc')
improved_data_wrangling(liver_res_gata3_vs_conu_4wpc, 'gata3', '4wpc')
improved_data_wrangling(liver_res_ivhd_vs_conu_4wpc, 'ivhd', '4wpc')
improved_data_wrangling(liver_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')

rm(list  = ls(pattern = '^liver'))

# downregulated
b <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC < 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC < 0, ]$ID,
  C = results_ivhd_4wpc[results_ivhd_4wpc$log2FC < 0, ]$ID,
  D = results_ivld_4wpc[results_ivld_4wpc$log2FC < 0, ]$ID
)

# add treatment names
names(b) <-
  c('DNA vaccine', 'EOMES', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
b_clean <- lapply(b, function(x) x[!is.na(x)])
b <- b_clean

png('venn_downregulated_liver.png', width = 800, height = 800, res = 100)
display_venn(
  b_clean,
  fill = c('#cdb4db', '#bde0fe', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.10, 0.08),
  cat.pos = c(340, 380, 330, 370)
)
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(b, length)), col.names = c('count'))

# common genes
EOMESvIVLD <- intersect(b$EOMES, b$IVLD)

unique_EOMESvIVLD <- setdiff(EOMESvIVLD,
                             c(b$`DNA vaccine`, b$IVHD))


downregulated_EOMESvIVLD <- gorth(
  query = unique_EOMESvIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs IV-LD')


IVLDvDNAvacc <- intersect(b$`DNA vaccine`, b$IVLD)

unique_IVLDvDNAvacc <- setdiff(IVLDvDNAvacc,
                               c(b$IVHD, b$EOMES))

downregulated_IVLDvDNAvacc <- gorth(
  query = unique_IVLDvDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'IV-LD vs DNA vaccine' )

downregulated_common_genes <- bind_rows(downregulated_EOMESvIVLD, downregulated_IVLDvDNAvacc)

table_md <- downregulated_common_genes %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/venn_diagrams/downregulated_common_genes.md')

# gene IDs
gene_ids <- ('ENSSSAG00000090325') ## unidentified geneIDs from downregulated liver venn diagram

# Query for InterPro domains associated with the gene
interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro', 'interpro_description'),
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)

# Display the results
print(interpro_results)



# upregulated
c <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC > 0, ]$ID,
  B = results_eomes_4wpc[results_eomes_4wpc$log2FC > 0, ]$ID,
  C = results_ivhd_4wpc[results_ivhd_4wpc$log2FC > 0, ]$ID,
  D = results_ivld_4wpc[results_ivld_4wpc$log2FC > 0, ]$ID
)

# add treatment names
names(c) <-
  c('DNA vaccine', 'EOMES', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
c_clean <- lapply(c, function(x) x[!is.na(x)])
c <- c_clean

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

png('venn_upregulated_liver.png', width = 800, height = 800, res = 100)
display_venn(
  c,
  fill = c('#cdb4db', '#bde0fe', '#d4a373', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20, 0.10, 0.08),
  cat.pos = c(340, 380, 330, 370)
)
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

# common genes
EOMESvIVLD <- intersect(c$EOMES, c$IVLD)

unique_EOMESvIVLD <- EOMESvIVLD  # no other overlaps, no need to setdiff

upregulated_EOMESvIVLD <- gorth(
  query = unique_EOMESvIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'EOMES vs IV-LD' )


IVLDvDNAvacc <- intersect(c$`DNA vaccine`, c$IVLD)

unique_IVLDvDNAvacc <- setdiff(IVLDvDNAvacc,
                               c(c$EOMES, c$IVHD))

upregulated_IVLDvDNAvacc <- gorth(
  query = unique_IVLDvDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'IV-LD vs DNA vaccine' )


upregulated_common_genes <- bind_rows(upregulated_EOMESvIVLD, upregulated_IVLDvDNAvacc)

table_md <- upregulated_common_genes %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc/venn_diagrams/upregulated_common_genes.md')


library(biomaRt)

# Select the Ensembl dataset (use your specific species if necessary)
ensembl <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

# gene IDs
gene_ids <- c('ENSSSAG00000003797', 'ENSSSAG00000119965') ## unidentified geneIDs from upregulated liver venn diagram

# Query for InterPro domains associated with the gene
interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro', 'interpro_description'),
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)

# Display the results
print(interpro_results)

### hkidney ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_4wpc')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(hkidney_res_dnavaccine_vs_conu_4wpc, 'dnavaccine', '4wpc')
improved_data_wrangling(hkidney_res_ivld_vs_conu_4wpc, 'ivld', '4wpc')

rm(list  = ls(pattern = '^hkidney'))

# upregulated
c <- list(
  A = results_dnavaccine_4wpc[results_dnavaccine_4wpc$log2FC > 0, ]$ID,
  B = results_ivld_4wpc[results_ivld_4wpc$log2FC > 0, ]$ID
)

# add treatment names
names(c) <-
  c('DNA vaccine', 'IVLD')

# remove NA values from each element in the list 'b'
c_clean <- lapply(c, function(x) x[!is.na(x)])
c <- c_clean

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

png('venn_upregulated_hkidney.png', width = 800, height = 800, res = 100)
display_venn(
  c,
  fill = c('#cdb4db', '#f08080'),
  lwd = 1,
  cex = 1,
  cat.cex = 1,
  cat.fontfamily = 'serif',
  # cat.fontface = 'bold',
  cat.default.pos = 'outer',
  cat.dist = c(0.20, 0.20),
  cat.pos = c(360, 360)
)
dev.off()

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

# common genes
GATA3vDNAvacc <- intersect(c$`DNA vaccine`, c$GATA3)

unique_GATA3vDNAvacc <- setdiff(GATA3vDNAvacc,
                                c(c$EOMES, c$IVHD, c$IVLD))


upregulated_GATA3vDNAvacc <- gorth(
  query = unique_GATA3vDNAvacc,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs DNA vaccine' )


GATA3vIVLD <- intersect(c$GATA3, c$IVLD)

unique_GATA3vIVLD <- setdiff(GATA3vIVLD,
                             c(c$`DNA vaccine`, c$EOMES, c$IVHD))

upregulated_GATA3vIVLD <- gorth(
  query = unique_GATA3vIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'GATA3 vs IV-LD' )


DNAvaccvIVLD <- intersect(c$`DNA vaccine`, c$IVLD)

unique_DNAvaccvIVLD <- setdiff(DNAvaccvIVLD,
                               c(c$GATA3, c$IVHD, c$EOMES))

upregulated_DNAvaccvIVLD <- gorth(
  query = unique_DNAvaccvIVLD,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = F
) %>% dplyr::select(.,
                    input,
                    ortholog_name,
                    ortholog_ensg,
                    description) %>% dplyr::rename(
                      .,
                      ssalar_ensembl = input,
                      hsapiens_ortholog = ortholog_name,
                      hsapiens_ensembl = ortholog_ensg,
                      description = description
                    ) %>% mutate(intersection = 'DNA vaccine vs IV-LD' )



upregulated_common_genes <- bind_rows(upregulated_GATA3vDNAvacc, upregulated_GATA3vIVLD, upregulated_DNAvaccvIVLD)

table_md <- upregulated_common_genes %>%
  kableExtra::kable(
    booktabs = TRUE,
    col.names = c(
      'Salmon ENSEMBL',
      'Human ortholog',
      'Human ENSEMBL',
      'Description',
      'Intersection'
    ),
    align = 'c'
  ) %>%
  kableExtra::row_spec(., row = 0, italic = TRUE) %>%
  kableExtra::kable_styling(font_size = 14) %>% 
  kable_styling(full_width = F)

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc/venn_diagrams/upregulated_common_genes.md')


