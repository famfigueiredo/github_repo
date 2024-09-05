library(knitr)
library(kableExtra)
library(tidyverse)
library(gprofiler2)

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

## Loading results files ----

### spleen ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi'
)

results_files <-
  list.files(pattern = '^spleen_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(spleen_res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')
improved_data_wrangling(spleen_res_eomes_vs_conu_10wpi, 'eomes', '10wpi')
improved_data_wrangling(spleen_res_gata3_vs_conu_10wpi, 'gata3', '10wpi')
improved_data_wrangling(spleen_res_ivhd_vs_conu_10wpi, 'ivhd', '10wpi')
improved_data_wrangling(spleen_res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

rm(list  = ls(pattern = '^spleen'))

# downregulated
b <- list(
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC < 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC < 0, ]$ID,
  C = results_gata3_10wpi[results_gata3_10wpi$log2FC < 0, ]$ID,
  D = results_ivhd_10wpi[results_ivhd_10wpi$log2FC < 0, ]$ID,
  E = results_ivld_10wpi[results_ivld_10wpi$log2FC < 0, ]$ID
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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/venn_diagrams/downregulated_common_genes.md')

# querying InterPro for genes with no human ortholog

# gene IDs
gene_ids <- c('ENSSSAG00000067877', 'ENSSSAG00000007048', 'ENSSSAG00000084311') ## unidentified geneIDs from downregulated heart venn diagram

# Query for InterPro domains associated with the gene
interpro_results <- getBM(
  attributes = c('ensembl_gene_id', 'interpro', 'interpro_description'),
  filters = 'ensembl_gene_id',
  values = gene_ids,
  mart = ensembl
)

print(interpro_results)



# upregulated
c <- list(
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC > 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC > 0, ]$ID,
  C = results_gata3_10wpi[results_gata3_10wpi$log2FC > 0, ]$ID,
  D = results_ivhd_10wpi[results_ivhd_10wpi$log2FC > 0, ]$ID,
  E = results_ivld_10wpi[results_ivld_10wpi$log2FC > 0, ]$ID
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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/venn_diagrams/upregulated_common_genes.md')



### heart ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')
improved_data_wrangling(res_eomes_vs_conu_10wpi, 'eomes', '10wpi')
improved_data_wrangling(res_gata3_vs_conu_10wpi, 'gata3', '10wpi')
improved_data_wrangling(res_ivhd_vs_conu_10wpi, 'ivhd', '10wpi')
improved_data_wrangling(res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

rm(list  = ls(pattern = '^res_'))

# downregulated
b <- list(
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC < 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC < 0, ]$ID,
  C = results_gata3_10wpi[results_gata3_10wpi$log2FC < 0, ]$ID,
  D = results_ivhd_10wpi[results_ivhd_10wpi$log2FC < 0, ]$ID,
  E = results_ivld_10wpi[results_ivld_10wpi$log2FC < 0, ]$ID
)

# add treatment names
names(b) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
b_clean <- lapply(b, function(x) x[!is.na(x)])
b <- b_clean

png('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/venn_diagrams/venn_downregulated_heart.png', width = 800, height = 800, res = 100)
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

table_md <- downregulated_DNAvaccvIVLD %>%
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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/venn_diagrams/downregulated_common_genes.md')

# Select the Ensembl dataset (use your specific species if necessary)
ensembl <- useMart("ensembl", dataset = "ssalar_gene_ensembl")

# gene IDs
gene_ids <- c('ENSSSAG00000003797') ## unidentified geneIDs from downregulated heart venn diagram

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
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC > 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC > 0, ]$ID,
  C = results_gata3_10wpi[results_gata3_10wpi$log2FC > 0, ]$ID,
  D = results_ivhd_10wpi[results_ivhd_10wpi$log2FC > 0, ]$ID,
  E = results_ivld_10wpi[results_ivld_10wpi$log2FC > 0, ]$ID
)

# add treatment names
names(c) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD', 'IVLD')

# remove NA values from each element in the list 'b'
c_clean <- lapply(c, function(x) x[!is.na(x)])
c <- c_clean

# check gene counts per treatment
kableExtra::kable((sapply(c, length)), col.names = c('count'))

png('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/venn_diagrams/venn_upregulated_heart.png', width = 800, height = 800, res = 100)
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

EOMES_GATA3_IVLD <- intersect(EOMES_GATA3, c$IVLD)
c$GATA3

unique_EOMES_GATA3_IVLD <- setdiff(EOMES_GATA3_IVLD, c$IVHD)

upregulated_EOMES_GATA3_IVLD <- gorth(
  query = unique_EOMES_GATA3_IVLD,
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
                    ) %>% mutate(intersection = 'EOMES vs GATA3 vs IV-LD' )


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


upregulated_common_genes <- bind_rows(upregulated_EOMES_GATA3_IVLD, upregulated_EOMES_IVLD)

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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/venn_diagrams/upregulated_common_genes.md')



### liver ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi'
)

results_files <-
  list.files(pattern = '^liver_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(liver_res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')
improved_data_wrangling(liver_res_eomes_vs_conu_10wpi, 'eomes', '10wpi')
improved_data_wrangling(liver_res_gata3_vs_conu_10wpi, 'gata3', '10wpi')
improved_data_wrangling(liver_res_ivhd_vs_conu_10wpi, 'ivhd', '10wpi')
improved_data_wrangling(liver_res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

rm(list  = ls(pattern = '^liver'))

# downregulated
b <- list(
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC < 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC < 0, ]$ID,
  C = results_ivhd_10wpi[results_ivhd_10wpi$log2FC < 0, ]$ID,
  D = results_ivld_10wpi[results_ivld_10wpi$log2FC < 0, ]$ID
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
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC > 0, ]$ID,
  B = results_eomes_10wpi[results_eomes_10wpi$log2FC > 0, ]$ID,
  C = results_ivhd_10wpi[results_ivhd_10wpi$log2FC > 0, ]$ID,
  D = results_ivld_10wpi[results_ivld_10wpi$log2FC > 0, ]$ID
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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/venn_diagrams/upregulated_common_genes.md')


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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_10wpi')  # regex matching results files
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

improved_data_wrangling(hkidney_res_dnavaccine_vs_conu_10wpi, 'dnavaccine', '10wpi')
improved_data_wrangling(hkidney_res_ivld_vs_conu_10wpi, 'ivld', '10wpi')

rm(list  = ls(pattern = '^hkidney'))

# upregulated
c <- list(
  A = results_dnavaccine_10wpi[results_dnavaccine_10wpi$log2FC > 0, ]$ID,
  B = results_ivld_10wpi[results_ivld_10wpi$log2FC > 0, ]$ID
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

writeLines(table_md, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/venn_diagrams/upregulated_common_genes.md')


