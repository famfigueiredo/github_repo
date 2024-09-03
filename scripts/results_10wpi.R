library('clusterProfiler')
library('tidyverse')
library('pathview')
library('ReactomePA')

# Helper function to create a Markdown table
markdown_table <- function(data) {
  # Get the header
  header <-
    paste("|", paste(names(data), collapse = " | "), "|")
  
  # Get the separator line
  separator <-
    paste("|", paste(rep("---", ncol(data)), collapse = " | "), "|")
  
  # Get the table rows
  rows <- apply(data, 1, function(row) {
    paste("|", paste(row, collapse = " | "), "|")
  })
  
  # Combine header, separator, and rows
  c(header, separator, rows)
}

# Helper function to select top and bottom 10
headtail <- function(x) {
  rbind(head(x, n=10), tail(x, n=10))
}

# Improved helper function to select top and bottom 10
get_top_bottom_pathways <- function(data, n = 10, database) {
  # Get the top N pathways with positive NES
  top_n_pathways <- data %>%
    as_tibble() %>%
    filter(NES > 0) %>%
    arrange(desc(NES)) %>%
    top_n(n, wt = NES) %>%
    dplyr::select(Description, NES)
  
  # Get the bottom N pathways with negative NES
  bottom_n_pathways <- data %>%
    as_tibble() %>%
    filter(NES < 0) %>%
    arrange(NES) %>%
    top_n(n, wt = NES) %>%
    dplyr::select(Description, NES)
  
  # Combine the top and bottom N pathways
  combined_pathways <- bind_rows(top_n_pathways, bottom_n_pathways)
  
  # Create the object name based on the database parameter
  object_name <- paste0(database, "_pathways")
  
  # Assign the combined pathways tibble to a global environment object
  assign(object_name, combined_pathways, envir = .GlobalEnv)
  
  return(combined_pathways)
}

# Loading functions
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

# spleen
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE
## dna vaccine ----
load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_dnavaccine_10wpi.RData'
)
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

kegg_dnavaccine_10wpi_spleen <- gseKEGG(
  spleen_entrez_gene_list_dnavaccine_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

as_tibble(kegg_results_dnavaccine_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)
pathview(
  gene.data = spleen_entrez_gene_list_dnavaccine_10wpi,
  species = 'hsa',
  pathway.id = '04610',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'complement_and_coagulation_cascades'
)

reactome_dnavaccine_10wpi_spleen <-
  gsePathway(
    spleen_entrez_gene_list_dnavaccine_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_dnavaccine_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/dnavaccine_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_dnavaccine_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/dnavaccine_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('TRAF6 mediated IRF7 activation',
            readable = T,
            foldChange = spleen_entrez_gene_list_dnavaccine_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/dnavaccine_traf6_mediated_irf7_activation.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Negative regulation of MAPK pathway',
            readable = T,
            foldChange = spleen_entrez_gene_list_dnavaccine_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/dnavaccine_negative_regulation_mapk_pathway.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)


kegg_pathways <-
  as_tibble(kegg_dnavaccine_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)
reactome_pathways <-
  as_tibble(reactome_dnavaccine_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_dnavaccine.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_dnavaccine.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")


## iv-ld ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivld_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

kegg_ivld_10wpi_spleen <- gseKEGG(
  spleen_entrez_gene_list_ivld_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

as_tibble(kegg_ivld_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = spleen_entrez_gene_list_ivld_10wpi,
  species = 'hsa',
  pathway.id = '04610',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'ivld_complement_and_coagulation_cascades'
)

reactome_ivld_10wpi_spleen <-
  gsePathway(
    spleen_entrez_gene_list_ivld_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )


viewPathway('Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_ivld_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/ivld_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_ivld_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/ivld_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Response to elevated platelet cytosolic Ca2+',
            readable = T,
            foldChange = spleen_entrez_gene_list_ivld_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/ivld_response_elevated_platelet_cytosolic_ca2.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)


kegg_pathways <-
  as_tibble(kegg_ivld_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES)
reactome_pathways <-
  as_tibble(reactome_ivld_10wpi_spleen) %>% dplyr::select(., Description, NES) %>% arrange(NES)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_ivld.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_ivld.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## iv-hd ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_ivhd_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

kegg_ivhd_10wpi_spleen <- gseKEGG(
  spleen_entrez_gene_list_ivhd_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

as_tibble(kegg_ivhd_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = spleen_entrez_gene_list_ivhd_10wpi,
  species = 'hsa',
  pathway.id = '04060',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'ivhd_cytokine-cytokine_receptor_interaction'
)

reactome_ivhd_10wpi_spleen <-
  gsePathway(
    spleen_entrez_gene_list_ivhd_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(reactome_ivhd_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)

viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_ivhd_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/ivhd_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

kegg_pathways <-
  as_tibble(kegg_ivhd_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES)
reactome_pathways <-
  as_tibble(reactome_ivhd_10wpi_spleen) %>% dplyr::select(., Description, NES) %>% arrange(NES)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_ivhd.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_ivhd.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## eomes ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_eomes_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

kegg_eomes_10wpi_spleen <- gseKEGG(
  spleen_entrez_gene_list_eomes_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  nPermSimple = 100000,
  eps = 0
)

as_tibble(kegg_eomes_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)
as_tibble(kegg_eomes_10wpi_spleen)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = spleen_entrez_gene_list_eomes_10wpi,
  species = 'hsa',
  pathway.id = '04610',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'eomes_complement_and_coagulation_cascades'
)

reactome_eomes_10wpi_spleen <-
  gsePathway(
    spleen_entrez_gene_list_eomes_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 1000000,
    verbose = F
  )

as_tibble(reactome_eomes_10wpi_spleen) %>% arrange(NES) %>% print(n = Inf)


viewPathway('Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_eomes_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/eomes_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_eomes_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/eomes_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)


viewPathway('Response to elevated platelet cytosolic Ca2+',
            readable = T,
            foldChange = spleen_entrez_gene_list_eomes_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/eomes_response_to_elevated_platelet_cytosolic_ca2.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('TNFR1-induced NF-kappa-B signaling pathway',
            readable = T,
            foldChange = spleen_entrez_gene_list_eomes_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/eomes_tnfr1_induced_nfkb_signaling_pathway.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

kegg_pathways <-
  as_tibble(kegg_eomes_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES)
reactome_pathways <-
  as_tibble(reactome_eomes_10wpi_spleen) %>% dplyr::select(., Description, NES) %>% arrange(NES)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_eomes.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_eomes.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")
## gata3 ----

rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/gsea_results_tables/spleen_entrez_gene_list_gata3_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

kegg_gata3_10wpi_spleen <- gseKEGG(
  spleen_entrez_gene_list_gata3_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  nPermSimple = 100000,
  eps = 0
)

kegg_pathways <-
  as_tibble(kegg_gata3_10wpi_spleen) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = spleen_entrez_gene_list_gata3_10wpi,
  species = 'hsa',
  pathway.id = '04610',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'gata3_complement_and_coagulation_cascades'
)

reactome_gata3_10wpi_spleen <-
  gsePathway(
    spleen_entrez_gene_list_gata3_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

reactome_pathways <-
  as_tibble(reactome_gata3_10wpi_spleen) %>% dplyr::select(., Description, NES) %>% arrange(NES) %>% headtail(.)


viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = spleen_entrez_gene_list_gata3_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/gata3_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Response to elevated platelet cytosolic Ca2+',
            readable = T,
            foldChange = spleen_entrez_gene_list_gata3_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/gata3_response_to_elevated_platelet_cytosolic_ca2.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('CD209 (DC-SIGN) signaling',
            readable = T,
            foldChange = spleen_entrez_gene_list_gata3_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/gata3_cd209_signaling.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)



write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_gata3.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_gata3.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/kegg/kegg_spleen_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi/pathways/reactome/reactome_spleen_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")



# liver
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE
## dna vaccine ----
load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/gsea_results_tables/liver_entrez_gene_list_dnavaccine_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

kegg_dnavaccine_10wpi_liver <- gseKEGG(
  liver_entrez_gene_list_dnavaccine_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

as_tibble(kegg_dnavaccine_10wpi_liver) %>% arrange(NES) %>% print(n = Inf)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = liver_entrez_gene_list_dnavaccine_10wpi,
  species = 'hsa',
  pathway.id = '04064',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'NF-kB_signaling_pathway_downregulated'
)

pathview(
  gene.data = liver_entrez_gene_list_dnavaccine_10wpi,
  species = 'hsa',
  pathway.id = '04062',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'chemokine_signaling_pathway_upregulated'
)

reactome_dnavaccine_10wpi_liver <-
  gsePathway(
    liver_entrez_gene_list_dnavaccine_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

viewPathway('RHO GTPase Effectors',
            readable = T,
            foldChange = liver_entrez_gene_list_dnavaccine_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/dnavaccine_rho_gtpase_effectors_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('RHOQ GTPase cycle',
            readable = T,
            foldChange = liver_entrez_gene_list_dnavaccine_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/dnavaccine_rhoq_gtpase_cycle_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('RHOB GTPase cycle',
            readable = T,
            foldChange = liver_entrez_gene_list_dnavaccine_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/dnavaccine_rhob_gtpase_cycle_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

kegg_pathways <-
  as_tibble(kegg_dnavaccine_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)
reactome_pathways <-
  as_tibble(reactome_dnavaccine_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_dnavaccine.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_dnavaccine.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")


## iv-ld ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/gsea_results_tables/liver_entrez_gene_list_ivld_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

kegg_ivld_10wpi_liver <- gseKEGG(
  liver_entrez_gene_list_ivld_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  eps = 0,
  nPermSimple = 10000,
  seed = FALSE
)

as_tibble(kegg_ivld_10wpi_liver) %>% arrange(NES) %>% print(n = Inf)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = liver_entrez_gene_list_ivld_10wpi,
  species = 'hsa',
  pathway.id = '04670',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'ivld_leukocyte_transendothelial_migration_upregulated'
)

pathview(
  gene.data = liver_entrez_gene_list_ivld_10wpi,
  species = 'hsa',
  pathway.id = '04664',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'ivld_fc_epsilon_RI_signaling_pathway_upregulated'
)

pathview(
  gene.data = liver_entrez_gene_list_ivld_10wpi,
  species = 'hsa',
  pathway.id = '04662',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'ivld_Bcell_receptor_signaling_pathway_upregulated'
)


reactome_ivld_10wpi_liver <-
  gsePathway(
    liver_entrez_gene_list_ivld_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(reactome_ivld_10wpi_liver) %>% arrange(NES) %>% headtail(.)

viewPathway('ROS and RNS production in phagocytes',
            readable = T,
            foldChange = liver_entrez_gene_list_ivld_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/ivld_ros_and_rns_production_in_phagocytes_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)


kegg_pathways <-
  as_tibble(kegg_ivld_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% arrange(NES)
reactome_pathways <-
  as_tibble(reactome_ivld_10wpi_liver) %>% dplyr::select(., Description, NES) %>% arrange(NES)


write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_ivld.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_ivld.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## iv-hd ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/gsea_results_tables/liver_entrez_gene_list_ivhd_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

kegg_ivhd_10wpi_liver <- gseKEGG(
  liver_entrez_gene_list_ivhd_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE
)

as_tibble(kegg_ivhd_10wpi_liver) %>% arrange(NES) %>% print(n = Inf)  # no immune related pathways

reactome_ivhd_10wpi_liver <-
  gsePathway(
    liver_entrez_gene_list_ivhd_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(reactome_ivhd_10wpi_liver) %>% arrange(NES) %>% print(n = Inf)  # no enriched pathways

viewPathway('Regulation of Complement cascade',
            readable = T,
            foldChange = liver_entrez_gene_list_ivhd_10wpi)  # downregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/ivhd_regulation_complement_cascade.png',
  width = 1000,
  height = 648,
  units = "px",
  dpi = 72,
  bg = 'white'
)

kegg_pathways <-
  as_tibble(kegg_ivhd_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_ivhd.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

## eomes ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/gsea_results_tables/liver_entrez_gene_list_eomes_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

kegg_eomes_10wpi_liver <- gseKEGG(
  liver_entrez_gene_list_eomes_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  nPermSimple = 100000,
  eps = 0
)

as_tibble(kegg_eomes_10wpi_liver) %>% arrange(NES) %>% headtail()
as_tibble(kegg_eomes_10wpi_liver)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

pathview(
  gene.data = liver_entrez_gene_list_eomes_10wpi,
  species = 'hsa',
  pathway.id = '04610',
  kegg.native = T,
  limit = list(gene = 1, cpd = 1),
  out.suffix = 'eomes_complement_and_coagulation_cascades'
)

reactome_eomes_10wpi_liver <-
  gsePathway(
    liver_entrez_gene_list_eomes_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

as_tibble(reactome_eomes_10wpi_liver) %>% arrange(NES) %>% headtail()


kegg_pathways <-
  as_tibble(kegg_eomes_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES)
reactome_pathways <-
  as_tibble(reactome_eomes_10wpi_liver) %>% dplyr::select(., Description, NES) %>% arrange(NES) %>% headtail(.)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_eomes.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_eomes.tsv')

kegg_pathways %>% print(n = Inf)

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")
## gata3 ----
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

load(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/gsea_results_tables/liver_entrez_gene_list_gata3_10wpi.RData'
)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg'
)

kegg_gata3_10wpi_liver <- gseKEGG(
  liver_entrez_gene_list_gata3_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  nPermSimple = 10000,
  eps = 0
)

kegg_pathways <-
  as_tibble(kegg_gata3_10wpi_liver) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% headtail(.)

reactome_gata3_10wpi_liver <-
  gsePathway(
    liver_entrez_gene_list_gata3_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

reactome_pathways <-
  as_tibble(reactome_gata3_10wpi_liver) %>% dplyr::select(., Description, NES) %>% arrange(NES)

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_gata3.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_gata3.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/kegg/kegg_liver_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi/pathways/reactome/reactome_liver_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")



# hkidney 
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE
## dna vaccine ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_dnavaccine_10wpi.RData')

kegg_dnavaccine_10wpi_hkidney <- gseKEGG(
  hkidney_entrez_gene_list_dnavaccine_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)

reactome_dnavaccine_10wpi_hkidney <-
  gsePathway(
    hkidney_entrez_gene_list_dnavaccine_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

## iv-ld ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivld_10wpi.RData')

kegg_ivld_10wpi_hkidney <- gseKEGG(
  hkidney_entrez_gene_list_ivld_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 10000
)

reactome_ivld_10wpi_hkidney <-
  gsePathway(
    hkidney_entrez_gene_list_ivld_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )


viewPathway('FCERI mediated NF-kB activation',
            readable = T,
            foldChange = hkidney_entrez_gene_list_ivld_10wpi)  # downregulated


ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/ivld_fceri_mediated_nfkb_activation_downregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Downstream signaling events of B Cell Receptor (BCR)',
            readable = T,
            foldChange = hkidney_entrez_gene_list_ivld_10wpi)  # downregulated


ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/ivld_downstream_signaling_events_of_bcr_downregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

## iv-hd ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_ivhd_10wpi.RData')

kegg_ivhd_10wpi_hkidney <- gseKEGG(
  hkidney_entrez_gene_list_ivhd_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 10000
)
as_tibble(kegg_ivhd_10wpi_hkidney)

reactome_ivhd_10wpi_hkidney <-
  gsePathway(
    hkidney_entrez_gene_list_ivhd_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )
as_tibble(reactome_ivhd_10wpi_hkidney)

reactome_pathways <-
  as_tibble(reactome_ivhd_10wpi_hkidney) %>% arrange(NES) %>% dplyr::select(., Description, NES) %>% tail(n = 20)


viewPathway('NIK-->noncanonical NF-kB signaling',
            readable = T,
            foldChange = hkidney_entrez_gene_list_ivhd_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/ivhd_nik-noncanonical_nfkb_signaling_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)


viewPathway('Dectin-1 mediated noncanonical NF-kB signaling',
            readable = T,
            foldChange = hkidney_entrez_gene_list_ivhd_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/ivhd_dectin-1_mediated_noncanonical_nfkb_signaling_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

viewPathway('Vpu mediated degradation of CD4',
            readable = T,
            foldChange = hkidney_entrez_gene_list_ivhd_10wpi)  # upregulated

ggsave(
  filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/ivhd_vpu_mediated_degradation_of_cd4_upregulated.png',
  width = 1500,
  height = 1000,
  units = "px",
  dpi = 72,
  bg = 'white'
)

write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/reactome_hkidney_10wpi_ivhd.tsv')

# Convert to a Markdown table ---
reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/pathways/reactome/reactome_hkidney_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )

cat(markdown_table(reactome_data), sep = "\n")


## eomes ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_eomes_10wpi.RData')

kegg_eomes_10wpi_hkidney <- gseKEGG(
  hkidney_entrez_gene_list_eomes_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)
as_tibble(kegg_eomes_10wpi_hkidney)

reactome_eomes_10wpi_hkidney <-
  gsePathway(
    hkidney_entrez_gene_list_eomes_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

as_tibble(reactome_eomes_10wpi_hkidney)

## gata3 ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi/gsea_results_tables/hkidney_entrez_gene_list_gata3_10wpi.RData')

kegg_gata3_10wpi_hkidney <- gseKEGG(
  hkidney_entrez_gene_list_gata3_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)
as_tibble(kegg_gata3_10wpi_hkidney)

reactome_gata3_10wpi_hkidney <-
  gsePathway(
    hkidney_entrez_gene_list_gata3_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

as_tibble(reactome_gata3_10wpi_hkidney)



# heart
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

## dna vaccine ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_dnavaccine_10wpi.RData')

kegg_dnavaccine_10wpi_heart <- gseKEGG(
  heart_entrez_gene_list_dnavaccine_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)


as_tibble(kegg_dnavaccine_10wpi_heart) %>% arrange(NES) %>% print(n = Inf)


reactome_dnavaccine_10wpi_heart <-
  gsePathway(
    heart_entrez_gene_list_dnavaccine_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .2,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )

get_top_bottom_pathways(kegg_dnavaccine_10wpi_heart, n = 10, 'kegg')
get_top_bottom_pathways(reactome_dnavaccine_10wpi_heart, n = 10, 'reactome')

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_dnavaccine.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_dnavaccine.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_dnavaccine.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## iv-ld ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_ivld_10wpi.RData')

kegg_ivld_10wpi_heart <- gseKEGG(
  heart_entrez_gene_list_ivld_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 10000
)

reactome_ivld_10wpi_heart <-
  gsePathway(
    heart_entrez_gene_list_ivld_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

get_top_bottom_pathways(kegg_ivld_10wpi_heart, n = 10, 'kegg')
get_top_bottom_pathways(reactome_ivld_10wpi_heart, n = 10, 'reactome')

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_ivld.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_ivld.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_ivld.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## iv-hd ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_ivhd_10wpi.RData')

kegg_ivhd_10wpi_heart <- gseKEGG(
  heart_entrez_gene_list_ivhd_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 10000
)
as_tibble(kegg_ivhd_10wpi_heart)

reactome_ivhd_10wpi_heart <-
  gsePathway(
    heart_entrez_gene_list_ivhd_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

get_top_bottom_pathways(kegg_ivhd_10wpi_heart, n = 10, 'kegg')
get_top_bottom_pathways(reactome_ivhd_10wpi_heart, n = 10, 'reactome')

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_ivhd.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_ivhd.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_ivhd.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## eomes ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_eomes_10wpi.RData')

kegg_eomes_10wpi_heart <- gseKEGG(
  heart_entrez_gene_list_eomes_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)


reactome_eomes_10wpi_heart <-
  gsePathway(
    heart_entrez_gene_list_eomes_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 10000,
    verbose = F
  )

get_top_bottom_pathways(kegg_eomes_10wpi_heart, n = 10, 'kegg')
get_top_bottom_pathways(reactome_eomes_10wpi_heart, n = 10, 'reactome')

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_eomes.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_eomes.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_eomes.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")

## gata3 ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/gsea_results_tables/heart_entrez_gene_list_gata3_10wpi.RData')

kegg_gata3_10wpi_heart <- gseKEGG(
  heart_entrez_gene_list_gata3_10wpi,
  organism = "hsa",
  keyType = "ncbi-geneid",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  seed = FALSE,
  eps = 0,
  nPermSimple = 100000
)

reactome_gata3_10wpi_heart <-
  gsePathway(
    heart_entrez_gene_list_gata3_10wpi,
    # the gsea_formatting function removes the duplicates from this object
    pvalueCutoff = .05,
    pAdjustMethod = 'BH',
    eps = 1e-300,
    nPermSimple = 100000,
    verbose = F
  )


get_top_bottom_pathways(kegg_gata3_10wpi_heart, n = 10, 'kegg')
get_top_bottom_pathways(reactome_gata3_10wpi_heart, n = 10, 'reactome')

write_tsv(kegg_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_gata3.tsv')
write_tsv(reactome_pathways, file = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_gata3.tsv')

# Convert to a Markdown table ---
kegg_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/kegg/kegg_heart_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(kegg_data), sep = "\n")

reactome_data <-
  read.delim(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi/pathways/reactome/reactome_heart_10wpi_gata3.tsv',
    header = TRUE,
    sep = "\t"
  )
cat(markdown_table(reactome_data), sep = "\n")
