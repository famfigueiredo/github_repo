## Reset Global Environment ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/functions_data-wrangling_march24.R'
)

# Loading packages ####
suppressPackageStartupMessages({
  library('tidyverse')
  library('BiocParallel')
  library('patchwork')
  library('ggrepel')
  library('ggVennDiagram')
  library('ggpmisc')
  library('hrbrthemes')
  library('clusterProfiler')
  library('gprofiler2')
  library('org.Hs.eg.db')
  library('openxlsx')
  register(MulticoreParam(10))
})
# 1wpc ----
### Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_1wpc'
)

results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*1wpc')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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
  c(
    'DNA vaccine vs IV-LD',
    'DNA vaccine vs pTagRFP',
    'EOMES vs IV-LD',
    'EOMES vs pTagRFP',
    'GATA3 vs IV-LD',
    'GATA3 vs pTagRFP',
    'IV-HD vs IV-LD',
    'IV-HD vs pTagRFP',
    'IV-LD vs pTagRFP'
  )

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'DNA vaccine vs IV-LD',
      'DNA vaccine vs pTagRFP',
      'EOMES vs IV-LD',
      'EOMES vs pTagRFP',
      'GATA3 vs IV-LD',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'IV-HD vs pTagRFP',
      'IV-LD vs pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Contrast'))

deg_regulation_summary <-
  deg_regulation_summary %>% add_row(Contrast = 'DNA vaccine vs IV-LD',
                                     Regulation = 'downregulated',
                                     n = 0) %>%
  add_row(Contrast = 'IV-HD vs IV-LD',
          Regulation = 'upregulated',
          n = 0)

### Reordering factors to plot ----
deg_regulation_summary$Contrast <-
  factor(
    deg_regulation_summary$Contrast,
    levels = c(
      'IV-LD vs pTagRFP',
      'IV-HD vs pTagRFP',
      'DNA vaccine vs pTagRFP',
      'EOMES vs pTagRFP',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'DNA vaccine vs IV-LD',
      'EOMES vs IV-LD',
      'GATA3 vs IV-LD'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Contrast, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Pastel1') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 1 WPC') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_vline(xintercept = 5.5,
             linetype = 'solid',
             linewidth = 0.3) +
  geom_hline(yintercept = 0, size = 0.2) +
  annotate(
    'text',
    x = 3,
    y = 400,
    label = 'vs pTagRFP',
    size = 4,
    family = 'serif'
  ) +
  annotate(
    'text',
    x = 7.5,
    y = 400,
    label = 'vs IV-LD',
    size = 4,
    family = 'serif'
  )

ggsave(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within sampling point_across treatments/1wpc/DGE_1wpc_treatments.png',
  plot = last_plot()
)

#### With conu as a reference ----
##### Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_1wpc'
)

results_files <- list.files(pattern = '_conu')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)
##### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*1wpc')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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
  c('DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD')

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(c(
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD'
  )), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

### Reordering factors to plot ----
deg_regulation_summary$Treatment <-
  factor(
    deg_regulation_summary$Treatment,
    levels = c('IV-LD',
               'IV-HD',
               'DNA vaccine',
               'EOMES',
               'GATA3')
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
  scale_fill_manual(values = c('#FC7380',
                               '#FCEF73')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 1 WPC, CONU as reference') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2)

ggsave(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within sampling point_across treatments/1wpc/DGE_1wpc_treatments_vs_conu.png',
  plot = last_plot()
)









# 4wpc ----
### Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc'
)

results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*4wpc')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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
  c(
    'DNA vaccine vs IV-LD',
    'DNA vaccine vs pTagRFP',
    'EOMES vs IV-LD',
    'EOMES vs pTagRFP',
    'GATA3 vs IV-LD',
    'GATA3 vs pTagRFP',
    'IV-HD vs IV-LD',
    'IV-HD vs pTagRFP',
    'IV-LD vs pTagRFP'
  )

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'DNA vaccine vs IV-LD',
      'DNA vaccine vs pTagRFP',
      'EOMES vs IV-LD',
      'EOMES vs pTagRFP',
      'GATA3 vs IV-LD',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'IV-HD vs pTagRFP',
      'IV-LD vs pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Contrast'))

deg_regulation_summary <-
  deg_regulation_summary %>% add_row(Contrast = 'DNA vaccine vs IV-LD',
                                     Regulation = 'downregulated',
                                     n = 0)

### Reordering factors to plot ----
deg_regulation_summary$Contrast <-
  factor(
    deg_regulation_summary$Contrast,
    levels = c(
      'IV-LD vs pTagRFP',
      'IV-HD vs pTagRFP',
      'DNA vaccine vs pTagRFP',
      'EOMES vs pTagRFP',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'DNA vaccine vs IV-LD',
      'EOMES vs IV-LD',
      'GATA3 vs IV-LD'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Contrast, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_brewer(palette = 'Dark2') +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 4 WPC') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_vline(xintercept = 5.5,
             linetype = 'solid',
             linewidth = 0.3) +
  geom_hline(yintercept = 0, size = 0.2) +
  annotate(
    'text',
    x = 3,
    y = 550,
    label = 'vs pTagRFP',
    size = 4,
    family = 'serif'
  ) +
  annotate(
    'text',
    x = 7.5,
    y = 550,
    label = 'vs IV-LD',
    size = 4,
    family = 'serif'
  )

ggsave(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within sampling point_across treatments/4wpc/DGE_4wpc_treatments.png',
  plot = last_plot()
)



# 10wpc ----
### Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc'
)

results_files <- list.files(pattern = 'res_')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)

### For loop, significant genes ----
results_files <-
  ls(pattern = 'res_.*10wpc')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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
  c(
    'DNA vaccine vs IV-LD',
    'DNA vaccine vs pTagRFP',
    'EOMES vs IV-LD',
    'EOMES vs pTagRFP',
    'GATA3 vs IV-LD',
    'GATA3 vs pTagRFP',
    'IV-HD vs IV-LD',
    'IV-HD vs pTagRFP',
    'IV-LD vs pTagRFP'
  )

### Delisting dataframes ----
list2env(all_significant_genes, envir = .GlobalEnv)  # splitting the dataframes to summarize their contents

### Summarizing differentially expressed genes ----
deg_regulation_summary <-
  lapply(mget(
    c(
      'DNA vaccine vs IV-LD',
      'DNA vaccine vs pTagRFP',
      'EOMES vs IV-LD',
      'EOMES vs pTagRFP',
      'GATA3 vs IV-LD',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'IV-HD vs pTagRFP',
      'IV-LD vs pTagRFP'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Contrast'))

### Reordering factors to plot ----
deg_regulation_summary$Contrast <-
  factor(
    deg_regulation_summary$Contrast,
    levels = c(
      'IV-LD vs pTagRFP',
      'IV-HD vs pTagRFP',
      'DNA vaccine vs pTagRFP',
      'EOMES vs pTagRFP',
      'GATA3 vs pTagRFP',
      'IV-HD vs IV-LD',
      'DNA vaccine vs IV-LD',
      'EOMES vs IV-LD',
      'GATA3 vs IV-LD'
    )
  )

### Plotting ----
deg_regulation_summary %>%
  ggplot(aes(x = Contrast, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    size = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 10 WPC') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_vline(xintercept = 5.5,
             linetype = 'solid',
             linewidth = 0.3) +
  geom_hline(yintercept = 0, size = 0.2) +
  annotate(
    'text',
    x = 3,
    y = 200,
    label = 'vs pTagRFP',
    size = 4,
    family = 'serif'
  ) +
  annotate(
    'text',
    x = 7.5,
    y = 200,
    label = 'vs IV-LD',
    size = 4,
    family = 'serif'
  )

ggsave(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/plots/DGE within sampling point_across treatments/10wpc/DGE_10wpc_treatments.png',
  plot = last_plot()
)


# 4wpc - Testing Venn diagrams ----
## Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc'
)

results_files <-
  list.files(pattern = 'res_.*_ptag.')  # regex matching results files using ptag as reference
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)

### Venn diagram ----
dnavaccine_4wpc <- significant_genes(res_dnavaccine_vs_ptag_4wpc)
eomes_4wpc <- significant_genes(res_eomes_vs_ptag_4wpc)
gata3_4wpc <- significant_genes(res_gata3_vs_ptag_4wpc)
ivhd_4wpc <- significant_genes(res_ivhd_vs_ptag_4wpc)

#### Upregulated ----
z <- list(
  A = dnavaccine_4wpc[dnavaccine_4wpc$log2FC > 0,]$ID,
  B = eomes_4wpc[eomes_4wpc$log2FC > 0,]$ID,
  C = gata3_4wpc[gata3_4wpc$log2FC > 0,]$ID,
  D = ivhd_4wpc[ivhd_4wpc$log2FC > 0,]$ID
)

names(z) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD')

ggVennDiagram(
  z,
  label_alpha = 0,
  color = 'black',
  lwd = 0.5,
  edge_size = 0.5,
  set_size = 3,
  label_geom = 'label',
  show_intersect = F
) +
  theme(legend.position = 'none') +
  ggtitle('Upregulated genes, heart, 4WPC') +
  theme(plot.title = element_text(hjust = .5)) +
  # ggtitle('Upregulated genes, heart, 4WPC') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_distiller(palette = 'Spectral', direction = -1) +
  theme(text = element_text(family = 'serif', size = 12))

#### Common genes
common_genes_upregulated_4WPC <-
  Reduce(
    intersect,
    list(
      dnavaccine_4wpc[dnavaccine_4wpc$log2FC > 0,]$ID,
      eomes_4wpc[eomes_4wpc$log2FC > 0,]$ID,
      gata3_4wpc[gata3_4wpc$log2FC > 0,]$ID,
      ivhd_4wpc[ivhd_4wpc$log2FC > 0,]$ID
    )
  )

dnavaccine_4wpc_up <- dnavaccine_4wpc[dnavaccine_4wpc$log2FC > 0,]$ID
eomes_4wpc_up <- eomes_4wpc[eomes_4wpc$log2FC > 0,]$ID
gata3_4wpc_up <- gata3_4wpc[gata3_4wpc$log2FC > 0,]$ID
ivhd_4wpc_up <- ivhd_4wpc[ivhd_4wpc$log2FC > 0,]$ID


common_upregulated_orthologs_4wpc <- gorth(
  query = common_genes_upregulated_4WPC,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

#### GATA3 exclusive genes - 182 exclusive genes
gata3_exclusive_genes <-
  setdiff(gata3_4wpc_up,
          c(dnavaccine_4wpc_up, eomes_4wpc_up, ivhd_4wpc_up))

gata3_up_orthologs_4wpc <- gorth(
  query = gata3_exclusive_genes,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

head(gata3_up_orthologs_4wpc)
nrow(gata3_up_orthologs_4wpc)


ora_gata3_up <-
  enrichGO(
    gene = gata3_up_orthologs_4wpc$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )


as_tibble(gofilter(ora_gata3_up, level = 6)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value', ) +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  theme(legend.position = 'right') +
  ggtitle('Upregulated GATA3 exclusive genes, heart, 4WPC') +
  theme(text = element_text(family = 'serif', size = 12))


ora_common_up <-
  enrichGO(
    gene = common_upregulated_orthologs_4wpc$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

as_tibble(gofilter(ora_common_up, level = 5)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(
    stat = "identity",
    colour = 'black',
    linewidth = 0.3,
    show.legend = T
  ) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value') +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  theme(legend.position = 'right') +
  ggtitle('Upregulated common genes, heart, 4WPC') +
  theme(text = element_text(family = 'serif', size = 12))

#### Downregulated ----
y <- list(
  A = dnavaccine_4wpc[dnavaccine_4wpc$log2FC < 0,]$ID,
  B = eomes_4wpc[eomes_4wpc$log2FC < 0,]$ID,
  C = gata3_4wpc[gata3_4wpc$log2FC < 0,]$ID,
  D = ivhd_4wpc[ivhd_4wpc$log2FC < 0,]$ID
)

names(y) <-
  c('DNA vaccine', 'EOMES', 'GATA3', 'IVHD')

ggVennDiagram(
  y,
  label_alpha = 0,
  color = 'black',
  lwd = 0.5,
  edge_size = 0.5,
  set_size = 3,
  label_geom = 'label',
  show_intersect = F
) +
  theme(legend.position = 'none') +
  ggtitle('Downregulated common genes, heart, 4WPC') +
  theme(plot.title = element_text(hjust = .5)) +
  # ggtitle('Upregulated genes, heart, 4WPC') +
  scale_x_continuous(expand = expansion(mult = .2)) +
  scale_fill_distiller(palette = 'Spectral', direction = -1) +
  theme(text = element_text(family = 'serif', size = 12))


common_genes_downregulated_4WPC <-
  Reduce(
    intersect,
    list(
      dnavaccine_4wpc[dnavaccine_4wpc$log2FC < 0,]$ID,
      eomes_4wpc[eomes_4wpc$log2FC < 0,]$ID,
      gata3_4wpc[gata3_4wpc$log2FC < 0,]$ID,
      ivhd_4wpc[ivhd_4wpc$log2FC < 0,]$ID
    )
  )

dnavaccine_4wpc_down <- dnavaccine_4wpc[dnavaccine_4wpc$log2FC < 0,]$ID
eomes_4wpc_down <- eomes_4wpc[eomes_4wpc$log2FC < 0,]$ID
gata3_4wpc_down <- gata3_4wpc[gata3_4wpc$log2FC < 0,]$ID
ivhd_4wpc_down <- ivhd_4wpc[ivhd_4wpc$log2FC < 0,]$ID


common_downregulated_orthologs_4wpc <- gorth(
  query = common_genes_downregulated_4WPC,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

ora_common_down <-
  enrichGO(
    gene = common_downregulated_orthologs_4wpc$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )

as_tibble(ora_common_down)  # only one term downregulated. GO:0140467 - integrated stress response signaling, with 3 genes


#### GATA3 exclusive genes - 272 exclusive downregulated genes
gata3_exclusive_genes_down <-
  setdiff(gata3_4wpc_down,
          c(dnavaccine_4wpc_down, eomes_4wpc_down, ivhd_4wpc_down))


gata3_down_orthologs_4wpc <- gorth(
  query = gata3_exclusive_genes_down,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

head(gata3_down_orthologs_4wpc)

ora_gata3_down <-
  enrichGO(
    gene = gata3_down_orthologs_4wpc$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )


as_tibble(gofilter(ora_gata3_down, level = 6)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value') +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  theme(legend.position = 'right') +
  ggtitle('Upregulated GATA3 exclusive genes, heart, 4WPC') +
  theme(text = element_text(family = 'serif', size = 12))


#### EOMES exclusive genes - 196 exclusive downregulated genes
eomes_exclusive_genes_down <-
  setdiff(eomes_4wpc_down,
          c(dnavaccine_4wpc_down, gata3_4wpc_down, ivhd_4wpc_down))


eomes_down_orthologs_4wpc <- gorth(
  query = eomes_exclusive_genes_down,
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

head(eomes_down_orthologs_4wpc)

ora_eomes_down <-
  enrichGO(
    gene = eomes_down_orthologs_4wpc$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )


as_tibble(gofilter(ora_eomes_down, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value') +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  theme(legend.position = 'right') +
  ggtitle('Downregulated EOMES exclusive genes, heart, 4WPC') +
  theme(text = element_text(family = 'serif', size = 12))


intersection <- intersect(gata3_4wpc_down, eomes_4wpc_down)

gata3EOMES_exclusive_genes_down <- setdiff(intersection,
        c(dnavaccine_4wpc_down, ivhd_4wpc_down))


ora_gata3EOMES_down <-
  enrichGO(
    gene = gata3EOMES_exclusive_genes_down$ortholog_ensg,
    # significantly upregulated genes
    OrgDb = 'org.Hs.eg.db',
    keyType = 'ENSEMBL',
    ont = 'BP',
    # biological processes
    pAdjustMethod = 'BH',
    pvalueCutoff = 0.05,
    readable = T
  )


as_tibble(gofilter(ora_eomes_down, level = 4)) %>%  # filtering for level x GO Terms
  mutate(Description = fct_reorder(Description, p.adjust)) %>%  # sorting 'Description' by Adjusted p-value
  ggplot(aes(x = Description, y = Count, fill = p.adjust)) +
  geom_bar(stat = "identity",
           colour = 'black',
           linewidth = 0.3) +
  coord_flip() +
  scale_fill_distiller(palette = 'Spectral',
                       name = 'Adjusted \n p-value') +
  scale_y_continuous() +
  guides(fill = guide_colourbar(barwidth = 2, barheight = 10)) +
  labs(x = '', y = 'Gene count', fill = 'Adjusted p-value') +
  theme_classic() +
  theme(legend.position = 'right') +
  ggtitle('Downregulated EOMES and GATA3 exclusive genes, heart, 4WPC') +
  theme(text = element_text(family = 'serif', size = 12))






