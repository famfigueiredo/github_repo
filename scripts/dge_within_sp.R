## Reset Global Environment ----
rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
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
  MulticoreParam(10)
})

## NEED TO TRY STRING FOR PATHWAY ENRICHMENT AND GSEA ##
# 10 wpi ----
### Loading results files ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi'
)

results_files <- list.files(pattern = '^res_.*10wpi')
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

rm(list.data, i, results_files)

### For loop, significant genes ----
results_files <-
  ls(pattern = '^res_.*(ivld|ptagrfp)_10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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

# deg_regulation_summary <-
#   deg_regulation_summary %>% add_row(Contrast = 'DNA vaccine vs IV-LD',
#                                      Regulation = 'downregulated',
#                                      n = 0) %>%
#   add_row(Contrast = 'IV-HD vs IV-LD',
#           Regulation = 'upregulated',
#           n = 0)

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
  scale_fill_manual(values = c('#C9C9FF', '#E1F7D5')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    linewidth = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 10 WPI') +
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
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/plots/DGE within sampling point_across treatments/DGE_10wpi_treatments.png',
  plot = last_plot()
)


### CONU as reference ----
results_files <-
  ls(pattern = '^res_.*(conu)_10wpi')  # listing Global Environment files matching regex pattern, to be used in the for loop.

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

names(all_significant_genes) <-
  c(
    'DNA vaccine',
    'EOMES',
    'GATA3',
    'IV-HD',
    'IV-LD'
  )

list2env(all_significant_genes, envir = .GlobalEnv)


deg_regulation_summary <-
  lapply(mget(
    c(
      'DNA vaccine',
      'EOMES',
      'GATA3',
      'IV-HD',
      'IV-LD'
    )
  ), sig_genes_metrics)

### Transforming summary into a tibble ----
deg_regulation_summary <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Contrast'))

deg_regulation_summary <-
  deg_regulation_summary %>% add_row(Contrast = 'GATA3',
                                     Regulation = 'downregulated',
                                     n = 0) %>% 
  add_row(Contrast = 'IV-HD',
          Regulation = 'downregulated',
          n = 0)


### Reordering factors to plot ----
deg_regulation_summary$Contrast <-
  factor(
    deg_regulation_summary$Contrast,
    levels = c(
      'IV-LD',
      'IV-HD',
      'DNA vaccine',
      'EOMES',
      'GATA3'
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
  scale_fill_manual(values = c('#C9C9FF', '#E1F7D5')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5,
    linewidth = 4
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 10 WPI, CONU as reference') +
  theme_light() +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 270, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  )

ggsave(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/plots/DGE within sampling point_across treatments/DGE_10wpi_CONUref.png',
  plot = last_plot()
)


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


