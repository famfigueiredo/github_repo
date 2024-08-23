## Loading functions ----
source(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/github_repo/scripts/functions_data-wrangling_march24.R'
)

# Heart ----
# 10 WPI ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_10wpi'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_10wpi')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Renaming GE objects to match filename
objects <- ls(pattern = '^res_.*_conu_10wpi')
for (obj in objects) {
  assign(paste0('heart_', obj), get(obj))
}
rm(list = objects)


# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_10wpi_heart <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment')) %>% mutate(tissue = 'heart')

# Reordering factors to plot
deg_regulation_summary_10wpi$Treatment <-
  factor(
    deg_regulation_summary_10wpi$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_10wpi %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_manual(values = c('#ffffba', '#baffc9')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Heart tissue at 10 WPI') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 


# 1WPC ----
## Loading results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_1wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Renaming GE objects to match filename
objects <- ls(pattern = '^res_.*_conu_1wpc')
for (obj in objects) {
  assign(paste0('heart_', obj), get(obj))
}
rm(list = objects)

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_1wpc <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

# Reordering factors to plot
deg_regulation_summary_1wpc$Treatment <-
  factor(
    deg_regulation_summary_1wpc$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_1wpc %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_manual(values = c('#bae1ff', '#ffb3ba')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 1 WPC') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 


# 4WPC ----
## Loading results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc'
)

results_files <-
  list.files(pattern = '^heart_.*_conu_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Renaming GE objects to match filename
# objects <- ls(pattern = '^res_.*_conu_4wpc')
# for (obj in objects) {
#   assign(paste0('heart_', obj), get(obj))
# }
# rm(list = objects)


# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_4wpc <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

# Reordering factors to plot
deg_regulation_summary_4wpc$Treatment <-
  factor(
    deg_regulation_summary_4wpc$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

theme_update(plot.title = element_text(hjust = 0.5))
### Plotting
deg_regulation_summary_4wpc %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
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
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Cardiac tissue at 4 WPC') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, linewidth = 0.2) 


# Spleen ----
# 10 WPI ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi'
)

results_files <-
  list.files(pattern = '^spleen_.*_conu_10wpi')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_10wpi_spleen <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment')) %>% mutate(tissue = 'spleen')

# Reordering factors to plot
deg_regulation_summary_10wpi$Treatment <-
  factor(
    deg_regulation_summary_10wpi$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_10wpi %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_manual(values = c('#ffffba', '#baffc9')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Spleen tissue at 10 WPI') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 

# 4 WPC ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_4wpc'
)

results_files <-
  list.files(pattern = '^spleen_.*_conu_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_4wpc <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

# Reordering factors to plot
deg_regulation_summary_4wpc$Treatment <-
  factor(
    deg_regulation_summary_4wpc$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_4wpc %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
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
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Spleen tissue at 4 WPC') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 


# Liver ----
# 10 WPI ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi'
)

results_files <-
  list.files(pattern = '^liver_.*_conu_10wpi')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_10wpi_liver <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment')) %>% mutate(tissue = 'liver')

# Reordering factors to plot
deg_regulation_summary_10wpi$Treatment <-
  factor(
    deg_regulation_summary_10wpi$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_10wpi %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_manual(values = c('#ffffba', '#baffc9')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Liver tissue at 10 WPI') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 

# 4 WPC ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc'
)

results_files <-
  list.files(pattern = '^liver_.*_conu_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_4wpc <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

# Reordering factors to plot
deg_regulation_summary_4wpc$Treatment <-
  factor(
    deg_regulation_summary_4wpc$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_4wpc %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
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
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Liver tissue at 4 WPC') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 


# Head-kidney ----
# 10 WPI ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_10wpi')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[-c(2, 3, 4)]  # removing treatments without significantly differentially regulated genes

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_10wpi_hkidney <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment')) %>% mutate(tissue = 'hkidney')

# Reordering factors to plot
deg_regulation_summary_10wpi$Treatment <-
  factor(
    deg_regulation_summary_10wpi$Treatment,
    levels = c(
      'IV-LD',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_10wpi %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
  geom_bar(
    stat = 'identity',
    position = 'dodge',
    colour = 'black',
    linewidth = .2
  ) +
  scale_fill_manual(values = c('#ffffba', '#baffc9')) +
  ylab('Number of genes') +
  geom_text(
    aes(label = n),
    position = position_dodge(width = 0.9),
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Head-kidney tissue at 10 WPI') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 

# 4 WPC ----
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc'
)

results_files <-
  list.files(pattern = '^hkidney_.*_conu_4wpc')  # regex matching results files' filename
list.data <- list()
for (i in 1:length(results_files)) {
  list.data[[i]] <- load(results_files[i])
}

# Create empty list to store results
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

# Reordering dataframes
all_significant_genes <- all_significant_genes[c(5, 4, 3, 2, 1)]

# Renaming dataframes
names(all_significant_genes) <-
  c(
    'IV-LD',
    'IV-HD',
    'GATA3',
    'EOMES',
    'DNA vaccine'
  )

# Delisting dataframes
list2env(all_significant_genes, envir = .GlobalEnv)

# Summarizing differentially expressed genes
deg_regulation_summary <-
  lapply(mget(
    c('IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine')
  ), sig_genes_metrics)

deg_regulation_summary_4wpc <-
  as_tibble(bind_rows(deg_regulation_summary, .id = 'Treatment'))

# Reordering factors to plot
deg_regulation_summary_4wpc$Treatment <-
  factor(
    deg_regulation_summary_4wpc$Treatment,
    levels = c(
      'IV-LD',
      'IV-HD',
      'GATA3',
      'EOMES',
      'DNA vaccine'
    )
  )

### Plotting
deg_regulation_summary_4wpc %>%
  ggplot(aes(x = Treatment, y = n, fill = Regulation)) +
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
    vjust = -0.5
  ) +
  ggtitle(label = 'Differentially Expressed Genes',
          subtitle = 'Head-kidney tissue at 4 WPC') +
  theme_light() +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5)) +
  theme(panel.grid = element_blank()) +
  theme(text = element_text(size = 14, family = 'serif')) +
  theme(axis.text.x = element_text(angle = 320, vjust = 0.5)) +
  geom_vline(
    xintercept = c(1.5, 2.5, 3.5, 4.5, 6.5, 7.5, 8.5),
    linetype = 'dotted',
    linewidth = 0.2
  ) +
  geom_hline(yintercept = 0, size = 0.2) 



# Volcano plots ----
library(EnhancedVolcano)
rm(list = ls()[sapply(ls(), function(x) !is.function(get(x)))])  # delete values, keep functions in GE

# Heart ----

## DNA vaccine 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_dnavaccine_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_dnavaccine_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_dnavaccine_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_dnavaccine_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_dnavaccine_1wpc)

EnhancedVolcano(results_dnavaccine_1wpc,
                lab = results_dnavaccine_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'DNA vaccine',
                subtitle = bquote(italic('1WPC'))) + 
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_dnavaccine.png', 
       width = 1000, 
       height = 632, 
       units = "px", 
       dpi = 72)

## IV-LD 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_ivld_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivld_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivld_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivld_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivld_1wpc)

EnhancedVolcano(results_ivld_1wpc,
                lab = results_ivld_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-LD',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_ivld.png', width = 1000, height = 632, units = "px", dpi = 72)
## IV-HD 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_ivhd_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivhd_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivhd_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivhd_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivhd_1wpc)

EnhancedVolcano(results_ivhd_1wpc,
                lab = results_ivhd_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-HD',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_ivhd.png', width = 1000, height = 632, units = "px", dpi = 72)
## EOMES 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_eomes_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_eomes_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_eomes_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_eomes_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_eomes_1wpc)

EnhancedVolcano(results_eomes_1wpc,
                lab = results_eomes_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'EOMES',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_eomes.png', width = 1000, height = 632, units = "px", dpi = 72)

## GATA3 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_gata3_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_gata3_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_gata3_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_gata3_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_gata3_1wpc)

EnhancedVolcano(results_gata3_1wpc,
                lab = results_gata3_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'GATA3',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_gata3.png', width = 1000, height = 632, units = "px", dpi = 72)

### DNA vaccine 4WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
orth_hs <- gorth(
  query = rownames(res_dnavaccine_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(res_dnavaccine_vs_conu_4wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_dnavaccine_4wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_dnavaccine_4wpc)

EnhancedVolcano(results_dnavaccine_4wpc,
                lab = results_dnavaccine_4wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'DNA vaccine',
                subtitle = bquote(italic('4WPC'))) + 
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))


### IV-LD 4WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_ivld_vs_conu_4wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivld_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivld_vs_conu_4wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivld_4wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivld_4wpc)

EnhancedVolcano(results_ivld_4wpc,
                lab = results_ivld_4wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-LD',
                subtitle = bquote(italic('4WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/volcanoplot_shrunk_ivld.png', width = 1000, height = 632, units = "px", dpi = 72)
### IV-HD 4WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_ivhd_vs_conu_4wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivhd_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivhd_vs_conu_4wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivhd_4wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivhd_4wpc)

EnhancedVolcano(results_ivhd_4wpc,
                lab = results_ivhd_4wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-HD',
                subtitle = bquote(italic('4WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/volcanoplot_shrunk_ivhd.png', width = 1000, height = 632, units = "px", dpi = 72)
### EOMES 4WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_eomes_vs_conu_4wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_eomes_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_eomes_vs_conu_4wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_eomes_4wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_eomes_4wpc)

EnhancedVolcano(results_eomes_4wpc,
                lab = results_eomes_4wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'EOMES',
                subtitle = bquote(italic('4WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/volcanoplot_shrunk_eomes.png', width = 1000, height = 632, units = "px", dpi = 72)

### GATA3 4WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/heart_res_gata3_vs_conu_4wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_gata3_vs_conu_4wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_gata3_vs_conu_4wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_gata3_4wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_gata3_4wpc)

EnhancedVolcano(results_gata3_4wpc,
                lab = results_gata3_4wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'GATA3',
                subtitle = bquote(italic('4WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc/volcanoplot_shrunk_gata3.png', width = 1000, height = 632, units = "px", dpi = 72)





## DNA vaccine 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_dnavaccine_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_dnavaccine_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_dnavaccine_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_dnavaccine_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_dnavaccine_1wpc)

EnhancedVolcano(results_dnavaccine_1wpc,
                lab = results_dnavaccine_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'DNA vaccine',
                subtitle = bquote(italic('1WPC'))) + 
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_dnavaccine.png', 
       width = 1000, 
       height = 632, 
       units = "px", 
       dpi = 72)

## IV-LD 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_ivld_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivld_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivld_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivld_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivld_1wpc)

EnhancedVolcano(results_ivld_1wpc,
                lab = results_ivld_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-LD',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_ivld.png', width = 1000, height = 632, units = "px", dpi = 72)
## IV-HD 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_ivhd_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_ivhd_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_ivhd_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_ivhd_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_ivhd_1wpc)

EnhancedVolcano(results_ivhd_1wpc,
                lab = results_ivhd_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'IV-HD',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_ivhd.png', width = 1000, height = 632, units = "px", dpi = 72)
## EOMES 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_eomes_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_eomes_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_eomes_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_eomes_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_eomes_1wpc)

EnhancedVolcano(results_eomes_1wpc,
                lab = results_eomes_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'EOMES',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_eomes.png', width = 1000, height = 632, units = "px", dpi = 72)

## GATA3 1WPC volcano plot ----
# converting DESeq2 result table rownames to human orthologs
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/heart_res_gata3_vs_conu_1wpc.RData')

orth_hs <- gorth(
  query = rownames(heart_res_gata3_vs_conu_1wpc),
  source_organism = 'ssalar',
  target_organism = 'hsapiens',
  mthreshold = 1,
  filter_na = T
)

# rownames to column to be able to use left_join
df <- as.data.frame(heart_res_gata3_vs_conu_1wpc) %>% rownames_to_column(var = 'ensembl')

# join results table with human ortholog names
results_gata3_1wpc <- df %>% left_join(orth_hs, by = c('ensembl' = 'input')) %>% na.omit() %>% 
  dplyr::select(.,
                ensembl,
                ortholog_name,
                log2FoldChange,
                padj,
                pvalue,
                ortholog_ensg)

head(results_gata3_1wpc)

EnhancedVolcano(results_gata3_1wpc,
                lab = results_gata3_1wpc$ortholog_name,
                x = 'log2FoldChange',
                y = 'pvalue',
                FCcutoff = 1,
                title = 'GATA3',
                subtitle = bquote(italic('1WPC'))) +
  theme(text = element_text(size = 8, family = 'Times New Roman')) +
  theme(plot.title = element_text(hjust = .5),
        plot.subtitle = element_text(hjust = .5))

ggsave(filename = '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc/volcanoplot_shrunk_gata3.png', width = 1000, height = 632, units = "px", dpi = 72)
