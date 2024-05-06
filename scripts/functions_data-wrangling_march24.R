# Loading functions ####
# Significant genes grabs a DESeq2 result table, subsets the genes with padj < 0.05, and selects the ID, log2FoldChange, and padj columns. Also arranges log2FC in a descending manner.
significant_genes <- function(results_files) {
  b <- as.data.frame(subset(results_files, padj < 0.05)) %>%
    rownames_to_column(var = 'ID') %>%
    as_tibble()

  sig_genes <- b %>%
    dplyr::select(ID, log2FC = log2FoldChange, adjusted_p.val = padj, pvalue) %>%
    dplyr::arrange(desc(log2FC))

  return(sig_genes)
}

## Testing GSEA without selecting by p-value
# significant_genes <- function(results_files) {
#   b <- as.data.frame(results_files) %>%
#     rownames_to_column(var = 'ID') %>%
#     as_tibble()
#   
#   sig_genes <- b %>%
#     dplyr::select(ID, log2FC = log2FoldChange, adjusted_p.val = padj, pvalue) %>%
#     dplyr::arrange(desc(log2FC))
#   
#   return(sig_genes)
# }

# Significant genes metrics creates a table with information about number of up/down significantly regulated genes
sig_genes_metrics <- function(significant_genes) {
  as.data.frame(significant_genes) %>% filter(adjusted_p.val < 0.05) %>%
    dplyr::mutate(Regulation = ifelse(
      log2FC < 0,
      'downregulated',
      ifelse(log2FC > 0, 'upregulated', NA)
    )) %>%
    filter(!is.na(Regulation)) %>%
    dplyr::count(Regulation)
}

# improved data wrangling function selects significantly regulated genes, converts ssalar IDs to hsapiens orthologs, and joins both
improved_data_wrangling <-
  function(results_table, treatment, sampling_point) {
    if (!requireNamespace("dplyr", quietly = TRUE)) {
      install.packages("dplyr")
    }
    
    if (!requireNamespace("gprofiler2", quietly = TRUE)) {
      install.packages("gprofiler2")
    }
    
    # Load the required packages
    library(dplyr)
    library(gprofiler2)
    
    # select significant genes from DESeq2 results table
    a <-
      significant_genes(results_table)
    
    # convert ssalar gene IDs to human orthologs
    orth_hs <- gorth(
      query = a$ID,
      source_organism = 'ssalar',
      target_organism = 'hsapiens',
      mthreshold = 1,
      filter_na = T
    )
    
    # join significant genes table with human ortholog names
    results <- a %>% left_join(orth_hs, by = c('ID' = 'input')) %>%
      dplyr::select(.,
                    ID,
                    ortholog_name,
                    log2FC,
                    adjusted_p.val,
                    pvalue,
                    ortholog_ensg,
                    description)
    
    # getting gene lists for ORA
    ora_up <<-
      results %>% drop_na() %>% dplyr::filter(log2FC > 0) %>%  pull(ortholog_ensg)
    
    ora_down <<-
      results %>% drop_na() %>% dplyr::filter(log2FC < 0) %>%  pull(ortholog_ensg)
    
    # create results name
    results_name <-
      paste('results', treatment, sampling_point, sep = '_')
    
    # Add treatment column
    results <- results %>%
      mutate(treatment = treatment)
    
    # assign the results data frame to the dynamic name
    assign(results_name, results, envir = .GlobalEnv)
    
  }

# combine and label pathways merges both downregulated and upregulated dataframes into one
combine_and_label_pathways <-
  function(downreg_tibble, upreg_tibble) {
    # Add a "regulation" column to downregulated tibble
    df_downregulated <- as_tibble(downreg_tibble) %>%
      mutate(regulation = "downregulated")
    
    # Add a "regulation" column to upregulated tibble
    df_upregulated <- as_tibble(upreg_tibble) %>%
      mutate(regulation = "upregulated")
    
    # Concatenate the two tibbles
    final_tibble <- bind_rows(df_downregulated, df_upregulated)
    
    return(final_tibble)
  }

# filter rows by GO term intersects the above merged dataframe with the immune related GO terms and returns a dataframe containing those
filter_rows_by_GO_term <- function(df1, df2, id_column_name) {
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    install.packages("dplyr")
  }
  # Perform the intersection
  common_GO_terms <- intersect(df1$ID, df2[[id_column_name]])
  
  # Filter rows based on common GO terms
  filtered_rows <- df1 %>%
    dplyr::filter(ID %in% common_GO_terms)
  
  return(filtered_rows)
}


# Helper function to display Venn diagram
display_venn <- function(x, ...){
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
