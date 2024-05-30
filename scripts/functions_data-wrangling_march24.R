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

# running gsea starting from a DESeq results table
gsea_formatting <-
  function(results_table, treatment, sampling_point) {
    # Install and load required packages
    required_packages <-
      c('dplyr', 'gprofiler2', 'clusterProfiler', 'org.Hs.eg.db')
    installed_packages <- rownames(installed.packages())
    
    for (pkg in required_packages) {
      if (!(pkg %in% installed_packages)) {
        install.packages(pkg)
      }
      library(pkg, character.only = TRUE)
    }
    
    # Convert rownames to column 'ensembl'
    results_df <-
      tibble::rownames_to_column(as.data.frame(results_table), var = 'ensembl')
    
    # Convert salmon genes to human orthologs
    orthologs <- gorth(
      query = rownames(results_table),
      source_organism = 'ssalar',
      target_organism = 'hsapiens',
      mthreshold = 1,
      filter_na = TRUE
    )
    
    # Select relevant variables and join with ortholog data
    merged_df <- results_df %>%
      left_join(orthologs, by = c('ensembl' = 'input')) %>%
      dplyr::select(ensembl,
                    ortholog_name,
                    ortholog_ensg,
                    log2FoldChange,
                    padj,
                    description) %>%
      na.omit()
    
    # Order genes by fold change
    ordered_df <- merged_df[order(-merged_df$log2FoldChange),]
    
    # Prepare matrix for GSEA
    gene_list <- ordered_df$log2FoldChange
    names(gene_list) <- ordered_df$ortholog_name
    
    # Prepare matrix for gsePathway
    ordered_entrez <- bitr(ordered_df$ortholog_name, 'SYMBOL', 'ENTREZID', OrgDb = org.Hs.eg.db)
    entrez_genes <- ordered_df %>% left_join(ordered_entrez, by = c('ortholog_name' = 'SYMBOL'), relationship = 'many-to-many') %>% dplyr::select(ENTREZID, log2FoldChange)
    distinct_genes <- entrez_genes %>% distinct(ENTREZID, .keep_all = T)
    entrez_gene_list <<- distinct_genes$log2FoldChange
    names(entrez_gene_list) <<- distinct_genes$ENTREZID
    
    # Run GSEA
    gsea_results <- gseGO(
      gene_list,
      keyType = 'SYMBOL',
      OrgDb = org.Hs.eg.db,
      ont = 'BP',
      pvalueCutoff = 0.05,
      pAdjustMethod = 'BH',
      verbose = T,
      eps = 1e-300
    )
    
    # Assign the results to a variable including treatment and sampling_point in the name
    results_name <-
      paste0('gsea_results_', treatment, '_', sampling_point)
    assign(results_name, gsea_results, envir = .GlobalEnv)
    
    return(gsea_results)
  }

# gsea formatting starting from a DESeq results table and using only significantly differentially regulated genes
gsea_formatting_significant <-
  function(results_table, treatment, sampling_point) {
    # Wrangle the data using the improved_data_wrangling function
    wrangled_data <-
      improved_data_wrangling(results_table, treatment, sampling_point)
    
    # Select ortholog_name and log2FC, and remove NA values
    ortholog_fc <- wrangled_data %>%
      dplyr::select(ortholog_name, log2FC) %>%
      na.omit()
    
    # Convert gene symbols to Entrez IDs
    entrez_ids <-
      bitr(
        ortholog_fc$ortholog_name,
        fromType = 'SYMBOL',
        toType = 'ENTREZID',
        OrgDb = org.Hs.eg.db
      )
    
    # Merge ortholog_fc with entrez_ids and select relevant columns
    enrichment <- ortholog_fc %>%
      left_join(entrez_ids, by = c('ortholog_name' = 'SYMBOL')) %>%
      dplyr::select(ENTREZID, log2FC)

    # Prepare the enrichment_gsea matrix
    enrichment_gsea <- enrichment$log2FC
    names(enrichment_gsea) <- enrichment$ENTREZID
    
    # Prepare the enrichment_gsea matrix
    enrichment_gsea_symbol <<- ortholog_fc$log2FC
    names(enrichment_gsea_symbol) <<- enrichment$ortholog_name
    
    # Create a dynamic name for the enrichment_gsea object
    results_name <-
      paste0("enrichment_gsea_", treatment, "_", sampling_point)
    assign(results_name, enrichment_gsea, envir = .GlobalEnv)
    
    return(enrichment_gsea)
  }

# Helper function to display Venn diagram
display_venn <- function(x, ...) {
  library(VennDiagram)
  grid.newpage()
  venn_object <- venn.diagram(x, filename = NULL, ...)
  grid.draw(venn_object)
}
