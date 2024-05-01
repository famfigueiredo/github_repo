library(httr)
library(jsonlite)
library(xml2)

# Loading QuickGO API scraping function ----
get_GO_data <- function(search.term) {
  requestURL <-
    paste0(
      "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=",
      search.term,
      "&limit=600&page=1"
    )
  
  r <- GET(requestURL, accept("application/json"))
  response_content <- toJSON(content(r))
  json_data <- fromJSON(response_content)
  data_frame <- as.data.frame(json_data)
  tibble_data <- as_tibble(apply(data_frame, 2, as.character))
  
  tbl <- tibble_data %>%
    dplyr::select(results.id,
                  results.isObsolete,
                  results.name,
                  results.definition,
                  results.aspect) %>%
    dplyr::filter(results.isObsolete == FALSE &
                    results.aspect == 'biological_process') %>%
    dplyr::mutate(., search_term = search.term) %>%
    dplyr::select(., results.id, results.name, results.aspect, search_term) %>%
    dplyr::rename(.,
                  goterms = results.id,
                  description = results.name,
                  ontology = results.aspect)
  
  assign(paste0(search.term, '_tbl'), tbl, envir = .GlobalEnv)
}


immune_terms <- c('antigen', 'response', 'environmental%20stimulus', 'JAK%20STAT', 'CD4', 'CD8', 'apoptosis', 'acquired%20immune%20response', 'defense%20response', 'MDA-5', 'innate%20immunity', 'adaptive%20immunity', 'pathogen', 'vaccine%20response', 'toll-like', 't-helper', 'inflammation', 'immune', 'cytokine', 'cytotoxic', 'immunoglobulin', 'viral', 'viral%20load', 'lymphocyte', 'leukocyte', 'interleukin', 'interferon', 'dendritic', 'effector', 'chemokine', 'humoral', 'histocompatibility', 'b%20cell', 't%20cell', 'phagocytosis', 'innate')

immune_terms <- sort(immune_terms)

lapply(immune_terms, get_GO_data)

# Concatenate tibbles
mget(list.files(pattern = '_tbl'), envir = .GlobalEnv)
list.files(pattern = '_tbl')

immune_GOterms_tibble <- bind_rows(mget(ls(pattern = 'tbl')))
nrow(immune_GOterms_tibble)

immune_related_GOterms <-
  immune_GOterms_tibble[!duplicated(immune_GOterms_tibble$goterms), ]  # removing duplicated GO terms
nrow(immune_related_GOterms)

rows_with_quotes <- grep("'", immune_related_GOterms$description)
print(immune_related_GOterms[rows_with_quotes, ])

immune_related_GOterms$description <- gsub("'", "", immune_related_GOterms$description)  # removing quotes from the description text, so it doesn't fuck with read.table


write_tsv(immune_related_GOterms, '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/immune_related_GOterms.tsv')


immune_related_GOterms <-
  read.table(
    '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/immune_related_GOterms.tsv',
    header = T,
    sep = '\t'
  )




