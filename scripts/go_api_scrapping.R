library(httr)
library(jsonlite)
library(xml2)

# Improved get GO data function
get_GO_data <- function(search_term) {
  requestURL <- paste0(
    "https://www.ebi.ac.uk/QuickGO/services/ontology/go/search?query=",
    URLencode(search_term),  # encodes special characters in the search term to avoid API errors
    "&limit=600&page=1"
  )
  
  # Make API request and handle potential errors
  result <- tryCatch({
    r <- GET(requestURL, accept("application/json"))
    if (status_code(r) != 200) {
      stop("API request failed with status code: ", status_code(r))
    }
    
    response_content <- content(r, as = "text", encoding = "UTF-8")
    json_data <- fromJSON(response_content, flatten = TRUE)
    data_frame <- as.data.frame(json_data$results)
    
    if (nrow(data_frame) == 0) {
      return(tibble(
        goterms = NA_character_,
        description = NA_character_,
        ontology = NA_character_,
        search_term = search_term
      ))
    }
    
    tibble_data <- as_tibble(data_frame)
    
    # Filter and select relevant columns
    tibble_data %>%
      filter(isObsolete == FALSE, aspect == "biological_process") %>%
      transmute(
        goterms = id,
        description = name,
        ontology = aspect,
        search_term = search_term
      )
    
  }, error = function(e) {
    # Return empty tibble in case of an error
    tibble(
      goterms = NA_character_,
      description = NA_character_,
      ontology = NA_character_,
      search_term = search_term
    )
  })
  
  return(result)
}


# immune_terms <- c('antigen', 'response', 'environmental%20stimulus', 'JAK%20STAT', 'CD4', 'CD8', 'apoptosis', 'acquired%20immune%20response', 'defense%20response', 'MDA-5', 'innate%20immunity', 'adaptive%20immunity', 'pathogen', 'vaccine%20response', 'toll-like', 't-helper', 'inflammation', 'immune', 'cytokine', 'cytotoxic', 'immunoglobulin', 'viral', 'viral%20load', 'lymphocyte', 'leukocyte', 'interleukin', 'interferon', 'dendritic', 'effector', 'chemokine', 'humoral', 'histocompatibility', 'b%20cell', 't%20cell', 'phagocytosis', 'innate')

## Updated 07/12/2024 with input from ChatGPT
immune_terms <- c('antibody', 'autoimmunity', 'CD28', 'checkpoint inhibition', 'complement', 'cross-presentation', 'fc receptor', 'granulocyte', 'helper', 'inflammasome', 'lymphoid', 'histocompatibility', 'memory', 'natural killer', 'pathogen-associated', 'pattern recognition', 'plasma cell', 'antigen', 'response', 'environmental stimulus', 'JAK STAT', 'CD4', 'CD8', 'apoptosis', 'acquired immune response', 'defense response', 'MDA-5', 'innate immunity', 'adaptive immunity', 'pathogen', 'vaccine response', 'toll-like', 't-helper', 'inflammation', 'immune', 'cytokine', 'cytotoxic', 'immunoglobulin', 'viral', 'viral load', 'lymphocyte', 'leukocyte', 'interleukin', 'interferon', 'dendritic', 'effector', 'chemokine', 'humoral', 'histocompatibility', 'b cell', 't cell', 'phagocytosis', 'innate')
##

immune_terms <- sort(immune_terms)

results <- lapply(immune_terms, get_GO_data)

immune_related_GOterms <- bind_rows(results)

nrow(immune_related_GOterms)

write_tsv(immune_related_GOterms, '~/Documents/PhD/Papers/Paper III/data/GOterms/immune_related_GOterms.tsv')

immune_related_GOterms <-
  read.table(
    '~/Documents/PhD/Papers/Paper III/data/GOterms/immune_related_GOterms.tsv',
    header = T,
    sep = '\t'
  )



