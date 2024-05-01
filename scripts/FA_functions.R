b <-
  as.data.frame(subset(res_dnavaccine_vs_ptag, padj < 0.05)) %>% rownames_to_column(., var = 'id') %>% as_tibble()

# converting entrez accession numbers to entrez gene id, to make sure every gene name has the same format.
gene_id <-
  getBM(
    attributes = c('entrezgene_accession', 'entrezgene_id'),
    filters = 'entrezgene_accession',
    values = b$id,
    uniqueRows = T,
    mart = ensembl_ssalar
  )

gene_id$entrezgene_id <-
  as.character(gene_id$entrezgene_id)  # from int to chr so they can be 'mutated'


sig_genes <- dplyr::select(b,
              id,
              log2FoldChange,
              padj) %>%
  dplyr::arrange(., desc(log2FoldChange))

enrichment_gsea <-
  sig_genes$log2FoldChange  # converting to matrix for GSEA

names(enrichment_gsea) <-
  sig_genes$id  # adding gene ids for GSEA

enrichment_ora <<-
  sig_genes$id  # converting to matrix for ORA

enrichment_ora_upreg <<-
  sig_genes %>% dplyr::filter(log2FoldChange > 0) %>% pull(id)  # splitting between upregulated

enrichment_ora_downreg <<-
  sig_genes %>% dplyr::filter(log2FoldChange < 0) %>% pull(id)
