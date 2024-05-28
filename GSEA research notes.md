There are no explicitly immune related pathways highly ranked in GSEA at 1WPC. However, they show up prominently at 4WPC.

Can this be related to the fact that the heart is not an 'immune organ'? What is the case in the liver, head-kidney, and spleen?

However, in EOMES, more general GO terms like *negative regulation of response to external stimulus*, *positive regulation of response to external stimulus*, and *regulation of response to external stimulus* are present, with relatively large gene sets, generally upregulated. Is this evidence of an early antiviral response?

Responses to external stimuli show up in the DNA vaccine, EOMES, and GATA3, but not in IV-HD? Evidence of viral infection recognition?







From ChatGPT when prompted to explain a negative *Normalised Enrichment Score* (NES):

*In summary, a negative NES in GSEA indicates that the gene set is  significantly downregulated in the condition of interest compared to the control, reflecting a negative correlation with the phenotype being  studied. This information can provide valuable insights into the  biological processes that are less active or suppressed in the  experimental condition.*



### Interpretation of GSEA Results:

- **Positive NES:**
  - Indicates that the gene set is enriched at the top of the ranked list, suggesting that genes in the set are upregulated in the condition of interest.
- **Negative NES:**
  - Indicates that the gene set is enriched at the bottom of the ranked list, suggesting that genes in the set are downregulated in the condition of interest.







Most studies related to viral infection use qPCR to 'measure' gene expression of target genes in relation to a *housekeeping gene*. The approach of using RNAseq as a tool to study immune responses is inherently different, in the sense that it does not compare gene expression against a housekeeping gene. Instead, it requires using a *reference*, such as a control group. This control group provides the baseline gene expression against which all other groups are compared to.

Because with RNAseq one looks at thousands of genes simultaneously, it makes more sense to have a holistic view of gene expression. Any emerging patterns in gene expression tend to stand out during pathway analysis, since pathways, by definition, incorporate hundreds, or even thousands, of differentially regulated genes.



WHAT GENES ARE EXPECTED TO BE DIFFERENTIALLY REGULATED BY EOMES, GATA3, AND DNA VACCINE?

CREATE A TIMELINE PLOT OF GENE EXPRESSION, I.E. 1WPC VS 10WPI, 4WPC VS 1WPC, 6WPC VS 4WPC, 10WPC VS 6WPC FOR CONTROL



Immune-related genes

Housekeeping: β-actin (ACTB): AF012125, ENSSSAG00000116649. A lot of reads everywhere.

SasaIFN-α: AY216594, no Ensembl

Μx: U66475, ENSSSAG00070007970 (mxa). No hits in read files.

ISG15: AY926456, ENSSSAG00000105254 (UBIL). Mapped but mostly no reads.

IFN-γ: AY795563, ENSSSAG00000105299. Mapped but no reads.

ΤΝF-α (1&2): NM_001123589 & NM_001123590, ENSSSAG00000065312 (mapped but no reads) & ENSSSAG00000053783 (no hits)

IL12-β: BT049114, ENSSSAG00000009655. Mapped but with few reads.

IL-10: EF165029, ENSSSAG00000107544. Mapped but with no reads.

IL-8: NM_001140710, ENSSSAG00000006498. Mapped but with no reads.

CD3-ε: NM_001123622, ENSSSAG00000076824. A lot of reads in s and hk. Fewer reads in h.

CD4: NM_001124539, rainbow trout gene

CD8-α: NM_001123583, ENSSSAG00000065860. Mapped but with few reads.

TCR-α: BT050114, no Ensembl. ENSSSAG00070030723 for T-cell receptor alpha/delta variable 31.0. No hits.

MHC I: AF508864, no Ensembl. ENSSSAG00000077407 for a novel gene. A lot of hits in h, hk, and s. 

MHC II: BT049430, no Ensembl. ENSSSAG00000004635 for a novel gene. A LOT of hits everywhere.



How did you look at specific genes in the RNAseq data (genes in Table 1, for example)? Are these manually curated, and you get their respective raw counts from gene count files?

So these were manually curated and then extracted into a table where we included the FC values from Deseq2 and then used total counts and the house keeping gene EF1a to make a normalised expression. It's maybe not the best way of doing it, but was a useful proxy of the genes expression. (Shahmir Naseer)



About normalising expression

Expression values were obtained by normalising gene count against total raw counts and then the ratio to ELF1a multiplied by 1000, trout (ENSOMYG00000038328) and salmon (ENSSSAG00000077892) counts,

blank cells = no counts for the gene.
