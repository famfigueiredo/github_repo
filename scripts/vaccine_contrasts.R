rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions\

# Loading packages ----
library('DESeq2')
library('limma')
library('tidyverse')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('BiocParallel')
MulticoreParam(10)


# Loading data ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsGroup_ensembl.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_group_ensembl.Rda')


summary(sampleTable_group_ensembl)
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)
results(ddsGroup_ensembl)

rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions. It uses setdiff to find the subset of objects in the global environment
# (as returned by ls()) that don't have mode function (as returned by lsf.str())

# 1wpc vs 10wpi within treatment ----
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

ptagrfp_1wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.1wpc - groupptagrfp.10wpi, levels = colnames(design))

res_ptagrfp_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ptagrfp_1wpc_vs_10wpi, parallel = T)

conu_1wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.1wpc - groupconu.10wpi, levels = colnames(design))

res_conu_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = conu_1wpc_vs_10wpi, parallel = T)

dnavaccine_1wpc_vs_10wpi <-
  makeContrasts(groupdnavaccine.1wpc - groupdnavaccine.10wpi, levels = colnames(design))

res_dnavaccine_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = dnavaccine_1wpc_vs_10wpi, parallel = T)

eomes_1wpc_vs_10wpi <-
  makeContrasts(groupeomes.1wpc - groupeomes.10wpi, levels = colnames(design))

res_eomes_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = eomes_1wpc_vs_10wpi, parallel = T)

gata3_1wpc_vs_10wpi <-
  makeContrasts(groupgata3.1wpc - groupgata3.10wpi, levels = colnames(design))

res_gata3_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = gata3_1wpc_vs_10wpi, parallel = T)

ivhd_1wpc_vs_10wpi <-
  makeContrasts(groupivhd.1wpc - groupivhd.10wpi, levels = colnames(design))

res_ivhd_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivhd_1wpc_vs_10wpi, parallel = T)

ivld_1wpc_vs_10wpi <-
  makeContrasts(groupivld.1wpc - groupivld.10wpi, levels = colnames(design))

res_ivld_1wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivld_1wpc_vs_10wpi, parallel = T)


## saving results files ##
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_1wpc_vs_10wpi'
)
obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}
##

# 4wpc vs 10wpi within treatment ----
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

ptagrfp_4wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.4wpc - groupptagrfp.10wpi, levels = colnames(design))

res_ptagrfp_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ptagrfp_4wpc_vs_10wpi, parallel = T)

conu_4wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.4wpc - groupconu.10wpi, levels = colnames(design))

res_conu_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = conu_4wpc_vs_10wpi, parallel = T)

dnavaccine_4wpc_vs_10wpi <-
  makeContrasts(groupdnavaccine.4wpc - groupdnavaccine.10wpi, levels = colnames(design))

res_dnavaccine_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = dnavaccine_4wpc_vs_10wpi, parallel = T)

eomes_4wpc_vs_10wpi <-
  makeContrasts(groupeomes.4wpc - groupeomes.10wpi, levels = colnames(design))

res_eomes_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = eomes_4wpc_vs_10wpi, parallel = T)

gata3_4wpc_vs_10wpi <-
  makeContrasts(groupgata3.4wpc - groupgata3.10wpi, levels = colnames(design))

res_gata3_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = gata3_4wpc_vs_10wpi, parallel = T)

ivhd_4wpc_vs_10wpi <-
  makeContrasts(groupivhd.4wpc - groupivhd.10wpi, levels = colnames(design))

res_ivhd_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivhd_4wpc_vs_10wpi, parallel = T)

ivld_4wpc_vs_10wpi <-
  makeContrasts(groupivld.4wpc - groupivld.10wpi, levels = colnames(design))

res_ivld_4wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivld_4wpc_vs_10wpi, parallel = T)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc_vs_10wpi'
)
obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

# 6wpc vs 10wpi within treatment ----

results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

ptagrfp_6wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.6wpc - groupptagrfp.10wpi, levels = colnames(design))

res_ptagrfp_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ptagrfp_6wpc_vs_10wpi, parallel = T)

conu_6wpc_vs_10wpi <-
  makeContrasts(groupptagrfp.6wpc - groupconu.10wpi, levels = colnames(design))

res_conu_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = conu_6wpc_vs_10wpi, parallel = T)

dnavaccine_6wpc_vs_10wpi <-
  makeContrasts(groupdnavaccine.6wpc - groupdnavaccine.10wpi, levels = colnames(design))

res_dnavaccine_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = dnavaccine_6wpc_vs_10wpi, parallel = T)

eomes_6wpc_vs_10wpi <-
  makeContrasts(groupeomes.6wpc - groupeomes.10wpi, levels = colnames(design))

res_eomes_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = eomes_6wpc_vs_10wpi, parallel = T)

gata3_6wpc_vs_10wpi <-
  makeContrasts(groupgata3.6wpc - groupgata3.10wpi, levels = colnames(design))

res_gata3_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = gata3_6wpc_vs_10wpi, parallel = T)

ivhd_6wpc_vs_10wpi <-
  makeContrasts(groupivhd.6wpc - groupivhd.10wpi, levels = colnames(design))

res_ivhd_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivhd_6wpc_vs_10wpi, parallel = T)

ivld_6wpc_vs_10wpi <-
  makeContrasts(groupivld.6wpc - groupivld.10wpi, levels = colnames(design))

res_ivld_6wpc_vs_10wpi <-
  results(ddsGroup_ensembl, contrast = ivld_6wpc_vs_10wpi, parallel = T)

setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_6wpc_vs_10wpi'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


# 10wpc vs 10wpi within treatment ----

results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

conu_10wpc_vs_10wpi <-
  makeContrasts(groupconu.10wpc - groupconu.10wpi, levels = colnames(design))  # changed this formula on 23/04/2024. It was groupptagrfp.10wpc - groupconu.10wpi

res_ptagrfp_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","ptagrfp.10wpc","ptagrfp.10wpi"), parallel = T)

res_conu_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","conu.10wpc","conu.10wpi"), parallel = T)

res_dnavaccine_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","dnavaccine.10wpc","conu.10wpi"), parallel = T)

res_eomes_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","eomes.10wpc","eomes.10wpi"), parallel = T)

res_gata3_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","gata3.10wpc","gata3.10wpi"), parallel = T)

res_ivhd_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","ivhd.10wpc","ivhd.10wpi"), parallel = T)

res_ivld_10wpc_vs_10wpi <- results(ddsGroup_ensembl, contrast=c("group","ivld.10wpc","ivld.10wpi"), parallel = T)


setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc_vs_10wpi'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

# At 1, 4, and 10wpc across treatment
## GATA3 vs pTag and/or CONU
## EOMES vs pTag and/or CONU
# Establish which vaccines worked best in terms of growth performance, antibody production and histoscore.
# After, look for genetic markers of vaccine efficacy.
# Create a timeline of virus response through the sampling points, for all vaccine groups + controls
# Compare relevant vaccine groups within sampling point. Use IVLD as a sort of within sampling point control?

# Treatment contrasts within sampling points - 1wpc ----
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

## Using CONU as the reference ----
### IVLD
res_ivld_vs_conu_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'ivld.1wpc', 'conu.1wpc'), parallel = T)

### IVHD
res_ivhd_vs_conu_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'ivhd.1wpc', 'conu.1wpc'), parallel = T)

### GATA3 
res_gata3_vs_conu_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'gata3.1wpc', 'conu.1wpc'), parallel = T)

### EOMES
res_eomes_vs_conu_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'eomes.1wpc', 'conu.1wpc'), parallel = T)

### DNA vaccine
res_dnavaccine_vs_conu_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'dnavaccine.1wpc', 'conu.1wpc'), parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

## Using pTagRFP as the reference ----
### IVLD
res_ivld_vs_ptag_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'ivld.1wpc', 'ptagrfp.1wpc'), parallel = T)
  
### IVHD
res_ivhd_vs_ptag_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'ivhd.1wpc', 'ptagrfp.1wpc'), parallel = T)

### GATA3 
res_gata3_vs_ptag_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'gata3.1wpc', 'ptagrfp.1wpc'), parallel = T)

### EOMES
res_eomes_vs_ptag_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'eomes.1wpc', 'ptagrfp.1wpc'), parallel = T)

### DNA vaccine
res_dnavaccine_vs_ptag_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'dnavaccine.1wpc', 'ptagrfp.1wpc'), parallel = T)


## Using IV-LD as the reference ----
### IVHD
res_ivhd_vs_ivld_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'ivhd.1wpc', 'ivld.1wpc'), parallel = T)

### GATA3 
res_gata3_vs_ivld_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'gata3.1wpc', 'ivld.1wpc'), parallel = T)

### EOMES
res_eomes_vs_ivld_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'eomes.1wpc', 'ivld.1wpc'), parallel = T)

### DNA vaccine
res_dnavaccine_vs_ivld_1wpc <- results(ddsGroup_ensembl, contrast = c('group', 'dnavaccine.1wpc', 'ivld.1wpc'), parallel = T)


# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_1wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


# Treatment contrasts within sampling points - 4wpc ----
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

## Using CONU as the reference ----
### IVLD
ivld_vs_conu_4wpc <-
  makeContrasts(groupivld.4wpc - groupconu.4wpc, levels = colnames(design))
res_ivld_vs_conu_4wpc <-
  results(ddsGroup_ensembl, contrast = ivld_vs_conu_4wpc, parallel = T)

### IVHD
ivhd_vs_conu_4wpc <-
  makeContrasts(groupivhd.4wpc - groupconu.4wpc, levels = colnames(design))
res_ivhd_vs_conu_4wpc <-
  results(ddsGroup_ensembl, contrast = ivhd_vs_conu_4wpc, parallel = T)

### GATA3 
gata3_vs_conu_4wpc <-
  makeContrasts(groupgata3.4wpc - groupconu.4wpc, levels = colnames(design))
res_gata3_vs_conu_4wpc <-
  results(ddsGroup_ensembl, contrast = gata3_vs_conu_4wpc, parallel = T)

### EOMES
eomes_vs_conu_4wpc <-
  makeContrasts(groupeomes.4wpc - groupconu.4wpc, levels = colnames(design))
res_eomes_vs_conu_4wpc <-
  results(ddsGroup_ensembl, contrast = eomes_vs_conu_4wpc, parallel = T)

### DNA vaccine
dnavaccine_vs_conu_4wpc <-
  makeContrasts(groupdnavaccine.4wpc - groupconu.4wpc, levels = colnames(design))
res_dnavaccine_vs_conu_4wpc <-
  results(ddsGroup_ensembl, contrast = dnavaccine_vs_conu_4wpc, parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

## Using pTagRFP as the reference ----
### IVLD
ivld_vs_ptag_4wpc <-
  makeContrasts(groupivld.4wpc - groupptagrfp.4wpc, levels = colnames(design))
res_ivld_vs_ptag_4wpc <-
  results(ddsGroup_ensembl, contrast = ivld_vs_ptag_4wpc, parallel = T)

### IVHD
ivhd_vs_ptag_4wpc <-
  makeContrasts(groupivhd.4wpc - groupptagrfp.4wpc, levels = colnames(design))
res_ivhd_vs_ptag_4wpc <-
  results(ddsGroup_ensembl, contrast = ivhd_vs_ptag_4wpc, parallel = T)

### GATA3 
gata3_vs_ptag_4wpc <-
  makeContrasts(groupgata3.4wpc - groupptagrfp.4wpc, levels = colnames(design))
res_gata3_vs_ptag_4wpc <-
  results(ddsGroup_ensembl, contrast = gata3_vs_ptag_4wpc, parallel = T)

### EOMES
eomes_vs_ptag_4wpc <-
  makeContrasts(groupeomes.4wpc - groupptagrfp.4wpc, levels = colnames(design))
res_eomes_vs_ptag_4wpc <-
  results(ddsGroup_ensembl, contrast = eomes_vs_ptag_4wpc, parallel = T)

### DNA vaccine
dnavaccine_vs_ptag_4wpc <-
  makeContrasts(groupdnavaccine.4wpc - groupptagrfp.4wpc, levels = colnames(design))
res_dnavaccine_vs_ptag_4wpc <-
  results(ddsGroup_ensembl, contrast = dnavaccine_vs_ptag_4wpc, parallel = T)

## Using IV-LD as the reference ----
### IVHD
ivhd_vs_ivld_4wpc <-
  makeContrasts(groupivhd.4wpc - groupivld.4wpc, levels = colnames(design))
res_ivhd_vs_ivld_4wpc <-
  results(ddsGroup_ensembl, contrast = ivhd_vs_ivld_4wpc, parallel = T)

### GATA3 
gata3_vs_ivld_4wpc <-
  makeContrasts(groupgata3.4wpc - groupivld.4wpc, levels = colnames(design))
res_gata3_vs_ivld_4wpc <-
  results(ddsGroup_ensembl, contrast = gata3_vs_ivld_4wpc, parallel = T)

### EOMES
eomes_vs_ivld_4wpc <-
  makeContrasts(groupeomes.4wpc - groupivld.4wpc, levels = colnames(design))
res_eomes_vs_ivld_4wpc <-
  results(ddsGroup_ensembl, contrast = eomes_vs_ivld_4wpc, parallel = T)

### DNA vaccine
dnavaccine_vs_ivld_4wpc <-
  makeContrasts(groupdnavaccine.4wpc - groupivld.4wpc, levels = colnames(design))
res_dnavaccine_vs_ivld_4wpc <-
  results(ddsGroup_ensembl, contrast = dnavaccine_vs_ivld_4wpc, parallel = T) 

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_4wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


# Treatment contrasts within sampling points - 10wpc ----
results_names <- resultsNames(ddsGroup_ensembl)
design <- model.matrix( ~ group, sampleTable_group_ensembl)

## Using CONU as the reference ---- RE-WROTE THESE RESULTS EXTRACTIONS ON 23/04/2024 after figuring out makeContrasts was not picking up group.conu10wpc

## IVLD
res_ivld_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","ivld.10wpc","conu.10wpc"))

## IVHD
res_ivhd_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","ivhd.10wpc","conu.10wpc"))

### GATA3 
res_gata3_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","gata3.10wpc","conu.10wpc"))

### EOMES
res_eomes_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","eomes.10wpc","conu.10wpc"))

### DNA vaccine
res_dnavaccine_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","dnavaccine.10wpc","conu.10wpc"))
  
  
# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

## Using pTagRFP as the reference ----
### IVLD
ivld_vs_ptag_10wpc <-
  makeContrasts(groupivld.10wpc - groupptagrfp.10wpc, levels = colnames(design))
res_ivld_vs_ptag_10wpc <-
  results(ddsGroup_ensembl, contrast = ivld_vs_ptag_10wpc, parallel = T)

### IVHD
ivhd_vs_ptag_10wpc <-
  makeContrasts(groupivhd.10wpc - groupptagrfp.10wpc, levels = colnames(design))
res_ivhd_vs_ptag_10wpc <-
  results(ddsGroup_ensembl, contrast = ivhd_vs_ptag_10wpc, parallel = T)

### GATA3 
gata3_vs_ptag_10wpc <-
  makeContrasts(groupgata3.10wpc - groupptagrfp.10wpc, levels = colnames(design))
res_gata3_vs_ptag_10wpc <-
  results(ddsGroup_ensembl, contrast = gata3_vs_ptag_10wpc, parallel = T)

### EOMES
eomes_vs_ptag_10wpc <-
  makeContrasts(groupeomes.10wpc - groupptagrfp.10wpc, levels = colnames(design))
res_eomes_vs_ptag_10wpc <-
  results(ddsGroup_ensembl, contrast = eomes_vs_ptag_10wpc, parallel = T)

### DNA vaccine
dnavaccine_vs_ptag_10wpc <-
  makeContrasts(groupdnavaccine.10wpc - groupptagrfp.10wpc, levels = colnames(design))
res_dnavaccine_vs_ptag_10wpc <-
  results(ddsGroup_ensembl, contrast = dnavaccine_vs_ptag_10wpc, parallel = T)

## Using IV-LD as the reference ----
### IVHD
ivhd_vs_ivld_10wpc <-
  makeContrasts(groupivhd.10wpc - groupivld.10wpc, levels = colnames(design))
res_ivhd_vs_ivld_10wpc <-
  results(ddsGroup_ensembl, contrast = ivhd_vs_ivld_10wpc, parallel = T)

### GATA3 
gata3_vs_ivld_10wpc <-
  makeContrasts(groupgata3.10wpc - groupivld.10wpc, levels = colnames(design))
res_gata3_vs_ivld_10wpc <-
  results(ddsGroup_ensembl, contrast = gata3_vs_ivld_10wpc, parallel = T)

### EOMES
eomes_vs_ivld_10wpc <-
  makeContrasts(groupeomes.10wpc - groupivld.10wpc, levels = colnames(design))
res_eomes_vs_ivld_10wpc <-
  results(ddsGroup_ensembl, contrast = eomes_vs_ivld_10wpc, parallel = T)

### DNA vaccine
dnavaccine_vs_ivld_10wpc <-
  makeContrasts(groupdnavaccine.10wpc - groupivld.10wpc, levels = colnames(design))
res_dnavaccine_vs_ivld_10wpc <-
  results(ddsGroup_ensembl, contrast = dnavaccine_vs_ivld_10wpc, parallel = T) 

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/quantseq_dataAnalysis/deseq2_january2024/scripts/objects/results_10wpc'
)

obj <- ls(pattern = 'res_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

























# Treatment contrasts within sampling points - 10wpi