rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions

# Loading packages ----
library('DESeq2')
library('limma')
library('tidyverse')
library('clusterProfiler')
library('gprofiler2')
library('org.Hs.eg.db')
library('BiocParallel')
MulticoreParam(8) # setting number of available cores for 'parallel'

## Loading data - HEART ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_heart.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_group_heart.Rda')

rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions. It uses setdiff to find the subset of objects in the global environment
# (as returned by ls()) that don't have mode function (as returned by lsf.str())

## Treatment contrasts within sampling points - HEART - 1wpc ----

### IVLD
heart_res_ivld_vs_conu_1wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'ivld.1wpc', 'conu.1wpc'), type = 'ashr', parallel = T)

### IVHD
heart_res_ivhd_vs_conu_1wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'ivhd.1wpc', 'conu.1wpc'), type = 'ashr', parallel = T)

### GATA3 
heart_res_gata3_vs_conu_1wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'gata3.1wpc', 'conu.1wpc'), type = 'ashr', parallel = T)

### EOMES
heart_res_eomes_vs_conu_1wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'eomes.1wpc', 'conu.1wpc'), type = 'ashr', parallel = T)

### DNA vaccine
heart_res_dnavaccine_vs_conu_1wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'dnavaccine.1wpc', 'conu.1wpc'), type = 'ashr', parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_1wpc'
)

obj <- ls(pattern = '^heart_')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


## Treatment contrasts within sampling points - HEART - 4wpc ----
### IVLD
res_ivld_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### IVHD
res_ivhd_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### GATA3 
res_gata3_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### EOMES
res_eomes_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

### DNA vaccine
res_dnavaccine_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_heart, contrast = c('group', 'dnavaccine.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)


objects <- ls(pattern = '^res_.*_conu_4wpc')
for (obj in objects) {
  assign(paste0('heart_', obj), get(obj))
}
rm(list = objects)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/heart/results_4wpc'
)

obj <- ls(pattern = '^heart_.*4wpc')  # regex pattern matching files starting with heart and containing 4wpc
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


## Treatment contrasts within sampling points - HEART - 10wpc ----
results_names <- resultsNames(ddsGroup_ensembl)

## IVLD
res_ivld_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","ivld.10wpc","conu.10wpc"), parallel = T)

## IVHD
res_ivhd_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","ivhd.10wpc","conu.10wpc"), parallel = T)

### GATA3 
res_gata3_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","gata3.10wpc","conu.10wpc"), parallel = T)

### EOMES
res_eomes_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","eomes.10wpc","conu.10wpc"), parallel = T)

### DNA vaccine
res_dnavaccine_vs_conu_10wpc <- results(ddsGroup_ensembl, contrast=c("group","dnavaccine.10wpc","conu.10wpc"), parallel = T)
  
  
# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpc'
)

obj <- ls(pattern = '^res_.*10wpc')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}



## Treatment contrasts within sampling points - HEART - 10wpi ----

## IVLD
res_ivld_vs_conu_10wpi <- results(ddsGroup_ensembl, contrast=c("group","ivld.10wpi","conu.10wpi"), parallel = T)

## IVHD
res_ivhd_vs_conu_10wpi <- results(ddsGroup_ensembl, contrast=c("group","ivhd.10wpi","conu.10wpi"), parallel = T)

### GATA3 
res_gata3_vs_conu_10wpi <- results(ddsGroup_ensembl, contrast=c("group","gata3.10wpi","conu.10wpi"), parallel = T)

### EOMES
res_eomes_vs_conu_10wpi <- results(ddsGroup_ensembl, contrast=c("group","eomes.10wpi","conu.10wpi"), parallel = T)

### DNA vaccine
res_dnavaccine_vs_conu_10wpi <- results(ddsGroup_ensembl, contrast=c("group","dnavaccine.10wpi","conu.10wpi"), parallel = T)


# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/results_10wpi'
)

obj <- ls(pattern = '^res_.*10wpi')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}




###########################

# Loading data - SPLEEN ----

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_spleen.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_grouped_spleen.Rda')

summary(sampleTable_grouped_spleen)
resultsNames(ddsDGE_grouped_spleen)

rm(list = setdiff(ls(), lsf.str()))  # remove everything from Global except functions. It uses setdiff to find the subset of objects in the global environment
# (as returned by ls()) that don't have mode function (as returned by lsf.str())
## Treatment contrasts within sampling points - 10wpi ----
### Using CONU as the reference ----

## IVLD
spleen_res_ivld_vs_conu_10wpi <- results(ddsDGE_grouped_spleen, contrast=c("group","ivld.10wpi","conu.10wpi"), parallel = T)

## IVHD
spleen_res_ivhd_vs_conu_10wpi <- results(ddsDGE_grouped_spleen, contrast=c("group","ivhd.10wpi","conu.10wpi"), parallel = T)

### GATA3 
spleen_res_gata3_vs_conu_10wpi <- results(ddsDGE_grouped_spleen, contrast=c("group","gata3.10wpi","conu.10wpi"), parallel = T)

### EOMES
spleen_res_eomes_vs_conu_10wpi <- results(ddsDGE_grouped_spleen, contrast=c("group","eomes.10wpi","conu.10wpi"), parallel = T)

### DNA vaccine
spleen_res_dnavaccine_vs_conu_10wpi <- results(ddsDGE_grouped_spleen, contrast=c("group","dnavaccine.10wpi","conu.10wpi"), parallel = T)


# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_10wpi'
)

obj <- ls(pattern = '^spleen_.*10wpi')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}







## Treatment contrasts within sampling points - 4wpc ----
### Using CONU as the reference ----
### IVLD
spleen_res_ivld_vs_conu_4wpc <- results(ddsDGE_grouped_spleen, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), parallel = T)

### IVHD
spleen_res_ivhd_vs_conu_4wpc <- results(ddsDGE_grouped_spleen, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), parallel = T)

### GATA3 
spleen_res_gata3_vs_conu_4wpc <- results(ddsDGE_grouped_spleen, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), parallel = T)

### EOMES
spleen_res_eomes_vs_conu_4wpc <- results(ddsDGE_grouped_spleen, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), parallel = T)

### DNA vaccine
spleen_res_dnavaccine_vs_conu_4wpc <- results(ddsDGE_grouped_spleen, contrast = c('group', 'dnavaccine.4wpc', 'conu.4wpc'), parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/spleen/results_4wpc'
)

obj <- ls(pattern = '^spleen_.*4wpc')  # regex pattern matching files starting with spleen and containing 4wpc
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}




##########################

# LIVER ----

load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_liver.RData')
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/sampleTable_grouped_liver.Rda')

summary(sampleTable_grouped_liver)
resultsNames(ddsDGE_grouped_liver)

## Treatment contrasts within sampling points - 10wpi ----
### Using CONU as the reference ----

## IVLD
liver_res_ivld_vs_conu_10wpi <- results(ddsDGE_grouped_liver, contrast=c("group","ivld.10wpi","conu.10wpi"), parallel = T)

## IVHD
liver_res_ivhd_vs_conu_10wpi <- results(ddsDGE_grouped_liver, contrast=c("group","ivhd.10wpi","conu.10wpi"), parallel = T)

### GATA3 
liver_res_gata3_vs_conu_10wpi <- results(ddsDGE_grouped_liver, contrast=c("group","gata3.10wpi","conu.10wpi"), parallel = T)

### EOMES
liver_res_eomes_vs_conu_10wpi <- results(ddsDGE_grouped_liver, contrast=c("group","eomes.10wpi","conu.10wpi"), parallel = T)

### DNA vaccine
liver_res_dnavaccine_vs_conu_10wpi <- results(ddsDGE_grouped_liver, contrast=c("group","dnavaccine.10wpi","conu.10wpi"), parallel = T)


# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_10wpi'
)

obj <- ls(pattern = '^liver_.*10wpi')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


## Treatment contrasts within sampling points - 4wpc ----
### Using CONU as the reference ----
### IVLD
liver_res_ivld_vs_conu_4wpc <- results(ddsDGE_grouped_liver, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), parallel = T)

### IVHD
liver_res_ivhd_vs_conu_4wpc <- results(ddsDGE_grouped_liver, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), parallel = T)

### GATA3 
liver_res_gata3_vs_conu_4wpc <- results(ddsDGE_grouped_liver, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), parallel = T)

### EOMES
liver_res_eomes_vs_conu_4wpc <- results(ddsDGE_grouped_liver, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), parallel = T)

### DNA vaccine
liver_res_dnavaccine_vs_conu_4wpc <- results(ddsDGE_grouped_liver, contrast = c('group', 'dnavaccine.4wpc', 'conu.4wpc'), parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/liver/results_4wpc'
)

obj <- ls(pattern = '^liver_.*4wpc')  # regex pattern matching files starting with liver and containing 4wpc
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}

















#########################

# HEAD-KIDNEY ----
load('~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/RData/ddsDGE_grouped_hkidney.RData')

## Treatment contrasts within sampling points - 10wpi ----
### Using CONU as the reference ----

## IVLD
hkidney_res_ivld_vs_conu_10wpi <- results(ddsDGE_grouped_hkidney, contrast=c("group","ivld.10wpi","conu.10wpi"), parallel = T)

## IVHD
hkidney_res_ivhd_vs_conu_10wpi <- results(ddsDGE_grouped_hkidney, contrast=c("group","ivhd.10wpi","conu.10wpi"), parallel = T)

### GATA3 
hkidney_res_gata3_vs_conu_10wpi <- results(ddsDGE_grouped_hkidney, contrast=c("group","gata3.10wpi","conu.10wpi"), parallel = T)

### EOMES
hkidney_res_eomes_vs_conu_10wpi <- results(ddsDGE_grouped_hkidney, contrast=c("group","eomes.10wpi","conu.10wpi"), parallel = T)

### DNA vaccine
hkidney_res_dnavaccine_vs_conu_10wpi <- results(ddsDGE_grouped_hkidney, contrast=c("group","dnavaccine.10wpi","conu.10wpi"), parallel = T)


# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_10wpi'
)

obj <- ls(pattern = '^hkidney_.*10wpi')
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


## Treatment contrasts within sampling points - 4wpc ----
### Using CONU as the reference ----
### IVLD
hkidney_res_ivld_vs_conu_4wpc <- results(ddsDGE_grouped_hkidney, contrast = c('group', 'ivld.4wpc', 'conu.4wpc'), parallel = T)

### IVHD
hkidney_res_ivhd_vs_conu_4wpc <- results(ddsDGE_grouped_hkidney, contrast = c('group', 'ivhd.4wpc', 'conu.4wpc'), parallel = T)

### GATA3 
hkidney_res_gata3_vs_conu_4wpc <- results(ddsDGE_grouped_hkidney, contrast = c('group', 'gata3.4wpc', 'conu.4wpc'), parallel = T)

### EOMES
hkidney_res_eomes_vs_conu_4wpc <- results(ddsDGE_grouped_hkidney, contrast = c('group', 'eomes.4wpc', 'conu.4wpc'), parallel = T)

### DNA vaccine
hkidney_res_dnavaccine_vs_conu_4wpc <- results(ddsDGE_grouped_hkidney, contrast = c('group', 'dnavaccine.4wpc', 'conu.4wpc'), parallel = T)

hkidney_shrunk_dnavaccine_vs_conu_4wpc <- lfcShrink(ddsDGE_grouped_hkidney, contrast = c('group', 'dnavaccine.4wpc', 'conu.4wpc'), type = 'ashr', parallel = T)

# Saving results files
setwd(
  '~/Documents/PhD/Thesis/quantseq_dataAnalysis/deseq2_dataAnalysis_2024/results/hkidney/results_4wpc'
)

obj <- ls(pattern = '^hkidney_.*4wpc')  # regex pattern matching files starting with hkidney and containing 4wpc
for (i in 1:length(obj)) {
  save(list = (obj[i]),
       file = paste(obj[i], ".RData", sep = ""))
}


DESeq2::plotMA(hkidney_shrunk_dnavaccine_vs_conu_4wpc)
