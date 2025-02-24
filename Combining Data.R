# source("C:/Users/Owner/Desktop/Merging mutation data.R", echo=TRUE)
# mutation_df = read.csv(paste0('Data/LIHC/Mutation Raw Data/',filename))
# source("C:/Users/Owner/Desktop/Merging mutation data.R", echo=TRUE)
# mutation_df = read.csv(paste0('Data/LIHC/Mutation Raw DATA/',filename))
# filename = list.files('Data/LIHC Mutation Raw DATA')[1]
# filename

#setwd('~/Research')
#personsalpc
#setwd("C:/Users/Owner/Desktop/Research")

#borrowedpc
setwd("C:/Users/jlopezco/Desktop/Research")

library(dplyr)
library(waldo)
library(sem)
library(tidyr)
library(tidyverse)
library(readr)
library(tidymodels)
library(ggplot2)

# Parameters --------------------------------------------------------------
RECOMBINE_ALL_MUTATIONS = FALSE
RECOMBINE_CLINICAL_DATA = FALSE

# Reading in clinical data ---------------------------------------------------------
df = read.csv("Data/Liver Cancer/LIHC.clin.merged.csv")
BLCA_df = read.csv("Data/Bladder Cancer/BLCA.clin.merged.csv")
colon_df = read.csv("Data/Colon/COAD.clin.merged.csv")
ESCA_df=(read.csv("Data/Esophageal/ESCA.clin.merged.csv"))

# Merging similar data ----------------------------------------------------
# Create indicator column for each cancer type
df$CancerType = "Liver"
BLCA_df$CancerType = "Bladder"
colon_df$CancerType = "Colon"
ESCA_df$CancerType = "Esophageal"

# Find common column names between the three data frames
common_columns <- Reduce(intersect, list(colnames(df), colnames(BLCA_df), colnames(colon_df), colnames(ESCA_df)))

# Subset each data frame to keep only the common columns
df1_common <- df[, common_columns, drop = FALSE]
df2_common <- BLCA_df[, common_columns, drop = FALSE]
df3_common <- colon_df[, common_columns, drop = FALSE]
df4_common <- ESCA_df[, common_columns, drop = FALSE]

# Bind the data frames together
combined_df <- rbind(df1_common, df2_common, df3_common, df4_common)


# Assigning BMI classes ---------------------------------------------------
df = combined_df %>%
  filter(!is.na(patient.height)) %>%
  filter(!is.na(patient.weight))
df$bmi = df$patient.weight/((df$patient.height)/100)^2
df$bmi_class = NA

df[df$bmi < 18.5, "bmi_class"] = "underweight"
df[(18.5 <= df$bmi) & (df$bmi < 25), "bmi_class"] = "normal"
df[(25 <= df$bmi) & (df$bmi < 30), "bmi_class"] = "overweight"
df[(df$bmi >= 30), "bmi_class"] = "obese"

# Create Clinical Data Frame ------------------------------------------------------

# Names of patient IDs
patient_id_cols = "patient.patient_id"

# Isolate only the columns we want --> BMI, BMI Class, and Patient ID
clinical_id_df = df[ ,c("bmi","bmi_class","CancerType", patient_id_cols,"patient.gender")]

# Clean up tumor sample barcodes
clinical_id_df = clinical_id_df %>%
  mutate(across(all_of(patient_id_cols), ~toupper(.)))

# Create Mutation Data Frame -------------------------------------------------------------
combine_all_mutations = function(folder, fname) {
  # List all the files in the folder
  all_mutation_files = list.files(folder)
  
  # Read in the first file so I have a data frame to add the other patients to
  f.temp = all_mutation_files[1]
  all_mutations_df = read.csv(paste0(folder, f.temp), sep="\t")
  
  # Read in the other patients and append their data
  i = 1
  for (filename in all_mutation_files[2:length(all_mutation_files)]) {
    mutation_df = read.csv(paste0(folder, filename), sep="\t")
    
    all_mutations_df = rbind(all_mutations_df, mutation_df)
    i = i + 1
    if (i %% 10 == 0) {
      print(paste0(i, '/', length(all_mutation_files)))
    }
  }
  
  # Save this file so that we don't have to recompute it next time
  write.csv(all_mutations_df, paste0('Data/', fname, '.csv'), row.names=FALSE)
  
  return(all_mutations_df)
}

# Combine files or load already combined files
if (RECOMBINE_ALL_MUTATIONS) {
  mutations_LIHC_with_clinical_df = combine_all_mutations('Data/Liver Cancer/LIHC Mutation Raw DATA/', fname='all_mutations_lihc')
  mutations_BLCA_with_clinical_df = combine_all_mutations('Data/Bladder Cancer/BLCA Mutation Raw Data/', fname='all_mutations_blca')
  mutations_COAD_with_clinical_df = combine_all_mutations('Data/Colon/COAD Mutation Raw Data/', fname='all_mutations_coad')
  mutations_ESCA_with_clinical_df = combine_all_mutations('Data/Esophageal/ESCA Mutation Raw Data/', fname='all_mutations_ESCA')
  
  #this will read the file that was already created for all mutation information ^^
} else {
  mutations_LIHC_with_clinical_df = read.csv('Data/all_mutations_lihc.csv')
  mutations_BLCA_with_clinical_df = read.csv('Data/all_mutations_blca.csv')
  mutations_COAD_with_clinical_df = read.csv('Data/all_mutations_coad.csv')
  mutations_ESCA_with_clinical_df = read.csv('Data/all_mutations_ESCA.csv')
}

#this will read the file that was already created for clinical information
if (RECOMBINE_CLINICAL_DATA) {
  
  common_columns <- Reduce(intersect, list(colnames( mutations_LIHC_with_clinical_df), colnames(mutations_BLCA_with_clinical_df), colnames(mutations_COAD_with_clinical_df),colnames(mutations_ESCA_with_clinical_df)))
  
  # Subset the data frames and keep the common columns
  lihc_common <- mutations_LIHC_with_clinical_df[, common_columns, drop = FALSE]
  blca_common <- mutations_BLCA_with_clinical_df[, common_columns, drop = FALSE]
  coad_common <- mutations_COAD_with_clinical_df[, common_columns, drop = FALSE]
  esca_common <- mutations_ESCA_with_clinical_df[, common_columns, drop = FALSE]
  
  # Combine mutation files
  all_mutations_df = rbind(lihc_common, blca_common,coad_common,esca_common)
  
  
  # Fix Tumor Sample Barcode
  all_mutations_df$Truncated_Barcode <- ifelse(
    grepl("^[^-]*-[^-]*-([^-]*).*", all_mutations_df$Tumor_Sample_Barcode),
    sub("^[^-]*-[^-]*-([^-]*).*", "\\1", all_mutations_df$Tumor_Sample_Barcode),
    NA  # Assign NA if it doesn't match
  )
  
  write.csv(all_mutations_df,'Data/all_mutations_df.csv',row.names = FALSE)
  
} else {
  all_mutations_df = read.csv('Data/all_mutations_df.csv')
}

# Combining Clinical and Mutation Data Frames  ----------------------------
clinical_mutation_df = merge(clinical_id_df, all_mutations_df,
                             by.x=patient_id_cols,
                             by.y='Truncated_Barcode',
                             all.x=TRUE)

# Isolates specific kinds of mutations
mutations_with_clinical_df = clinical_mutation_df %>%
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site'))

# DATA WITH COMBINED INFO -------------------------------------------------
# write.csv(all_mutations_df, paste0('Data/', fname, 'mutations_with_clinical.csv'), row.names=FALSE
mutations_clinical_df = read.csv('Data/mutations_with_clinical.csv')




