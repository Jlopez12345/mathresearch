# source("C:/Users/Owner/Desktop/Merging mutation data.R", echo=TRUE)
# mutation_df = read.csv(paste0('Data/LIHC/Mutation Raw Data/',filename))
# source("C:/Users/Owner/Desktop/Merging mutation data.R", echo=TRUE)
# mutation_df = read.csv(paste0('Data/LIHC/Mutation Raw DATA/',filename))
# filename = list.files('Data/LIHC Mutation Raw DATA')[1]
# filename

# setwd('~/Research')
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
library(scales)

# Parameters --------------------------------------------------------------
RECOMBINE_ALL_MUTATIONS = FALSE
RECOMBINE_CLINICAL_DATA = FALSE

# Reading in clinical data ---------------------------------------------------------
df = read.csv("Data/Liver Cancer/LIHC.clin.merged.csv")
BLCA_df = read.csv("Data/Bladder Cancer/BLCA.clin.merged.csv")
colon_df = read.csv("Data/Colon/COAD.clin.merged.csv")

# Create indicator column for each df
df$CancerType = "Liver"
BLCA_df$CancerType = "Bladder"
colon_df$CancerType = "Colon"


# Find common column names between the three data frames
common_columns <- Reduce(intersect, list(colnames(df), colnames(BLCA_df), colnames(colon_df)))

# Subset each data frame to keep only the common columns
df1_common <- df[, common_columns, drop = FALSE]
df2_common <- BLCA_df[, common_columns, drop = FALSE]
df3_common <- colon_df[, common_columns, drop = FALSE]

# Bind the data frames together
combined_df <- rbind(df1_common, df2_common, df3_common)


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
  
  
  #this will read the file that was already created for all mutation information ^^
} else {
  mutations_LIHC_with_clinical_df = read.csv('Data/all_mutations_lihc.csv')
  mutations_BLCA_with_clinical_df = read.csv('Data/all_mutations_blca.csv')
  mutations_COAD_with_clinical_df = read.csv('Data/all_mutations_coad.csv')
}


#this will read the file that was already created for clinical information
if (RECOMBINE_CLINICAL_DATA) {
  
  common_columns <- Reduce(intersect, list(colnames( mutations_LIHC_with_clinical_df), colnames(mutations_BLCA_with_clinical_df), colnames(mutations_COAD_with_clinical_df)))
  
  # Subset the data frames and keep the common columns
  lihc_common <- mutations_LIHC_with_clinical_df[, common_columns, drop = FALSE]
  blca_common <- mutations_BLCA_with_clinical_df[, common_columns, drop = FALSE]
  coad_common <- mutations_COAD_with_clinical_df[, common_columns, drop = FALSE]
  
  # Keep only relevant columns from BLCA
  #mutations_LIHC_with_clinical_df = mutations_LIHC_with_clinical_df[, intersect(names(mutations_LIHC_with_clinical_df), names(mutations_BLCA_with_clinical_df))]
  #mutations_BLCA_with_clinical_df = mutations_BLCA_with_clinical_df[, intersect(names(mutations_LIHC_with_clinical_df), names(mutations_BLCA_with_clinical_df))]
  
  # Combine mutation files
  all_mutations_df = rbind(lihc_common, blca_common,coad_common)
  
  
  #fix Tumor Sample Barcode
  all_mutations_df$Truncated_Barcode <- ifelse(
    grepl("^[^-]*-[^-]*-([^-]*).*", all_mutations_df$Tumor_Sample_Barcode),
    sub("^[^-]*-[^-]*-([^-]*).*", "\\1", all_mutations_df$Tumor_Sample_Barcode),
    NA  # Assign NA if it doesn't match
  )
  
  write.csv(all_mutations_df,'Data/all_mutations.csv',row.names = FALSE)
  } else { 
    all_mutations_df = read.csv('Data/all_mutations.csv')
  }
  
# Combining Clinical and Mutation Data Frames  ----------------------------
  
clinical_id_df = merge(clinical_id_df, all_mutations_df,
                       by.x=patient_id_cols,
                       by.y='Truncated_Barcode',
                       all.x=TRUE)

  
# Isolates specific kinds of mutations
mutations_with_clinical_df = clinical_id_df %>%
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site'))

clinical_id_df = clinical_id_df %>%
  filter(Variant_Classification %in% c('Missense_Mutation', 'Nonsense_Mutation', 'Splice_Site'))

#Creation of New Column to determine whether they HAVE liver cancer or NOT (for 0 and 1) 
#clinical_id_df = clinical_id_df %>%
 # mutate(is_liver = as.numeric(CancerType=="Liver"))

#Creation of New Column to determine whether they HAVE FAT4 mutation or NOT
clinical_id_df = clinical_id_df %>%
  mutate(has_FAT4 = as.numeric(Hugo_Symbol =="FAT4"))

FAT4_BMI_df = clinical_id_df %>%
  select(Hugo_Symbol,has_FAT4,bmi,bmi_class,CancerType,patient.gender, all_of(patient_id_cols)) %>%
  filter(has_FAT4==1)

# Top 10 Mutated Genes ----------------------------------------------------
clinical_id_df %>%
  select(Hugo_Symbol) %>%       # Grab the Hugo_Symbol column
  table %>%                     # Count up how many time each occurs
  as.data.frame %>%             # Make it into a data frame (prettier)
  arrange(-Freq) %>%            # Sort by frequency
  head(10)


# #NUMBERS
# clinical_id_df %>%
#   group_by(CancerType, Hugo_Symbol) %>%
#   summarise(num_patients = n()) %>% 
#   group_by(CancerType) %>%
#   arrange(-num_patients) %>%    # Sort (descending) by num_patients
#   mutate(row_num = row_number()) %>% 
#   filter(row_num <= 10) %>% View

#PERCENTAGES
No_filter <- clinical_id_df %>%
  group_by(CancerType) %>%
  mutate(num_patients_cancer_type = n_distinct(patient.patient_id)) %>%
  group_by(CancerType, Hugo_Symbol) %>%
  summarise(pct_patients = n_distinct(patient.patient_id) / first(num_patients_cancer_type)) %>% 
  group_by(CancerType) %>%
  arrange(-pct_patients) %>%    # Sort (descending) by pct_patients
  mutate(row_num = row_number()) %>% 
  filter(row_num <= 5) %>% View

#OBESE TOP 10
filtered_data <- clinical_id_df %>%
  filter(bmi_class %in% c("obese", "overweight")) %>%  # Filter for BMI classes
  group_by(CancerType) %>%
  mutate(num_patients_cancer_type = n_distinct(patient.patient_id)) %>%
  group_by(CancerType, Hugo_Symbol) %>%
  summarise(pct_patients = n_distinct(patient.patient_id) / first(num_patients_cancer_type)) %>%
  group_by(CancerType) %>%
  arrange(-pct_patients) %>%    # Sort (descending) by pct_patients
  mutate(row_num = row_number()) %>%
   filter(row_num <= 5) %>% 
   ggplot() +
  geom_bar(aes(x=Hugo_Symbol, y=pct_patients, fill=CancerType), 
            stat='identity', position='dodge', width=0.5) +
   scale_y_continuous(labels=percent, limits=c(0, 1)) +
   coord_flip()



# Data manipulation
filtered_data <- mutations_with_clinical_df %>%
  filter(bmi_class %in% c("obese", "overweight")) %>%  # Filter for BMI classes
  group_by(CancerType) %>%
  mutate(num_patients_cancer_type = n_distinct(patient.patient_id)) %>%
  group_by(CancerType, Hugo_Symbol) %>%
  summarise(pct_patients = n_distinct(patient.patient_id) / first(num_patients_cancer_type), .groups = 'drop') %>%
  group_by(CancerType) %>%
  arrange(-pct_patients) %>%    # Sort (descending) by pct_patients
  mutate(row_num = row_number()) %>%
  filter(row_num <= 5)

# Create the plot
ggplot(filtered_data) +
  geom_bar(aes(x = Hugo_Symbol, y = pct_patients, fill = CancerType), 
           stat = 'identity', position = 'dodge', width = 0.5) +
  geom_text(aes(x = Hugo_Symbol, 
                y = pct_patients, 
                label = scales::percent(pct_patients)), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5, # Adjust vertical position
            size = 4) +  # Size of the text
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  coord_flip() +
  labs(title = "Top 10 Genes by Cancer Type for Obese and Overweight Patients",
       x = "Gene Symbol (Hugo)",
       y = "Percentage of Patients")


clinical_id_df %>%
  group_by(CancerType) %>%
  mutate(num_patients_cancer_type = n_distinct(patient.patient_id)) %>%
  group_by(CancerType, Hugo_Symbol) %>%
  summarise(pct_patients = n_distinct(patient.patient_id) / first(num_patients_cancer_type)) %>% 
  group_by(CancerType) %>%
  arrange(-pct_patients) %>%    # Sort (descending) by pct_patients
  mutate(row_num = row_number()) %>% 
  filter(row_num <= 5) %>%
  ggplot() +
  geom_bar(aes(x=Hugo_Symbol, y=pct_patients, fill=CancerType), 
           stat='identity', position='dodge', width=0.5) +
  scale_y_continuous(labels=percent, limits=c(0, 1)) +
  geom_text(aes(x = Hugo_Symbol, 
                y = pct_patients, 
                label = scales::percent(pct_patients)), 
            position = position_dodge(width = 0.5), 
            vjust = -0.5, # Adjust vertical position
            size = 4) +  # Size of the text
  scale_y_continuous(labels = percent, limits = c(0, 1)) +
  coord_flip()


# BMI Classes for Gene '__'   --------------------------------------------------------
plot_Hugo_Symbol = 'FAT4'
gene_BMI_df = clinical_id_df %>%
  filter(Hugo_Symbol == plot_Hugo_Symbol) %>%
  group_by(Hugo_Symbol, bmi_class,bmi,CancerType) %>%
  summarise(num_patients=n()) %>%
  mutate(total_num_patients=sum(num_patients)) %>%
  arrange(-total_num_patients)



# Bar Graph -------------------------------------------------------------------
# all_patient_clinical_FAT4 = FAT4_BMI_df %>%
#   group_by(patient.samples.sample.portions.portion.analytes.analyte.aliquots.aliquot.2.bcr_aliquot_barcode) %>%
#   summarize(bmi=first(bmi),has_FAT4=any(has_FAT4==1),bmi_class=first(bmi_class),)
# 
# #Plot FAT4 class against the target variable (has TTN y/n)
# ggplot(all_patient_clinical_FAT4, aes(bmi_class, fill = has_FAT4)) +
#   geom_bar() +
#   coord_flip()






# Chi-Square TEST---------------------------------------------

#Creation of New Column to determine whether they HAVE FAT4 mutation or NOT

#clinical_id_df = clinical_id_df %>%
 # mutate(has_FAT4 = as.numeric(Hugo_Symbol =="FAT4"))

#new_df <- clinical_id_df %>%
  #mutate(is_bladder=(is_liver==0)) %>%
  #mutate(group = case_when(
  # `has_FAT4` == 1 & `bmi` > 25 ~ "has FAT4 + obese",
  # `has_FAT4` == 1 & `bmi` < 25  ~ "has FAT4 + not obese",
  # `has_FAT4` == 0 & `bmi` > 25 ~ "no FAT4 + obese",
  # `has_FAT4` == 0 & `bmi` < 25 ~ "no FAT4 + not obese"
 # )) %>%
  #select(group, `is_bladder`)


#new_df %>%
 # group_by(group) %>%
  #summarise(bladder_patients=sum(is_bladder),
   #         non_bladder_patients=n()-sum(is_bladder))

#contingency_table <- table(new_df$group, new_df$is_bladder)

#clinical_id_df_obese <- clinical_id_df %>%
  #filter(bmi>25)

#table(clinical_id_df_obese$is_liver, clinical_id_df_obese$has_FAT4) %>% chisq.test
#contingency_table <- new_df %>%
  #filter('bmi'>25)
#group_by(group) %>%
  #summarise(bladder_patients=sum(is_bladder),
            #non_bladder_patients=n()-sum(is_bladder))

#print(contingency_table)

#gender vs bmi class table
#contingency_table <- table(clinical_id_df$patient.gender, clinical_id_df$bmi_class)




# Create new dataframe categorizing FAT4 status and BMI
#new_df <- clinical_id_df %>%
 # mutate(is_bladder = (is_liver == 0)) %>%
 # mutate(group = case_when(
#    `has_FAT4` == 1 & `bmi` > 25 ~ "has FAT4 + obese",
#    `has_FAT4` == 1 & `bmi` <= 25 ~ "has FAT4 + not obese",
#    `has_FAT4` == 0 & `bmi` > 25 ~ "no FAT4 + obese",
#    `has_FAT4` == 0 & `bmi` <= 25 ~ "no FAT4 + not obese"
#  )) %>%
#  select(group, `is_bladder`)

# Summarise the data for the chi-square test
#summary_df <- new_df %>%
#  group_by(group) %>%
#  summarise(
#    bladder_patients = sum(is_bladder),
#    non_bladder_patients = n() - sum(is_bladder)
# )

# Create contingency matrix based on bladder vs non-bladder patients
#contingency_matrix <- as.matrix(summary_df %>%
                                  #select(bladder_patients, non_bladder_patients) %>%
                                  #t()) # Transpose to fit the chi-square test

# Add appropriate row and column names
#rownames(contingency_matrix) <- c("Bladder", "Non-Bladder")
#colnames(contingency_matrix) <- summary_df$group

# Perform the Chi-Square test
#chi_square_result <- chisq.test(contingency_matrix)

# Print the result
#print(chi_square_result)



#obese_fat4_df = new_df %>%
#  group_by(group) %>%
#  filter('bmi'>25)
#  summarise(bladder_patients=sum(is_bladder),
#          non_bladder_patients=n()-sum(is_bladder))

# contingency_table <- new_df %>%
#   count(group) %>%
#   pivot_wider(names_from = group, values_from = n, values_fill = 0)

# contingency_table <- table(new_df$group, new_df$is_bladder)































# Cochran-Mantel-Haenszel Test --------------------------------------------
#partial_tables <- margin.table(clinical_id_df)
#partial_tables






