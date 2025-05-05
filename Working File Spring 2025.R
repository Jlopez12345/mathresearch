setwd("C:/Users/jlopezco/Desktop/Research")

library(tidyverse)
library(samplesizeCMH)
library(graphics)

wdf = read.csv('Data/clinical_mutation_df.csv')

wdf$focuscancer <- ifelse(wdf$CancerType == "ESCA","Yes", "No")

wdf$bi_bmi <- ifelse(wdf$bmi_class %in% c("obese", "overweight"), 
                                "High", "Low")
#HIGH
wdf_highbmi = wdf %>% 
  filter(bi_bmi=="High") %>% 
  select(Hugo_Symbol, focuscancer)



#HIGH BMI & GENES CHI-SQUARE 
common_hugo_symbols = wdf_highbmi %>%
  group_by(Hugo_Symbol) %>%
  summarise(num_people = n()) %>%
  filter(num_people >= 100) %>%
  select(Hugo_Symbol)

wdf_highbmi = wdf_highbmi %>%
  filter(Hugo_Symbol %in% common_hugo_symbols$Hugo_Symbol)

#chi-square test on high bmi and ESCA cancer
wdf_highbmi_table = table(wdf_highbmi$Hugo_Symbol, wdf_highbmi$focuscancer)
chisq_results <- chisq.test(wdf_highbmi_table)


#check the residuals of each gene
observed <- wdf_highbmi_table
expected <- chisq_results$expected
standardized_residuals <- (observed - expected) / sqrt(expected)

#make into data frame
residuals_df <- as.data.frame(as.table(standardized_residuals))
colnames(residuals_df) <- c("Gene", "ESCA_Status", "Standardized_Residual")

#filter for genes with residuals >2
significant_residuals <- residuals_df %>% 
  filter((Standardized_Residual) > 1)


top_significant_residuals <- significant_residuals %>%
  arrange(desc(abs(Standardized_Residual))) %>%
  head(10)

top_significant_residuals


# Create a contingency table of genes and cancer status
gene_cancer_table <- table(wdf_highbmi$Hugo_Symbol, wdf_highbmi$focuscancer)

# Run the mosaic plot
mosaicplot(gene_cancer_table, shade = TRUE, las = 2, main = "Gene vs Cancer Status")

#LOW BMI
wdf_lowbmi = wdf %>%
  filter(bi_bmi=="Low")%>%
  select(Hugo_Symbol, focuscancer)

#LOW BMI & GENES CHI-SQUARE
common_hugo_symbols = wdf_lowbmi %>%
  group_by(Hugo_Symbol) %>%
  summarise(num_people = n()) %>%
  filter(num_people >= 100) %>%
  select(Hugo_Symbol)

wdf_lowbmi = wdf_lowbmi %>%
  filter(Hugo_Symbol %in% common_hugo_symbols$Hugo_Symbol)

wdf_lowbmi_table = table(wdf_lowbmi$Hugo_Symbol, wdf_lowbmi$focuscancer)

chisq.test(wdf_lowbmi_table)

wdf_lowbmi = wdf_lowbmi %>%
  filter(Hugo_Symbol %in% common_hugo_symbols$Hugo_Symbol)

#chi-square test on high bmi and ESCA cancer
wdf_lowbmi_table = table(wdf_lowbmi$Hugo_Symbol, wdf_lowbmi$focuscancer)
chisq_results <- chisq.test(wdf_lowbmi_table)

#check the residuals of each gene
observed <- wdf_lowbmi_table
expected <- chisq_results$expected
standardized_residuals <- (observed - expected) / sqrt(expected)

#make into data frame
residuals_df <- as.data.frame(as.table(standardized_residuals))
colnames(residuals_df) <- c("Gene", "ESCA_Status", "Standardized_Residual")

#filter for genes with residuals >2
significant_residuals <- residuals_df %>% 
  filter((Standardized_Residual) >1)

top_significant_residuals <- significant_residuals %>%
  arrange(desc(abs(Standardized_Residual))) %>%
  head(10)


print(head(top_significant_residuals))


