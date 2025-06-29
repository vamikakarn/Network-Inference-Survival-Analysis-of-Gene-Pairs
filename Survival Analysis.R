#install libraries
install.packages("survival")
install.packages("survminer")

#Load the libraries
library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(survival)
library(survminer)
library(tidyverse)

#Download and prepare gene expression data
query <- GDCquery(
  project = "TCGA-KIRC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  access = "open"
)
GDCdownload(query = query)
counts <- GDCprepare(query = query)
count_matrix <- assay(counts, 'unstranded')
sample_metadata <- colData(counts)
sample_metadata$condition <- ifelse(sample_metadata$sample_type == "Primary Tumor", "Tumor", "Normal")
sample_metadata$condition <- as.factor(sample_metadata$condition)

#Create DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_metadata,
  design = ~ condition
)

#Apply variance-stabilizing transformation (VST)
vsd <- vst(dds, blind = TRUE)
vst_counts <- assay(vsd)

#Getting clinical data for TCGA-KIRC cohort 
clinical_kirc <- GDCquery_clinic("TCGA-KIRC")
any(colnames(clinical_kirc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
which(colnames(clinical_kirc) %in% c("vital_status", "days_to_last_follow_up", "days_to_death"))
clinical_kirc[,c(9,40,46)]
#Identify the variables associated with survival 
table(clinical_kirc$vital_status)
#Change certain values the way they are encoded
clinical_kirc$deceased <- ifelse(clinical_kirc$vital_status == "Alive", FALSE, TRUE)
clinical_kirc$overall_survival <- ifelse(clinical_kirc$vital_status == "Alive",
                                         clinical_kirc$days_to_last_follow_up,
                                         clinical_kirc$days_to_death)
gene_metadata <- as.data.frame(rowData(counts))

#Get data for respective gene and add gene metadata information to it
gene1 <- vst_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "CORRESPONDING GENE NAME")
gene2 <- vst_counts %>% 
  as.data.frame() %>% 
  rownames_to_column(var = 'gene_id') %>% 
  gather(key = 'case_id', value = 'counts', -gene_id) %>% 
  left_join(., gene_metadata, by = "gene_id") %>% 
  filter(gene_name == "CORRESPONDING GENE NAME")

# Rename counts columns to distinguish gene 1 and gene 2
gene1 <- gene1 %>% select(case_id, gene1_counts = counts)
gene2 <- gene2 %>% select(case_id, gene2_counts = counts)

#Differentiate the cases with low or high counts
gene1$strata <- ifelse(gene1$gene1_counts >= median(gene1$gene1_counts), "HIGH", "LOW")
gene2$strata <- ifelse(gene2$gene2_counts >= median(gene2$gene2_counts), "HIGH", "LOW")
gene1$case_id <- gsub('-01.*', '', gene1$case_id)
gene1$case_id <- gsub('-11.*', '', gene1$case_id)
gene2$case_id <- gsub('-01.*', '', gene2$case_id)
gene2$case_id <- gsub('-11.*', '', gene2$case_id)

#Merge with clinical data
gene1 <- merge(gene1, clinical_kirc, by.x = 'case_id', by.y = 'submitter_id')
gene2 <- merge(gene2, clinical_kirc, by.x = 'case_id', by.y = 'submitter_id')

#Combine the data for the two genes
kirc_combined <- merge(gene1, gene2, by = "case_id")
kirc_combined$combined_strata <- paste(kirc_combined$strata.x, kirc_combined$strata.y, sep = "_")

#Plot Kaplan-Meier
fit <- survfit(Surv(time = kirc_combined$overall_survival.x, event = kirc_combined$deceased.x) ~ combined_strata, data = kirc_combined)
ggsurvplot(fit, 
           data = kirc_combined, 
           pval = TRUE,                     # Add p-value
           risk.table = TRUE,                # Add risk table
           conf.int = TRUE,                  # Show confidence intervals
           palette = c("#E7B800", "#2E9FDF", "#D95F02", "#7570B3"), # Custom colors
           xlab = "Time (in months)",        # X-axis label
           ggtheme = theme_minimal())

#Calculate Hazard Ratio
cox_model <- coxph(Surv(time = kirc_combined$overall_survival.x, event = kirc_combined$deceased.x) ~ combined_strata, data = kirc_combined)
summary(cox_model)

