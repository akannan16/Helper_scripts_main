#APractice

setwd("Desktop/folder")

library(limma)
library(edgeR)
library(ggplot2)
library(pheatmap)
library(readr)
library(dplyr)
library(corrplot)
library(FactoMineR)
library(factoextra)
library(gplots)
library(RColorBrewer)
library(survival)
library(survminer)

# Load the data
rsem_data <- read.csv("tcga.brca.rsem.csv")
#View(rsem_data)

# Transposing the gene expression data
gene_expression_counts = t(rsem_data[, 4:ncol(rsem_data)])
colnames(gene_expression_counts) = rsem_data$bcr_patient_barcode


# Creating group factor based on the sample_type
group <- factor(rsem_data$sample_type)
group <- relevel(group, ref = "Adjacent Normal")
#print(str(group))
#print(table(group))

# Creating the design matrix
design <- model.matrix(~ group)
#print(head(design))

# Creating DGEList object
dge <- DGEList(counts = gene_expression_counts)
print(str(dge))
print(summary(dge))

# Filtering
keep = filterByExpr(dge, design)
dge = dge[keep, , keep.lib.sizes = FALSE]
print(summary(dge))

# Normalization (TMM)
dge = calcNormFactors(dge)
v = voom(dge, design, plot = TRUE)

# Fit model to data given design
fit = lmFit(v, design)
fit = eBayes(fit)

# Extract top100 genes
topGenes = topTable(fit, coef="groupPrimary Tumor", n=100, sort.by="p")
print(topGenes)

#PCA visualization for normalized data
pca <- prcomp(t(v$E))
df_pca <- data.frame(Sample = colnames(v$E), PC1 = pca$x[,1], PC2 = pca$x[,2], Group = group)
ggplot(df_pca, aes(x=PC1, y=PC2, color=Group)) + geom_point() + ggtitle("PCA Plot")

#PCA visualization for topGene data
top_genes_names <- rownames(topGenes)
expression_data_subset <- v$E[top_genes_names, ]
pca_top <- prcomp(t(expression_data_subset))
df_pca_top <- data.frame(Sample = colnames(expression_data_subset), PC1 = pca_top$x[,1], PC2 = pca_top$x[,2], Group = group)
ggplot(df_pca_top, aes(x=PC1, y=PC2, color=Group)) + geom_point() + ggtitle("topGenes Plot")

# Load clinical data
clinical_data <- read.csv("brca_tcga_clinical_data.csv")
names(clinical_data) <- gsub(" ", "_", names(clinical_data))
names(clinical_data) <- gsub("-", "_", names(clinical_data))
#View(clinical_data)

# Creating a consistent sample identifier since SampleIDs and PatientIDs differ in all datasets to merge the tables together
# Only considering the Top 100 Genes for this analysis
df_pca_top$ID_Match <- substr(df_pca_top$Sample, 1, 12)
clinical_data$ID_Match <- substr(clinical_data$Patient.ID, 1, 12)
merged_data <- merge(df_pca_top, clinical_data, by.x="ID_Match", by.y="ID_Match")
View(merged_data)
# Plotting the PCA with IHC.HER2 status
ggplot(merged_data, aes(x = PC1, y = PC2, color = IHC.HER2)) +
  geom_point() +
  theme_minimal() +
  ggtitle("PCA Plot with IHC.HER2 Status") +
  labs(x = "PC1", y = "PC2") +
  scale_color_manual(values = c("Negative" = "red", "Positive" = "green", "Equivocal" = "darkblue"))

#heatmap
# Extract the sample IDs from the column names of 'expression_data_subset'
# Load necessary library
library(gplots)

# Assuming 'expression_data_subset' and 'merged_data' have already been loaded and preprocessed

# Extract the sample IDs from the column names of 'expression_data_subset'
sample_ids <- colnames(expression_data_subset)

# Create a named vector of HER2 status using 'Sample' as names
her2_status_vector <- setNames(merged_data$IHC.HER2, merged_data$Sample)

# Align HER2 status with the order of 'expression_data_subset'
aligned_her2_status <- her2_status_vector[sample_ids]

# Order the expression_data_subset matrix by HER2 status
ordered_her2_status <- sort(aligned_her2_status, decreasing = TRUE)  # Sort so that "Positive" comes first
ordered_sample_ids <- names(ordered_her2_status)  # Get the sample IDs in the order of HER2 status
ordered_expression_data <- expression_data_subset[, ordered_sample_ids]  # Order the columns of expression data

# Create a color vector for HER2 status
her2_colors <- c("Negative" = "red", "Positive" = "green", "Equivocal" = "orange")
col_annotation_her2 <- her2_colors[as.character(ordered_her2_status)]

# Create a color palette for the heatmap
hmcol <- colorRampPalette(c("blue", "white", "yellow"))(100)

# Perform clustering on the ordered expression data
dist_genes <- dist(ordered_expression_data)
hc_genes <- hclust(dist_genes)
dist_samples <- dist(t(expression_data_subset))
hc_samples <- hclust(dist_samples)
# Create the heatmap with ordered samples by HER2 status
heatmap.2(as.matrix(ordered_expression_data),
          scale="row",
          col=hmcol,
          ColSideColors=col_annotation_her2,
          Rowv=as.dendrogram(hc_genes),  # Cluster genes
          Colv=FALSE,  # Do not cluster samples, keep the HER2 status order
          labRow=rownames(ordered_expression_data),
          labCol=ordered_sample_ids,
          dendrogram="row",
          trace="none",
          cexRow=0.5,
          cexCol=0.5,
          key=TRUE,
          key.title="Expression",
          main="Differentially Expressed Genes Heatmap by HER2 Status")



#Identifying HER2+ patients based on IHC-HER2 status
her2_positive_patients <- clinical_data %>% 
  filter(`IHC.HER2` == "Positive")
summary_table1 <- her2_positive_patients %>% 
  count(`Cancer.Type.Detailed`, `IHC.HER2`)

# Filtering HER2+ patients based on HER2 FISH status
her2_fish_analysis <- clinical_data %>% 
  filter(`HER2.fish.status` == "Positive")
summary_table2 <- her2_fish_analysis %>% 
  count(`Cancer.Type.Detailed`, `HER2.fish.status`)

# Filtering HER2+ patients based on ICH and FISH
her2_combined <- her2_positive_patients %>% 
  filter(`HER2.fish.status` == "Positive")
summary_table3 <- her2_combined %>% 
  count(`Cancer.Type.Detailed`, `IHC.HER2`, `HER2.fish.status`)

#Exploring sub-types
# Create LumA table and add subtype
LumA <- clinical_data %>% 
  filter(`IHC.HER2` == "Negative", ER.Status.By.IHC == "Positive", PR.status.by.ihc == "Positive") %>%
  count(`IHC.HER2`, `ER.Status.By.IHC`, `PR.status.by.ihc`) %>%
  mutate(subtype = "LumA")

# Create LumB table and add subtype
LumB <- clinical_data %>% 
  filter(`IHC.HER2` == "Positive", ER.Status.By.IHC == "Positive", PR.status.by.ihc == "Positive") %>%
  count(`IHC.HER2`, `ER.Status.By.IHC`, `PR.status.by.ihc`) %>%
  mutate(subtype = "LumB")

# Create HER2_high table and add subtype
HER2_high <- clinical_data %>% 
  filter(`IHC.HER2` == "Positive", ER.Status.By.IHC == "Negative", PR.status.by.ihc == "Negative") %>%
  count(`IHC.HER2`, `ER.Status.By.IHC`, `PR.status.by.ihc`) %>%
  mutate(subtype = "HER2_high")

# Create TNBC table and add subtype
TNBC <- clinical_data %>% 
  filter(`IHC.HER2` == "Negative", ER.Status.By.IHC == "Negative", PR.status.by.ihc == "Negative") %>%
  count(`IHC.HER2`, `ER.Status.By.IHC`, `PR.status.by.ihc`) %>%
  mutate(subtype = "TNBC")

# Combine all tables into one
combined_summary_table <- bind_rows(LumA, LumB, HER2_high, TNBC)

# View the combined summary table
View(combined_summary_table)

# Visualize the data
#IHC/HER2 status
ggplot(data = summary_table1, aes(x = `Cancer.Type.Detailed`, y = n, fill = `Cancer.Type.Detailed`)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "HER2+ Patient Distribution by Cancer Sub Types ",
       x = "Cancer Type Detailed",
       y = "Number of Patients") + theme(axis.text.x = element_text(angle = 90, hjust = 1), plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette = "Set1")
#HER2/FISH status
ggplot(data = summary_table2, aes(x = `Cancer.Type.Detailed`, y = n, fill = `Cancer.Type.Detailed`)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "HER2+ Patient Distribution by FISH Status",
       x = "Cancer Type Detailed",
       y = "Number of Patients") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5) )+ scale_fill_brewer(palette = "Set2")

#HER2/IHC-FISH status
ggplot(data = summary_table3, aes(x = `Cancer.Type.Detailed`, y = n, fill = `Cancer.Type.Detailed`)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "HER2+ Patient Distribution by Combined IHC/FISH Status",
       x = "Cancer Type Detailed",
       y = "Number of Patients") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette = "Set3")

#Sub-type exploration
ggplot(data = combined_summary_table, aes(x = `subtype`, y = n, fill = `subtype`)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(title = "Cancer Sub-types in clinical data",
       x = "Cancer Sub-types",
       y = "Number of Patients") + theme(axis.text.x = element_text(angle = 45, hjust = 1), plot.title = element_text(hjust = 0.5)) + scale_fill_brewer(palette = "RdYlGn")


# Survival Analysis
install.packages("dplyr")
library(dplyr)
survival_data2 <- clinical_data %>%
  select(Sample.Type, Sample.ID, Diagnosis.Age,
         Overall.Survival..Months., Overall.Survival.Status, IHC.HER2, 
         ER.Status.By.IHC, PR.status.by.ihc, Sex) %>%
  mutate(
    Subtype = case_when(
      ER.Status.By.IHC == "Positive" & PR.status.by.ihc == "Positive" & IHC.HER2 == "Negative" ~ "LumA",
      ER.Status.By.IHC == "Positive" & PR.status.by.ihc == "Positive" & IHC.HER2 == "Positive" ~ "LumB",
      ER.Status.By.IHC == "Negative" & PR.status.by.ihc == "Negative" & IHC.HER2 == "Positive" ~ "Her2-high",
      ER.Status.By.IHC == "Negative" & PR.status.by.ihc == "Negative" & IHC.HER2 == "Negative" ~ "Triple Neg BC",
      TRUE ~ NA_character_  # Exclude other subtypes
    )
  ) %>%
  filter(!is.na(Subtype))  

# Convert Overall Survival Status to a binary format
survival_data2$Overall_Survival_Status_Binary <- ifelse(survival_data2$Overall.Survival.Status == "1:DECEASED", 1, 0)

# Create a Surv object
surv_obj <- Surv(time = survival_data2$`Overall.Survival..Months.`, event = survival_data2$Overall_Survival_Status_Binary)

# Fit a Kaplan-Meier survival curve stratified by Subtype
km_fit <- survfit(surv_obj ~ Subtype, data = survival_data2)

# Plot the Kaplan-Meier curve for LumA vs LumB
ggsurvplot(km_fit, data = survival_data2, risk.table = TRUE, risk.table.height = 0.30, xlab = "Months", ylab = "Survival Probability", 
           title = "Kaplan-Meier Survival Curve: Subtype")


#Copy Number Analysis
# Load the CNV data
cnv_data <- read.csv("brca_tcga_erbb2_copy_number.csv")
cnv_data$ID_Match <- substr(cnv_data$sample_id, 1, 12)
# Now merge the two data frames based on the new ID_Match column
cnv_merged_data <- merge(merged_data, cnv_data, by.x = "ID_Match", by.y = "ID_Match")
cnv_merged_data_selected <- select(cnv_merged_data, Sample, PC1, PC2, Group, IHC.HER2, erbb2_copy_number)
View(cnv_merged_data_selected)
sum(is.na(cnv_merged_data_selected))

# Impute missing values or remove rows with missing values
miss_IHC_HER2 <- names(which.max(table(cnv_merged_data_selected$IHC.HER2, useNA = "no")))
cnv_merged_data_selected$IHC.HER2[is.na(cnv_merged_data_selected$IHC.HER2)] <- miss_IHC_HER2

# Now the data is ready for comparative analysis
# Save the cleaned data if necessary
write.csv(cnv_merged_data_selected, "cleaned_data.csv", row.names = FALSE)

# Print out the cleaned data to verify
head(cnv_merged_data_selected)

# Step a. Descriptive Statistics
# Calculate descriptive statistics for numerical variables
summary(cnv_merged_data_selected$PC1)
summary(cnv_merged_data_selected$PC2)
summary(cnv_merged_data_selected$erbb2_copy_number)

# Calculate frequencies for categorical variables
table(cnv_merged_data_selected$Group)
table(cnv_merged_data_selected$IHC.HER2)

# Step b. Group Comparison
anova_pc1_group <- aov(PC1 ~ Group, data = cnv_merged_data_selected)
summary(anova_pc1_group)

anova_pc1_ihc <- aov(PC1 ~ IHC.HER2, data = cnv_merged_data_selected)
summary(anova_pc1_ihc)


anova_pc2_group <- aov(PC2 ~ Group, data = cnv_merged_data_selected)
summary(anova_pc2_group)

anova_pc2_ihc <- aov(PC2 ~ IHC.HER2, data = cnv_merged_data_selected)
summary(anova_pc2_ihc)

anova_erbb2_group <- aov(erbb2_copy_number ~ Group, data = cnv_merged_data_selected)
summary(anova_erbb2_group)

anova_erbb2_ihc <- aov(erbb2_copy_number ~ IHC.HER2, data = cnv_merged_data_selected)
summary(anova_erbb2_ihc)

#Check distribution
hist(cnv_merged_data_selected$PC1)
hist(cnv_merged_data_selected$PC2)
hist(cnv_merged_data_selected$erbb2_copy_number)

shapiro.test(cnv_merged_data_selected$PC1)
shapiro.test(cnv_merged_data_selected$PC2)
shapiro.test(cnv_merged_data_selected$erbb2_copy_number)

# Step c. Correlation Analysis
# Calculate Spearman's rank correlation matrix
spearman_correlation <- cor(cnv_merged_data_selected[, c("PC1", "PC2", "erbb2_copy_number")], 
                            method = "spearman", 
                            use = "complete.obs")

# Print the correlation matrix
print(spearman_correlation)

# Install and load the ggplot2 package for visualization (if not already installed)
# install.packages("ggplot2")
library(ggplot2)
library(ggcorrplot)

# Visualize the correlation matrix
ggcorrplot(spearman_correlation, hc.order = TRUE, type = "lower", lab = TRUE)


cnv_merged_data_selected$IHC.HER2_binary <- ifelse(cnv_merged_data_selected$IHC.HER2 == "Positive", 1, 0)

library(pROC)
# ROC analysis for ERBB2 copy number
roc_dna <- roc(cnv_merged_data_selected$IHC.HER2_binary, cnv_merged_data_selected$erbb2_copy_number)
auc_dna <- auc(roc_dna)

# ROC analysis for RNA (expression levels)
# Assuming that 'rna_expression' is the column with RNA expression data
# ROC analysis for PC1
roc_pc1 <- roc(cnv_merged_data_selected$IHC.HER2_binary, cnv_merged_data_selected$PC1, levels = c("0", "1"), direction = "<")
auc_pc1 <- auc(roc_pc1)
print(paste("AUC for PC1: ", auc_pc1))

# ROC analysis for PC2
roc_pc2 <- roc(cnv_merged_data_selected$IHC.HER2_binary, cnv_merged_data_selected$PC2, levels = c("0", "1"), direction = "<")
auc_pc2 <- auc(roc_pc2)
print(paste("AUC for PC2: ", auc_pc2))

# Plot ROC curves for both PC1 and PC2
plot(roc_pc1, col = "red", main="ROC Curves for PC1 and PC2")
plot(roc_pc2, col = "blue", add = TRUE)
legend("bottomright", legend=c("PC1", "PC2"), col=c("red", "blue"), lwd=2)


# Compare the AUC values
print(paste("AUC for DNA: ", auc_dna))
# Print AUC values for PC1 and PC2 together
print(paste("AUC for PC1: ", auc_pc1, "; AUC for PC2: ", auc_pc2))


# Plot ROC curves for both DNA and RNA
plot(roc_dna, col = "red", main="ROC Curves for DNA and RNA")
plot(roc_pc1, col = "blue", add = TRUE)
plot(roc_pc2, col = "green", add = TRUE)
legend("bottomright", legend=c("DNA", "PC1", "PC2"), col=c("red", "blue", "green"), lwd=5)

# Statistical comparison of the AUCs
# Comparing DNA ROC with PC1 ROC
test_dna_pc1 <- roc.test(roc_dna, roc_pc1)
print(test_dna_pc1)

# Comparing DNA ROC with PC2 ROC
test_dna_pc2 <- roc.test(roc_dna, roc_pc2)
print(test_dna_pc2)

