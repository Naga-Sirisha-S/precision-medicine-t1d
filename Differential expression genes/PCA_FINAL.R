renv::status()
renv::activate()
renv::restore()

renv::install(c("dplyr", "ggplot2", "patchwork", "pheatmap", "sva", "mclust"))
renv::install("bioc::sva")
renv::install("bioc::DESeq2")
renv::snapshot()
library(dplyr)
library(ggplot2)
library(DESeq2)
library(sva)

rm(list = ls())
############################################################
# STEP 1: Read Expression Matrix
expr <- read.csv("D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/read_count1_updated.csv",
                 row.names = 1, check.names = FALSE)

# Check dimensions
dim(expr)      # Should be 16785 x 82
head(expr)

# Convert to matrix
expr_matrix <- as.matrix(expr)


############################################################
# STEP 2: Remove Low Expression Genes
# Keep genes expressed >1 in at least 3 samples
expr_matrix <- expr_matrix[rowSums(expr_matrix > 1) >= 3, ]
dim(expr_matrix)

####################################################################################################
################### NAMING THE COLUMNS MANUALLY -- REMOVE STAGE
## name the columns

# Check total columns
ncol(expr_matrix)

# Create new names
healthy_names <- paste0("healthy.", 0:42)   # 43 samples
t1d_names     <- paste0("T1D.", 0:38)       # 39 samples

# Combine
new_names <- c(healthy_names, t1d_names)

# Assign to expression matrix
colnames(expr_matrix) <- new_names


colnames(expr_matrix)

####################################################################
####### REMOVE THE PARTCULAR SAMPLE OVERLAPPING MANUALLY ###########
##### removal of 6 h  and 3 t1d 

# List of samples to remove
remove_samples <- c("healthy.8","healthy.17","healthy.26","healthy.32",
                    "healthy.34","healthy.37","healthy.22","healthy.24","healthy.40",
                    "T1D.28","T1D.12","T1D.2","T1D.38","T1D.31","T1D.32", "T1D.19")

# Check if names exist
remove_samples %in% colnames(expr_matrix)

# Remove columns
expr_matrix_clean <- expr_matrix[, !colnames(expr_matrix) %in% remove_samples]

# Check new dimension
dim(expr_matrix_clean)

#### Update Sample Labels (IMPORTANT)
##You must remove labels too, otherwise mismatch error will occur.

# Get cleaned sample names
clean_samples <- colnames(expr_matrix_clean)

# Create condition vector based on names
condition_clean <- ifelse(grepl("^healthy", clean_samples), "healthy", "T1D")

# Convert to factor
condition_clean <- factor(condition_clean)

# Check
table(condition_clean)
length(condition_clean)
ncol(expr_matrix_clean)   # MUST be same

##########################################################################################
## save expr_matrix_clean -- for doing differential expression analysis

write.csv(expr_matrix_clean, 
          file = "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/expr_matrix_clean.csv",
          row.names = TRUE)


dim(expr_matrix_clean)      # Should be 16785 x 82
head(expr_matrix_clean)

##########################################

expr_matrix_log_clean <- log2(expr_matrix_clean + 1)

pca_result_clean <- prcomp(t(expr_matrix_log_clean), scale. = TRUE)
summary(pca_result_clean)

#################### Create PCA DataFrame
pca_df_clean <- data.frame(
  PC1 = pca_result_clean$x[,1],
  PC2 = pca_result_clean$x[,2],
  Sample = clean_samples,
  Condition = condition_clean
)

######################## Plot Clean PCA

library(ggplot2)

ggplot(pca_df_clean, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4, alpha = 0.8) +
  theme_minimal() +
  labs(title = "PCA After Manual Outlier Removal",
       x = paste0("PC1 (", round(summary(pca_result_clean)$importance[2,1]*100,1), "%)"),
       y = paste0("PC2 (", round(summary(pca_result_clean)$importance[2,2]*100,1), "%)"))




### (Optional) Label Samples
ggplot(pca_df_clean, aes(PC1, PC2, color = Condition)) +
  geom_point(size = 4) +
  geom_text(aes(label = Sample), size = 3, vjust = -1) +
  theme_minimal()
