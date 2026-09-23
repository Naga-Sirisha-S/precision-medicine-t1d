# Load
library(DESeq2)
library(ggplot2)
library(pheatmap)
library(dplyr)
library(ggrepel)
library(plotly)
library(org.Hs.eg.db)
library(AnnotationDbi)

### read count
counts <- read.csv("D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/expr_matrix_clean.csv", row.names = 1)

# Convert to integer
counts <- round(as.matrix(counts))

dim(counts)    ### 16067  66
head(counts)

#### Metadata Creation

samples <- colnames(counts)

# Detect conditions
condition <- ifelse(grepl("^healthy", samples, ignore.case = TRUE), "Healthy", "T1D")
condition <- factor(condition, levels = c("Healthy","T1D"))

colData <- data.frame(Condition = condition)
rownames(colData) <- samples

table(colData$Condition)

######### Create DESeq2 Dataset

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = colData,
                              design = ~ Condition)


#### Low-Count Filtering
keep <- rowSums(counts(dds) >= 10) >= 3    ## Keep genes that are expressed (≥10 counts) in at least 3 samples.
dds <- dds[keep, ]
dim(dds)

##########  Run DESeq2

dds <- DESeq(dds)
res <- results(dds)
res_df <- as.data.frame(res)
res_df$gene_id <- rownames(res_df)

############# Add Gene Symbols

res_df$symbol <- mapIds(org.Hs.eg.db,
                        keys = res_df$gene_id,
                        column = "SYMBOL",
                        keytype = "ENSEMBL",
                        multiVals = "first")

# Replace missing symbols with Ensembl ID
res_df$symbol[is.na(res_df$symbol)] <- res_df$gene_id

#################################################################################################################
################################################################################################################
################################ save both significant and unsignificant data
dim(res_df)   #### 14182     8
head(res_df)

sum(grepl("^ENSG", res_df$symbol))   #### 431



###########3 Add Significant / Not Significant in res_df #############

library(dplyr)

# Add significance label column
res_df <- res_df %>%
  mutate(Significance = ifelse(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1,
                               "Significant",
                               "Not_Significant"))

head(res_df)

table(res_df$Significance)  ### not_significant (14104)  ,  significant (78)

###### remove the rows with contains symbols start with ENSG

res_df <- res_df[!grepl("^ENSG00", res_df$symbol), ]
head(res_df)
dim(res_df)                       #### 13751  8
table(res_df$Significance)         ##### not significant(13692), significant(59)

##################### divide up and down regulated in the significant column


library(dplyr)

res_df <- res_df %>%
  mutate(Significance = case_when(
    !is.na(padj) & padj < 0.05 & log2FoldChange > 1  ~ "Upregulated",
    !is.na(padj) & padj < 0.05 & log2FoldChange < -1 ~ "Downregulated",
    TRUE ~ "Not_Significant"
  ))

head(res_df)

table(res_df$Significance)    #### downregulated (19) , upregulated (40), not_significant (13692)

write.csv(res_df, "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/TOTAL_DFA.csv", row.names = FALSE)

#### volcano plot 

library(ggrepel)

library(dplyr)

# Create -log10(padj)
res_df <- res_df %>%
  mutate(neg_log10_padj = -log10(padj))


ggplot(res_df, aes(x = log2FoldChange, y = neg_log10_padj, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not_Significant" = "yellow")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "Volcano Plot",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")


#### with labels volcano plot

# Select top genes
top_up   <- res_df %>% filter(Significance == "Upregulated") %>% arrange(padj) %>% head(10)
top_down <- res_df %>% filter(Significance == "Downregulated") %>% arrange(padj) %>% head(10)

top_genes <- bind_rows(top_up, top_down)

# Plot with labels
ggplot(res_df, aes(log2FoldChange, neg_log10_padj, color = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_genes, aes(label = symbol),
                  size = 3, max.overlaps = 50) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not_Significant" = "yellow")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme_classic() +
  labs(title = "Volcano Plot with Gene Labels",
       x = "Log2 Fold Change",
       y = "-Log10 Adjusted P-value")

##################
Warning message:
  Removed 236 rows containing missing values or values outside the scale range (`geom_point()`).



########################### biomarker find out
### TOP 10 BIOMARKERS

biomarkers <- res_df %>%
  filter(Significance != "Not_Significant") %>%
  arrange(padj) %>%
  head(10) %>%
  pull(symbol)

biomarkers     ### "AGAP9"   "MRPL12"  "MIF"     "MZT2B"   "NPIPB3"  "NICOL1"  "CCDC85B" "NDUFB7"  "MTFP1"   "SEMA6B" 


###### Top Upregulated + Downregulated genes

top_up   <- res_df %>% filter(Significance == "Upregulated") %>% arrange(padj) %>% head(5)
top_down <- res_df %>% filter(Significance == "Downregulated") %>% arrange(padj) %>% head(5)
biomarkers <- c(top_up$symbol, top_down$symbol)

biomarkers    #### "MRPL12"  "MIF"     "MZT2B"   "NICOL1"  "CCDC85B" "AGAP9"   "NPIPB3"  "NHERF4"  "AMY1B"   "PRR16" 

#################################################################################################################################################
###### Strong Fold Change Biomarkers

biomarkers <- res_df %>%
  filter(Significance != "Not_Significant") %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(10) %>%
  pull(symbol)

###### Manually choose known disease genes (Literature-baseD)
#### NA ###########
##########################################################################################################################
##############################################################################################################

########## Significant DE Genes
res_sig <- res_df %>%
  filter(!is.na(padj) & padj < 0.05 & abs(log2FoldChange) > 1)

write.csv(res_sig, "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/differential_exp_genes_final2.csv", row.names = FALSE)


####################  Up & Down Genes

upregulated   <- res_sig %>% filter(log2FoldChange > 1)
downregulated <- res_sig %>% filter(log2FoldChange < -1)

write.csv(upregulated,   "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/Upregulated_genes_final2.csv", row.names = FALSE)
write.csv(downregulated, "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/Downregulated_genes_final2.csv", row.names = FALSE)

################ Volcano Plot

# Replace NA padj
res_df$padj[is.na(res_df$padj)] <- 1

# Classify genes
res_df$Significant <- "Not Sig"
res_df$Significant[res_df$padj < 0.05 & res_df$log2FoldChange > 1]  <- "Up"
res_df$Significant[res_df$padj < 0.05 & res_df$log2FoldChange < -1] <- "Down"

ggplot(res_df, aes(log2FoldChange, -log10(padj), color=Significant)) +
  geom_point(alpha=0.6) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_minimal() +
  labs(title="Volcano Plot: T1D vs Healthy")

##################### Heatmap Top 50 Genes

top_genes <- res_df %>%
  filter(!is.na(padj)) %>%
  arrange(padj) %>%
  head(50) %>%
  pull(gene_id)

norm_counts <- counts(dds, normalized=TRUE)
top_genes <- intersect(top_genes, rownames(norm_counts))

pheatmap(norm_counts[top_genes, ],
         scale = "row",
         show_rownames = FALSE)


################## High Expression Biomarkers

biomarkers <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 1 & baseMean > 50)

write.csv(biomarkers, "D:/BIOINFORMATICS/OMICS LOGIC/PROJECTS/combined_data/PCA2/High_expression_genes_final2.csv", row.names = FALSE)

################### Volcano with Labels + Arrows

top_biomarkers <- res_df %>%
  filter(padj < 0.05 & abs(log2FoldChange) > 2) %>%
  arrange(desc(abs(log2FoldChange))) %>%
  head(15)

ggplot(res_df, aes(log2FoldChange, -log10(padj))) +
  geom_point(aes(color = Significant), alpha=0.6) +
  geom_point(data = top_biomarkers, color="black", size=3) +
  geom_text_repel(data = top_biomarkers,
                  aes(label = symbol),
                  size = 4,
                  fontface = "bold",
                  arrow = arrow(length = unit(0.02, "npc")),
                  box.padding = 0.5) +
  geom_vline(xintercept=c(-1,1), linetype="dashed") +
  geom_hline(yintercept=-log10(0.05), linetype="dashed") +
  theme_classic() +
  labs(title="Volcano Plot: T1D vs Healthy",
       x="Log2 Fold Change",
       y="-Log10 Adjusted P-value")

################ Interactive Volcano (Plotly)

p <- ggplot(res_df, aes(log2FoldChange, -log10(padj),
                        text = paste("Gene:", symbol,
                                     "<br>log2FC:", log2FoldChange,
                                     "<br>padj:", padj))) +
  geom_point(aes(color=Significant), alpha=0.6) +
  theme_minimal()

ggplotly(p, tooltip="text")

#############################################

colSums(is.na(res_df))
sum(is.na(res_df))
dim(res_df)  ## 1885 removed
dim(counts)  

#############################################################################################################

######## Differential expression analysis was performed using DESeq2 (vX.X), with p-values adjusted for multiple testing using the Benjamini–Hochberg method 
## to control the false discovery rate (FDR). Genes were considered significantly differentially expressed at FDR < 0.05 and |log₂ fold change| > 1.


