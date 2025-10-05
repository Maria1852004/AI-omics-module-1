
# ==========================================
# 1️⃣ Install and load required packages
# ==========================================
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c("GEOquery", "Biobase", "limma", "arrayQualityMetrics", "EnhancedVolcano", "pheatmap"))
install.packages("ggplot2")

library(GEOquery)
library(Biobase)
library(arrayQualityMetrics)
library(limma)
library(ggplot2)
library(EnhancedVolcano)
library(pheatmap)

# ==========================================
# 2️⃣ Download GEO dataset (GSE13911)
# ==========================================
gse_data <- getGEO("GSE13911", GSEMatrix = TRUE)
show(gse_data)

# Extract expression values and phenotype data
exprs_data <- exprs(gse_data[[1]])
pheno_data <- pData(gse_data[[1]])

# ==========================================
# 3️⃣ Create ExpressionSet object
# ==========================================
eset <- ExpressionSet(assayData = exprs_data,
                      phenoData = AnnotatedDataFrame(pheno_data))

# ==========================================
# 4️⃣ Quality control before normalization
# ==========================================
arrayQualityMetrics(expressionset = eset, 
                    outdir = "QC_before_norm", 
                    force = TRUE)
# Check the report for flagged arrays

# ==========================================
# 5️⃣ Normalize data (quantile normalization)
# ==========================================
exprs_data_norm <- normalizeBetweenArrays(exprs_data, method = "quantile")

eset_norm <- ExpressionSet(assayData = exprs_data_norm,
                           phenoData = AnnotatedDataFrame(pheno_data))

arrayQualityMetrics(expressionset = eset_norm, 
                    outdir = "QC_after_norm", 
                    force = TRUE)

# ==========================================
# 6️⃣ Filter low-expression probes
# ==========================================
threshold <- 5
keep <- rowSums(exprs_data_norm > threshold) >= (ncol(exprs_data_norm)/2)
exprs_data_filtered <- exprs_data_norm[keep, ]
cat("Number of transcripts after filtering:", nrow(exprs_data_filtered), "\n")

# ==========================================
# 7️⃣ Relabel phenotype groups
# ==========================================
pheno_data$group <- ifelse(pheno_data$disease_status == "healthy", "normal", "cancer")
table(pheno_data$group)

# ==========================================
# 8️⃣ Differential Expression Analysis using limma
# ==========================================
design <- model.matrix(~ group, data = pheno_data)
fit <- lmFit(exprs_data_filtered, design)
fit <- eBayes(fit)

results <- topTable(fit, adjust = "fdr", number = Inf)
write.csv(results, "DEG_results.csv")
head(results)

# ==========================================
# 9️⃣ Boxplot after normalization
# =====================
colnames(pheno_data)
pheno_data$characteristics_ch1.1
# Extract disease information
pheno_data$disease_status <- ifelse(grepl("cancer", pheno_data$characteristics_ch1.1, ignore.case = TRUE), 
                                    "cancer", 
                                    "normal")

# Add grouping column
pheno_data$group <- pheno_data$disease_status
table(pheno_data$group)
unique(pheno_data$characteristics_ch1)
unique(pheno_data$source_name_ch1)
# Create disease status and group columns correctly
pheno_data$disease_status <- ifelse(grepl("tumor", pheno_data$source_name_ch1, ignore.case = TRUE),
                                    "cancer",
                                    "normal")

pheno_data$group <- pheno_data$disease_status

# Check number of samples in each group
table(pheno_data$group)

# ===============================
# Differential expression + plots
# (Assumes exprs_data_filtered and pheno_data already exist)
# ===============================

# Load required packages
library(limma)
library(EnhancedVolcano)
library(pheatmap)
library(ggplot2)

# -------------------------------
# 1) Define group factor from metadata
# -------------------------------
# Use the column that contains "tumor biopsy" / "adjacent normal tissue"
group <- factor(pheno_data$source_name_ch1)
table(group)   # check counts (should show both groups)

# -------------------------------
# 2) Design matrix (use syntactic names)
# -------------------------------
design <- model.matrix(~0 + group)
levels_names <- levels(group)
colnames(design) <- make.names(levels_names)  # make names safe for contrasts
design

# -------------------------------
# 3) Fit linear model with limma
# -------------------------------
fit <- lmFit(exprs_data_filtered, design)

# -------------------------------
# 4) Create contrast: Tumor vs Normal
# -------------------------------
# Build a contrast string using safe names
contrast_string <- paste0(make.names("tumor biopsy"), "-", make.names("adjacent normal tissue"))
contrast_matrix <- makeContrasts(contrasts = contrast_string, levels = colnames(design))
contrast_matrix

# Apply contrasts and empirical Bayes moderation
fit2 <- contrasts.fit(fit, contrast_matrix)
fit2 <- eBayes(fit2)

# -------------------------------
# 5) Extract differential expression results
# -------------------------------
deg_results <- topTable(fit2, number = Inf, adjust.method = "fdr", sort.by = "P")
head(deg_results)
# Save full results
write.csv(deg_results, file = "DEG_results_tumor_vs_normal.csv", row.names = TRUE)

# Quick counts of significant genes (FDR < 0.05, |logFC| >= 1)
sig_idx <- deg_results$adj.P.Val < 0.05 & abs(deg_results$logFC) >= 1
sig_count <- sum(sig_idx, na.rm = TRUE)
sig_up <- sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC > 0, na.rm = TRUE)
sig_down <- sum(deg_results$adj.P.Val < 0.05 & deg_results$logFC < 0, na.rm = TRUE)
cat("Significant (FDR<0.05 & |logFC|>=1):", sig_count, "  Up:", sig_up, "  Down:", sig_down, "\n")

# -------------------------------
# 6) Volcano plot (EnhancedVolcano)
# -------------------------------
# Use P.Value for raw p-value and adj.P.Val for FDR; adjust thresholds as needed
volcano_plot <- EnhancedVolcano(deg_results,
                                lab = rownames(deg_results),
                                x = 'logFC',
                                y = 'P.Value',
                                pCutoff = 0.05,
                                FCcutoff = 1,
                                title = 'Volcano Plot: Tumor vs Adjacent Normal',
                                subtitle = 'GSE13911',
                                caption = paste0('n=', ncol(exprs_data_filtered)),
                                pointSize = 2.0,
                                labSize = 3.0)
print(volcano_plot)
ggsave("Volcano_Tumor_vs_Normal.pdf", plot = volcano_plot, width = 8, height = 6)

# -------------------------------
# 7) Heatmap of top DEGs
# -------------------------------
# Choose top N genes (by P-value or adjust to your preference)
n_top <- min(50, nrow(deg_results))
top_probes <- rownames(deg_results)[1:n_top]
exprs_top <- exprs_data_filtered[top_probes, , drop = FALSE]

# Scale rows (z-score) for visualization
exprs_top_scaled <- t(scale(t(exprs_top)))

# Column annotation (samples)
annotation_col <- data.frame(Group = group)
rownames(annotation_col) <- colnames(exprs_top)

# Save heatmap to PDF
pdf("Heatmap_Top_DEGs_Tumor_vs_Normal.pdf", width = 8, height = 10)
pheatmap(exprs_top_scaled,
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = FALSE,
         main = paste0("Top ", n_top, " DEGs (Tumor vs Normal)"))
dev.off()

# -------------------------------
# 8) Optional: simple base volcano (if EnhancedVolcano isn't available)
# -------------------------------
with(deg_results, {
  plot(logFC, -log10(P.Value), pch = 20,
       main = "Volcano Plot (base R)",
       xlab = "log2 Fold Change", ylab = "-log10(P-value)")
  points(logFC[adj.P.Val < 0.05 & abs(logFC) >= 1],
         -log10(P.Value[adj.P.Val < 0.05 & abs(logFC) >= 1]),
         pch = 20, col = "red")
})
# Save base volcano
dev.copy(pdf, "Volcano_base_Tumor_vs_Normal.pdf", width = 8, height = 6)
dev.off()

# ===============================
# End of script — outputs:
# - DEG_results_tumor_vs_normal.csv
# - Volcano_Tumor_vs_Normal.pdf
# - Heatmap_Top_DEGs_Tumor_vs_Normal.pdf
# - Volcano_base_Tumor_vs_Normal.pdf
# ===============================
head(deg_results)
sum(deg_results$adj.P.Val < 0.05)
library(arrayQualityMetrics)
arrayQualityMetrics(expressionset = exprs_data_filtered, 
                    outdir = "QC_Report", 
                    force = TRUE)
# Load the required package
library(Biobase)
library(arrayQualityMetrics)

# Convert matrix to ExpressionSet
exprs_set <- ExpressionSet(assayData = as.matrix(exprs_data_filtered))

# Run the quality control
arrayQualityMetrics(expressionset = exprs_set,
                    outdir = "QC_Report",
                    force = TRUE)
# ==========================================
# 🔹 PCA Plot After Normalization
# ==========================================

# Load required package
library(ggplot2)

# Transpose normalized expression data (samples as rows)
pca_data <- prcomp(t(exprs_data_filtered), scale. = TRUE)

# Extract principal components
pca_df <- data.frame(pca_data$x)

# Add phenotype (group) info
pca_df$group <- pheno_data$group

# Save PCA plot to PDF
pdf("PCA_After_Normalization.pdf")
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot After Normalization",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )
dev.off()

# Display PCA plot on screen
ggplot(pca_df, aes(x = PC1, y = PC2, color = group)) +
  geom_point(size = 3) +
  labs(title = "PCA Plot After Normalization (Displayed)",
       x = "Principal Component 1",
       y = "Principal Component 2") +
  theme_minimal() +
  theme(
    text = element_text(size = 14),
    plot.title = element_text(hjust = 0.5)
  )
# ==========================================
# 🔹 Boxplot After Normalization
# ==========================================

# Save boxplot as PDF
pdf("Boxplot_After_Normalization.pdf")

boxplot(exprs_data_filtered,
        main = "Boxplot After Normalization",
        las = 2,           # make sample names vertical
        outline = FALSE,   # remove outliers
        col = "lightblue") # color for boxes

dev.off()

# Display boxplot on screen
boxplot(exprs_data_filtered,
        main = "Boxplot After Normalization",
        las = 2,
        outline = FALSE,
        col = "lightblue")