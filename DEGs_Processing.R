
# Read CSV files
DEGs1 <- read.csv("DEGs_Data_1.csv")
DEGs2 <- read.csv("DEGs_Data_2.csv")

# Handle missing padj values
DEGs1$padj[is.na(DEGs1$padj)] <- 1
DEGs2$padj[is.na(DEGs2$padj)] <- 1

# Function to classify genes
classify_gene <- function(logFC, padj) {
  if (padj < 0.05) {
    if (logFC > 1) {
      return("Upregulated")
    } else if (logFC < -1) {
      return("Downregulated")
    }
  }
  return("Not_Significant")
}

# Apply classification to DEGs_Data_1
status1 <- c()
for (i in 1:nrow(DEGs1)) {
  status1[i] <- classify_gene(DEGs1$logFC[i], DEGs1$padj[i])
}
DEGs1$status <- status1

# Apply classification to DEGs_Data_2
status2 <- c()
for (i in 1:nrow(DEGs2)) {
  status2[i] <- classify_gene(DEGs2$logFC[i], DEGs2$padj[i])
}
DEGs2$status <- status2

# Save results to Results folder
write.csv(DEGs1, "Results/DEGs_Data_1_processed.csv", row.names = FALSE)
write.csv(DEGs2, "Results/DEGs_Data_2_processed.csv", row.names = FALSE)

# Show summary tables
table(DEGs1$status)
table(DEGs2$status)
