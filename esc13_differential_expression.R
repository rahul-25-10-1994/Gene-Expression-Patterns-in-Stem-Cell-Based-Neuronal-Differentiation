# Load necessary libraries
library(Seurat)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

# Load the esc cell Seurat object
esc_seurat <- readRDS("C:/Users/panka/OneDrive/Desktop/ESC/escmerged13_seurat_object.rds")

# Load the neuronal progenitor cell count matrix
npc_counts <- read.table("C:/Users/panka/OneDrive/Desktop/adipose/pjr/npc1_raw_counts.DataMatrix.txt", 
                         header = TRUE, 
                         row.names = 1, 
                         sep = "\t")

# Extract the counts from the Seurat object
esc_counts <- GetAssayData(esc_seurat, assay = "RNA", layer = "counts")

npc_gene_names <- sub("^ENSG[0-9]+\\.[0-9]+>(.+)", "\\1", rownames(npc_counts))

# Create a data frame to map Ensembl IDs to gene names
npc_mapping <- data.frame(ensembl_id = rownames(npc_counts), gene_name = npc_gene_names)

# Remove duplicates, keeping the first occurrence
npc_mapping <- npc_mapping[!duplicated(npc_mapping$gene_name), ]

esc_gene_names <- rownames(esc_counts)
common_genes <- intersect(esc_gene_names, npc_mapping$gene_name)

# Filter both datasets to keep only the common genes
filtered_mesen_counts <- esc_counts[esc_gene_names %in% common_genes, ]
filtered_npc_counts <- npc_counts[npc_gene_names %in% common_genes, ]
common_ensembl_ids <- npc_mapping$ensembl_id[npc_mapping$gene_name %in% common_genes]

# Combine the counts
combined_counts <- cbind(filtered_mesen_counts, filtered_npc_counts[rownames(filtered_npc_counts) %in% common_ensembl_ids, ])

# Prepare metadata
esc_metadata <- esc_seurat@meta.data
esc_metadata$celltype <- "esc"
esc_metadata$orig.ident <- NULL
esc_metadata$nCount_RNA <- NULL
esc_metadata$nFeature_RNA <- NULL
esc_metadata$percent.mt <- NULL
esc_metadata$RNA_snn_res.0.5 <- NULL
esc_metadata$seurat_clusters <- NULL

npc_metadata <- data.frame(celltype = rep("neuronalprogenitorcells", ncol(filtered_npc_counts)), 
                           row.names = colnames(filtered_npc_counts), 
                           stringsAsFactors = FALSE)

# Combine metadata
combined_metadata <- rbind(esc_metadata, npc_metadata)
combined_metadata$celltype <- as.factor(combined_metadata$celltype)

# Ensure numeric columns
is_numeric <- sapply(combined_counts, is.numeric)
if (any(!is_numeric)) {
  print(colnames(combined_counts)[!is_numeric])
}

# Convert columns to numeric where possible
for (col in colnames(combined_counts)[!is_numeric]) {
  combined_counts[[col]] <- as.numeric(as.character(combined_counts[[col]]))
}

# Check for NAs
na_columns <- colnames(combined_counts)[sapply(combined_counts, function(x) any(is.na(x)))]
print(na_columns)

# Add pseudo-count to avoid log issues
combined_counts <- combined_counts + 1

# Create DESeq dataset
dds <- DESeqDataSetFromMatrix(countData = combined_counts, colData = combined_metadata, design = ~ celltype)

# Filter low-expression genes
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

# Run DESeq2 analysis
dds <- DESeq(dds)
results <- results(dds)

# Save results
write.csv(as.data.frame(significant_results), file = "C:/Users/panka/OneDrive/Desktop/esc13(ref)RES.csv")

# Filter significant results
significant_results <- results[which(results$padj < 0.05), ]

# MA PLOT (log2 fold-change vs mean expression)
plotMA(results, ylim = c(-2, 2), main = "MA Plot")
plotMA(results)

# VOLCANO PLOT using ggplot2
results_df <- as.data.frame(results)
results_df$Gene <- rownames(results_df)
results_df$Significant <- ifelse(results_df$padj < 0.05 & abs(results_df$log2FoldChange) > 1, "Yes", "No")

ggplot(results_df, aes(x = log2FoldChange, y = -log10(padj), color = Significant)) +
  geom_point(alpha = 0.6) +
  scale_color_manual(values = c("No" = "grey", "Yes" = "red")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")

# VOLCANO PLOT using EnhancedVolcano
EnhancedVolcano(results_df,
                lab = results_df$Gene,
                x = "log2FoldChange",
                y = "padj",
                pCutoff = 0.1,
                FCcutoff = 1,
                col = c("grey", "grey", "blue", "red"))
results_shrunk <- lfcShrink(dds, coef="celltype_neuronalprogenitorcells_vs_esc", type="apeglm")
plotMA(results_shrunk, ylim=c(-1,1))

