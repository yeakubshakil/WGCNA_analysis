



#==========Install Packages
install.packages("BiocManager")
BiocManager::install("WGCNA")
install.packages("fastcluster")
BiocManager::install("biomaRt")
install.packages("corrplot")
install.packages("fastcluster")
install.packages("doParallel")
install.packages("flashClust")


#============Load packages
library(tidyverse)
library(readr)
library(dplyr)
library(ggplot2)
library(corrplot)
library(doParallel)
library(fastcluster)
library(flashClust)
library(biomaRt)
library(WGCNA)

#========Enable multithreading 
nCores <- max(2, parallel::detectCores() - 1)
registerDoParallel(cores = nCores)
enableWGCNAThreads(nThreads = nCores)
cat("Using", getDoParWorkers(), "cores for parallel processing.\n")

#===== Step 1: Load and preprocess data
# Load counts (genes as rows, samples as columns) and metadata
counts <- read_csv("/content/GSE159984_counts_data.csv", show_col_types = FALSE)

# Handle duplicate gene symbols
if (any(duplicated(counts$gene_symbol))) {
  cat("Duplicate gene symbols found. Aggregating by mean expression.\n")
  counts <- counts %>%
    group_by(gene_symbol) %>%
    summarise(across(where(is.numeric), \(x) mean(x, na.rm = TRUE))) %>%
    ungroup()
} else {
  cat("No duplicate gene symbols found.\n")
}

# ======Set gene symbols as row names
counts <- counts %>% column_to_rownames("gene_symbol")

# ======Load metadata
metadata <- read_csv("/content/GSE159984_filtered_metadata.csv", show_col_types = FALSE) %>%
  mutate(Sample = as.character(Sample))


# =======Ensure sample ID consistency
common_samples <- intersect(colnames(counts), metadata$Sample)
if (length(common_samples) == 0) {
  stop("No matching sample IDs between counts and metadata!")
}
counts <- counts %>% dplyr::select(all_of(common_samples))
metadata <- metadata %>% filter(Sample %in% common_samples)
cat("Number of matching samples:", length(common_samples), "\n")


# ======Transpose counts for WGCNA (samples as rows, genes as columns)
datExpr <- t(counts)
cat("Expression data dimensions:", dim(datExpr), "\n")

#========== Step 2: Clean Data
# Remove genes/samples with excessive missing values or low variance
gsg <- WGCNA::goodSamplesGenes(datExpr, verbose = 3)
if (!gsg$allOK) {
  cat("Removing samples/genes with excessive missing values or low variance.\n")
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  metadata <- metadata %>% filter(Sample %in% rownames(datExpr))
}

#============ Step 3: Detect outliers
sampleTree <- flashClust(dist(datExpr), method = "average")
pdf("sample_clustering.pdf", width = 10, height = 6)
par(mar = c(0, 4, 2, 0))
plot(sampleTree, main = "Sample Clustering to Detect Outliers", sub = "", xlab = "",
     cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
dev.off()


# ================Step 4: Pick soft threshold
powers <- c(1:10, seq(12, 30, by = 2))
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5,
                         networkType = "unsigned", corFnc = "bicor")
if (is.na(sft$powerEstimate)) {
  warning("No suitable soft threshold found. Selecting power with highest R^2.")
  optimal_power <- sft$fitIndices$Power[which.max(sft$fitIndices$SFT.R.sq)]
} else {
  optimal_power <- sft$fitIndices$Power[min(which(sft$fitIndices$SFT.R.sq >= 0.80))]
  if (length(optimal_power) == 0) {
    optimal_power <- sft$powerEstimate
  }
}
cat("Selected soft power:", optimal_power, "\n")

# Plot soft threshold diagnostics
pdf("soft_threshold_plots.pdf", width = 10, height = 5)
par(mfrow = c(1, 2))
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free R^2", type = "n",
     main = "Scale Independence")
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3]) * sft$fitIndices[,2],
     labels = powers, cex = 0.9, col = "red")
abline(h = 0.80, col = "blue")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n",
     main = "Mean Connectivity")
text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, cex = 0.9, col = "red")
dev.off()


# ===========Step 5: Construct co-expression network
adjacency <- adjacency(datExpr, power = optimal_power, type = "unsigned", corFnc = "bicor")
TOM <- TOMsimilarity(adjacency, TOMType = "unsigned")
dissTOM <- 1 - TOM


# ==========Step 6: Hierarchical clustering and module identification
geneTree <- flashClust(as.dist(dissTOM), method = "average")
minModuleSize <- 30
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                             deepSplit = 2, pamRespectsDendro = FALSE,
                             minClusterSize = minModuleSize)
dynamicColors <- labels2colors(dynamicMods)
cat("Number of initial modules:", length(unique(dynamicColors)), "\n")


# Plot initial modules
pdf("gene_clustering.pdf", width = 12, height = 8)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Modules",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Initial Modules")
dev.off()


# ==============Step 7: Merge similar modules
MEs <- moduleEigengenes(datExpr, dynamicColors)$eigengenes
METree <- flashClust(as.dist(1 - cor(MEs)), method = "average")
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 3)
mergedColors <- merge$colors
mergedMEs <- merge$newMEs
cat("Number of merged modules:", length(unique(mergedColors)), "\n")


# Plot merged modules
pdf("merged_modules.pdf", width = 12, height = 8)
plotDendroAndColors(geneTree, mergedColors, "Merged Modules",
                    dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05,
                    main = "Gene Dendrogram and Merged Modules")
dev.off()

# ================Step 8: Correlate modules with traits
if ("condition" %in% colnames(metadata)) {
  traitData <- metadata %>%
    select(Sample, condition) %>%
    mutate(condition = as.numeric(factor(condition)) - 1) %>%
    column_to_rownames("Sample")
  traitVector <- traitData$condition[match(rownames(datExpr), rownames(traitData))]
  
  moduleTraitCor <- cor(mergedMEs, traitVector, use = "p", method = "pearson")
  moduleTraitP <- corPvalueStudent(moduleTraitCor, nrow(datExpr))
  
  
  # Plot module-trait correlation heatmap
  pdf("module_trait_heatmap.pdf", width = 8, height = 8)
  textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitP, 1), ")", sep = "")
  dim(textMatrix) <- dim(moduleTraitCor)
  labeledHeatmap(Matrix = moduleTraitCor, xLabels = "Condition", yLabels = colnames(mergedMEs),
                 colorLabels = FALSE, colors = blueWhiteRed(50), textMatrix = textMatrix,
                 setStdMargins = FALSE, cex.text = 0.7, main = "Module-Trait Correlations")
  dev.off()
} else {
  warning("No 'condition' column in metadata. Skipping trait correlation.")
  moduleTraitCor <- NULL
  moduleTraitP <- NULL
}



# ================Step 9: Identify hub genes in hub gene module
# Find hub gene module
hub_module <- mergedColors[colnames(datExpr) == "hub gene"]
if (length(hub gene_module) == 0 || is.na(hub gene_module)) {
  warning("hub gene not found or not assigned to a module. Selecting module with highest trait correlation.")
  if (!is.null(moduleTraitCor)) {
    ccna2_module <- names(moduleTraitCor)[which.max(abs(moduleTraitCor))]
  } else {
    ccna2_module <- names(table(mergedColors))[1] # Fallback to largest module
  }
}
cat("hub gene module:", ccna2_module, "\n")


# Calculate module membership (MM) and gene significance (GS)
geneModuleMembership <- cor(datExpr, mergedMEs, use = "p", method = "bicor")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))
geneTraitSignificance <- cor(datExpr, traitVector, use = "p", method = "bicor")
GSPvalue <- corPvalueStudent(geneTraitSignificance, nrow(datExpr))


# Select genes in hub gene module
module_genes <- names(mergedColors)[mergedColors == hub gene_module]


# Identify hub genes (MM >= 0.85, GS >= 0.85)
MM_threshold <- 0.85
GS_threshold <- 0.85
important_MM_genes <- module_genes[abs(geneModuleMembership[module_genes, paste0("ME", hub gene_module)]) >= MM_threshold]
important_GS_genes <- module_genes[abs(geneTraitSignificance[module_genes, 1]) >= GS_threshold]
hub_genes <- intersect(important_MM_genes, important_GS_genes)


# Annotate hub genes using biomaRt
if (length(hub_genes) > 0) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  geneOverall: The warning in the `summarise()` step has been resolved by updating the `across()` syntax to use an anonymous function, ensuring compatibility with `dplyr` 1.1.0+. The pipeline remains robust, efficient, and focused on identifying hub genes in the CCNA2 module, with all previous enhancements preserved.
  
  
# =========Calculate module membership (MM) and gene significance (GS)
geneModuleMembership <- cor(datExpr, mergedMEs, use = "p", method = "bicor")
MMPvalue <- corPvalueStudent(geneModuleMembership, nrow(datExpr))
geneTraitSignificance <- cor(datExpr, traitVector, use = "p", method = "bicor")
GSPvalue <- corPvalueStudent(geneTraitSignificance, nrow(datExpr))
  

# Select genes in hub gene module
module_genes <- names(mergedColors)[mergedColors == hub gene_module]


# Identify hub genes (MM >= 0.85, GS >= 0.85)
MM_threshold <- 0.85
GS_threshold <- 0.85
important_MM_genes <- module_genes[abs(geneModuleMembership[module_genes, paste0("ME", ccna2_module)]) >= MM_threshold]
important_GS_genes <- module_genes[abs(geneTraitSignificance[module_genes, 1]) >= GS_threshold]
hub_genes <- intersect(important_MM_genes, important_GS_genes)


# Annotate hub genes using biomaRt
if (length(hub_genes) > 0) {
  ensembl <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  gene_annotations <- getBM(attributes = c("hgnc_symbol", "description", "entrezgene_id"),
                            filters = "hgnc_symbol", values = hub_genes, mart = ensembl)
} else {
  gene_annotations <- data.frame(hgnc_symbol = character(), description = character(), entrezgene_id = numeric())
}



# Save hub genes with annotations
hub_genes_df <- data.frame(
  Gene = hub_genes,
  Module = hub gene_module,
  MM = geneModuleMembership[hub_genes, paste0("ME", ccna2_module)],
  GS = geneTraitSignificance[hub_genes, 1],
  MMPvalue = MMPvalue[hub_genes, paste0("ME", ccna2_module)],
  GSPvalue = GSPvalue[hub_genes, 1]
) %>% left_join(gene_annotations, by = c("Gene" = "hgnc_symbol"))
write_csv(hub_genes_df, "hub gene_hub_genes.csv")


# ===============Step 10: Visualize MM vs GS for hub gene module
if (length(module_genes) > 0) {
  pdf("MM_vs_GS_plot.pdf", width = 8, height = 6)
  plot(geneModuleMembership[module_genes, paste0("ME", ccna2_module)],
       geneTraitSignificance[module_genes, 1],
       xlab = "Module Membership", ylab = "Gene Significance",
       main = paste("MM vs GS in", ccna2_module, "Module"),
       pch = 19, col = "blue", cex = 0.5)
  abline(h = GS_threshold, v = MM_threshold, col = "red", lty = 2)
  text(geneModuleMembership[hub_genes, paste0("ME", ccna2_module)],
       geneTraitSignificance[hub_genes, 1], labels = hub_genes, pos = 3, cex = 0.7, col = "black")
  dev.off()
}


# ============Step 11: Save comprehensive results
# Module assignments
module_assignments <- data.frame(Gene = colnames(datExpr), Module = mergedColors)
write_csv(module_assignments, "wgcna_module_assignments.csv")

# Module eigengenes
write_csv(as.data.frame(mergedMEs), "wgcna_module_eigengenes.csv")


# Results summary
results_summary <- list(
  optimal_power = optimal_power,
  n_initial_modules = length(unique(dynamicColors)),
  n_merged_modules = length(unique(mergedColors)),
  hub gene_module = hub gene_module,
  n_hub_genes = length(hub_genes)
)
write_csv(as.data.frame(results_summary), "wgcna_results_summary.csv")

# Save workspace for reproducibility
save.image("wgcna_analysis.RData")

cat("WGCNA pipeline completed successfully. Results saved in CSV and PDF files.\n")