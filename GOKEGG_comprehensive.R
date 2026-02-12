# Load required libraries
library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(ReactomePA)
library(enrichplot)
library(openxlsx)
library(dplyr)
library(stringr)
# =============================================================================
rm(list = ls())  # Clear the environment
gc()  # Run garbage collection




# Read the gene list
print("Reading gene list from Excel file...")
gene_data <- read_excel("48_intersection.xlsx")
gene_list <- gene_data[[1]]  # First column
gene_list <- gene_list[!is.na(gene_list)]  # Remove NA values
print(paste("Loaded", length(gene_list), "genes"))

# Display first few genes to verify
print("First 10 genes:")
print(head(gene_list, 10))

# =============================================================================
# Gene ID Conversion
# =============================================================================

# Convert gene symbols to Entrez IDs (required for most analyses)
print("Converting gene symbols to Entrez IDs...")
gene_entrez <- bitr(gene_list, 
                    fromType = "SYMBOL",
                    toType = c("ENTREZID", "ENSEMBL", "UNIPROT"),
                    OrgDb = org.Hs.eg.db)

print(paste("Successfully converted", nrow(gene_entrez), "genes"))
print(paste("Conversion rate:", round(nrow(gene_entrez)/length(gene_list)*100, 1), "%"))

# Extract Entrez IDs for analysis
entrez_ids <- gene_entrez$ENTREZID

# =============================================================================
# LAYER 2: Standard Enrichment Analysis
# =============================================================================

print("Starting Layer 2 analysis...")

# 1. GO Biological Process
print("Running GO Biological Process analysis...")
go_bp <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "BP",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# 2. GO Molecular Function
print("Running GO Molecular Function analysis...")
go_mf <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "MF",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# 3. GO Cellular Component
print("Running GO Cellular Component analysis...")
go_cc <- enrichGO(gene = entrez_ids,
                  OrgDb = org.Hs.eg.db,
                  ont = "CC",
                  pAdjustMethod = "BH",
                  pvalueCutoff = 0.05,
                  qvalueCutoff = 0.2,
                  readable = TRUE)

# 4. KEGG Pathway Analysis
print("Running KEGG pathway analysis...")
kegg_result <- enrichKEGG(gene = entrez_ids,
                          organism = 'hsa',
                          pvalueCutoff = 0.05,
                          pAdjustMethod = "BH",
                          qvalueCutoff = 0.2)

# 5. Reactome Pathway Analysis
print("Running Reactome pathway analysis...")
reactome_result <- enrichPathway(gene = entrez_ids,
                                 pvalueCutoff = 0.05,
                                 pAdjustMethod = "BH",
                                 qvalueCutoff = 0.2,
                                 readable = TRUE)




# =============================================================================
# LAYER 3: Advanced Analysis
# =============================================================================

print("Starting Layer 3 analysis...")

# 1. Gene Set Enrichment Analysis (GSEA) preparation
# First, we need a ranked gene list (usually fold changes)
# Since we don't have expression data, we'll use degree centrality as proxy
print("Preparing for advanced analysis...")

# 2. GO Semantic Similarity Analysis
library(GOSemSim)
hsGO <- godata('org.Hs.eg.db', ont="BP", computeIC=TRUE)

# Calculate semantic similarity for top GO terms
if(nrow(go_bp@result) > 0) {
  top_go_terms <- go_bp@result$ID[1:min(20, nrow(go_bp@result))]
  go_similarity <- mgoSim(top_go_terms, top_go_terms, 
                          semData=hsGO, measure="Wang", combine=NULL)
}

# 3. Pathway-Pathway Interaction Analysis
# Compare KEGG and Reactome overlaps
if(nrow(kegg_result@result) > 0 & nrow(reactome_result@result) > 0) {
  kegg_genes <- strsplit(kegg_result@result$geneID, "/")
  reactome_genes <- strsplit(reactome_result@result$geneID, "/")
  
  # Calculate pathway overlap matrix
  pathway_overlap <- data.frame()
  for(i in 1:min(10, length(kegg_genes))) {
    for(j in 1:min(10, length(reactome_genes))) {
      overlap <- length(intersect(kegg_genes[[i]], reactome_genes[[j]]))
      if(overlap > 0) {
        pathway_overlap <- rbind(pathway_overlap, data.frame(
          KEGG_pathway = kegg_result@result$Description[i],
          Reactome_pathway = reactome_result@result$Description[j],
          overlap_genes = overlap,
          jaccard_index = overlap / length(union(kegg_genes[[i]], reactome_genes[[j]]))
        ))
      }
    }
  }
}

# 4. Multi-level GO Analysis
# Analyze GO terms at different hierarchy levels
if(nrow(go_bp@result) > 0) {
  library(ontologyIndex)
  # This would require additional GO hierarchy analysis
  go_levels <- data.frame(
    GO_ID = go_bp@result$ID,
    Description = go_bp@result$Description,
    Level = "Multiple",  # Would calculate actual levels
    p.adjust = go_bp@result$p.adjust
  )
}

# =============================================================================
# Results Compilation and Export
# =============================================================================

print("Compiling results for export...")

# Create a comprehensive results list
results_list <- list()

# Add input data
results_list[["Input_Genes"]] <- data.frame(
  Original_Symbol = gene_list,
  stringsAsFactors = FALSE
)

results_list[["Gene_ID_Mapping"]] <- gene_entrez

# Add Layer 2 results
if(nrow(go_bp@result) > 0) {
  results_list[["GO_Biological_Process"]] <- go_bp@result
}

if(nrow(go_mf@result) > 0) {
  results_list[["GO_Molecular_Function"]] <- go_mf@result
}

if(nrow(go_cc@result) > 0) {
  results_list[["GO_Cellular_Component"]] <- go_cc@result
}

if(nrow(kegg_result@result) > 0) {
  results_list[["KEGG_Pathways"]] <- kegg_result@result
}

if(nrow(reactome_result@result) > 0) {
  results_list[["Reactome_Pathways"]] <- reactome_result@result
}

# if(nrow(do_result@result) > 0) {
#   results_list[["Disease_Ontology"]] <- do_result@result
# }

# Add Layer 3 results
if(exists("pathway_overlap") && nrow(pathway_overlap) > 0) {
  results_list[["Pathway_Crosstalk"]] <- pathway_overlap
}

if(exists("go_levels")) {
  results_list[["GO_Hierarchy_Analysis"]] <- go_levels
}

# Create summary statistics
summary_stats <- data.frame(
  Database = c("GO_BP", "GO_MF", "GO_CC", "KEGG", "Reactome"),
  Total_Terms = c(
    ifelse(exists("go_bp"), nrow(go_bp@result), 0),
    ifelse(exists("go_mf"), nrow(go_mf@result), 0),
    ifelse(exists("go_cc"), nrow(go_cc@result), 0),
    ifelse(exists("kegg_result"), nrow(kegg_result@result), 0),
    ifelse(exists("reactome_result"), nrow(reactome_result@result), 0)
  ),
  Significant_Terms = c(
    ifelse(exists("go_bp"), sum(go_bp@result$p.adjust < 0.05), 0),
    ifelse(exists("go_mf"), sum(go_mf@result$p.adjust < 0.05), 0),
    ifelse(exists("go_cc"), sum(go_cc@result$p.adjust < 0.05), 0),
    ifelse(exists("kegg_result"), sum(kegg_result@result$p.adjust < 0.05), 0),
    ifelse(exists("reactome_result"), sum(reactome_result@result$p.adjust < 0.05), 0)
  )
)

results_list[["Analysis_Summary"]] <- summary_stats

# Export to Excel
output_file <- "Comprehensive_GO_KEGG_Analysis_48.xlsx"
write.xlsx(results_list, file = output_file, rowNames = FALSE)

print(paste("Analysis complete! Results exported to:", output_file))
print(paste("Total sheets created:", length(results_list)))

# Print summary
print("\n=== ANALYSIS SUMMARY ===")
for(i in 1:nrow(summary_stats)) {
  cat(paste(summary_stats$Database[i], ":", 
            summary_stats$Significant_Terms[i], 
            "significant terms\n"))
}

# Show top results from each analysis
print("\n=== TOP RESULTS PREVIEW ===")
if(exists("go_bp") && nrow(go_bp@result) > 0) {
  print("Top 3 GO Biological Processes:")
  print(go_bp@result[1:min(3, nrow(go_bp@result)), c("Description", "p.adjust", "Count")])
}

if(exists("kegg_result") && nrow(kegg_result@result) > 0) {
  print("Top 3 KEGG Pathways:")
  print(kegg_result@result[1:min(3, nrow(kegg_result@result)), c("Description", "p.adjust", "Count")])
}

print("\nAnalysis complete!")
