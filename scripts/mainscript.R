# Load required libraries 
library(dplyr)  
library(ggplot2)
library(ggrepel)
library(eulerr)
library(reshape2)
library(amap)
library(enrichplot) 
library(clusterProfiler)
library(org.Hs.eg.db)
library(STRINGdb)

# Source the function script to load parse_de_table function
source("/Users/claire/Downloads/3053994_BIOL5379_functionscript.R")

# Load files into RStudio using read.table()
sample_sheet=read.table("/Users/claire/Downloads/sample_sheet.csv",header=TRUE,row.names=1,sep="\t")
annotations=read.table("/Users/claire/Downloads/Human_Background_GRCh38.p13.csv",header=TRUE,row.names=1,sep="\t")
em=read.table("/Users/claire/Downloads/EM-2.csv",header=TRUE,row.names=1,sep="\t")
de_senes_vs_prolif=read.table("/Users/claire/Downloads/DE_Senes_vs_Prolif.csv",header=TRUE,row.names=1,sep="\t")
de_senesMtD_vs_senes=read.table("/Users/claire/Downloads/DE_Senes_MtD_vs_Senes.csv",header=TRUE,row.names=1,sep="\t")
de_senesMtD_vs_prolif=read.table("/Users/claire/Downloads/DE_Senes_MtD_vs_Prolif.csv",header=TRUE,row.names=1,sep="\t")

# Parse each differential expression table using the function
# The function produces 6 df and 1 vector: master, em_symbols, em_scaled, master_sig, sig_genes, em_symbols_sig, em_scaled_sig
result_senes_vs_prolif = parse_de_table(de_senes_vs_prolif, em, annotations)
result_senesMtD_vs_senes = parse_de_table(de_senesMtD_vs_senes, em, annotations)
result_senesMtD_vs_prolif = parse_de_table(de_senesMtD_vs_prolif, em, annotations)

---------------------------------------------------------------------------------------
# MY PLOTS:
# VOLCANO PLOT
# MA PLOT
# VENN DIAGRAMS
# FOLD VS FOLD PLOTS  
# PCA PLOTS  
# DENSITY PLOTS 
# HEAT MAPS  
# ORA
# GSEA
# STRING

---------------------------------------------------------------------------------------  
# VOLCANO PLOTS
    
# Variables for customization (Volcano plot)
p_threshold = 0.01
fold_threshold = 2
title = "Senes vs Prolif Volcano Plot"
x_label = "Log2 Fold Change"
y_label = "-log10 P-Value"
point_size = 2
label_size = 4

# Generate the volcano plot
volcano_plot_1 = plot_volcano(result_senes_vs_prolif$master, 
                            p_threshold = p_threshold, 
                            fold_threshold = fold_threshold,
                            title = title, 
                            x_label = x_label, 
                            y_label = y_label, 
                            point_size = point_size, 
                            label_size = label_size)

png("/Users/claire/Downloads/volcano_plot_1.png", height = 500, width = 500)
print(volcano_plot_1)  
dev.off()

volcano_plot_2 = plot_volcano(result_senesMtD_vs_prolif$master, 
                            p_threshold = p_threshold, 
                            fold_threshold = fold_threshold,
                            title = title, 
                            x_label = x_label, 
                            y_label = y_label, 
                            point_size = point_size, 
                            label_size = label_size)

png("/Users/claire/Downloads/volcano_plot_2.png", height = 500, width = 500)
print(volcano_plot_2)  
dev.off()

volcano_plot_3 = plot_volcano(result_senesMtD_vs_senes$master, 
                            p_threshold = p_threshold, 
                            fold_threshold = fold_threshold,
                            title = title, 
                            x_label = x_label, 
                            y_label = y_label, 
                            point_size = point_size, 
                            label_size = label_size)

png("/Users/claire/Downloads/volcano_plot_3.png", height = 500, width = 500)
print(volcano_plot_3)  
dev.off()
----------------------------------------------------------------------------------
# MA PLOTS

# Plot customization variables
up_color = "red"
down_color = "blue"
non_sig_color = "gray"
pt_size = 2
lbl_size = 4
fold_thresh = 2
plot_title = "SenesMtD vs Prolif MA Plot"
x_label = "Log10 Mean Expression"
y_label = "Log2 Fold Change"

# MA Plot for Prolif vs Senes
ma_plot_prolif_vs_senes = plot_MA(result_senes_vs_prolif$master, 
                                  mean_column = "mean", 
                                  fold_threshold = fold_thresh, 
                                  upregulated_color = up_color, 
                                  downregulated_color = down_color, 
                                  nonsignificant_color = non_sig_color, 
                                  point_size = pt_size, 
                                  label_size = lbl_size, 
                                  plot_title = "MA Plot: Prolif vs Senes", 
                                  x_label = x_label, 
                                  y_label = y_label)

png("/Users/claire/Downloads/MA_plot_1.png", height = 500, width = 500)
print(ma_plot_prolif_vs_senes)
dev.off()

# MA Plot for SenesMtD vs Senes
ma_plot_senesMtD_vs_senes = plot_MA(result_senesMtD_vs_senes$master, 
                                    mean_column = "mean", 
                                    fold_threshold = fold_thresh, 
                                    upregulated_color = up_color, 
                                    downregulated_color = down_color, 
                                    nonsignificant_color = non_sig_color, 
                                    point_size = pt_size, 
                                    label_size = lbl_size, 
                                    plot_title = "MA Plot: SenesMtD vs Senes", 
                                    x_label = x_label, 
                                    y_label = y_label)

png("/Users/claire/Downloads/MA_plot_2.png", height = 500, width = 500)
print(ma_plot_senesMtD_vs_senes)
dev.off()

# MA Plot for SenesMtD vs Prolif
ma_plot_senesMtD_vs_prolif = plot_MA(result_senesMtD_vs_prolif$master, 
                                     mean_column = "mean", 
                                     fold_threshold = fold_thresh, 
                                     upregulated_color = up_color, 
                                     downregulated_color = down_color, 
                                     nonsignificant_color = non_sig_color, 
                                     point_size = pt_size, 
                                     label_size = lbl_size, 
                                     plot_title = "MA Plot: SenesMtD vs Prolif", 
                                     x_label = x_label, 
                                     y_label = y_label)

png("/Users/claire/Downloads/MA_plot_3.png", height = 500, width = 500)
print(ma_plot_senesMtD_vs_prolif)
dev.off()
----------------------------------------------------------------------------------
# Venn Diagrams
  
# Venn diagram customization variables
venn_labels = c("Senes vs Prolif", "SenesMtD vs Prolif", "SenesMtD vs Senes")

# Prepare the list of master_sig tables
master_sig_list = list(result_senes_vs_prolif$master_sig,
                       result_senesMtD_vs_prolif$master_sig,
                       result_senesMtD_vs_senes$master_sig)

# Generate the Venn diagram and get gene lists
venn_results = venn_and_overlap(master_sig_list, labels = venn_labels)

# Print the Venn diagram

png("/Users/claire/Downloads/VennDiagram1.png", height = 500, width = 500)
print(venn_results$plot)
dev.off()

# Extract gene lists
group1_only = venn_results$group1_only
group2_only = venn_results$group2_only
group3_only = venn_results$group3_only
group2_and_group3 = venn_results$group2_and_group3
overlap_all = venn_results$overlap_all

# Print the number of shared genes between Group 2 and Group 3
print(paste("Number of shared genes between Group 2 and Group 3:", length(group2_and_group3)))

# Hypergeometric test
group1_size = length(row.names(result_senes_vs_prolif$master_sig))
group2_size = length(row.names(result_senesMtD_vs_prolif$master_sig))
total_genes = nrow(result_senes_vs_prolif$master)  # Assuming total gene count is the same across all comparisons
overlap_size = length(overlap_all)

# Calculate hypergeometric p-value
p_value = phyper(overlap_size - 1, group2_size, total_genes - group2_size, group1_size, lower.tail = FALSE)
print(paste("Hypergeometric test p-value:", p_value))
-----------------------------------------------------------------------------------
# FOLD VS FOLD PLOTS  
  
# Fold vs Fold plot customization variables
point_size = 2
x_label = "Log2 Fold Change (Senes vs Prolif)"
y_label = "Log2 Fold Change (SenesMtD vs Prolif)"

# 1. Fold vs Fold Plot: Senes vs Prolif (X) vs SenesMtD vs Prolif (Y)
fold_vs_fold_result1 = plot_fold_vs_fold(result_senes_vs_prolif$master, 
                                         result_senesMtD_vs_prolif$master, 
                                         common_column = "SYMBOL", 
                                         fold_col1 = "log2fold", 
                                         fold_col2 = "log2fold", 
                                         point_size = point_size, 
                                         title = "Fold vs Fold Plot: Senes vs Prolif vs SenesMtD vs Prolif", 
                                         x_label = x_label, 
                                         y_label = y_label)

print(fold_vs_fold_result1$correlation)
png("/Users/claire/Downloads/Fold_plot_1.png", height = 500, width = 500)
print(fold_vs_fold_result1$plot)
dev.off()

# 2. Fold vs Fold Plot: SenesMtD vs Prolif (X) vs SenesMtD vs Senes (Y)
fold_vs_fold_result2 = plot_fold_vs_fold(result_senesMtD_vs_prolif$master, 
                                         result_senesMtD_vs_senes$master, 
                                         common_column = "SYMBOL", 
                                         fold_col1 = "log2fold", 
                                         fold_col2 = "log2fold", 
                                         point_size = point_size, 
                                         title = "Fold vs Fold Plot: SenesMtD vs Prolif vs SenesMtD vs Senes", 
                                         x_label = x_label, 
                                         y_label = y_label)

print(fold_vs_fold_result2$correlation)
png("/Users/claire/Downloads/Fold_plot_2.png", height = 500, width = 500)
print(fold_vs_fold_result2$plot)
dev.off()
---------------------------------------------------------------------------------
# PCA PLOTS  
  
# PCA plot customization variables
pc_x = "PC1"
pc_y = "PC2"
point_size = 3
pca_title = "PCA Plot: Expression Data"

# Generate PCA plot for the scaled expression matrix (em_scaled)
pca_plot = plot_pca(result_senes_vs_prolif$em_scaled, sample_info = sample_sheet, sample_col = "row.names", group_col = "SAMPLE_GROUP", 
                    pc_x = pc_x, pc_y = pc_y, 
                    point_size = point_size, title = pca_title)

# Print PCA plot
png("/Users/claire/Downloads/PCA_plot.png", height = 500, width = 500)
print(pca_plot)
dev.off()
----------------------------------------------------------------------------------
# DENSITY PLOTs  
  
# Density plot customization variables
facet_columns = 3  # Number of columns for the facets
alpha_value = 0.6  # Transparency for density plots
line_size = 1      # Line size for density curves

# Generate the faceted density plot for em_symbols
density_plot = plot_exp_density(result_senes_vs_prolif$em_symbols, 
                                facet_columns = facet_columns, 
                                alpha_value = alpha_value, 
                                line_size = line_size)

png("/Users/claire/Downloads/den_plot.png", height = 500, width = 500)
print(density_plot) # I printed and screenshoted the plot instead of saving 
dev.off()
---------------------------------------------------------------------------------
# HEAT MAPS
  
# Heatmap customization variables
cluster_x = TRUE
cluster_y = TRUE
color_palette = c("blue", "white", "red")  # Customize as needed

# Venn Diagram function returns list of unique sig genes for senes_vs_prolif as vector: group1_only
# Subset em_scaled_sig for group1_only genes (group1: unique sig genes in senes vs prolif in Venn Diagram)
em_scaled_group1 = result_senes_vs_prolif$em_scaled_sig[row.names(result_senes_vs_prolif$em_scaled_sig) %in% group1_only, ]

# Subset em_scaled_sig for group3_only genes (group3: unique sig genes in senesMtD vs senes in Venn Diagram)
em_scaled_group3 = result_senesMtD_vs_senes$em_scaled_sig[row.names(result_senesMtD_vs_senes$em_scaled_sig) %in% group3_only, ]

# Subset em_scaled_sig for group2 and group3 shared genes (shared genes between senesMtD vs senes & senesMtD vs prolif)
em_scaled_shared_group2_group3 = result_senesMtD_vs_prolif$em_scaled_sig[row.names(result_senesMtD_vs_prolif$em_scaled_sig) %in% group2_and_group3, ]

# Subset em_scaled_sig for shared genes in all 3 groups/comparisons (82 shared sig genes appeared in all 3 comparisons)
em_scaled_overlap_all = result_senesMtD_vs_prolif$em_scaled_sig[row.names(result_senesMtD_vs_prolif$em_scaled_sig) %in% overlap_all, ]

# Generate the heatmap (group1: unique sig genes in senes vs prolif in Venn Diagram) with a rug
heatmap_plot_group1 = plot_heatmap(em_scaled_group1, 
                            sample_info = sample_sheet, 
                            group_col = "SAMPLE_GROUP", 
                            cluster_x = cluster_x, 
                            cluster_y = cluster_y, 
                            color_palette = color_palette)

# Print the heatmap
print(heatmap_plot_group1)

# Generate heatmap for group3_only genes (group3: unique sig genes in senesMtD vs senes in Venn Diagram)
heatmap_plot_group3 = plot_heatmap(em_scaled_group3, 
                                   sample_info = sample_sheet, 
                                   group_col = "SAMPLE_GROUP", 
                                   cluster_x = cluster_x, 
                                   cluster_y = cluster_y, 
                                   color_palette = color_palette)

print(heatmap_plot_group3)

# Generate heatmap for group2 and group3 shared sig genes (shared genes between senesMtD vs senes & senesMtD vs prolif)
heatmap_plot_group2_and_group3 = plot_heatmap(em_scaled_shared_group2_group3, 
                                   sample_info = sample_sheet, 
                                   group_col = "SAMPLE_GROUP", 
                                   cluster_x = cluster_x, 
                                   cluster_y = cluster_y, 
                                   color_palette = color_palette)

print(heatmap_plot_group2_and_group3)

# Generate heatmap for shared genes in all 3 groups/comparisons (82 shared sig genes appeared in all 3 comparisons)
heatmap_plot_overlap_all = plot_heatmap(em_scaled_overlap_all, 
                                              sample_info = sample_sheet, 
                                              group_col = "SAMPLE_GROUP", 
                                              cluster_x = cluster_x, 
                                              cluster_y = cluster_y, 
                                              color_palette = color_palette)

print(heatmap_plot_overlap_all)
---------------------------------------------------------------------------------
  
# Pathway analysis
#ORA
#run ORA for upregulated and downregulated genes
# Run ORA for Group 1
ora_results_group1 = run_ora(result_senes_vs_prolif$master_sig, group1_only)

# Run ORA for Group 3
ora_results_group3 = run_ora(result_senesMtD_vs_senes$master_sig, group3_only)

# Run ORA for Group 2 and Group 3 shared genes
ora_results_group2_and_3 = run_ora(result_senesMtD_vs_prolif$master_sig, group2_and_group3)

# Visualize ORA for Group 1 Upregulated  and Downregulated Genes
barplot(ora_results_group1$upregulated, showCategory = 10) + ggtitle("Group 1 Upregulated ORA")
barplot(ora_results_group1$downregulated, showCategory = 10) + ggtitle("Group 1 Downregulated ORA")

# Visualize ORA for Group 3 Upregulated  and Downregulated Genes
barplot(ora_results_group3$upregulated, showCategory = 10) + ggtitle("Group 3 Upregulated ORA")
barplot(ora_results_group3$downregulated, showCategory = 10) + ggtitle("Group 3 Downregulated ORA")

# Visualize ORA for Group 2 and Group 3 shared Upregulated  and Downregulated Genes
barplot(ora_results_group2_and_3$upregulated, showCategory = 10) + ggtitle("Group 2 and 3 Upregulated ORA")
barplot(ora_results_group2_and_3$downregulated, showCategory = 10) + ggtitle("Group 2 and 3 Downregulated ORA")
---------------------------------------------------------------------------------
# Group 1 (senes vs prolif) upregulated genes   
  
# Extract ORA results as a data frame for group 1 Upregulated genes
ora_results_table_group1_up = extract_ora_results(ora_results_group1$upregulated)

# Prepare candidate genes: genes from most enriched gene set for group 1 Upregulated genes
candidate_genes_group1_up = get_candidate_genes(ora_results_table_group1_up, rank = 1) 

# Filter expression matrices for candidate genes: genes from most enriched gene set for group 1 Upregulated genes
filtered_em_group1_up = em_scaled_group1[row.names(em_scaled_group1) %in% candidate_genes_group1_up, ] #use this one

# Genes in most enriched gene set in group 1 Upregulated: Plot multi-gene boxplots using plot_gene_boxplot function
boxplot_group1_up = plot_gene_boxplot_facet(filtered_em_group1_up)

# Print the boxplots
print(boxplot_group1_up)
----------------------------------------------------------------------------------------------------
# Group 1 (senes vs prolif) downregulated genes 
  
# Extract ORA results as a data frame for group 1 downregulated genes
ora_results_table_group1_down = extract_ora_results(ora_results_group1$downregulated)

# Prepare candidate genes: genes from most enriched gene set for group 1 downregulated genes
candidate_genes_group1_down = get_candidate_genes(ora_results_table_group1_down, rank = 1) 

# Filter expression matrices for candidate genes: genes from most enriched gene set for group 1 downregulated genes
filtered_em_group1_down = em_scaled_group1[row.names(em_scaled_group1) %in% candidate_genes_group1_down, ] 

# Genes in most enriched gene set in group 1 downregulated: Plot multi-gene boxplots using plot_gene_boxplot function
boxplot_group1_down = plot_gene_boxplot_facet(filtered_em_group1_down)

# Print the boxplots
print(boxplot_group1_down)
---------------------------------------------------------------------------------
# Group 3 (senesMtD vs senes) upregulated genes   
  
# Extract ORA results as a data frame for group 3 Upregulated genes
ora_results_table_group3_up = extract_ora_results(ora_results_group3$upregulated)

# Prepare candidate genes: genes from most enriched gene set for group 3 Upregulated genes
candidate_genes_group3_up = get_candidate_genes(ora_results_table_group3_up, rank = 1) 

# Filter expression matrices for candidate genes: genes from most enriched gene set for group 3 Upregulated genes
filtered_em_group3_up = em_scaled_group3[row.names(em_scaled_group3) %in% candidate_genes_group3_up, ] 

# Genes in most enriched gene set in group 3 Upregulated: Plot multi-gene boxplots using plot_gene_boxplot function
boxplot_group3_up = plot_gene_boxplot_facet(filtered_em_group3_up)

# Print the boxplots
print(boxplot_group3_up)
---------------------------------------------------------------------------------
# Group 3 (senesMtD vs senes) downregulated genes   
  
# Extract ORA results as a data frame for group 3 downregulated genes
ora_results_table_group3_down = extract_ora_results(ora_results_group3$downregulated)

# Prepare candidate genes: genes from most enriched gene set for group 3 downregulated genes
candidate_genes_group3_down = get_candidate_genes(ora_results_table_group3_down, rank = 1) 

# Filter expression matrices for candidate genes: genes from most enriched gene set for group 3 downregulated genes
filtered_em_group3_down = em_scaled_group3[row.names(em_scaled_group3) %in% candidate_genes_group3_down, ] 

# Genes in most enriched gene set in group 3 downregulated: Plot multi-gene boxplots using plot_gene_boxplot function
boxplot_group3_down = plot_gene_boxplot_facet(filtered_em_group3_down)

# Print the boxplots
print(boxplot_group3_down)
  
----------------------------------------------------------------------------------
# GSEA customization variables
ont = "BP"
keyType = "SYMBOL"
nPerm = 10000
minGSSize = 3
maxGSSize = 800
pvalueCutoff = 0.05
OrgDb = org.Hs.eg.db
pAdjustMethod = "none"

# List of results and their labels
result_list <- list(
  "Senes_vs_Prolif" = result_senes_vs_prolif$master,
  "SenesMtD_vs_Prolif" = result_senesMtD_vs_prolif$master,
  "SenesMtD_vs_Senes" = result_senesMtD_vs_senes$master
)

# Function to prepare the input vector and run GSEA for each comparison
gsea_results <- lapply(names(result_list), function(label) {
  cat("Running GSEA for", label, "...\n")
  
  # Prepare the log2fold change vector
  gsea_input <- result_list[[label]]$log2fold
  names(gsea_input) <- row.names(result_list[[label]])
  gsea_input <- na.omit(gsea_input)
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  # Run GSEA using the run_gsea function
  gse_results <- run_gsea(
    log2fc_vector = gsea_input,
    ont = ont,
    keyType = keyType,
    nPerm = nPerm,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pvalueCutoff = pvalueCutoff,
    OrgDb = OrgDb,
    pAdjustMethod = pAdjustMethod
  )
  
  # Return the results with a label
  list(label = label, gse_results = gse_results)
})

# Step 3: Visualize the results for each comparison
for (res in gsea_results) {
  cat("Visualizing GSEA results for", res$label, "...\n")
  plot <- ridgeplot(res$gse_results, showCategory = 10) + 
    ggtitle(paste("GSEA Ridgeplot for Top Enriched Terms -", res$label))
  print(plot)  # Ensure plot is printed in the loop
}
---------------------------------------------------------------------------------

# STRING  
# Variables
species_id = 9606  # Human species ID in STRING
score_threshold = 400  # Minimum interaction score threshold for STRING

# Extract candidate genes from the ORA results table
enriched_pathway = 1  # Change this to select a different pathway
candidate_genes = get_candidate_genes(ora_results_table_group1, rank = enriched_pathway)

# Run STRING analysis for the candidate genes
string_interactions = get_string_network(
  gene_list = candidate_genes, 
  species_id = species_id, 
  score_threshold = score_threshold
)

# Check and print the first few rows of the interaction network
if (!is.null(string_interactions)) {
  head(string_interactions)
}
#-----------------------------------END OF SCRIPT----------------------------------#








