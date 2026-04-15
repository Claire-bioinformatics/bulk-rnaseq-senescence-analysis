# Function to parse de tables, em, annotations to produce various master and em tables
parse_de_table = function(de_table, em, annotations) {
  # Step 1: Merge expression matrix with annotations
  master_temp = merge(em, annotations, by.x = 0, by.y = 0)
  
  # Step 2: Merge the differential expression table with master_temp
  master = merge(de_table, master_temp, by.x = 0, by.y = 1)
  
  # Step 3: Set row names to the SYMBOL column and rename the first column to "gene ID"
  row.names(master) = master$SYMBOL
  names(master)[1] = "gene ID"
  
  # Step 4: Extract expression columns 
  em_symbols = master[, c(5,6,7,8,9,10,11,12,13)]

  # Step 5: Remove rows with missing values
  master = na.omit(master)
  
  # Step 6: Sort by absolute log2 fold change (descending) and then by adjusted p-value (ascending)
  sorted_order = order(-abs(master$log2fold), master$p.adj)
  master = master[sorted_order, ]
  
  # Step 7: Calculate mean expression 
  master$mean = rowMeans(master[, c(5,6,7,8,9,10,11,12,13)])

  # Step 8: Ensure p-values are greater than zero to avoid log10(0) issues
  master$p[master$p == 0] = min(master$p[master$p > 0], na.rm = TRUE) / 10
  
  # Step 9: Add -log10(p-value) for easier visualization
  master$mlog10p = -log10(master$p)
  
  # Step 10: Flag significant genes (p.adj < 0.01 and abs(log2fold) > 2)
  master$sig = as.factor(master$p.adj < 0.01 & abs(master$log2fold) > 2.0)
  
  # Step 11: Scale the expression matrix by rows and remove NAs
  em_scaled = na.omit(data.frame(t(scale(t(em_symbols)))))
  
  # Step 12: Create a subset of significant genes
  master_sig = master[master$sig == TRUE, ]
  
  # Step 13: Extract significant gene names as a vector
  sig_genes = rownames(master_sig)
  
  # Step 14: Subset original and scaled expression matrices for significant genes
  em_symbols_sig = em_symbols[rownames(em_symbols) %in% sig_genes, ]
  em_scaled_sig = em_scaled[rownames(em_scaled) %in% sig_genes, ]
  
  # Return all the generated tables as a list
  return(list(
    master = master,
    em_symbols = em_symbols,
    em_scaled = em_scaled,
    master_sig = master_sig,
    sig_genes = sig_genes,
    em_symbols_sig = em_symbols_sig,
    em_scaled_sig = em_scaled_sig
  ))
}



# Volcano plot function
plot_volcano = function(data, 
                        x_col = "log2fold", 
                        y_col = "mlog10p", 
                        symbol_col = "SYMBOL", 
                        p_adj_col = "p.adj", 
                        significance_col = "Significance",
                        
                        p_threshold = 0.01, 
                        fold_threshold = 2,
                        
                        upregulated_color = "red", 
                        downregulated_color = "blue", 
                        nonsignificant_color = "gray",
                        
                        title = "Volcano Plot", 
                        x_label = "Log2 Fold Change", 
                        y_label = "-log10 Adjusted P-Value",
                        
                        point_size = 2, 
                        label_size = 4) 
{
  # Ensure necessary columns exist
  required_cols = c(x_col, y_col, symbol_col, p_adj_col)
  if (!all(required_cols %in% colnames(data))) {
    stop(paste("The data frame must contain the following columns:", paste(required_cols, collapse = ", ")))
  }
  
  # Remove rows with NA values in key columns
  data = na.omit(data)
  
  # Add significance column
  data[[significance_col]] = "Not Significant"
  data[[significance_col]][data[[p_adj_col]] < p_threshold & data[[x_col]] > fold_threshold] = "Upregulated"
  data[[significance_col]][data[[p_adj_col]] < p_threshold & data[[x_col]] < -fold_threshold] = "Downregulated"
  
  # Subset significantly upregulated and downregulated genes
  data_sig = subset(data, data[[p_adj_col]] < p_threshold & abs(data[[x_col]]) > fold_threshold)
  data_sig_up = subset(data_sig, data[[x_col]] > fold_threshold)
  data_sig_down = subset(data_sig, data[[x_col]] < -fold_threshold)
  
  # Select top 10 up/downregulated genes for labeling (EXCLUDING NA SYMBOLS)
  data_sig_up = data_sig_up[!is.na(data_sig_up[[symbol_col]]), ]
  data_sig_down = data_sig_down[!is.na(data_sig_down[[symbol_col]]), ]
  
  top10_up = head(data_sig_up[order(-data_sig_up[[x_col]]), ], 10)
  top10_down = head(data_sig_down[order(data_sig_down[[x_col]]), ], 10)
  top10_genes = rbind(top10_up, top10_down)
  
  # Set axis limits dynamically
  x_min = min(data[[x_col]], na.rm = TRUE) - 1
  x_max = max(data[[x_col]], na.rm = TRUE) + 1
  y_max = max(data[[y_col]], na.rm = TRUE) + 1
  
  # Generate Volcano plot
  plot = ggplot(data, aes_string(x = x_col, y = y_col, color = significance_col)) +
    geom_point(alpha = 0.6, size = point_size) +
    scale_color_manual(values = c("Upregulated" = upregulated_color, 
                                  "Downregulated" = downregulated_color, 
                                  "Not Significant" = nonsignificant_color)) +
    geom_vline(xintercept = c(-fold_threshold, fold_threshold), linetype = "dashed") +
    geom_hline(yintercept = -log10(p_threshold), linetype = "dashed") +
    geom_text_repel(
      data = top10_genes, aes_string(label = symbol_col),
      size = label_size, color = "black",
      box.padding = 1, 
      point.padding = 0.5,
      nudge_y = 1,
      max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(0, y_max)) + 
    theme_minimal(base_size = 14) +
    labs(title = title, x = x_label, y = y_label, color = "Significance")
  
  return(plot)
}



# MA Plot Function
plot_MA = function(de_table, mean_column = "mean", fold_threshold = 2,
                   upregulated_color = "red", downregulated_color = "blue", nonsignificant_color = "gray",
                   point_size = 2, label_size = 4,
                   plot_title = "MA Plot", x_label = "Log10 Mean Expression", y_label = "Log2 Fold Change") {
  
  # Ensure necessary columns exist
  if (!all(c(mean_column, "log2fold") %in% colnames(de_table))) {
    stop(paste("Input data frame must contain '", mean_column, "' and 'log2fold' columns.", sep = ""))
  }
  
  # Calculate log10 of the mean expression (avoiding log10(0))
  de_table$log10_mean = log10(de_table[[mean_column]] + 1)  # Add 1 to avoid issues with log10(0)
  
  # Remove rows with NA values in key columns
  de_table = na.omit(de_table)
  
  # Add significance column
  de_table$Significance = "Not Significant"
  de_table$Significance[de_table$log2fold > fold_threshold] = "Upregulated"
  de_table$Significance[de_table$log2fold < -fold_threshold] = "Downregulated"
  
  # Select top 10 upregulated and downregulated genes for labeling
  de_table_sig_up = subset(de_table, log2fold > fold_threshold)
  de_table_sig_down = subset(de_table, log2fold < -fold_threshold)
  
  top10_up = head(de_table_sig_up[order(-de_table_sig_up$log2fold), ], 10)
  top10_down = head(de_table_sig_down[order(de_table_sig_down$log2fold), ], 10)
  top10_genes = rbind(top10_up, top10_down)
  
  # Set axis limits dynamically
  x_min = min(de_table$log10_mean, na.rm = TRUE) - 0.5
  x_max = max(de_table$log10_mean, na.rm = TRUE) + 0.5
  y_min = min(de_table$log2fold, na.rm = TRUE) - 1
  y_max = max(de_table$log2fold, na.rm = TRUE) + 1
  
  # Generate MA Plot
  ggp = ggplot(de_table, aes(x = log10_mean, y = log2fold, color = Significance)) +
    geom_point(alpha = 0.6, size = point_size) +
    scale_color_manual(values = c("Upregulated" = upregulated_color,
                                  "Downregulated" = downregulated_color,
                                  "Not Significant" = nonsignificant_color)) +
    geom_hline(yintercept = c(-fold_threshold, fold_threshold), linetype = "dashed") +
    geom_text_repel(
      data = top10_genes, aes(label = rownames(top10_genes)),
      size = label_size, color = "black",
      box.padding = 1, point.padding = 0.5, nudge_y = 0.5, max.overlaps = Inf
    ) +
    coord_cartesian(xlim = c(x_min, x_max), ylim = c(y_min, y_max)) +
    theme_minimal() +
    labs(title = plot_title, x = x_label, y = y_label, color = "Significance")
  
  return(ggp)
}



# Function to generate a Venn diagram from significant gene lists
plot_venn = function(master_sig_list, labels = c("Group 1", "Group 2", "Group 3"),
                     title = "Venn Diagram of Significant Genes") {
  
  # Extract row names (significant gene SYMBOLs) from each master_sig table
  sig_list = lapply(master_sig_list, row.names)
  
  # Create a named list for Venn diagram
  venn_data = setNames(sig_list, labels)
  
  # Generate the Venn diagram using eulerr
  venn_plot = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, main = title)
  
  # Return the plot
  return(venn_plot)
}

# Function to generate Venn diagram and return unique/overlapping genes
venn_and_overlap = function(master_sig_list, labels = c("Group 1", "Group 2", "Group 3")) {
  # Extract significant gene SYMBOLs from each master_sig table
  sig_list = lapply(master_sig_list, row.names)
  # Extract the first column ("gene ID") from each master_sig table
  #sig_list = lapply(master_sig_list, function(x) x[, 1])
  
  
  # Create a named list for Venn diagram
  venn_data = setNames(sig_list, labels)
  
  # Generate the Venn diagram
  venn_plot = plot(euler(venn_data, shape = "ellipse"), quantities = TRUE, main = "Venn Diagram")
  
  # Identify unique and overlapping genes
  group1_only = setdiff(sig_list[[1]], union(sig_list[[2]], sig_list[[3]]))
  group2_only = setdiff(sig_list[[2]], union(sig_list[[1]], sig_list[[3]]))
  group3_only = setdiff(sig_list[[3]], union(sig_list[[1]], sig_list[[2]]))
  
  # Shared genes between specific groups
  #group2_and_group3 = intersect(sig_list[[2]], sig_list[[3]])
  #group2_and_group3 = intersect(unique(sig_list[[2]]), unique(sig_list[[3]]))
  group2_and_group3 = setdiff(intersect(sig_list[[2]], sig_list[[3]]), sig_list[[1]])
  
  overlap_all = Reduce(intersect, sig_list)
  
  # Return the plot and gene lists
  return(list(
    plot = venn_plot,
    group1_only = group1_only,
    group2_only = group2_only,
    group3_only = group3_only,
    group2_and_group3 = group2_and_group3,
    overlap_all = overlap_all
  ))
}



# Function to generate Fold vs Fold plot and calculate Pearson correlation
plot_fold_vs_fold = function(master_table1, master_table2, common_column = "SYMBOL", 
                             fold_col1 = "log2fold", fold_col2 = "log2fold", 
                             point_size = 2, title = "Fold vs Fold Plot", 
                             x_label = "Log2 Fold Change (Comparison 1)", 
                             y_label = "Log2 Fold Change (Comparison 2)") {
  
  # Merge the two tables by the common column (e.g., SYMBOL)
  combined_table = merge(master_table1, master_table2, by = common_column, suffixes = c(".comp1", ".comp2"))
  
  # Ensure the fold change columns exist
  if (!all(c(paste0(fold_col1, ".comp1"), paste0(fold_col2, ".comp2")) %in% colnames(combined_table))) {
    stop("Fold change columns not found in the combined table.")
  }
  
  # Generate Fold vs Fold scatter plot
  fold_plot = ggplot(combined_table, aes_string(x = paste0(fold_col1, ".comp1"), y = paste0(fold_col2, ".comp2"))) +
    geom_point(alpha = 0.6, size = point_size) +
    geom_smooth(method = "lm", se = FALSE, color = "red", linetype = "dashed") +
    theme_minimal() +
    labs(title = title, x = x_label, y = y_label)
  
  # Perform Pearson correlation test
  correlation_result = cor.test(combined_table[[paste0(fold_col1, ".comp1")]], combined_table[[paste0(fold_col2, ".comp2")]], method = "pearson")
  
  # Return both the plot and the correlation result
  return(list(plot = fold_plot, correlation = correlation_result))
}



# Function to generate PCA plot without sample labels
plot_pca = function(expression_matrix, sample_info, sample_col, group_col, 
                    pc_x = "PC1", pc_y = "PC2", 
                    point_size = 3, title = "PCA Plot") {
  
  # Convert expression matrix to numeric matrix and perform PCA
  numeric_matrix = as.matrix(sapply(expression_matrix, as.numeric))
  pca = prcomp(t(numeric_matrix), scale = TRUE)
  
  # Extract PCA coordinates and convert to data frame
  pca_coordinates = data.frame(pca$x)
  pca_coordinates$sample = rownames(pca_coordinates)
  
  # Merge with sample information for group/color labels
  pca_coordinates = merge(pca_coordinates, sample_info, by.x = "sample", by.y = sample_col)
  
  # Calculate % variance for axis labels
  vars = apply(pca$x, 2, var)
  prop_x = round(vars[pc_x] / sum(vars), 4) * 100
  prop_y = round(vars[pc_y] / sum(vars), 4) * 100
  x_axis_label = paste(pc_x, " (", prop_x, "%)", sep = "")
  y_axis_label = paste(pc_y, " (", prop_y, "%)", sep = "")
  
  # Generate PCA plot without sample labels
  pca_plot = ggplot(pca_coordinates, aes_string(x = pc_x, y = pc_y, color = group_col)) +
    geom_point(size = point_size, alpha = 0.8) +
    scale_color_manual(values = c("red", "blue", "green")) +  # Customize colors as needed
    theme_minimal() +
    labs(title = title, x = x_axis_label, y = y_axis_label, color = "Group") +
    theme(axis.text = element_text(size = 12), axis.title = element_text(size = 14))
  
  return(pca_plot)
}



# Function to create faceted density plots for each sample
plot_exp_density = function(expression_matrix, facet_columns = 3, alpha_value = 0.6, line_size = 1) {
  # Melt the expression matrix
  expression_matrix_melted = melt(expression_matrix, variable.name = "Sample", value.name = "Expression")
  
  # Generate the density plot
  density_plot = ggplot(expression_matrix_melted, aes(x = log10(Expression + 1), color = Sample, fill = Sample)) +
    geom_density(alpha = alpha_value, size = line_size) +
    facet_wrap(~ Sample, ncol = facet_columns, scales = "free_y") +
    theme_minimal() +
    labs(title = "Expression Density Plot for Each Sample", 
         x = "Log10(Expression Value)", y = "Density") +
    theme(
      panel.grid = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 1),
      strip.text.x = element_text(size = 12, family = "Arial", face = "bold"),
      panel.spacing = unit(1, "lines"),
      legend.position = "none"
    )
  
  return(density_plot)
}



# Function to generate a heatmap with clustering
plot_heatmap = function(expression_matrix, sample_info, group_col, cluster_x = TRUE, cluster_y = TRUE, color_palette = c("blue", "white", "red")) {
  # Convert the expression matrix to a numeric matrix
  hm_matrix = as.matrix(expression_matrix)
  
  # Cluster rows (genes) if cluster_y is TRUE
  if (cluster_y) {
    y_dist = Dist(hm_matrix, method = "spearman")
    y_cluster = hclust(y_dist, method = "average")
    y_dd = as.dendrogram(y_cluster)
    y_order = order.dendrogram(y_dd)
    hm_matrix = hm_matrix[y_order, ]  # Reorder rows
  }
  
  # Cluster columns (samples) if cluster_x is TRUE
  if (cluster_x) {
    x_dist = Dist(t(hm_matrix), method = "spearman")
    x_cluster = hclust(x_dist, method = "average")
    x_dd = as.dendrogram(x_cluster)
    x_order = order.dendrogram(x_dd)
    hm_matrix = hm_matrix[, x_order]  # Reorder columns
  }
  
  # Melt the clustered matrix for plotting
  hm_matrix_melted = melt(hm_matrix)
  
  # Generate the heatmap
  heatmap_plot = ggplot(hm_matrix_melted, aes(x = Var2, y = Var1, fill = value)) +
    geom_tile() +
    scale_fill_gradientn(colours = colorRampPalette(color_palette)(100)) +
    theme_minimal() +
    labs(title = "Heatmap of Scaled Expression", x = "", y = "") +
    theme(
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      legend.title = element_blank(),
      legend.spacing.x = unit(0.25, 'cm')
    )
  
  # Generate the rug (sample group color bar)
  group_data = as.matrix(as.numeric(as.factor(sample_info[[group_col]])))
  group_data_melted = melt(group_data)
  rug_plot = ggplot(group_data_melted, aes(x = Var1, y = 1, fill = value)) +
    geom_tile(linetype = "blank") +
    scale_fill_gradientn(colours = gg_color_hue(length(unique(sample_info[[group_col]])))) +
    labs(x = "", y = "") +
    theme(
      legend.position = "none",
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      panel.background = element_blank()
    )
  
  # Combine heatmap and rug using patchwork
  library(patchwork)
  combined_plot = heatmap_plot / rug_plot
  
  return(combined_plot)
}

# Function to get default ggplot colors
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}



# Function to perform GSEA using gseGO from clusterProfiler
run_gsea <- function(log2fc_vector, 
                     ont = "BP", 
                     keyType = "SYMBOL", 
                     nPerm = 10000, 
                     minGSSize = 3, 
                     maxGSSize = 800, 
                     pvalueCutoff = 0.05, 
                     OrgDb = org.Hs.eg.db, 
                     pAdjustMethod = "none") {
  
  # Ensure the input vector is sorted and has no NA values
  gsea_input <- na.omit(log2fc_vector)
  gsea_input <- sort(gsea_input, decreasing = TRUE)
  
  # Run the GSEA
  gse_results <- gseGO(geneList = gsea_input,
                       ont = ont,
                       keyType = keyType,
                       nPerm = nPerm,
                       minGSSize = minGSSize,
                       maxGSSize = maxGSSize,
                       pvalueCutoff = pvalueCutoff,
                       verbose = TRUE,
                       OrgDb = OrgDb,
                       pAdjustMethod = pAdjustMethod)
  
  return(gse_results)
}

# Function to extract GSEA results
extract_gsea_results = function(gsea_results) {
  # Extract relevant columns from the GSEA result
  gene_sets = gsea_results@result$core_enrichment
  description = gsea_results@result$Description
  p_adj = gsea_results@result$p.adjust
  
  # Create a data frame with gene sets and p-values, using the description as row names
  gsea_results_table = data.frame(
    gene_sets = gene_sets,
    p_adj = p_adj,
    row.names = description
  )
  
  return(gsea_results_table)
}

# Function to get candidate genes from a GSEA result
get_candidate_genes_gsea = function(gsea_results_table, rank = 1) {
  # Extract the gene set from the specified rank
  enriched_gene_set = as.character(gsea_results_table[rank, "gene_sets"])
  
  # Split the gene set into a vector of individual genes
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  
  return(candidate_genes)
}

# Function to visualize GSEA results (ridgeplot and dotplot)
plot_gsea_results = function(gsea_results, plot_type = "ridgeplot", showCategory = 10) {
  library(clusterProfiler)
  library(ggplot2)
  
  if (plot_type == "ridgeplot") {
    plot = ridgeplot(gsea_results, showCategory = showCategory) +
      ggtitle("GSEA Ridgeplot for Top Enriched Terms")
  } else if (plot_type == "dotplot") {
    plot = dotplot(gsea_results, showCategory = showCategory) +
      ggtitle("GSEA Dotplot for Top Enriched Terms")
  } else {
    stop("Invalid plot_type. Choose 'ridgeplot' or 'dotplot'.")
  }
  
  return(plot)
}



# Function to split genes, convert to ENTREZ, and run ORA
run_ora = function(master_sig, gene_list, ont = "ALL", OrgDb = org.Hs.eg.db, 
                   pvalueCutoff = 0.05, qvalueCutoff = 0.10) {
  # Filter the master_sig table for the gene list
  filtered_master_sig = master_sig[row.names(master_sig) %in% gene_list, ]
  
  # Split into upregulated and downregulated
  upregulated = row.names(subset(filtered_master_sig, log2fold > 0))
  downregulated = row.names(subset(filtered_master_sig, log2fold < 0))
  
  # Convert SYMBOLs to ENTREZ IDs
  upregulated_entrez = bitr(upregulated, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  downregulated_entrez = bitr(downregulated, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
  
  # Run ORA for upregulated genes
  ora_results_up = enrichGO(gene = upregulated_entrez$ENTREZID, OrgDb = OrgDb, readable = TRUE, 
                            ont = ont, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
  
  # Run ORA for downregulated genes
  ora_results_down = enrichGO(gene = downregulated_entrez$ENTREZID, OrgDb = OrgDb, readable = TRUE, 
                              ont = ont, pvalueCutoff = pvalueCutoff, qvalueCutoff = qvalueCutoff)
  
  # Return the results
  return(list(upregulated = ora_results_up, downregulated = ora_results_down))
}

extract_ora_results = function(ora_results) {
  # Extract relevant columns from the ORA result
  gene_sets = ora_results@result$geneID
  description = ora_results@result$Description
  p_adj = ora_results@result$p.adjust
  
  # Create a data frame with gene sets and p-values, using the description as row names
  ora_results_table = data.frame(
    gene_sets = gene_sets,
    p_adj = p_adj,
    row.names = description
  )
  
  return(ora_results_table)
}

get_candidate_genes = function(ora_results_table, rank = 1) {
  # Extract the gene set from the specified rank
  enriched_gene_set = as.character(ora_results_table[rank, "gene_sets"])
  
  # Split the gene set into a vector of individual genes
  candidate_genes = unlist(strsplit(enriched_gene_set, "/"))
  
  return(candidate_genes)
}

plot_gene_boxplot = function(expression_matrix) {
  library(reshape2)
  library(ggplot2)
  
  # Ensure the matrix is converted to a data frame and retains row names as a column
  expression_df = as.data.frame(expression_matrix)
  expression_df$Gene = row.names(expression_df)
  
  # Melt the data frame into long format for plotting
  expression_melted = melt(expression_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  
  # Generate the box plot
  boxplot_plot = ggplot(expression_melted, aes(x = Sample, y = Expression, fill = Gene)) +
    geom_boxplot(outlier.size = 1, alpha = 0.7) +
    theme_minimal() +
    labs(title = "Multi-Gene Boxplot for Candidate Genes", x = "Sample", y = "Expression Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(boxplot_plot)
}

plot_combined_violin_boxplot = function(expression_matrix) {
  library(reshape2)
  library(ggplot2)
  
  # Convert the expression matrix to long format
  expression_df = as.data.frame(expression_matrix)
  expression_df$Gene = row.names(expression_df)
  expression_melted = melt(expression_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  
  # Generate a combined violin and box plot
  combined_plot = ggplot(expression_melted, aes(x = Gene, y = Expression, fill = Sample)) +
    geom_violin(trim = TRUE, alpha = 0.6) +
    geom_boxplot(width = 0.2, outlier.size = 1, alpha = 0.7, position = position_dodge(width = 0.75)) +
    theme_minimal() +
    labs(title = "Expression of Candidate Genes", x = "Gene", y = "Expression (Z-score)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(combined_plot)
}

plot_gene_boxplot_facet = function(expression_matrix) {
  library(reshape2)
  library(ggplot2)
  
  # Convert the expression matrix to long format
  expression_df = as.data.frame(expression_matrix)
  expression_df$Gene = row.names(expression_df)
  expression_melted = melt(expression_df, id.vars = "Gene", variable.name = "Sample", value.name = "Expression")
  
  # Generate a faceted box plot
  boxplot_plot = ggplot(expression_melted, aes(x = Sample, y = Expression, fill = Gene)) +
    geom_boxplot(outlier.size = 1, alpha = 0.7) +
    theme_minimal() +
    facet_wrap(~ Gene, scales = "free_y", ncol = 4) +
    labs(title = "Multi-Gene Faceted Boxplot for Candidate Genes", x = "Sample", y = "Expression Value") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  
  return(boxplot_plot)
}



# Function to get and visualize STRING network for candidate genes
get_string_network = function(gene_list, species_id = 9606, score_threshold = 400) {
  # Initialize STRINGdb object for human (species_id = 9606)
  string_db = STRINGdb$new(version = "11.5", species = species_id, score_threshold = score_threshold)
  
  # Map the genes to STRING IDs
  mapped_genes = string_db$map(data.frame(SYMBOL = gene_list), "SYMBOL", removeUnmappedRows = TRUE)
  
  if (nrow(mapped_genes) == 0) {
    cat("No genes mapped to STRING IDs.\n")
    return(NULL)
  }
  
  # Get the interaction network for the mapped genes
  interactions = string_db$get_interactions(mapped_genes$STRING_id)
  
  if (nrow(interactions) == 0) {
    cat("No interactions found for the mapped genes.\n")
    return(NULL)
  }
  
  # Plot the interaction network
  string_db$plot_network(mapped_genes$STRING_id)
  
  # Return the interactions data frame
  return(interactions)
}