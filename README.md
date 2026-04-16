# Bulk RNA-seq Analysis of Cellular Senescence

## Overview
This project investigates transcriptomic changes in human fibroblasts during replicative senescence and mitochondrial depletion using bulk RNA sequencing data.

## Biological Question
- What happens to the fibroblast transcriptome during senescence?
- How does mitochondrial depletion affect gene expression in senescent cells?

## Dataset
- Human IMR90 fibroblasts
- 3 conditions:
  - Proliferating (prolif)
  - Senescent (senes)
  - Senescent with mitochondrial depletion (senes_MtD)
- 3 replicates per condition (9 samples total)

## Methods
- Differential gene expression analysis in R. There are 2 R scripts. 
- One R script includes all the analysis, the other R script contains the functions. The analysis script should source the functions script. 
- Visualization:
  - PCA
  - Volcano plots
  - Heatmaps
- Functional focus:
  - Cell cycle genes
  - Immune-related genes
  - Extracellular matrix genes

## Key Findings
- Senescent cells show reduced expression of cell cycle genes
- Mitochondrial depletion alters immune-related pathways
- Clear transcriptomic differences between all three conditions

## Reproducibility Note
This analysis was originally performed in a university environment using specific datasets and configurations. The script and report are provided to demonstrate the workflow and interpretation.

## Skills Demonstrated
- RNA-seq data analysis
- R programming
- Differential expression analysis
- Data visualization
- Biological interpretation of transcriptomic data
