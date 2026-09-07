# Github Repo for the Article : 'Distinct Temporal Responses to Xanthomonas translucens Strains Shape Bacterial Leaf Streak Development in Triticale'

## Description

Repositories with all the codes for the analysis and figures compiled for this article

## Repository Contents

### 1. RNASeq Pipeline

**Folder:** [RNASeq Pipeline](https://github.com/zlianglab/Triticale_BLS/tree/main/1.%20RNASeq%20Pipeline) \
**Description:** A complete RNAseq pipeline with fastp, hisat2, htseq-count. Can be used with standard fastq files. (This paper used kallisto for read quantiification) \
**Contents:** \
- [rnaseq.pipeline.py](https://github.com/zlianglab/Triticale_BLS/blob/main/1.%20RNASeq%20Pipeline/rnaseq-pipeline.py)

### 2. DESeq2 Analysis

**Folder:** [DESeq2 Analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/2.%20DESeq2%20Analysis) \
**Description:** DESeq2 analysis between mock inoculation as control as LB10- and P3-inoculated leaves. DEGs were identified if padj < 0.05 and |log2FC| > 1 \
**Contents:**
- [DESeq2_run.R](https://github.com/zlianglab/Triticale_BLS/blob/main/2.%20DESeq2%20Analysis/DESeq2_run.R)
- [PCA.R](https://github.com/zlianglab/Triticale_BLS/blob/main/2.%20DESeq2%20Analysis/PCA.R)
- [gene_level_annotation_building.R](https://github.com/zlianglab/Triticale_BLS/blob/main/2.%20DESeq2%20Analysis/gene_level_annotation_building.R)
- [triad_bar_plot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/2.%20DESeq2%20Analysis/triad_bar_plot.R)

### 3. Syntenic Homolog Analysis

**Folder:** [Syntenic homolog analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/3.%20Syntenic%20homoelog%20analysis) \
**Description:** Code for the analysis and figures for syntenic traid analysis in triticale subgenome \
**Contents:** 
- [DEG-ratio-in-syntenic-table.R](https://github.com/zlianglab/Triticale_BLS/blob/main/3.%20Syntenic%20homoelog%20analysis/DEG-ratio-in-syntenic-table.R)
- [Donut-plot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/3.%20Syntenic%20homoelog%20analysis/Donut-plot.R)
- [sankeyplot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/3.%20Syntenic%20homoelog%20analysis/sankeyplot.R)
- [combined_sankey_plot_LB10_P3.R](https://github.com/zlianglab/Triticale_BLS/blob/main/3.%20Syntenic%20homoelog%20analysis/combined_sankey_plot_LB10_P3.R)
- [ternary-plot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/3.%20Syntenic%20homoelog%20analysis/ternary-plot.R)


### 4. K-means Clustering

**Folder:** [K-means clustering](https://github.com/zlianglab/Triticale_BLS/tree/main/4.%20K-means%20clustering) \
**Description:** K-Means clustering, wss, and second-derivative calculation for combined DEGs between LB10 v control and P3 v control \
**Contents:**
- [kmeans-clustering.R](https://github.com/zlianglab/Triticale_BLS/blob/main/4.%20K-means%20clustering/kmeans-clustering.R)


### 5. Gene Ontology Enrichment Analysis

**Folder:** [Gene ontology enrichment analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/5.%20Gene%20ontology%20enrichment%20analysis) \
**Description:** Gene Ontology (GO) analysis between each clusters from K-means \
**Contents:**
- [CPM_file_processing_for_background_gene_list.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/CPM_file_processing_for_background_gene_list.R)
- [gene_seperation_by_cluster.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/gene_seperation_by_cluster.R)
- [go-collapse.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/go-collapse.R)
- [heatmap.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/heatmap.R)
- [piechart-GO.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/piechart-GO.R)
- [piechart.R](https://github.com/zlianglab/Triticale_BLS/blob/main/5.%20Gene%20ontology%20enrichment%20analysis/piechart.R)

### 6. Sequence Similarity Analysis

**Folder:** [Sequence similarity analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/6.%20Sequence%20similarity%20analysis) \
**Description:** Analysis and figure for comparing sequence similarities between CDS sequence and 1000bp upstream sequence from transcription start site (TSS) \
**Contents:**
- [align_homeolog.sh](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/align_homeolog.sh)
- [boxplot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/boxplot.R)
- [ext_seq_from_annotation.py](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/ext_seq_from_annotation.py)
- [mafft_data_processing.R](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/mafft_data_processing.R)
- [pairwise_similarity_check.py](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/pairwise_similarity_check.py)
- [subgenome_expression_bias_between_clusters.R](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/subgenome_expression_bias_between_clusters.R)
- [subgenome_expression_bias_file_processing.R](https://github.com/zlianglab/Triticale_BLS/blob/main/6.%20Sequence%20similarity%20analysis/subgenome_expression_bias_file_processing.R)

### 7. Motif Enrichment Analysis

**Folder:** [Motif enrichment analysis](https://github.com/zlianglab/Triticale_BLS/tree/main#:~:text=7.%20Motif%20enrichment%20analysis) \
**Description:** Analysis for motif enrichment analysis \
**Contents:**
- [arabidopsis_export.R](https://github.com/zlianglab/Triticale_BLS/blob/main/7.%20Motif%20enrichment%20analysis/arabidopsis_export.R)
- [background_gene.R](https://github.com/zlianglab/Triticale_BLS/blob/main/7.%20Motif%20enrichment%20analysis/background_gene.R)
- [interpret_MP.py](https://github.com/zlianglab/Triticale_BLS/blob/main/7.%20Motif%20enrichment%20analysis/interpret_MP.py)

### 8. AnnoTAL Prediction

**Folder:** [AnnoTAL prediction](https://github.com/zlianglab/Triticale_BLS/tree/main/8.%20AnnoTAL%20prediction) \
**Description:** Analysis and figures for the AnnoTALE predicted genes in LB10 and P3 \
**Contents:**
- [heatmap.R](https://github.com/zlianglab/Triticale_BLS/blob/main/8.%20AnnoTAL%20prediction/heatmap.R)
- [lb10.R](https://github.com/zlianglab/Triticale_BLS/blob/main/8.%20AnnoTAL%20prediction/lb10.R)
- [p3.R](https://github.com/zlianglab/Triticale_BLS/blob/main/8.%20AnnoTAL%20prediction/p3.R)

### 9. Chlorophyll Fluorescence Phenotype Analysis

**Folder:** [Chlorophyll fluorescence phenotype analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/9.%20Chlorophyll%20fluoresence%20phenotype%20analysis) \
**Description:** Analysis for chorolophyll fluorescence analysis
**Contents:**
- [fvfm.R](https://github.com/zlianglab/Triticale_BLS/blob/main/9.%20Chlorophyll%20fluoresence%20phenotype%20analysis/fvfm.R)
- [linear_model.R](https://github.com/zlianglab/Triticale_BLS/blob/main/9.%20Chlorophyll%20fluoresence%20phenotype%20analysis/linear_model.R)

### 10. Knockout Validation Analysis

**Folder:** [Knockout validation analysis](https://github.com/zlianglab/Triticale_BLS/tree/main/10.%20Knockout%20validation%20analysis)
**Description:** Analysis for DEG identification between P3 and P3TAL5KO
**Contents:**
- [DEGs-alluvial plot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/10.%20Knockout%20validation%20analysis/DEGs-alluvial%20plot.R)
- [boxplot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/10.%20Knockout%20validation%20analysis/boxplot.R)
- [degs.R](https://github.com/zlianglab/Triticale_BLS/blob/main/10.%20Knockout%20validation%20analysis/degs.R)
- [gene-boxplot.R](https://github.com/zlianglab/Triticale_BLS/blob/main/10.%20Knockout%20validation%20analysis/gene-boxplot.R)

## Contributor

- **Name:** Fahad Hasan
- **Department:** Department of Plant Sciences, Genomics, Phenomics & Bioinformaitcs Program
- **Institution:** North Dakota State Universtiy, Fargo, ND
- **Email:** mdfahad.hasan@ndsu.edu
- **ORCID:** [0009-0006-6209-5289](https://orcid.org/0009-0006-6209-5289)
- **GitHub:** [Fahad Hasan](https://github.com/hasanfahad)

## Contact

For questions about this project, contact **Zhikai Liang** at **[zhikai.liang@ndsu.edu]**.
