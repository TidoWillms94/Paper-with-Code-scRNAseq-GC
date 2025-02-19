# Paper-with-Code-scRNAseq-GC
Analyzing the expression level of LGALS9 (Galectin 9) and HAVCR2 (TIM3) in single-cell RNA sequencing data from gastric cancer patients (blood, paratumor, and tumor tissue) across different cell types (Epithelial cells, Tumor cells, CD8+ T cells, CD4+ T cells, regulatory T cells, NK cells )
# Data avalibility
Single-cell sequencing data (.fastq files) were accessed through the Genome Sequencing Archive from the China National Center for Bioinformation (CNCB) under the accession number [HRA000704](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000704).  
The sequencing data were made public through the publication by K. Sun et al. in Nat. Commun 13, 4943 (2022):  
[scRNA-seq of gastric tumor shows complex intercellular interaction with an alternative T cell exhaustion trajectory]( https://doi.org/10.1038/s41467-022-32627-z) Information for sample collection, processing, and library preparation can be accessed through the linked paper in the Methods section.
Reference genome GRCh38-2024A
VDJ-GRCh38-alts-ensembl-7.1.0

# Method
## Single cell RNA-seq data processing
The FASTQ reads from the gene expression (GEX) and V(D)J sequences were aligned to the respective reference genomes (GRCh38-2024A, VDJ-GRCh38-alts-ensembl-7.1.0) using Cell Ranger 9.0.0 and the corresponding alignment pipelines for counting and V(D)J analysis. Further preprocessing and analysis were performed using [Scanpy](https://doi.org/10.1186/s13059-017-1382-0) 1.10.4, which included quality control steps by filtering for cells with a low UMI count (<400 UMIs), a low gene number (<200 genes), and a high mitochondrial transcript count (>30%). Moreover, transcripts that were detected in fewer than 3 cells were excluded. Gene expression counts were normalized to the median of the total transcript count per cell in each batch to correct for differences in sequencing depth per cell. Highly expressed genes, accounting for more than 10% of the transcript count per cell, were excluded. Gene expression counts were log₂(x+1) transformed and concatenated into an AnnData object for downstream analysis. First, doublets were removed using [scDblFinder](10.12688/f1000research.73600.2) 1.16.0, and additional doublets were identified during iterative subclustering by detecting contradictory cell marker expression. Further doublet removal was performed by utilizing the V(D)J sequence to exclude T cells with more than 2 copies of TRA and TRB subunits. Cell types of interest were identified by subclustering and analyzing the expression of specific marker genes (Supplement information with marker genes). Tumor cells were identified through clustering, copy number variation analysis, and comparison with the results from [K. Sun et al. 2022](https://doi.org/10.1038/s41467-022-32627-z).  
Finally, after annotating the cell clusters as epithelial cells, tumor cells, CD8+ T cells, CD4+ T cells, regulatory T cells, and NK cells, the expression of LGALS9 (Galectin 9) and HAVCR2 (TIM3) was retrieved for each cell type.

