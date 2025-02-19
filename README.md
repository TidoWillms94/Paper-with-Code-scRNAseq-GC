# Paper-with-Code-scRNAseq-GC
Analyzing the expression level of LGALS9 (Galectin 9) and HAVCR2 (TIM3) in single-cell RNA sequencing data from gastric cancer patients (blood, paratumor, and tumor tissue).  
# Data avalibility
Single-cell sequencing data (.fastq files) were accessed through the Genome Sequencing Archive from the China National Center for Bioinformation (CNCB) under the accession number [HRA000704](https://ngdc.cncb.ac.cn/gsa-human/browse/HRA000704).  
The sequencing data were made public through the publication by K. Sun et al. in Nat. Commun 13, 4943 (2022):  
[scRNA-seq of gastric tumor shows complex intercellular interaction with an alternative T cell exhaustion trajectory]( https://doi.org/10.1038/s41467-022-32627-z) Information for sample collection, processing, and library preparation can be accessed through the linked paper in the Methods section.
Reference genome GRCh38-2024A
VDJ-GRCh38-alts-ensembl-7.1.0

# Method
## Single cell RNA-seq data processing
The FASTQ reads from the gene expression (GEX) and V(D)J sequences were aligned to the respective reference genomes (GRCh38-2024A, VDJ-GRCh38-alts-ensembl-7.1.0) using Cell Ranger 9.0.0 and the corresponding alignment pipelines for counting and V(D)J analysis. Further preprocessing and analyse were perfomed using [scanpy](https://doi.org/10.1186/s13059-017-1382-0) 1.10.4 which included quality controll steps by filtering for cells with low UMI count (>400 UMIs), low gene number (>200 genes) and high mirochondrial transcript count (<30%). Moreover, transcripts which were detected in at least  less tahn 3 cells were neglected. Gene expression counts were normalized for divergent sequencing depth per cell to the median of the total transcript count per batch. Highly expressed genes which are accounting for more than 10% of the transcript count per cell were excluded. Gene expression counts were log2(x+1) tranformed and concatinated into an AnnData object for downstream analyse. [scDblFinder](10.12688/f1000research.73600.2) 1.16.0 was used to identify doublet  During inital clustering of all cell types from gastric cancer tissue and blood samples Doublets were detected and removed via [scDblFinder](10.12688/f1000research.73600.2) 1.16.0


doublet removel. transcript detected in at least 3 cells 
