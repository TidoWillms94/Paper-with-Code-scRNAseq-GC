import os
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
import scanpy.external as sce
import itertools
import infercnvpy as cnv
import matplotlib.pyplot as plt 
import seaborn as sns
from tqdm import tqdm


sc.settings.set_figure_params( frameon=True,color_map='Spectral_r') 
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False



cell_type_marker_lvl1 = {
    "B Cells": ["MS4A1", "CD19", "BANK1", "VPREB3", "CD79A"],
    "Plasma Cells": ["JCHAIN", "IGKC", "MZB1", "DERL3", "CD79A"],
    "Endocrine Cells": ["CHGA", "PCSK1N", "TTR", "SCG3", "SCG5", "EPCAM", "KRT18", "KRT19"],
    "Epithelial Cells": ["MUC5AC", "MUC1", "S100P", "LIPF", "TFF1", "TFF2", "PSCA", "EPCAM", "KRT18", "KRT19"],
    "T & NK Cells": ["CD2", "CD3E", "CD3D", "CD3G", "CD7"],
    "Erythrocytes": ["HBB", "HBA1", "ALAS2", "HBA2", "CA1"],
    "Mast Cells": ["TPSAB1", "TPSB2", "CPA3", "TPSD1", "GATA2"],
    "Myeloid Cells": ["AIF1", "CD68", "CD14", "FCN1", "S100A9", "MS4A7"],
    "Endothelial Cells": ["PLVAP", "VWF", "PECAM1", "ACKR1", "CLDN5", "CD34"],
    "Fibroblasts": ["PDGFRA", "PDPN", "DCN", "DPT", "TRPA1"],
    "Smooth Muscle Cells": ["ACTA2", "ACTG2", "MYH11", "RGS5", "NDUFA4L2"]
}



marker_genes_epithilial_cells = {
    "Chief cell" : ["PGA4", "PGA3", "LIPF"],
    "GMC": ["MUC6", "FUT9"],
    "PMC" : ["MUC5AC", "TFF1", "TFF2", "GKN1"],
    "Parietal cell" :  ["ATP4A", "ATP4B"],
    "Goblet cell" : ["MUC2", "ATOH1", "TFF3", "SPINK4", "CLCA1", "FCGBP"],
    "Enterocyte" : ["FABP1", "VIL1", "CDX1", "CDX2", "REG4", "KRT20"]
       } 


def get_Tcell_doublets_and_subtypes(VDJ_samples):
    list_doublets_t_cells = []
    list_t_ab = []
    list_t_gd = []
    for sample in list(VDJ_samples.keys()):
        #First filtering for doublets 
        df = pd.read_csv(VDJ_samples[str(sample)])
        df1 = pd.DataFrame(data=df.groupby(['barcode','chain']).size(), columns=["chain_count"])
        df_filtered = df1.loc[df1['chain_count'] <= 2]
        df_doublets = df1.loc[df1['chain_count'] > 2]
        barcode_doublets = df_doublets.index.get_level_values('barcode')
        barcode_doublets = barcode_doublets +"-"+ sample
        list_doublets_t_cells.append(barcode_doublets)
        #get alpha/beta and delta/gamma cells
        
        df_alpha_beta = df_filtered.loc[df_filtered.index.get_level_values('chain').isin(['TRA', 'TRB'])]
        df_gamma_delta = df_filtered.loc[df_filtered.index.get_level_values('chain').isin(['TRG', 'TRD'])]
        
        barcode_ab = df_alpha_beta.index.get_level_values('barcode')
        barcode_ab = barcode_ab +"-"+ sample
        list_t_ab.append(barcode_ab)
        
        barcode_gd = df_gamma_delta.index.get_level_values('barcode')
        barcode_gd = barcode_gd +"-"+ sample
        list_t_gd.append(barcode_gd)
    
    t_cell_doublets = list(itertools.chain.from_iterable(list_doublets_t_cells))
    t_cells_ab = list(itertools.chain.from_iterable(list_t_ab))
    t_cells_gd = list(itertools.chain.from_iterable(list_t_gd))

    return t_cell_doublets, t_cells_ab, t_cells_gd




"""
#########################################################################################################
Workflow scRNAseq celltype annotation

1.:  Read AnnData Object from 1_scRNAseq_preprocessing
2.:  Remove Doublets from VDJ-seqeuencing
3.:  Tumor Score Calculation
4.:  Copy Numver Variation Analyse 
5.:  
6.:
7.:
8.:
9.: 
    
#########################################################################################################
"""





"""
#########################################################################################################
---------------------------1----Read-AnnData-Onject------------------------------------------------------
#########################################################################################################
"""

path_save_adata_combined = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scRNA_GC_David_Tim3_Gal9/AnnData_objects/adata_combined.h5ad"
adata_all = sc.read_h5ad(path_save_adata_combined)
adata_all = adata_all[adata_all.obs.scDblFinder_class ==1]  #Remove doublets 



"""
#########################################################################################################
---------------------------2----Remove-T-Cells-Doublets-Identifyed-by-VDJ-seq----------------------------
#########################################################################################################
"""

VDJ_samples = {
    "GC01T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207251/outs/all_contig_annotations.csv",
    "GC02T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207252/outs/all_contig_annotations.csv",
    "GC04P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207254/outs/all_contig_annotations.csv",
    "GC04T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207255/outs/all_contig_annotations.csv",
    "GC05P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207256/outs/all_contig_annotations.csv",
    "GC05T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207257/outs/all_contig_annotations.csv",
    "GC06B" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207258/outs/all_contig_annotations.csv",
    "GC06P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207259/outs/all_contig_annotations.csv",
    "GC06T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207260/outs/all_contig_annotations.csv",
    "GC07B" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207261/outs/all_contig_annotations.csv",
    #"GC07P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207262/outs/all_contig_annotations.csv",  #No alignment posssible
    "GC07T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207263/outs/all_contig_annotations.csv",
    "GC08B" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207264/outs/all_contig_annotations.csv",
    "GC08P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207265/outs/all_contig_annotations.csv",
    "GC08T1" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207266/outs/all_contig_annotations.csv",
    "GC08T2" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207267/outs/all_contig_annotations.csv",
    "GC09P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207268/outs/all_contig_annotations.csv",
    "GC09T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207269/outs/all_contig_annotations.csv",
    "GC10B" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207270/outs/all_contig_annotations.csv",
    "GC10P1" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207271/outs/all_contig_annotations.csv",
    "GC10P2" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207272/outs/all_contig_annotations.csv",
    "GC10T" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207273/outs/all_contig_annotations.csv",
    }


t_cell_doublets, t_cells_ab, t_cells_gd = get_Tcell_doublets_and_subtypes(VDJ_samples)


adata_all =  adata_all[~adata_all.obs.index.isin(t_cell_doublets)]




"""
#########################################################################################################
---------------------------3----Tumor-Score-Calculation------------------------------------------------
#########################################################################################################
"""


malignant_genes = [
  "REG4", "TFF3", "PHGR1", "FABP1", "PRAP1", "APOA1", "CLDN4", "LEFTY1", 
  "SPINK4", "ANPEP", "CLDN7", "ALDOB", "S100A10", "IL32", "MALAT1", "LGALS4", 
  "CLDN3", "MDK", "OLFM4", "PRSS3", "LGALS3", "DMBT1", "CES2", "HLA-A", 
  "CDHR5", "FABP2", "CCT2", "PPP1R1B", "AOC1", "TSPAN7", "KRT7", "NFKBIA", 
  "ISG15", "CCL20", "TCIM", "FTL", "AGR3", "CXCL1", "RBP2", "TMC5", "PCK1", 
  "TRIM31", "CDH17", "PI3", "HES1", "LCN2", "TMEM17", "TM4SF4", "SERINC2", 
  "MISP"
]

nonmalignant_genes = [
  "REP15", "TM9SF3", "SLC12A2", "IL33", "RPS4Y1", "BRI3", "GHRL", "S100P", 
  "CXCL17", "PAPSS1", "FOXQ1", "FABP5", "CD63", "CLDN18", "MUC1", "KLF2", 
  "CTSE", "SOSTDC1", "EEF1A1", "ANXA10", "KCNE2", "TAGLN2", "TSPO", "VSIG1", 
  "AGR2", "METTL7B", "ALDH1A1", "VSIG2", "PSCA", "ALDH3A1", "PGA4", "TMSB4X", 
  "MUCL3", "SPINK1", "IGFBP2", "CA9", "GLUL", "PGA3", "TFF1", "MSMB", "CA2", 
  "CYSTM1", "RNASE1", "MUC6", "MUC5AC", "TFF2", "GKN2", "GKN1", "PGC", "LIPF"
]

#no layer option in scanpy in older versions
#sc.tl.score_genes(adata_combined, gene_list=malignant_genes, score_name="MalignantScore", random_state=1, layer='counts_normalized')
#sc.tl.score_genes(adata_combined, gene_list=nonmalignant_genes, score_name="NonMalignantScore", random_state=1, layers='counts_normalized')
#adata_combined.obs["TumorScore"] = adata_combined.obs["MalignantScore"] - adata_combined.obs["NonMalignantScore"]

counts_df = pd.DataFrame(adata_all.layers["counts_normalized"].toarray(), 
                         index=adata_all.obs_names, 
                         columns=adata_all.var_names)


adata_all.obs["MalignantScore"] = counts_df[malignant_genes].mean(axis=1)
adata_all.obs["NonMalignantScore"] = counts_df[nonmalignant_genes].mean(axis=1)
adata_all.obs["TumorScore"] = (adata_all.obs["MalignantScore"] - adata_all.obs["NonMalignantScore"])

sc.pl.umap(adata_all, color=['cluster_lvl1'],size=1, legend_loc='on data', legend_fontsize=10)

sc.pl.umap(adata_all, color=["TumorScore"])




"""
#########################################################################################################
---------------------------4----Copy-Number-Variation----------------------------------------------------
#########################################################################################################
"""
"""
4.1 Get Cell Annotation for potential tumro cells and non-tumor cell (here potetnial Tumor cells and Epithelial cells vs Stroma Cells and Lymphocytes)
"""

def get_potential_tumor_non_tumor_cells(adata, cluster_non_tumor, cluster_tumor):

    barcodes_potential_tumor_cells = adata[adata.obs['cluster_lvl1'].astype(str).isin(potential_tumor_cells), :].obs.index
    barcodes_non_tumor_cells = adata[adata.obs['cluster_lvl1'].astype(str).isin(non_tumor_cells), :].obs.index
    
    barcodes = [s.rsplit("-", 1)[0] for s in barcodes_potential_tumor_cells]
    batch = [s.rsplit("-", 1)[1] for s in barcodes_potential_tumor_cells]
    pot_tumor = pd.DataFrame(list(zip(barcodes, batch)), columns=["CellName", "Batch"])
    pot_tumor["Group"] = "tumor"
    
    barcodes = [s.rsplit("-", 1)[0] for s in barcodes_non_tumor_cells]
    batch = [s.rsplit("-", 1)[1] for s in barcodes_non_tumor_cells]
    non_tumor = pd.DataFrame(list(zip(barcodes, batch)), columns=["CellName", "Batch"])
    non_tumor["Group"] = "normal"
    
    df_annotation_cnv = pd.concat([pot_tumor, non_tumor], axis=0).reset_index(drop=True)
    
    return df_annotation_cnv

"""
######################################################################################################################################
-------------------------------------run Infer CNV with infercnvpy----------------------------------------------------------------
######################################################################################################################################
"""

potential_tumor_cells   = ["5","8","14"]
non_tumor_cells         = ["0","1","2","3","4","6","7","9","10","11","12","13","15","16","17","18","19","20","21"]

df_annotation_cnv =  get_potential_tumor_non_tumor_cells(adata=adata_all, cluster_non_tumor=non_tumor_cells , cluster_tumor=potential_tumor_cells)


"""
4.2.Get Alignment Dataframe for gene position 
"""

def parse_attributes(attr_str):
    attr_dict = {}
    pairs = attr_str.split("; ")
    for pair in pairs:
        if " " in pair:
            key, value = pair.split(" ", 1)
            attr_dict[key] = value.strip('"')  # Remove quotes if present
    return attr_dict


file_gtf = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/references_transcriptome/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"

df_gft = pd.read_csv(file_gtf, sep="\t", comment="#", header=None, compression="gzip",
                 names=["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"])


df = df_gft[df_gft["feature"] == "gene"]

chr_list = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
       'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
       'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chr22']

df = df[df["chromosome"].isin(chr_list)]


# Function to parse attributes
parsed_df = df["attribute"].apply(parse_attributes).apply(pd.Series)

df = pd.merge(df,parsed_df, left_index=True, right_index=True )
df = df[df["gene_type"]== "protein_coding"]
df_alignment = df[["gene_name", "chromosome", "start", "end"]].set_index("gene_name", drop=True) #18479
 

#add alignment file to anndata object
common_genes = adata_all.var.index.intersection(df_alignment.index)    
df_alignment = df_alignment.loc[common_genes]
df_alignment = df_alignment[~df_alignment.index.duplicated(keep='first')] 
                            
adata_all.var = adata_all.var.merge(df_alignment, how="left", left_index=True, right_index=True) #31837


df_annotation_cnv["Barcodes"] = df_annotation_cnv["CellName"] + "-" + df_annotation_cnv["Batch"]
df_annotation_cnv = df_annotation_cnv.set_index("Barcodes", drop=True)


common_Barcodes = adata_all.obs.index.intersection(df_annotation_cnv.index)    
df_annotation_cnv = df_annotation_cnv.loc[common_Barcodes]
                                             
adata_all.obs = adata_all.obs.merge(df_annotation_cnv, how="left", left_index=True, right_index=True)


sc.pl.umap(adata_all, color="Group")


adata_tumor = adata_all[
                                (adata_all.obs.batch == "GC01T") |
                                (adata_all.obs.batch == "GC02T") |
                                (adata_all.obs.batch == "GC03T1") |
                                (adata_all.obs.batch == "GC03T2") |
                                (adata_all.obs.batch == "GC04T") |
                                (adata_all.obs.batch == "GC05T") |
                                (adata_all.obs.batch == "GC06T") |
                                (adata_all.obs.batch ==  "GC07T") |
                                (adata_all.obs.batch ==  "GC08T1") |
                                (adata_all.obs.batch ==  "GC08T2") |
                                (adata_all.obs.batch ==  "GC10T") 
                                ].copy()




"""
4.3 run INFERCNV
"""

cnv.tl.infercnv(
    adata_tumor,
    reference_key="Group",
    reference_cat="normal",
    window_size=250
)

"""
INFERCNV VISUALISATION
"""
cnv.pl.chromosome_heatmap(adata_tumor, groupby="Group")


#Clustering
cnv.tl.pca(adata_tumor)
cnv.pp.neighbors(adata_tumor)
cnv.tl.leiden(adata_tumor)

cnv.pl.chromosome_heatmap(adata_tumor, groupby="cnv_leiden", dendrogram=True)


cnv.tl.umap(adata_tumor)
cnv.tl.cnv_score(adata_tumor)

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(11, 11))
ax4.axis("off")
cnv.pl.umap(
    adata_tumor,
    color="cnv_leiden",
    legend_loc="on data",
    legend_fontoutline=2,
    ax=ax1,
    show=False,
)
cnv.pl.umap(adata_tumor, color="cnv_score", ax=ax2, show=False)
cnv.pl.umap(adata_tumor, color="Group", ax=ax3)

cnv.pl.umap(adata_tumor, color="cnv_score")
cnv.pl.umap(adata_tumor, color="Group")




def get_subcluster(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes,theta, cluster_lvl_parent, cluster_lvl_child, harmony=True, ):

    barcodes_cell_types = adata[adata.obs[cluster_lvl_parent].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]

    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')
    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
    
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    
    #adjust PC's for Batch influnce   
    if harmony:
        #calculate community graph with embedding from X_pca_harmony PC's
        sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=theta, max_iter_harmony=15)
        sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors)
    else:
        #calculate community graph with embedding from X_pca PC's
        sc.pp.neighbors(cluster_adata, use_rep='X_pca', n_pcs=n_pcs, n_neighbors=n_neighbors)
    
   
    sc.tl.umap(cluster_adata, init_pos="spectral", min_dist=min_dist, spread=spread)  
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = cluster_lvl_child )
    
    fig = sc.pl.umap(cluster_adata, color=[cluster_lvl_child], size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}",show=False, return_fig=True)
   
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size
    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
    
    return cluster_adata







cell_type = "Epi_tumor"
sample_type = "Tumor Tissue"
cell_type_cluster = ["5","8","14"]

path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/results/lvl2/cnv/epi_tumor"
os.makedirs(path_save, exist_ok=True)

adata_cluster_lvl2 = get_subcluster(adata_tumor, cell_type_cluster, path_save, sample_type=sample_type, 
                                           cell_type=cell_type, n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4,
                                           n_top_genes=3000,theta=0.5, cluster_lvl_parent="cluster_lvl1", cluster_lvl_child="cluster_lvl2", harmony=True )


sc.pl.umap(adata_cluster_lvl2, color="cluster_lvl2", size=3, legend_loc='on data', legend_fontsize=10) 
sc.pl.umap(adata_cluster_lvl2, color="batch", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl2, color="cnv_score", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl2, color="percent_mito", palette=sns.color_palette() ) 


sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Epithelial Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["T & NK Cells"])                                                   
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Erythrocytes"])                                                   
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Endothelial Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Fibroblasts"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Myeloid Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Fibroblasts"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Endocrine Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Myeloid Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Smooth Muscle Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["B Cells"]) 
sc.pl.umap(adata_cluster_lvl2, color=cell_type_marker_lvl1["Plasma Cells"]) 

clusters_to_remove = ["12"]
cluster_retain = list(set(adata_cluster_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))



path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/results/lvl3/cnv/epi_tumor"
os.makedirs(path_save, exist_ok=True)


#Remove GC02T sample because it was define as uncertain if tumor or non-tumor
adata_cluster_lvl2 = adata_cluster_lvl2[adata_cluster_lvl2.obs.batch != "GC02T"]



adata_cluster_lvl3 = get_subcluster(adata_cluster_lvl2, cluster_retain, path_save, sample_type=sample_type, 
                                           cell_type=cell_type, n_pcs=15, n_neighbors=20, min_dist=0.5, spread=2, resolution=0.3,
                                           n_top_genes=1000,theta=0.5, cluster_lvl_parent="cluster_lvl2", cluster_lvl_child="cluster_lvl3", harmony=False )



sc.pl.umap(adata_cluster_lvl3, color="cluster_lvl3", size=3, legend_loc='on data', legend_fontsize=10) 
sc.pl.umap(adata_cluster_lvl3, color="batch", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl3, color="cnv_score", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl3, color="percent_mito", palette=sns.color_palette() ) 


sc.tl.rank_genes_groups(adata_cluster_lvl3, groupby="cluster_lvl3", method="wilcoxon")

sc.pl.rank_genes_groups(adata_cluster_lvl3, n_genes=25, sharey=False)

sc.pl.rank_genes_groups_heatmap(adata_cluster_lvl3, n_genes=5,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True)

sc.pl.rank_genes_groups_dotplot(adata_cluster_lvl3, n_genes=4)

sc.pl.heatmap(adata_cluster_lvl3, marker_genes_epithilial_cells ,groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True, figsize=(11, 5))

sc.pl.dotplot(adata_cluster_lvl3, marker_genes_epithilial_cells, groupby='cluster_lvl3',
              dendrogram=True,standard_scale="group")



clusters_to_remove = ["10"]
cluster_retain = list(set(adata_cluster_lvl3.obs['cluster_lvl3'].cat.categories) - set(clusters_to_remove))

adata_cluster_lvl4 = get_subcluster(adata_cluster_lvl3, cluster_retain, path_save, sample_type=sample_type, 
                                           cell_type=cell_type, n_pcs=15, n_neighbors=20, min_dist=0.5, spread=1, resolution=0.3,
                                           n_top_genes=1000,theta=0.5, cluster_lvl_parent="cluster_lvl3", cluster_lvl_child="cluster_lvl4", harmony=False )



sc.tl.rank_genes_groups(adata_cluster_lvl4, groupby="cluster_lvl4", method="wilcoxon")

sc.pl.rank_genes_groups(adata_cluster_lvl4, n_genes=25, sharey=False)

sc.pl.rank_genes_groups_heatmap(adata_cluster_lvl4, n_genes=5,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True)


sc.pl.rank_genes_groups_dotplot(adata_cluster_lvl4, n_genes=4, standard_scale="group")

sc.pl.heatmap(adata_cluster_lvl4, marker_genes_epithilial_cells ,groupby="cluster_lvl4",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True, figsize=(11, 5))

sc.pl.dotplot(adata_cluster_lvl4, marker_genes_epithilial_cells, groupby='cluster_lvl4',
              dendrogram=True,standard_scale="group")


sc.pl.umap(adata_cluster_lvl4, color="cluster_lvl4", size=3, legend_loc='on data', legend_fontsize=10) 
sc.pl.umap(adata_cluster_lvl4, color="batch", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl4, color="cnv_score", palette=sns.color_palette() ) 
sc.pl.umap(adata_cluster_lvl4, color="percent_mito", palette=sns.color_palette() ) 


cell_type_clusters_lvl4 ={ 
    "Chief cells" : ["7"],
    "GMC" : ["2"] ,  
    "PMC" : ["0"],
    "Parietal cells": ["9"],
    "na" : ["6"] ,  
    "Enterocyte" : [""],
    "Tumor cells" : ["8","3","5","1","4"]
    #"GC02T":["11"]
    }



annotated_clusters = [item for sublist in cell_type_clusters_lvl4.values() for item in sublist]
adata_cluster_lvl4 = adata_cluster_lvl4[adata_cluster_lvl4.obs["cluster_lvl4"].isin(annotated_clusters)]


if "cell_type_by_cluster" not in adata_cluster_lvl4.obs.columns:
    adata_cluster_lvl4.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl4.items():
    adata_cluster_lvl4.obs.loc[adata_cluster_lvl4.obs["cluster_lvl4"].isin(clusters), "cell_type_by_cluster"] = cell_type



sc.pl.dotplot(adata_cluster_lvl4, marker_genes_epithilial_cells, groupby='cell_type_by_cluster', standard_scale="group",use_raw=False)

sc.pl.heatmap(adata_cluster_lvl4, marker_genes_epithilial_cells, groupby="cell_type_by_cluster", use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6))



tumor_clusters = adata_cluster_lvl4[adata_cluster_lvl4.obs["batch"].isin(["GC07T", "GC08T1", "GC08T2", "GC10T"])]
tumor_clusters.obs["cnv_status"] = "normal"
tumor_clusters.obs.loc[tumor_clusters.obs["cluster_lvl4"].isin(["8","3","5","1","4"]), "cnv_status"] = ("tumor")
tumor_cells = tumor_clusters.obs[tumor_clusters.obs["cnv_status"] == "tumor"].index


adata_all.obs["Tumor_cells"] = "no"
adata_all.obs.loc[adata_all.obs.index.isin(tumor_cells), "Tumor_cells"] = "yes"

sc.pl.umap(adata_all, color="Tumor_cells") 
"""
#########################################################################################################
---------------------------6----Add-CNV-Score-to-AnnData-object------------------------------------------
#########################################################################################################
"""

path_save_adata_cnv = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scRNA_GC_David_Tim3_Gal9/AnnData_objects/adata_all_cnv.h5ad"
adata_all.write(filename=path_save_adata_combined, compression=None, compression_opts=None, as_dense=())



