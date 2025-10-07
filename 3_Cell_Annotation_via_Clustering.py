import os
import numpy as np
import pandas as pd
import scanpy as sc
import harmonypy
import scanpy.external as sce
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


marker_genes_lymphocytes1 = {
"Tcells": ["CD3E", "CD3G", "CD3D", "TRAC", "TRBC1", "TRBC2", "CD2", "CD8A", "CD8B", "CD4"],
"NKcells": ["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB" ],
"NKTcells":["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB","CD3E","CD4", "CD8A", "CD8B"  ]
}



marker_genes_lympho = {
"T cell":["CD3E"]   ,
"TCD4+" : [ "CD4", "CD40LG", "IL7R", "CCR7"],# MAILFTH1
"TCD8+" : [ "CD8A", "CD8B"],
"Treg" : ["FOXP3","IL2RA", "CTLA4", "TNFRSF18"],
"NK" : ["FCGR3A","NKG7", "TYROBP", "NCAM1"]
}




"""
###################################################################################################################################################
-----------------------------------------------------------Functions-for-Subclusters--------------------------------------------------------------
###################################################################################################################################################
"""

def get_gene_expression_for_cell_type(adata, cell_type_by_cluster, sample_type):

    cell_barcodes = adata.obs["cell_type_by_cluster"] == cell_type_by_cluster
     
    counts_matrix = adata[cell_barcodes].layers["counts_normalized"]
    counts_df = pd.DataFrame(counts_matrix.toarray(), index=adata.obs_names[cell_barcodes], columns=adata.var_names)
    
    df_Gal_tim = counts_df[["HAVCR2", "LGALS9"]]
    df_Gal_tim["Sample"] = adata.obs["batch"]
   
    mean_expression = df_Gal_tim.groupby("Sample")[["HAVCR2", "LGALS9"]].mean().reset_index()
    
    mean_expression["Sample_type"] = sample_type
    mean_expression["Cell_type"] = cell_type_by_cluster
    
    mean_expression_melted = mean_expression.melt(id_vars=["Sample", "Sample_type", "Cell_type"], var_name="Gene", value_name="Mean Expression")
    
    return mean_expression_melted, cell_barcodes


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




"""
#######################################################################################################################
------------------------------------------Annotating cells------------------------------------------------------------ 
#######################################################################################################################
"""

#right
path_save_adata_cnv = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scRNA_GC_David_Tim3_Gal9/AnnData_objects/adata_all_cnv.h5ad"
adata_all = sc.read_h5ad(path_save_adata_cnv)
#wrong
path_save_adata_combined = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scRNA_GC_David_Tim3_Gal9/AnnData_objects/adata_combined.h5ad"
adata_all = sc.read_h5ad(path_save_adata_combined)


"""
#######################################################################################################################
------------------------------------------1. Lymphocytes in Tumor and Paratumor Tissue---------------------------------
#######################################################################################################################
"""


adata_tumor_tissue_nonTumorCells = adata_all[                                 
                                 (adata_all.obs.batch == "GC01T")   |
                                 (adata_all.obs.batch == "GC02T")   |
                                 (adata_all.obs.batch == "GC03T1")  |
                                 (adata_all.obs.batch == "GC03T2")  |
                                 (adata_all.obs.batch == "GC04T")   |
                                 (adata_all.obs.batch == "GC05T")   |
                                 (adata_all.obs.batch == "GC06T")   |
                                 (adata_all.obs.batch ==  "GC07T")  |
                                 (adata_all.obs.batch ==  "GC08T1") |
                                 (adata_all.obs.batch ==  "GC08T2") |
                                 (adata_all.obs.batch ==  "GC10T")                                   
                                 ]

adata_tumor_tissue_nonTumorCells = adata_tumor_tissue_nonTumorCells[adata_tumor_tissue_nonTumorCells.obs["Tumor_cells"] == "no"]


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/results/lvl2/TumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

cell_type = "Lymphocytes"
sample_type ="TumorTissue"
cell_type_cluster = ["0","1","2","3","4","6","7","13"]

adata_cluster_lvl2 = get_subcluster(adata_tumor_tissue_nonTumorCells, cell_type_cluster, path_save, sample_type=sample_type, 
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


                                                                                                                        #res1.5
#Look at general expreseeion markers
fig = sc.pl.umap(adata_cluster_lvl2, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type_plot = "Tcells"
fig = sc.pl.umap(adata_cluster_lvl2, color=marker_genes_lymphocytes1[cell_type_plot], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type_plot} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type_plot}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type_plot = "NKTcells"
fig = sc.pl.umap(adata_cluster_lvl2, color=marker_genes_lymphocytes1[cell_type_plot], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type_plot} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type_plot}_{sample_type}.png", dpi=300, bbox_inches='tight')


#Differential gene expression between markers 
sc.tl.rank_genes_groups(adata_cluster_lvl2, groupby="cluster_lvl2", method="wilcoxon")

fig = sc.pl.rank_genes_groups(adata_cluster_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(adata_cluster_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(adata_cluster_lvl2, marker_genes_lympho, groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["Bcells"], wspace=0.2,size=3, ncols=4)

"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""




clusters_to_remove = ["11","8"]
cluster_retain = list(set(adata_cluster_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))

path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/results/lvl3/cnv//TumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

"""
######################################################################################################################################
------------------------------------------------------Subcluster-Level-3--------------------------------------------------------------
######################################################################################################################################
"""


#Remove GC02T sample because it was define as uncertain if tumor or non-tumor
#adata_cluster_lvl2 = adata_cluster_lvl2[adata_cluster_lvl2.obs.batch != "GC02T"]


adata_cluster_lvl3 = get_subcluster(adata_cluster_lvl2, cluster_retain, path_save, sample_type=sample_type, 
                                           cell_type=cell_type, n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1,
                                           n_top_genes=2000,theta=0.5, cluster_lvl_parent="cluster_lvl2", cluster_lvl_child="cluster_lvl3", harmony=True )



sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Epithelial Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["T & NK Cells"])                                                   
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Erythrocytes"])                                                   
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Endothelial Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Fibroblasts"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Myeloid Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Fibroblasts"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Endocrine Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Myeloid Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Smooth Muscle Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["B Cells"]) 
sc.pl.umap(adata_cluster_lvl3, color=cell_type_marker_lvl1["Plasma Cells"]) 


fig = sc.pl.umap(adata_cluster_lvl3, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type_plot = "Tcells"
fig = sc.pl.umap(adata_cluster_lvl3, color=marker_genes_lymphocytes1[cell_type_plot], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type_plot} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type_plot}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type_plot = "NKTcells"
fig = sc.pl.umap(adata_cluster_lvl3, color=marker_genes_lymphocytes1[cell_type_plot], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type_plot} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type_plot}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type_plot = "TCD8_CD4_NK"
fig = sc.pl.dotplot(adata_cluster_lvl3, marker_genes_lympho, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type_plot} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type_plot}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(adata_cluster_lvl3, marker_genes_lympho, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')



cell_type_clusters_lvl3 ={ 
    "TCD8+" : ["","","","","","",""],
    "TCD4+" : ["","","", ""] ,  
    "Treg" : ["", ""],
    "NK": [""]     
    }



annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(adata_cluster_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

adata_cluster_lvl3 = adata_cluster_lvl3[adata_cluster_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]


if "cell_type_by_cluster" not in adata_cluster_lvl3.obs.columns:
    adata_cluster_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    adata_cluster_lvl3.obs.loc[adata_cluster_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type

cell_type = "TcellsNKcells"
fig = sc.pl.dotplot(adata_cluster_lvl3, marker_genes_lympho, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers for TILs", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(adata_cluster_lvl3, marker_genes_lympho, groupby="cell_type_by_cluster",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Expression of Cell Markers for TILs", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


"""
Get GAL9 and Tim3 Expression from lvl3 clusters
##########################################################################################################
""" 

list_gene_readout_tumor_tissue = []
list_barcodes_tumor_tissue = []
sample_type="Tumor_Tissue"
for cell_type in list(cell_type_clusters_lvl3.keys()):
     print(cell_type)
     gene_expression, barcodes = get_gene_expression_for_cell_type(adata=adata_cluster_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
     list_gene_readout_tumor_tissue.append(gene_expression)
     
     df_temp = pd.DataFrame({
         "barcode": barcodes,
         "cell_type": cell_type,  # this adds the cell type for every barcode
         "sample_type": sample_type}) 
     list_barcodes_tumor_tissue.append(df_temp)
     
df_gene_tumor_lymph = pd.concat(list_gene_readout_tumor_tissue,ignore_index=True)
df_barcode_tumor_lymph = pd.concat(list_barcodes_tumor_tissue)
   
#Dotplot for genes
#

adata_subset_1 = adata_cluster_lvl3.copy()

adata_subset_1.obs["cell_type_by_cluster"] = adata_subset_1.obs["cell_type_by_cluster"].replace({
    "TCD8+": "Tumor-TCD8+",
    "TCD4+": "Tumor-TCD4+",
    "Treg": "Tumor-Treg",
    "NK": "Tumor-NK"})
