import os
import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import itertools
import matplotlib.pyplot as plt 
import seaborn as sns


sc.settings.set_figure_params( frameon=True,color_map='Spectral_r') 
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False






"""
###################################################################################################
--------------------------------------Load-Preprocessed-AnnData-Object-----------------------------
###################################################################################################
"""
/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/python/save_AnnData/
save_path_adata_scaled_umap = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/python/save_AnnData/GC_dublets_corrected_scaled_umap_from_prepocessing_all_parralel_4.h5ad"
adata_all = sc.read_h5ad(save_path_adata_scaled_umap)
adata_all = adata_all[adata_all.obs.scDblFinder_class ==1]

#adata_combined contains all GCxxT, GCxxP and GCxxB samples 


#save_path = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/python/save_AnnData/adata_CNV_samples_P_T.h5ad"
#adata_tissue = sc.read_h5ad(save_path)
#adata_tumor contains all GCxxT samples and GCxxP samples with CNV anaylse 





"""
###################################################################################################
---------------------------------------Cell-type-marker-genes--------------------------------------
###################################################################################################
"""







marker_genes_lvl1 = {
    "Bcells": ["MS4A1", "CD19", "BANK1", "VPREB3", "CD79A"],
    "PlasmaCells": ["JCHAIN", "IGKC", "MZB1", "DERL3", "CD79A"],
    "EndocrineCells": ["CHGA", "PCSK1N", "TTR", "SCG3", "SCG5", "EPCAM", "KRT18", "KRT19"],
    "EpithelialCells": ["MUC5AC", "MUC1", "S100P", "LIPF", "TFF1", "TFF2", "PSCA", "EPCAM", "KRT18", "KRT19"],
    "T&NKcells": ["CD2", "CD3E", "CD3D", "CD3G", "CD7", "FCGR3A", "TYROBP"],
    "Erythrocytes": ["HBB", "HBA1", "ALAS2", "HBA2", "CA1"],
    "MastCells": ["TPSAB1", "TPSB2", "CPA3", "TPSD1", "GATA2"],
    "MyeloidCells": ["AIF1", "CD68", "CD14", "FCN1", "S100A9", "MS4A7"],
    "EndothelialCells": ["PLVAP", "VWF", "PECAM1", "ACKR1", "CLDN5", "CD34"],
    "Fibroblasts": ["PDGFRA", "PDPN", "DCN", "DPT", "TRPA1"],
    "SmoothMuscleCells": ["ACTA2", "ACTG2", "MYH11", "RGS5", "NDUFA4L2"]   
}

marker_genes_lvl2_epithilial = {
    "chief cell" : ["PGA4", "PGA3", "LIPF"],
    "GMC": ["MUC6", "FUT9"],
    "PMC" : ["MUC5AC", "TFF1", "TFF2", "GKN1"],
    "parietal cell" :  ["ATP4A", "ATP4B", "GIF"],
    "goblet cell" : ["MUC2", "ATOH1", "TFF3", "SPINK4", "CLCA1", "FCGBP"],
    "enterocyte" : ["FABP1", "VIL1", "CDX1", "CDX2", "REG4", "KRT20"]
       } 

marker_genes_lvl2_lymphocytes = {
"Tcells": ["CD3E", "CD3G", "CD3D", "TRAC", "TRBC1", "TRBC2", "CD2", "CD8A", "CD8B", "CD4"],
"NKcells": ["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB" ],
"NKTcells":["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB","CD3E","CD4", "CD8A", "CD8B"  ]
}

marker_genes_lvl3_T_NKcells = {
    "HAVCR2": ["HAVCR2"],
    "LGALS9": ["LGALS9"], 
    "Tcells": ["CD3E", "CD3G", "CD3D", "TRAC", "TRBC1", "TRBC2", "CD2", "CD8A", "CD8B", "CD4"],
    "NKTcells":["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB","CD3E","CD4", "CD8A", "CD8B"  ]
    }



a = {
"Tcell": ["CD3E","CD4", "CD8A", "CD8B",                 #Tcell + subsets 
          "CCR7", "TCF7", "PTPRC", "SELL",              # Naive
          "CD69", "IL2RA","TNFRSF9","CD3E",             #eary activation CD69, CD25 (IL2RA), CD137 (4-1BB,TNFRSF9)
          "GZMB", "PRF1", "IFNG", "IL2",                #Cytotoxix t cells
          "FOXP3", "IL2RA",                             #Treg
          "TBX21","GATA3",                              #
          "CCR7", "TCF7", "PTPRC", "SELL",              #PTPRC CD4
          "IFNG", "GZMB","PRF1","CD69","TNF", "IL2",    #
          "IL7R", "FAS",                                #Memory CD45RO
          ],

"TcellCD8": ["CD3E","CD4", "CD8A", "CD8B",              #subtype marker
             "CCR7", "TCF7", "LEF1", "SELL",            #Naive
             "CD69", "IL2RA","TNFRSF9","CD3E",          #early activation
             "GZMB", "PRF1", "IFNG", "NKG7",            #Cytotoxix t cells
             "GNLY" , "TNF", "CD3E", "CD3E",            #Cytotoxix t cells
             "KLRG1" , "CX3CR1", "CD3E", "CD3E",        #Late / Terminally Differentiated CD8 T Cells
             "PDCD1" , "CTLA4", "HAVCR2", "LAG3",       #Exhausted CD8 T Cells
             "TIGIT", "TOX", "EOMES"                    #Exhausted CD8 T Cells
             ],

"TcellCD4": ["CD3E","CD4", "CD8A", "CD8B",          #
            "CCR7","SELL", "TCF7", "LEF1",          #naive t cells
            "CCR7","SELL", "IL7R", "CD8B",          #central memory t cells
            "CCR7","SELL", "IFNG", "CD8B",          #effector memory t cells
            "TBX21","IFNG", "CXCR3", "IL12RB2",     #T Helper 1 (Th1)
            "GATA3","IL4", "IL5", "IL13",           #T Helper 2 (Th2)
            "RORC","IL17A", "IL17F", "CCR6",       #T Helper 17 (Th17)   RORC (RORγt) — master transcription factor.
            "BCL6","CXCR5", "PDCD1", "ICOS",           #T Follicular Helper (Tfh)
            "FOXP3","IL2RA", "CTLA4", "TNFRSF18"            #Regulatory T (Treg)
             ],
"NKTCells": ["FCGR3A", "KLRG1", "NCAM1","TYROBP", "NKG7","GNLY","GZMB","CD3E","CD4", "CD8A", "CD8B"  ]
}

#Invariant TCR: (if specifically iNKT) might see TRAV10–TRAJ18 in humans.

marker_genes_lvl3_NKT_cells = {
"TIM3": ["HAVCR2"],
"Gal9": ["LGALS9"],  
"Tcells": ["CD3E", "CD3G", "CD3D", "TRAC", "TRBC1", "TRBC2", "CD2", "CD8A", "CD8B", "CD4"],
"NKcells": ["FCGR3A", "KLRG1", "TYROBP", "NKG7","GNLY","GZMB" ]}



"""
########################################################################################
---------------------------------------FUNCTIONS----------------------------------------
########################################################################################
"""
#or check on cluster lvl2 for heterotypic doublets 
def doublets_check_lvl1_marker_genes():
    for x in marker_genes:
        fig = sc.pl.umap(adata_all, color=marker_genes[x] ,size=3, legend_loc='on data', legend_fontsize=10,show=False, return_fig=True)
        plt.suptitle(f"UMAP Marker Genes for {x}", y=0.96, fontsize=20, weight='bold')  
        fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
        fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
        fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size
        fig.savefig(f"{path_save}/Umap_Marker_genes_{x}.png", dpi=300, bbox_inches='tight')
        
        

def cluster_lvl_2(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes):
    #Retreive non scaled data for new HVG and run clustering pipline
    barcodes_cell_types = adata[adata.obs['cluster_lvl1_res_1'].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
        
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    
    sc.pp.neighbors(cluster_adata , n_pcs=n_pcs, n_neighbors=n_neighbors) #math.isqrt(cluster_adata.n_obs))
    sc.tl.umap(cluster_adata , init_pos="spectral", min_dist=min_dist, spread=spread)  #
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl2' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl2'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}", show=False, return_fig=True)
   
  
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size

    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

    return cluster_adata


def cluster_general(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes, theta):
    #Retreive non scaled data for new HVG and run clustering pipline
    barcodes_cell_types = adata[adata.obs['cluster_lvl1_res_1'].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
        
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    
    sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=theta)

    sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors) #math.isqrt(cluster_adata.n_obs))
    
    sc.tl.umap(cluster_adata , init_pos="spectral", min_dist=min_dist, spread=spread)  #
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl2' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl2'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}", show=False, return_fig=True)
   
  
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size

    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

    return cluster_adata






def cluster_lvl_2_with_harmony(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes, theta, harmony):
    #Retreive non scaled data for new HVG and run clustering pipline
    barcodes_cell_types = adata[adata.obs['cluster_lvl1_res_1'].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
        
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    
    if harmony:
        sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=theta, max_iter_harmony=15)
        sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors)
    else:
        #calculate community graph with embedding from X_pca_harmony PC's
        sc.pp.neighbors(cluster_adata, use_rep='X_pca', n_pcs=n_pcs, n_neighbors=n_neighbors)
        
    sc.tl.umap(cluster_adata , init_pos="spectral", min_dist=min_dist, spread=spread)  #
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl2' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl2'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}", show=False, return_fig=True)
   
  
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size

    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

    return cluster_adata



def cluster_lvl_2_with_harmony_and_diffusion(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes):
    #Retreive non scaled data for new HVG and run clustering pipline
    barcodes_cell_types = adata[adata.obs['cluster_lvl1_res_1'].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
    #get linare depenend variables with highest variation    
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    #adjust PC's for Batch influnce    
    sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=1)
    #calculate community graph with embedding from X_pca_harmony PC's
    sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors) #math.isqrt(cluster_adata.n_obs))
    # diffmap uses community graph .uns[‘neighbors’] for neighbors settings and .obsp[‘connectivities’], .obsp[‘distances’] for connectivities and distances respectively 
    sc.tl.diffmap(cluster_adata, n_comps=15)
    #Get new community graph according to X_diffmap embedding
    sc.pp.neighbors(cluster_adata, use_rep='X_diffmap')
    #specify embedding and neighbor graph for umap
    sc.tl.umap(cluster_adata, init_pos="spectral", min_dist=min_dist, spread=spread)  #
    #Get Cluster
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl2' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl2'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}" ) #,show=False, return_fig=True
   
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size
    
    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
    
    return cluster_adata


def cluster_lvl_3_with_harmony_and_diffusion(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes,theta, cluster_lvl='cluster_lvl2' ):

    barcodes_cell_types = adata[adata.obs[cluster_lvl].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
    #get linare depenend variables with highest variation    
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    #adjust PC's for Batch influnce    
    sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=1)
    #calculate community graph with embedding from X_pca_harmony PC's
    sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors) #math.isqrt(cluster_adata.n_obs))
    # diffmap uses community graph .uns[‘neighbors’] for neighbors settings and .obsp[‘connectivities’], .obsp[‘distances’] for connectivities and distances respectively 
    sc.tl.diffmap(cluster_adata, n_comps=15)
    #Get new community graph according to X_diffmap embedding
    sc.pp.neighbors(cluster_adata, use_rep='X_diffmap')
    #specify embedding and neighbor graph for umap
    sc.tl.umap(cluster_adata, init_pos="spectral", min_dist=min_dist, spread=spread)  #
    #Get Cluster
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl3_diff' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl3_diff'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}",show=False, return_fig=True)
   
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size
    
    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
    
    return cluster_adata


def cluster_lvl_3_with_harmony(adata, cell_type_cluster, path_save, sample_type, cell_type, n_pcs, n_neighbors, min_dist, spread, resolution, n_top_genes,theta, harmony=True, cluster_lvl='cluster_lvl2'):

    barcodes_cell_types = adata[adata.obs[cluster_lvl].astype(str).isin(cell_type_cluster), :].obs.index
    cluster_adata = adata[adata.obs.index.isin(barcodes_cell_types), :]
    #sc.pp.highly_variable_genes(cluster_adata, layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat')  
    sc.pp.highly_variable_genes(cluster_adata,layer="counts_normalized", n_top_genes = n_top_genes, flavor ='seurat', batch_key='batch')

    cluster_adata = sc.pp.scale(cluster_adata, max_value=10, copy=True, zero_center=True) 
    #get linare depenend variables with highest variation    
    sc.tl.pca(cluster_adata ,svd_solver='arpack', use_highly_variable = True)
    sc.pl.pca_variance_ratio(cluster_adata)
    #adjust PC's for Batch influnce   
    if harmony:
        sce.pp.harmony_integrate(cluster_adata, 'batch',  theta=theta, max_iter_harmony=15)
        sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony', n_pcs=n_pcs, n_neighbors=n_neighbors)
    else:
        #calculate community graph with embedding from X_pca_harmony PC's
        sc.pp.neighbors(cluster_adata, use_rep='X_pca', n_pcs=n_pcs, n_neighbors=n_neighbors)
    
    # diffmap uses community graph .uns[‘neighbors’] for neighbors settings and .obsp[‘connectivities’], .obsp[‘distances’] for connectivities and distances respectively 
    #sc.tl.diffmap(cluster_adata, n_comps=15)
    #Get new community graph according to X_diffmap embedding
    #sc.pp.neighbors(cluster_adata, use_rep='X_pca_harmony')
    #specify embedding and neighbor graph for umap
    sc.tl.umap(cluster_adata, init_pos="spectral", min_dist=min_dist, spread=spread)  #
    #Get Cluster
    sc.tl.leiden(cluster_adata , resolution=resolution, key_added = 'cluster_lvl3' )
    fig = sc.pl.umap(cluster_adata, color=['cluster_lvl3'],size=3, legend_loc='on data', legend_fontsize=10,
                     title=f"UMAP Clustering for {cell_type} in {sample_type}",show=False, return_fig=True)
   
    fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
    fig.axes[0].set_xlabel("UMAP 1", fontsize=10)  # X-axis Label Font Size
    fig.axes[0].set_ylabel("UMAP 2", fontsize=10)  # Y-axis Label Font Size
    fig.axes[0].tick_params(axis='both', labelsize=10)  # Tick Font Size
    
    fig.savefig(f"{path_save}/Umap_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
    
    return cluster_adata






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

"""
########################################################################################
---------------------------------------FUNCTIONS-END------------------------------------
########################################################################################
"""


"""
########################################################################################
---------------------------------------Filter-for-T-cell-Doublets-----------------------
########################################################################################
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
    #"GC07P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207262/outs/all_contig_annotations.csv",
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


"""
remove t_cell_doublets from anndataobject
"""
#adata_all contains all GCxxT, GCxxP and GCxxB samples
adata_all =  adata_all[~adata_all.obs.index.isin(t_cell_doublets)]
#adata_tissue contains all GCxxT samples and GCxxP samples with CNV anaylse 
#adata_tissue =  adata_tissue[~adata_tissue.obs.index.isin(t_cell_doublets)]


"""
GET GENE EXPRESSION FROM TUMOR CELLS
########################################################################################

last code in inferCNVpy
tumor_cells = cluster_adata_lvl3.obs[cluster_adata_lvl3.obs["cnv_status"] == "tumor"].index
Barcodes of Tumorcells 
"""





adata_tumor_tissue = adata_all[
                                (adata_all.obs.batch ==  "GC07T")  |
                                (adata_all.obs.batch ==  "GC08T1") |
                                (adata_all.obs.batch ==  "GC08T2") |
                                (adata_all.obs.batch ==  "GC10T")  
                                ] 


common_barcodes = adata_tumor_tissue.obs_names.intersection(tumor_cells)

subset = adata_tumor_tissue[adata_tumor_tissue.obs.index.isin(common_barcodes)]
counts_matrix = subset.layers["counts_normalized"]


counts_df = pd.DataFrame(counts_matrix.toarray(), 
                         index=common_barcodes, 
                         columns=adata_tumor_tissue.var_names)


df_Gal_tim = counts_df[["HAVCR2", "LGALS9"]]
df_Gal_tim["Sample"] = adata_tumor_tissue.obs["batch"]

mean_expression = df_Gal_tim.groupby("Sample")[["HAVCR2", "LGALS9"]].mean().reset_index()
mean_expression["Sample_type"] = "Tumor_Tissue"
mean_expression["Cell_type"] = "Tumor_cells"
mean_expression_melted_tumor_cells = mean_expression.melt(id_vars=["Sample", "Sample_type", "Cell_type"], var_name="Gene", value_name="Mean Expression")

df_tumor_cells= pd.DataFrame(tumor_cells).reset_index()
df_tumor_cells.to_csv("/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/df_tumto_cells_barcode_readout.csv")


"""
########################################################################################################
--------------------------------------------LYMPHOCYTE CELLS IN TUMOR-TISSUE----------------------------
########################################################################################################

WORKFLOW
1. Define Set of genes for relevant cell types
2. select relavent clusters
3. 

"""


marker_genes = {
"T cell":["CD3E"]   ,
"TCD4+" : [ "CD4", "CD40LG", "IL7R", "CCR7"],# MAILFTH1
"TCD8+" : [ "CD8A", "CD8B"],
"Treg" : ["FOXP3","IL2RA", "CTLA4", "TNFRSF18"],
"NK" : ["FCGR3A","NKG7", "TYROBP", "NCAM1"]

}



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

adata_tumor_tissue_nonTumorCells = adata_tumor_tissue_nonTumorCells[~adata_tumor_tissue_nonTumorCells.obs_names.isin(tumor_cells)]


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl2/TumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

cell_type = "Lymphocytes"
sample_type ="TumorTissue"
cell_type_cluster = ["0","1","2","3","4","6","7","13"]

cluster_adata_lymphocytes_lvl2 = cluster_lvl_2_with_harmony(adata_tumor_tissue_nonTumorCells, cell_type_cluster, path_save, sample_type="Tumor Tissue", 
                                           cell_type="Lymphocytes", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4,
                                           n_top_genes=3000,theta=0.5, harmony=True )
                                                                                                                            #res1.5
#Look at general expreseeion markers
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#Differential gene expression between markers 
sc.tl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, groupby="cluster_lvl2", method="wilcoxon")

fig = sc.pl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(cluster_adata_lymphocytes_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(cluster_adata_lymphocytes_lvl2,marker_genes,groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["Bcells"], wspace=0.2,size=3, ncols=4)

"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""

clusters_to_remove = ["12","5"]
selected_lvl2_clusters = list(set(cluster_adata_lymphocytes_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl3/TumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

cluster_adata_lvl3 = cluster_lvl_3_with_harmony(cluster_adata_lymphocytes_lvl2, selected_lvl2_clusters, path_save, sample_type="Tumor Tissue", 
                                           cell_type="TCD8_CD4_NK", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1, n_top_genes=2000, theta=0.5,
                                           cluster_lvl='cluster_lvl2', harmony=True )

fig = sc.pl.umap(cluster_adata_lvl3, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "TCD8_CD4_NK"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')



cell_type_clusters_lvl3 ={ 
    "TCD8+" : ["3","4","5","1","7","11","9"],
    "TCD4+" : ["0","5","10", "6"] ,  
    "Treg" : ["2", "13"],
    "NK": ["8","12"]     
    }


annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]


if "cell_type_by_cluster" not in cluster_adata_lvl3.obs.columns:
    cluster_adata_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    cluster_adata_lvl3.obs.loc[cluster_adata_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type

cell_type = "TcellsNKcells"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers for TILs", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cell_type_by_cluster",use_raw=False,
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
     gene_expression, barcodes = get_gene_expression_for_cell_type(adata=cluster_adata_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
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

adata_subset_1 = cluster_adata_lvl3.copy()

adata_subset_1.obs["cell_type_by_cluster"] = adata_subset_1.obs["cell_type_by_cluster"].replace({
    "TCD8+": "Tumor-TCD8+",
    "TCD4+": "Tumor-TCD4+",
    "Treg": "Tumor-Treg",
    "NK": "Tumor-NK"})
    

    
    
"""
##########################################################################################################################################
-------------------------------------------LYMPHOCYTES in PERITUMOR TISSUE----------------------------------------------------------------
##########################################################################################################################################
"""


adata_peritumor = adata_all[
                                (adata_all.obs.batch == "GC03P") |
                                (adata_all.obs.batch == "GC04P") |
                                (adata_all.obs.batch == "GC05P") |
                                (adata_all.obs.batch == "GC06P") |
                                (adata_all.obs.batch == "GC07P") |
                                (adata_all.obs.batch == "GC08P") |
                                (adata_all.obs.batch == "GC09P") |
                                (adata_all.obs.batch ==  "GC10P1") |
                                (adata_all.obs.batch ==  "GC10P2") 
                                ]


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl2/ParatumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

cell_type = "Lymphocytes"
sample_type ="Paratumor"
cell_type_cluster = ["0","1","2","3","4","6","7","13"]
cluster_adata_lymphocytes_lvl2 = cluster_lvl_2_with_harmony(adata_peritumor, cell_type_cluster, path_save, sample_type="Paratumor Tissue", 
                                           cell_type="Lymphocytes", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4,
                                           n_top_genes=3000,theta=0.5, harmony=True  )
                                                                                                                            #res1.5
#Look at general expreseeion markers
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#Differential gene expression between markers 
sc.tl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, groupby="cluster_lvl2", method="wilcoxon")

fig = sc.pl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(cluster_adata_lymphocytes_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(cluster_adata_lymphocytes_lvl2,marker_genes,groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl1["EpithelialCells"], wspace=0.2,size=3, ncols=4)
"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""

clusters_to_remove = ["11","10"]
selected_lvl2_clusters = list(set(cluster_adata_lymphocytes_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl3/ParaumorTissue/Lympo"
os.makedirs(path_save, exist_ok=True)

cluster_adata_lvl3 = cluster_lvl_3_with_harmony(cluster_adata_lymphocytes_lvl2, selected_lvl2_clusters, path_save, sample_type="Paratumor Tissue", 
                                           cell_type="TCD8_CD4_NK", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1, n_top_genes=2000, theta=0.5,
                                           cluster_lvl='cluster_lvl2', harmony=True )


fig = sc.pl.umap(cluster_adata_lvl3, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "TCD8_CD4_NK"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["EpithelialCells"], size=3,)

cell_type_clusters_lvl3 ={ 
    "TCD8+" : ["0","5","1","12","4","10","13"],
    "TCD4+" : ["3","6","2","7","",""] ,  
    "Treg" : ["9"],
    "NK": ["11","17"]     
    }

annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]



if "cell_type_by_cluster" not in cluster_adata_lvl3.obs.columns:
    cluster_adata_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    cluster_adata_lvl3.obs.loc[cluster_adata_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type

cell_type = "TcellsNKcells"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers in {sample_type} Tissue", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cell_type_by_cluster",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Expression of Cell Markers {sample_type} Tissue", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')




"""
Get GAL9 and Tim3 Expression from lvl3 clusters
##########################################################################################################
"""
list_gene_readout_paratumor_tissue = []
list_barcodes_paratumor_tissue = []
sample_type="Paratumor_Tissue"
for cell_type in list(cell_type_clusters_lvl3.keys()):
    print(cell_type)
    gene_expression, barcodes = get_gene_expression_for_cell_type(adata=cluster_adata_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
    list_gene_readout_paratumor_tissue.append(gene_expression)
    
    df_temp = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": cell_type,  # this adds the cell type for every barcode
        "sample_type": sample_type}) 
    list_barcodes_paratumor_tissue.append(df_temp)
    
df_gene_paratumor_lymph = pd.concat(list_gene_readout_paratumor_tissue,ignore_index=True)
df_barcode_paratumor_lymph = pd.concat(list_barcodes_paratumor_tissue)


adata_subset_2 = cluster_adata_lvl3.copy()

adata_subset_2.obs["cell_type_by_cluster"] = adata_subset_2.obs["cell_type_by_cluster"].replace({
    "TCD8+": "Paratumor-TCD8+",
    "TCD4+": "Paratumor-TCD4+",
    "Treg": "Paratumor-Treg",
    "NK": "Paratumor-NK"})



"""
##########################################################################################################################################
------------------------------------- TUMOR TISSUE (Stroma/Epithilial)-----------------------------------------------------------
##########################################################################################################################################
"""



marker_genes_lvl2_epithilial = {
    "Chief cell" : ["PGA4", "PGA3", "LIPF"],
    "GMC": ["MUC6", "FUT9"],
    "PMC" : ["MUC5AC", "TFF1", "TFF2", "GKN1"],
    "Parietal cell" :  ["ATP4A", "ATP4B"],
    "Goblet cell" : ["MUC2", "ATOH1", "TFF3", "SPINK4", "CLCA1", "FCGBP"],
    "Enterocyte" : ["FABP1", "VIL1", "CDX1", "CDX2", "REG4", "KRT20"]
       } 


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

adata_tumor_tissue_nonTumorCells = adata_tumor_tissue_nonTumorCells[~adata_tumor_tissue_nonTumorCells.obs_names.isin(tumor_cells)]


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl2/TumorTissue/Stroma"
os.makedirs(path_save, exist_ok=True)

cell_type = "EpithelialCells"
sample_type = "TumorTissue"
cell_type_cluster = ["5","8","14"]
cluster_adata_lymphocytes_lvl2 = cluster_lvl_2_with_harmony(adata_tumor_tissue_nonTumorCells, cell_type_cluster, path_save, sample_type="Tumor Tissue", 
                                           cell_type=cell_type, n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4,
                                           n_top_genes=3000,theta=0.5, harmony=False )


                                                                                                                            #res1.5
#Look at general expreseeion markers
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=['batch'],size=3,show=False, return_fig=True, palette=sns.color_palette())
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


cell_type = "Epithelial"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl1["EpithelialCells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl1["T&NKcells"], wspace=0.2,size=3, ncols=4)
sc.pl.umap(cluster_adata_lymphocytes_lvl2, color="cluster_lvl2", wspace=0.2,size=3, ncols=4)


#Differential gene expression between markers 
sc.tl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, groupby="cluster_lvl2", method="wilcoxon")

cell_type = "Epithelial"
fig = sc.pl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(cluster_adata_lymphocytes_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 24), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(cluster_adata_lymphocytes_lvl2,marker_genes_lvl2_epithilial,groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""

clusters_to_remove = ["5","15","12"]
selected_lvl2_clusters = list(set(cluster_adata_lymphocytes_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl3/TumorTissue/stroma"
os.makedirs(path_save, exist_ok=True)

cluster_adata_lvl3 = cluster_lvl_3_with_harmony(cluster_adata_lymphocytes_lvl2, selected_lvl2_clusters, path_save, sample_type="Paratumor Tissue", 
                                           cell_type="Epithelial", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1, n_top_genes=2000, theta=0.5,
                                           cluster_lvl='cluster_lvl2', harmony=True )


fig = sc.pl.umap(cluster_adata_lvl3, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


cell_type = "Epithelial"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["EpithelialCells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["T&NKcells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Epithelial"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

sc.pl.umap(cluster_adata_lvl3, color="cnv_score", wspace=0.2,size=3, ncols=4)
sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["EndothelialCells"])

annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]


cell_type_clusters_lvl3 ={ 
    "Chief cells" : ["11","15","1"],
    "GMC" : ["3"] ,  
    "PMC" : ["10","4","5","12","6","0","8","17"],
    "Parietal cells": ["13","19"],
    "Globlet cells" : [] ,  
    "Enterocyte" : ["14","18"]
    }


annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]




if "cell_type_by_cluster" not in cluster_adata_lvl3.obs.columns:
    cluster_adata_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    cluster_adata_lvl3.obs.loc[cluster_adata_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type


cell_type = "Epithelial"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers for {cell_type} Cells in {sample_type} Tissue", fontsize=10, weight='bold', y=1.11)  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby="cell_type_by_cluster",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6),show=False)
plt.suptitle(f"Expression of Cell Markers {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
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
    gene_expression, barcodes = get_gene_expression_for_cell_type(adata=cluster_adata_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
    list_gene_readout_tumor_tissue.append(gene_expression)
    
    df_temp = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": cell_type,  # this adds the cell type for every barcode
        "sample_type": sample_type}) 
    list_barcodes_tumor_tissue.append(df_temp)
    
df_gene_tumor_epi = pd.concat(list_gene_readout_tumor_tissue,ignore_index=True)
df_barcode_tumor_epi = pd.concat(list_barcodes_tumor_tissue)

adata_subset_7 = cluster_adata_lvl3.copy()
adata_subset_3 = cluster_adata_lvl3.copy()

adata_subset_3.obs["cell_type_by_cluster"] = adata_subset_3.obs["cell_type_by_cluster"].replace({
    "Chief cells": "Tumor-Epithelial-cells",
    "GMC": "Tumor-Epithelial-cells",
    "PMC": "Tumor-Epithelial-cells",
    "Parietal cells": "Tumor-Epithelial-cells",
    "Globlet cells": "Tumor-Epithelial-cells",
    "Enterocyte": "Tumor-Epithelial-cells"
    })


  

"""
##########################################################################################################################################
------------------------------------- ParaTUMOR TISSUE (Stroma/Epithilial)-----------------------------------------------------------
##########################################################################################################################################
"""


adata_peritumor = adata_all[
                                (adata_all.obs.batch == "GC03P") |
                                (adata_all.obs.batch == "GC04P") |
                                (adata_all.obs.batch == "GC05P") |
                                (adata_all.obs.batch == "GC06P") |
                                (adata_all.obs.batch == "GC07P") |
                                (adata_all.obs.batch == "GC08P") |
                                (adata_all.obs.batch == "GC09P") |
                                (adata_all.obs.batch ==  "GC10P1") |
                                (adata_all.obs.batch ==  "GC10P2") 
                                ]



path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl2/ParatumorTissue/Stroma"
os.makedirs(path_save, exist_ok=True)

cell_type = "Epithelial"
sample_type ="Paratumor"
cell_type_cluster = ["5","8","14"]
cluster_adata_lymphocytes_lvl2 = cluster_lvl_2_with_harmony(adata_peritumor, cell_type_cluster, path_save, sample_type="Paratumor Tissue", 
                                           cell_type=cell_type, n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4, 
                                           n_top_genes=3000,theta=0.5, harmony=True )
                                                                                                                            #res1.5
#Look at general expreseeion markers
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=['batch'],size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Epithelial"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl1["EpithelialCells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#Differential gene expression between markers 
sc.tl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, groupby="cluster_lvl2", method="wilcoxon")

cell_type = "Epithelial"
fig = sc.pl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(cluster_adata_lymphocytes_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 24), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(cluster_adata_lymphocytes_lvl2,marker_genes_lvl2_epithilial,groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""

clusters_to_remove = ["13"]
selected_lvl2_clusters = list(set(cluster_adata_lymphocytes_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl3/ParaumorTissue/stroma"
os.makedirs(path_save, exist_ok=True)

cluster_adata_lvl3 = cluster_lvl_3_with_harmony(cluster_adata_lymphocytes_lvl2, selected_lvl2_clusters, path_save, sample_type="Paratumor Tissue", 
                                           cell_type="Epithelial", n_pcs=20, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1, n_top_genes=2000, theta=0.5,
                                           cluster_lvl='cluster_lvl2', harmony=True )

cell_type = "Epithelial"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["EpithelialCells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["T&NKcells"], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Epithelial"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 4),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["EndothelialCells"])
sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl1["Erythrocytes"])





cell_type_clusters_lvl3 ={ 
    "Chief cells" : ["5","11","0","4",],
    "GMC" : ["2","16","7"] ,  
    "PMC" : ["8","6","14","9","1","13","3"],
    "Parietal cells": ["15"],
    "Globlet cells" : ["17"] ,  
    "Enterocyte" : ["12","10","18"]
    }


annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]


if "cell_type_by_cluster" not in cluster_adata_lvl3.obs.columns:
    cluster_adata_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    cluster_adata_lvl3.obs.loc[cluster_adata_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type

cell_type = "Epithelial"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers for {cell_type} in {sample_type}", fontsize=10, weight='bold', y=1.11)  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes_lvl2_epithilial, groupby="cell_type_by_cluster",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 6),show=False)
plt.suptitle(f"Expression of Cell Markers {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')




"""
Get GAL9 and Tim3 Expression from lvl3 clusters
##########################################################################################################
"""
list_gene_readout_paratumor_tissue = []
list_barcodes_paratumor_tissue = []
sample_type="Paratumor_Tissue"
for cell_type in list(cell_type_clusters_lvl3.keys()):
    print(cell_type)
    gene_expression, barcodes = get_gene_expression_for_cell_type(adata=cluster_adata_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
    list_gene_readout_paratumor_tissue.append(gene_expression)
    
    df_temp = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": cell_type,  # this adds the cell type for every barcode
        "sample_type": sample_type}) 
    list_barcodes_paratumor_tissue.append(df_temp)
    
df_gene_paratumor_epi = pd.concat(list_gene_readout_paratumor_tissue,ignore_index=True)
df_barcode_paratumor_epi = pd.concat(list_barcodes_paratumor_tissue)


adata_subset_6 = cluster_adata_lvl3.copy()
adata_subset_4 = cluster_adata_lvl3.copy()

adata_subset_4.obs["cell_type_by_cluster"] = adata_subset_4.obs["cell_type_by_cluster"].replace({
    "Chief cells": "Paratumor-Epithelial-cells",
    "GMC": "Paratumor-Epithelial-cells",
    "PMC": "Paratumor-Epithelial-cells",
    "Parietal cells": "Paratumor-Epithelial-cells",
    "Globlet cells": "Paratumor-Epithelial-cells",
    "Enterocyte": "Paratumor-Epithelial-cells"
    })








"""
##########################################################################################################
-----------------------------------------lymphocytes in Blood---------------------------------------------
##########################################################################################################
"""


adata_blood = adata_all[
                                (adata_all.obs.batch == "GC06B") |
                                (adata_all.obs.batch == "GC07B") |
                                (adata_all.obs.batch == "GC08B") |
                                (adata_all.obs.batch == "GC10B")
                                ]



path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl2/Blood/Lymphocytes"
os.makedirs(path_save, exist_ok=True)

cell_type = "Lymphocytes"
sample_type ="Blood"
cell_type_cluster = ["0","1","2","3","4","6","7","13"]
cluster_adata_lymphocytes_lvl2 = cluster_lvl_2_with_harmony(adata_blood, cell_type_cluster, path_save, sample_type="Paratumor Tissue", 
                                           cell_type="Lymphocytes", n_pcs=8, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=0.4, n_top_genes=3000,theta=0.1, 
                                           harmony=True)
                                                                                                                            #res1.5
#Look at general expreseeion markers
custom_palette = ['#FF5733', '#33FF57', '#3357FF', '#FF33A6']  # example colors
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=['batch'],palette=custom_palette,size=3,show=False, return_fig=True)
fig.axes[0].set_title(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_sample_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lymphocytes_lvl2, color=marker_genes_lvl2_lymphocytes[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Ćell Marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
fig.savefig(f"{path_save}/Umap_Marker_genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


#Differential gene expression between markers 
sc.tl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, groupby="cluster_lvl2", method="wilcoxon")
cell_type = "NKTcells"
fig = sc.pl.rank_genes_groups(cluster_adata_lymphocytes_lvl2, n_genes=25, sharey=False,show=False, return_fig=True)
plt.suptitle(f"Rank Genes for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Rank_genes_leiden_cluster_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.rank_genes_groups_heatmap(cluster_adata_lymphocytes_lvl2, n_genes=8,use_raw=False,
    swap_axes=True,vmin=-3,vmax=3, cmap="bwr",figsize=(10, 20), show_gene_labels=True, show=False, return_fig=True)
plt.suptitle(f"Differential Expressed Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Rank_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')
   
fig = sc.pl.heatmap(cluster_adata_lymphocytes_lvl2,marker_genes,groupby="cluster_lvl2",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5), show=False)
plt.suptitle(f"Clusters of Marker Genes for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


"""
Level 2 Clusters
grobes clustering nur für erste selektion und auswerfen von dead stressed and non-of-intereest cell-types
"""

clusters_to_remove = ["9","8"]
selected_lvl2_clusters = list(set(cluster_adata_lymphocytes_lvl2.obs['cluster_lvl2'].cat.categories) - set(clusters_to_remove))


path_save ="/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/lvl3/Blood/Lympo"
os.makedirs(path_save, exist_ok=True)

cluster_adata_lvl3 = cluster_lvl_3_with_harmony(cluster_adata_lymphocytes_lvl2, selected_lvl2_clusters, path_save, sample_type="Tumor Tissue", 
                                           cell_type="TCD8_CD4_NK", n_pcs=12, n_neighbors=10, min_dist=0.3, spread=1.5, resolution=1.2, n_top_genes=2500, theta=0.1,
                                           cluster_lvl='cluster_lvl2', harmony=True )

cell_type = "Tcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP cell marker for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_Tcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "NKTcells"
fig = sc.pl.umap(cluster_adata_lvl3, color=marker_genes_lvl3_T_NKcells[cell_type], wspace=0.2,size=3, ncols=4, show=False, return_fig=True)
plt.suptitle(f"UMAP Clustering for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Umap_Marker_genes_NKTcells_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

cell_type = "TCD8_CD4_NK"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cluster_lvl3', dendrogram=True,standard_scale="group",show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type} per Cluster", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Marker_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')

fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cluster_lvl3",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 4),show=False)
plt.suptitle(f"Cell Marker Expression for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Marker_Genes_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


cell_type_clusters_lvl3 ={ 
    "TCD8+" : ["8","14","16","11"],
    "TCD4+" : ["9","13","5","2","0","6"] ,  
    "Treg" : ["15"],
    "NK": ["7","1","10"]
    }
annotated_clusters = [item for sublist in cell_type_clusters_lvl3.values() for item in sublist]

na_clusters = list(set(list(cluster_adata_lvl3.obs.cluster_lvl3.unique())) - set(annotated_clusters))
#combine clusters

cluster_adata_lvl3 = cluster_adata_lvl3[cluster_adata_lvl3.obs["cluster_lvl3"].isin(annotated_clusters)]



if "cell_type_by_cluster" not in cluster_adata_lvl3.obs.columns:
    cluster_adata_lvl3.obs["cell_type_by_cluster"] = "na"  # Initialize with default
    
for cell_type, clusters in cell_type_clusters_lvl3.items():
    cluster_adata_lvl3.obs.loc[cluster_adata_lvl3.obs["cluster_lvl3"].isin(clusters), "cell_type_by_cluster"] = cell_type

cell_type = "TcellsNKcells"
fig = sc.pl.dotplot(cluster_adata_lvl3, marker_genes, groupby='cell_type_by_cluster', standard_scale="group",show=False)
plt.suptitle(f"Expression of Cell Markers for {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Dotplot_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')


fig = sc.pl.heatmap(cluster_adata_lvl3, marker_genes, groupby="cell_type_by_cluster",use_raw=False,
    vmin=-3,vmax=3,cmap="RdBu_r",dendrogram=True,swap_axes=True,figsize=(11, 5),show=False)
plt.suptitle(f"Expression of Cell Markers {cell_type} in {sample_type}", fontsize=10, weight='bold')  # Title Font Size
plt.savefig(f"{path_save}/Heatmap_Cell_Markers_{cell_type}_{sample_type}.png", dpi=300, bbox_inches='tight')




"""
Get GAL9 and Tim3 Expression from lvl3 clusters
##########################################################################################################
"""
list_gene_readout_blood = []
list_barcodes_blood = []
sample_type="Blood"
for cell_type in list(cell_type_clusters_lvl3.keys()):
    print(cell_type)
    gene_expression, barcodes = get_gene_expression_for_cell_type(adata=cluster_adata_lvl3, cell_type_by_cluster=cell_type, sample_type=sample_type)
    list_gene_readout_blood.append(gene_expression)
    
    df_temp = pd.DataFrame({
        "barcode": barcodes,
        "cell_type": cell_type,  # this adds the cell type for every barcode
        "sample_type": sample_type}) 
    list_barcodes_blood.append(df_temp)   
df_gene_blood_lymph = pd.concat(list_gene_readout_blood,ignore_index=True)
df_barcode_blood_lymph = pd.concat(list_barcodes_blood)


adata_subset_5 = cluster_adata_lvl3.copy()

adata_subset_5.obs["cell_type_by_cluster"] = adata_subset_5.obs["cell_type_by_cluster"].replace({
    "TCD8+": "Blood-TCD8+",
    "TCD4+": "Blood-TCD4+",
    "Treg": "Blood-Treg",
    "NK": "Blood-NK"})





adata_subset_T = adata_tumor_tissue[adata_tumor_tissue.obs_names.isin(tumor_cells)]
adata_subset_T.obs["cell_type_by_cluster"] = "Tumor-cells"







adata_subset_concate = sc.concat([adata_subset_1,adata_subset_2,adata_subset_3,adata_subset_4,adata_subset_5,adata_subset_T])

adata_subset_concate.write("/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/plots_final/adata_cellsubsets.h5ad")






"""
####################################################################################
Hier geht's weiter
###############################################################################
"""



adata_subset_concate = sc.read("/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/plots_final/adata_cellsubsets.h5ad")




sample_groups_dict = {
    "GC06" : ["GC06B", "GC06P", "GC06T"],
    "GC07" : ["GC07B", "GC07P", "GC07T"],
    "GC08" : ["GC08B", "GC08P", "GC08T1"],
    "GC10" : ["GC10B", "GC10P1", "GC10T"]
    }

#Create new column in adata to annotate samples on patient level
adata_list = []
for x in list(sample_groups_dict.keys()):  
   adata = adata_subset_concate[adata_subset_concate.obs["batch"].isin(sample_groups_dict[x])]
   adata.obs["Patient ID"] = x
   adata_list.append(adata)
   
adata_concate = sc.concat(adata_list)
 

adata_tumor_cd8 = adata_concate[adata_concate.obs["cell_type_by_cluster"] == "Tumor-TCD8+"]
adata_tumor_cd8.obs.batch.unique()
# Output ['GC06T', 'GC07T', 'GC08T1', 'GC10T']

df_adata = pd.DataFrame(data = adata_concate.layers["counts_normalized"].toarray(),
             index = adata_concate.obs.index,
             columns = adata_concate.var.index)



df_adata = df_adata.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")



















df_adata[df_adata["cell_type_by_cluster"] == "Tumor-TCD8+"]["batch"].unique()
df_adata.to_csv(os.path.join("/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/plots_final/","df_adata_cells_genes_annotations.csv"), sep=";")



"""
#####################################################################################################
-----------------------------------------  Correlation Plot  ----------------------------------------
#####################################################################################################
1.: Get all transcripts of interest 
    PDCD1, LAG3, TIGIT, ENTPD1, CTLA4, TOX, IL2, IFNG, TNF, GZMB


2. ADAM Cleaved TIM-3    
ADAM17 TACE CD156B
ADAM10 CD156C
CEACAM1 
"""


df_adata_ex = df_adata[["HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "IL2", "IFNG", "TNF", "GZMB", "ADAM10", "ADAM17", "CEACAM1"]]
df_adata_ex = df_adata_ex.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")


df_adata_ex.columns
"""
['HAVCR2', 'LGALS9', 'PDCD1', 'LAG3', 'TIGIT', 'ENTPD1', 'CTLA4', 'TOX',
       'IL2', 'IFNG', 'TNF', 'GZMB', 'Patient ID', 'batch',
       'cell_type_by_cluster']


batch: ['GC06T', 'GC06P', 'GC06B', 'GC07T', 'GC07P', ..., 'GC08P', 'GC08B', 'GC10T', 'GC10P1', 'GC10B']
cell_type_by_cluster: 'Tumor-TCD8+', 'Tumor-TCD4+', 'Tumor-NK', 'Tumor-Treg', 'Paratumor-TCD8+', ..., 'Blood-NK', 'Blood-TCD4+', 'Blood-TCD8+', 'Blood-Treg', 'Tumor-cells']

what do i want to plot? 
are the specific genes expressed? 
(Dotplot)


 
(heatmap)

1.: z-norm per variable over all sample and celltypes 
2.: cap to +-3
3.: heatmap mean expression per sample and cell type
4.: sort by clusters 
 
correlation map and cluster -> gives an exhaution profile 
"""

"""
#####################################################################################################
1.: 
    ADAM10/17 in tissue vs. sTim3 in blood. vgl. nicht möglich
    
2.: 
    Korrelation of all markers 
"""


"""
2.: 
"""
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy import stats

df_adata_ex.columns


df_adata_ex = df_adata[["HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "TOX", "IL2", "IFNG", "TNF", "GZMB", "ADAM10", "ADAM17", "CEACAM1"]]
df_adata_ex = df_adata_ex.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")


conditions = df_adata_ex.cell_type_by_cluster.unique()
variables = ["HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "TOX"]
path_save_rebuttel = "/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/plots_final/rebuttel"


corr_mats = {}   # condition -> correlation matrix (DataFrame)
pval_mats = {}   # condition -> p-value matrix (DataFrame)


for condition in conditions:
    df_condition = df_adata_ex[df_adata_ex["cell_type_by_cluster"] == condition]
    print(condition)
    
    corr_mat = pd.DataFrame(np.nan, index=variables, columns=variables, dtype=float)
    pval_mat = pd.DataFrame(np.nan, index=variables, columns=variables, dtype=float)
    np.fill_diagonal(corr_mat.values, 1.0)
    np.fill_diagonal(pval_mat.values, 0.0)

    
    
    for variable_1, variable_2 in itertools.combinations(variables, 2):
        print(variable_1)
        print(variable_2)

            #take non zero values for both variables 
        df_filtered = df_condition[~ ((df_condition[variable_1] == 0) &  (df_condition[variable_2] == 0)) ]
        
        if len(df_filtered) < 2:
            continue
        
        x = df_filtered[variable_1]
        y = df_filtered[variable_2]
       
        r, p = stats.pearsonr(x, y)
        
#####################################################################################
#pearson correlation heatmap       
                
        
        
        corr_mat.loc[variable_1, variable_2] = r
        corr_mat.loc[variable_2, variable_1] = r
        pval_mat.loc[variable_1, variable_2] = p
        pval_mat.loc[variable_2, variable_1] = p

    corr_mat.index = [f"{condition}_{v}" for v in corr_mat.index]
    corr_mat.columns = [f"{condition}_{v}" for v in corr_mat.columns]
    
    corr_mats[condition] = corr_mat
    
    plt.figure(figsize=(7, 6))
    sns.heatmap(
          corr_mat, vmin=-1, vmax=1, cmap="vlag", annot=True, fmt=".2f",
          square=True, cbar_kws={"label": "Pearson r"})
    plt.title(f"Correlation heatmap ({condition})")
    plt.tight_layout()
    out_png = os.path.join(path_save_rebuttel, f"corr_heatmap_{condition}.png")
    plt.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.close()

     
#####################################################################################
#pearson correlation plot         
        
        hue = df_filtered["batch"].cat.remove_unused_categories()
                 
        x = df_filtered[variable_1]
        y = df_filtered[variable_2]
       
        r, p = stats.pearsonr(x, y)
                
        plt.figure(figsize=(6, 6))
        # Scatter with hue
        sns.scatterplot(x=x, y=y, hue=hue, palette="tab10", alpha=0.7, edgecolor="none")
        
        # Regression line (all points together, no hue)
        sns.regplot(x=x, y=y, scatter=False, line_kws={'color':'red'})
        
        # Annotate correlation
        plt.text(0.05, 0.95, f"ρ = {r:.2f}, p = {p:.2e}", 
                 transform=plt.gca().transAxes, 
                 ha='left', va='top', fontsize=12, 
                 bbox=dict(facecolor='white', alpha=0.7, edgecolor='none'))
        
        plt.xlabel(variable_1)
        plt.ylabel(variable_2)
        plt.title(f"Scatterplot of {variable_1} vs {variable_2} for {condition}")
        plt.legend(title="Batch", bbox_to_anchor=(1.05, 1), loc='upper left')
        plt.tight_layout()
        outpath = os.path.join(path_save_rebuttel, f"{condition}_{variable_1}_{variable_2}.png")
        plt.savefig(outpath, dpi=300, bbox_inches="tight")
        plt.close() 


#####################################################################################
#pearson correlation heatmap all combinations          
        




                
            
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

df_adata_ex = df_adata[["HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "TOX", "IL2", "IFNG", "TNF", "GZMB" ]]
df_adata_ex = df_adata_ex.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")
df_adata_ex = df_adata_ex.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")

df_adata_ex = df_adata_ex[["batch","cell_type_by_cluster","HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "TOX", "IL2", "IFNG", "TNF", "GZMB" ]]


df_adata_ex.columns
df_adata_ex.cell_type_by_cluster

df = df_adata_ex.copy()

df[['sample_type','cell_type']] = df['cell_type_by_cluster'].str.split('-', n=1, expand=True)
mask = df["cell_type"]=="cells"
df.loc[mask, 'cell_type'] = 'Tumor-cells' 

path = "/home/willmstid/Desktop/New Folder/"
df.to_excel(os.path.join(path,"readout_scRNA.xlsx"), index=False)




df_adata_ex.columns
df_adata_ex.head()






# select only gene expression columns (exclude metadata)
gene_cols = ["HAVCR2","LGALS9", "PDCD1", "LAG3", "TIGIT", 
             "ENTPD1", "CTLA4", "TOX", "IL2", "IFNG", "TNF", "GZMB"]

# z-normalize per gene (subtract mean, divide by std)
df_z = df_adata_ex[gene_cols].apply(lambda x: (x - x.mean()) / x.std(), axis=0)

# calculating corrmatrix
corr_matrix = df_z.corr()

#

plt.figure(figsize=(10,8))
sns.heatmap(corr_matrix, 
            cmap="RdBu_r",    # blue = negative, red = positive
            center=0,         # 0 correlation at white
            annot=True,       # show correlation values
            fmt=".2f",        # format values
            cbar_kws={'label': 'Correlation'})

plt.title("Correlation Matrix of Immune Checkpoint & Cytokine Genes", fontsize=14)
plt.tight_layout()
plt.show()
################################################################################

# your genes of interest
df_adata_ex.cell_type_by_cluster.unique()



genes = ['HAVCR2','LGALS9','PDCD1','LAG3','TIGIT','ENTPD1','CTLA4','TOX','IL2','IFNG','TNF','GZMB']

df = df_adata_ex.copy()




# Split "Tumor-TCD8+" into sample_type="Tumor", cell_type_raw="TCD8+"
df[['sample_type','cell_type']] = df['cell_type_by_cluster'].str.split('-', n=1, expand=True)
mask = df["cell_type"]=="cells"
df.loc[mask, 'cell_type'] = 'Tumor-cells' 

df.columns

pb = (df.groupby(['batch','sample_type','cell_type'])[genes]
        .mean()
        .reset_index())

group_gene = pb.groupby(['sample_type','cell_type'])[genes].mean()  # rows = groups, cols = genes


corr_groups = group_gene.T.corr()   # index/columns are (sample_type, cell_type)

# Pretty labels like "Tumor • CD8"
nice = corr_groups.copy()
nice.index   = [f"{a} • {b}" for a,b in nice.index]
nice.columns = [f"{a} • {b}" for a,b in nice.columns]

# Heatmap with matplotlib
fig, ax = plt.subplots(figsize=(9,7))
im = ax.imshow(nice.values, vmin=-1, vmax=1)
ax.set_xticks(range(len(nice.columns))); ax.set_xticklabels(nice.columns, rotation=90)
ax.set_yticks(range(len(nice.index)));   ax.set_yticklabels(nice.index)
cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04); cbar.set_label('Pearson r')
fig.tight_layout()
plt.show()


"""
1.: get number of cells per sample and celltype 
2.: get number of pos cells per sample and celltype for specific gene 
3.: get relative number of cells, percentage of pos cells for specific gene 
4.: 
"""

genes = ['HAVCR2','LGALS9','PDCD1','LAG3','TIGIT','ENTPD1','CTLA4','','IL2','IFNG','TNF','GZMB']

genes_exhaustion = ["PDCD1", "LAG3", "TIGIT", "ENTPD1", "CTLA4", "TOX"]
genes_effector =["IL2", "IFNG", "TNF", "GZMB"]
gene_other = ["ADAM10", "ADAM17", "CEACAM1"]

path_save_rebuttel = "/home/willmstid/Desktop/New Folder/"




def get_dotplot(df_adata, gene_name):
    
    gene_name= "TOX"
    
    df_adata_tim_gal = df_adata[[gene_name]]
    df_adata_tim_gal = df_adata_tim_gal.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")
    
    #get number of cells
    df_count_per_batch_and_cell_type = df_adata_tim_gal.groupby(["Patient ID", "cell_type_by_cluster"]).count()
    #get number of non zero cells
    df_count_per_cell_type_tim_non_zero = df_adata_tim_gal[df_adata_tim_gal[gene_name] > 0.0].groupby(["Patient ID", "cell_type_by_cluster"]).count()[gene_name]
    #get relative number of cells
    df_relative_cell_number_for_tim = df_count_per_cell_type_tim_non_zero / df_count_per_batch_and_cell_type[gene_name]
    #get mean expression of non zero cells
    df_mean_per_cell_type_tim_non_zero = df_adata_tim_gal[["Patient ID", "cell_type_by_cluster", gene_name]].groupby(["Patient ID", "cell_type_by_cluster"]).mean()
    # [1,0] normalization 
    df_mean_per_cell_type_tim_non_zero[gene_name] = (df_mean_per_cell_type_tim_non_zero[gene_name] - df_mean_per_cell_type_tim_non_zero[gene_name].min()) / (df_mean_per_cell_type_tim_non_zero[gene_name].max() - df_mean_per_cell_type_tim_non_zero[gene_name].min())
    
        
    df_mean_melted = df_mean_per_cell_type_tim_non_zero.reset_index().melt(
        id_vars=["Patient ID", "cell_type_by_cluster"],
        var_name="Gene",
        value_name="Mean Expression")
    
    
    df_fraction_melted = df_relative_cell_number_for_tim.reset_index().melt(
        id_vars=["Patient ID", "cell_type_by_cluster"],
        var_name="Gene",
        value_name="Fraction")

    df_merged = pd.merge(
        df_mean_melted, df_fraction_melted,
        on=["Patient ID", "cell_type_by_cluster", "Gene"],
        how="inner")
        
    df_merged["Fraction of Cells"] = (df_merged["Fraction"] * 100).round(decimals=3)
    df_merged["Patient ID"] = df_merged["Patient ID"].replace({"GC06": "1", "GC07": "2", "GC08": "3", "GC10": "4"})
            
    plt.figure(figsize=(1.5, 4.7))
    sns.set_style("white")    
    ax = sns.scatterplot(
            data=df_merged,
            x="Patient ID",
            y="cell_type_by_cluster",
            size="Fraction of Cells",         
            hue="Mean Expression",  
            sizes=(1, 300),         
            #size_norm= (df_gene["Fraction of Cells"].min(),df_gene["Fraction of Cells"].max()),
            size_norm = (1,80),
            palette="Blues",  #Purples Reds
            #edgecolor="black",
            alpha=0.8
    
        )
    
    
   # ax.legend(bbox_to_anchor=(1.05, 0.9), loc='upper left')# borderaxespad=0.)
    
    ax.legend(
    bbox_to_anchor=(1.05, 0.9),
    loc="upper left",
    borderaxespad=0.,
    handleheight=0.5,   # increases vertical spacing
    labelspacing=0.8)
    
    ax.get_legend().remove()
    ax.set_yticklabels([])
    ax.margins(x=0.2)  #
    ax.margins(y=0.1)
    plt.xlabel("")
    plt.ylabel("")
    plt.title(f"{gene_name}")
    outpath = os.path.join(path_save_rebuttel, f"TOX_without_boarder_Dotplot.png")
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    plt.show() 
    plt.close() 
    
    return print("save dot plot for {gene_name}")

    gene_other
    
    
for gene in gene_other:
    
    get_dotplot(df_adata, gene)


"""
#####################################################################################################################
Dotplot for tim3 and gal-9

"""









df_adata_tim_gal = df_adata[["HAVCR2","LGALS9"]]
df_adata_tim_gal = df_adata_tim_gal.merge(adata_concate.obs[["Patient ID","batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")


df_adata_tim_gal.columns
"""
Count number of cell_type_by_cluster for each batch and cell type
Count number of cell_type_by_cluster for each batch and cell type which are non zero
"""

df_count_per_batch_and_cell_type = df_adata_tim_gal.groupby(["Patient ID", "cell_type_by_cluster"]).count()


df_count_per_cell_type_tim_non_zero = df_adata_tim_gal[df_adata_tim_gal["HAVCR2"] > 0.0].groupby(["Patient ID", "cell_type_by_cluster"]).count()["HAVCR2"]
df_count_per_cell_type_Gal_non_zero = df_adata_tim_gal[df_adata_tim_gal["LGALS9"] > 0.0].groupby(["Patient ID", "cell_type_by_cluster"]).count()["LGALS9"]

df_relative_cell_number_for_tim = df_count_per_cell_type_tim_non_zero / df_count_per_batch_and_cell_type["HAVCR2"]
df_relative_cell_number_for_gal = df_count_per_cell_type_Gal_non_zero / df_count_per_batch_and_cell_type["LGALS9"]

df_mean_per_cell_type_tim_non_zero = df_adata_tim_gal[["Patient ID", "cell_type_by_cluster", "HAVCR2"]].groupby(["Patient ID", "cell_type_by_cluster"]).mean()
df_mean_per_cell_type_Gal_non_zero = df_adata_tim_gal[["Patient ID", "cell_type_by_cluster", "LGALS9"]].groupby(["Patient ID", "cell_type_by_cluster"]).mean()

df_mean_per_cell_type_Gal_non_zero["LGALS9"] = (df_mean_per_cell_type_Gal_non_zero["LGALS9"] - df_mean_per_cell_type_Gal_non_zero["LGALS9"].min()) / (df_mean_per_cell_type_Gal_non_zero["LGALS9"].max() - df_mean_per_cell_type_Gal_non_zero["LGALS9"].min())
df_mean_per_cell_type_tim_non_zero["HAVCR2"] = (df_mean_per_cell_type_tim_non_zero["HAVCR2"] - df_mean_per_cell_type_tim_non_zero["HAVCR2"].min()) / (df_mean_per_cell_type_tim_non_zero["HAVCR2"].max() - df_mean_per_cell_type_tim_non_zero["HAVCR2"].min())


df_merge_genes = pd.merge(
    df_mean_per_cell_type_tim_non_zero, df_mean_per_cell_type_Gal_non_zero,
    on=["Patient ID", "cell_type_by_cluster"],
    how="inner"
)



df_merge_counts = pd.merge(
    df_relative_cell_number_for_tim, df_relative_cell_number_for_gal,
    on=["Patient ID", "cell_type_by_cluster"],
    how="inner"
)




df_mean_melted = df_merge_genes.reset_index().melt(
    id_vars=["Patient ID", "cell_type_by_cluster"],
    var_name="Gene",
    value_name="Mean Expression"
)


df_fraction_melted = df_merge_counts.reset_index().melt(
    id_vars=["Patient ID", "cell_type_by_cluster"],
    var_name="Gene",
    value_name="Fraction"
)



df_merged = pd.merge(
    df_mean_melted, df_fraction_melted,
    on=["Patient ID", "cell_type_by_cluster", "Gene"],
    how="inner"
)


df_merged["Fraction of Cells"] = (df_merged["Fraction"] * 100).round(decimals=3)


df_merged['Patient ID'] = df_merged['Patient ID'].replace({
    'GC06': '1',
    'GC07': '2',
    'GC08' : "3",
    'GC10' : "4"
})


df_gene = df_merged[df_merged["Gene"] == "HAVCR2"]
    
plt.figure(figsize=(4, 6))
sns.set_style("white")    
ax = sns.scatterplot(
        data=df_gene,
        x="Patient ID",
        y="cell_type_by_cluster",
        size="Fraction of Cells",         
        hue="Mean Expression",  
        sizes=(1, 300),         
        #size_norm= (df_gene["Fraction of Cells"].min(),df_gene["Fraction of Cells"].max()),
        size_norm = (1,60),
        palette="Oranges",
        edgecolor="black",
        alpha=0.8

    )

ax.legend(bbox_to_anchor=(1.05, 0.9), loc='upper left', borderaxespad=0.)
ax.margins(x=0.2)  #
ax.margins(y=0.1)
plt.ylabel("Cell Types")
plt.title("HAVCR2 Expression across Samples and Cell Types")
plt.show()


 


df_gene = df_merged[df_merged["Gene"] == "LGALS9"]

plt.figure(figsize=(2, 6))
sns.set_style("white")    
ax = sns.scatterplot(
        data=df_gene,
        x="Patient ID",
        y="cell_type_by_cluster",
        size="Fraction of Cells",        
        hue="Mean Expression",   
        sizes=(1, 300),        
        palette="Greens",
        #size_norm= (df_gene["Fraction of Cells"].min(),df_gene["Fraction of Cells"].max()),
        edgecolor="black",
        size_norm = (1,60),
        alpha=0.8

    )

ax.legend(bbox_to_anchor=(1.05, 0.9), loc='upper left', borderaxespad=0.)
ax.margins(x=0.2)  #
ax.margins(y=0.1)
plt.ylabel("Cell Types")
plt.title("LGALS9 Expression across Samples and Cell Types")
plt.tight_layout()
plt.show()







"""
sc.pl.dotplot(adata_concate, groupby="cell_type_by_cluster", var_names=["HAVCR2"] , layer="counts_normalized", standard_scale="var",use_raw=False,
              #mean_only_expressed=True,
              colorbar_title="Mean expression\nper gene",
              categories_order=[
              "Blood-NK",
              "Blood-TCD4+",
              "Blood-TCD8+",
              "Blood-Treg",
              "Paratumor-NK",
              "Paratumor-TCD4+",
              "Paratumor-TCD8+",
              "Paratumor-Treg",
              "Tumor-NK",
              "Tumor-TCD4+",
              "Tumor-TCD8+",
              "Tumor-Treg",
              "Tumor-Epithelial-cells",
              "Paratumor-Epithelial-cells",
              "Tumor-cells"],
              palette="Oranges",
              #figsize=(3, 6)
              )

sc.pl.dotplot(adata_concate, groupby="cell_type_by_cluster", var_names=["LGALS9"] , layer="counts_normalized", standard_scale="var",use_raw=False,
              #mean_only_expressed=True,
              colorbar_title="Mean expression\nper gene",
              categories_order=[
              "Blood-NK",
              "Blood-TCD4+",
              "Blood-TCD8+",
              "Blood-Treg",
              "Paratumor-NK",
              "Paratumor-TCD4+",
              "Paratumor-TCD8+",
              "Paratumor-Treg",
              "Tumor-NK",
              "Tumor-TCD4+",
              "Tumor-TCD8+",
              "Tumor-Treg",
              "Tumor-Epithelial-cells",
              "Paratumor-Epithelial-cells",
              "Tumor-cells"],
              palette="Greens",
              #figsize=(3, 6)
              )


sc.pl.matrixplot(adata_subset_concate, var_names=["HAVCR2", "LGALS9"], groupby="cell_type_by_cluster",layer="counts_normalized", cmap="Reds",
                 standard_scale="var", colorbar_title="Mean expression\nper gene",
                 categories_order=[
                 "Blood-NK",
                 "Blood-TCD4+",
                 "Blood-TCD8+",
                 "Blood-Treg",
                 "Paratumor-NK",
                 "Paratumor-TCD4+",
                 "Paratumor-TCD8+",
                 "Paratumor-Treg",
                 "Tumor-NK",
                 "Tumor-TCD4+",
                 "Tumor-TCD8+",
                 "Tumor-Treg",
                 "Tumor-Epithelial-cells",
                 "Paratumor-Epithelial-cells",
                 "Tumor-cells"
                 ])


"""



"""
############################################################################################################################################
--------------------------------------------------------------Normal Dotplot for all samples------------------------------------------------
############################################################################################################################################
"""








df_adata = pd.DataFrame(data = adata_concate.layers["counts_normalized"].toarray(),
             index = adata_concate.obs.index,
             columns = adata_concate.var.index)


df_adata_tim_gal = df_adata[["HAVCR2","LGALS9"]]
df_adata_tim_gal = df_adata_tim_gal.merge(adata_concate.obs[["batch", "cell_type_by_cluster"]], left_index=True, right_index=True, how="left")


"""
Count number of cell_type_by_cluster for each batch and cell type
Count number of cell_type_by_cluster for each batch and cell type which are non zero
"""

df_count_per_batch_and_cell_type = df_adata_tim_gal.groupby(["cell_type_by_cluster"]).count()


df_count_per_cell_type_tim_non_zero = df_adata_tim_gal[df_adata_tim_gal["HAVCR2"] > 0.0].groupby(["cell_type_by_cluster"]).count()["HAVCR2"]
df_count_per_cell_type_Gal_non_zero = df_adata_tim_gal[df_adata_tim_gal["LGALS9"] > 0.0].groupby(["cell_type_by_cluster"]).count()["LGALS9"]

df_relative_cell_number_for_tim = df_count_per_cell_type_tim_non_zero / df_count_per_batch_and_cell_type["HAVCR2"]
df_relative_cell_number_for_gal = df_count_per_cell_type_Gal_non_zero / df_count_per_batch_and_cell_type["LGALS9"]

df_mean_per_cell_type_tim_non_zero = df_adata_tim_gal[["cell_type_by_cluster", "HAVCR2"]].groupby(["cell_type_by_cluster"]).mean()
df_mean_per_cell_type_Gal_non_zero = df_adata_tim_gal[[ "cell_type_by_cluster", "LGALS9"]].groupby(["cell_type_by_cluster"]).mean()

df_mean_per_cell_type_Gal_non_zero["LGALS9"] = (df_mean_per_cell_type_Gal_non_zero["LGALS9"] - df_mean_per_cell_type_Gal_non_zero["LGALS9"].min()) / (df_mean_per_cell_type_Gal_non_zero["LGALS9"].max() - df_mean_per_cell_type_Gal_non_zero["LGALS9"].min())
df_mean_per_cell_type_tim_non_zero["HAVCR2"] = (df_mean_per_cell_type_tim_non_zero["HAVCR2"] - df_mean_per_cell_type_tim_non_zero["HAVCR2"].min()) / (df_mean_per_cell_type_tim_non_zero["HAVCR2"].max() - df_mean_per_cell_type_tim_non_zero["HAVCR2"].min())


df_merge_genes = pd.merge(
    df_mean_per_cell_type_tim_non_zero, df_mean_per_cell_type_Gal_non_zero,
    on=[ "cell_type_by_cluster"],
    how="inner"
)



df_merge_counts = pd.merge(
    df_relative_cell_number_for_tim, df_relative_cell_number_for_gal,
    on=["cell_type_by_cluster"],
    how="inner"
)




df_mean_melted = df_merge_genes.reset_index().melt(
    id_vars=["cell_type_by_cluster"],
    var_name="Gene",
    value_name="Mean Expression"
)


df_fraction_melted = df_merge_counts.reset_index().melt(
    id_vars=[ "cell_type_by_cluster"],
    var_name="Gene",
    value_name="Fraction"
)



df_merged_2 = pd.merge(
    df_mean_melted, df_fraction_melted,
    on=[ "cell_type_by_cluster", "Gene"],
    how="inner"
)


df_merged_2["Fraction of Cells"] = (df_merged_2["Fraction"] * 100).round(decimals=3)





df_gene = df_merged_2[df_merged_2["Gene"] == "HAVCR2"]
    
plt.figure(figsize=(0.5, 6))
sns.set_style("white")    
ax = sns.scatterplot(
        data=df_gene,
        x="Gene",
        y="cell_type_by_cluster",
        size="Fraction of Cells",         
        hue="Mean Expression",  
        sizes=(1, 300),         
        #size_norm= (df_gene["Fraction of Cells"].min(),df_gene["Fraction of Cells"].max()),
        size_norm = (1,60),
        palette="Oranges",
        edgecolor="black",
        alpha=0.8

    )

ax.legend(bbox_to_anchor=(1.05, 0.9), loc='upper left', borderaxespad=0.)

plt.ylabel("Cell Types")
plt.title("HAVCR2 Expression across Samples and Cell Types")
plt.tight_layout()
plt.show()


 


df_gene = df_merged_2[df_merged_2["Gene"] == "LGALS9"]

plt.figure(figsize=(0.5, 6))
sns.set_style("white")  
ax = sns.scatterplot(
        data=df_gene,
        x="Gene",
        y="cell_type_by_cluster",
        size="Fraction of Cells",        
        hue="Mean Expression",   
        sizes=(1, 300),        
        palette="Greens",
        #size_norm= (df_gene["Fraction of Cells"].min(),df_gene["Fraction of Cells"].max()),
        size_norm = (1,60),
        edgecolor="black",
        #hue_norm=(0, 1),
        alpha=0.8

    )

ax.legend(bbox_to_anchor=(1.05, 0.9), loc='upper left', borderaxespad=0.)

plt.ylabel("Cell Types")
plt.title("LGALS9 Expression across Samples and Cell Types")
plt.tight_layout()
plt.show()

