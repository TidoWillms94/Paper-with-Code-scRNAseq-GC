
import numpy as np
import scanpy as sc
import matplotlib.pyplot as plt 
from tqdm import tqdm
from joblib import Parallel, delayed
from scipy.sparse import csr_matrix


sc.settings.set_figure_params( frameon=True,color_map='Spectral_r') 
sc.settings.verbosity = 3    
plt.rcParams["axes.grid"] = False



dir_samples = [
#Tumortissue
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207224/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207225/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207228/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207229/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207231/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207233/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207236/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207239/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207242/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207243/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207245/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207249/outs/filtered_feature_bc_matrix",
#Paratumor Tissue
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207227/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207230/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207232/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207235/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207238/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207241/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207244/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207247/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207248/outs/filtered_feature_bc_matrix",
#Blood
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207234/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207237/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207240/outs/filtered_feature_bc_matrix",
"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/data/HRR207246/outs/filtered_feature_bc_matrix",
    ]

sample_names = [
    "GC01T","GC02T","GC03T1","GC03T2","GC04T", "GC05T","GC06T","GC07T","GC08T1","GC08T2","GC09T","GC10T",
    "GC03P","GC04P","GC05P","GC06P","GC07P","GC08P","GC09P","GC10P1","GC10P2",
    "GC06B","GC07B","GC08B","GC10B"             
        ]
   




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




"""
#######################################################################################################
scRNAseq Workflow - Preprocessing

1.  Read 10X files                
2.  Filter Cells             
3.  Normalizing Counts       
4.  Transfroming Counts      
5.  Doublets Detection via scDblFinder in R per batch
6.  Scaling, z-norm          
7.  Get HVG                  
8.  Run PCA                 
9.  Get Next Neighbour       
10. Get UMAP                 
11. Get Cluster, Leiden      
12. Cluster Annotation with Cell Marker
13. Save AnnData Object
#######################################################################################################
"""

"""
###########################################################################################################
-----------------------------------------------1-4---------------------------------------------------------
###########################################################################################################
"""

def data_filtering(data, min_genes, min_cells, min_counts, mt_pct, mt_pct2, ncount2, max_fraction, plot):

    data.var_names_make_unique()                                #check for unique gene names
                                       
	#filter cells and genes 
    sc.pp.filter_cells(data,  min_counts=min_counts)            #min number of transcripts per cell eg. 400
    sc.pp.filter_cells(data,  min_genes=min_genes)              #min feature count per cell eg. 200 
    sc.pp.filter_genes(data,  min_cells=min_cells)              # min cells eg.3

    mito_genes = data.var_names.str.startswith('MT-')
    data.obs['percent_mito'] = np.sum(data[:, mito_genes].X, axis=1).A1 / np.sum(data.X, axis=1).A1 
    data.obs['n_counts'] = data.X.sum(axis=1).A1 
    
    if plot:
        sc.pl.violin(data, ['n_genes', 'n_counts', 'percent_mito'],jitter=0.4, multi_panel=True)
        sc.pl.scatter(data, x='n_counts', y='percent_mito')
        sc.pl.scatter(data, x='n_counts', y='n_genes')
        
    data = data[data.obs['percent_mito'] < mt_pct, :]
    data = data[ ~((data.obs['percent_mito'] > mt_pct2) & (data.obs['n_counts']<ncount2))  , :]  # remove cells with low gene count and high mt
    sc.pp.normalize_total(data, target_sum=None, exclude_highly_expressed=True, max_fraction=max_fraction) 
    sc.pp.log1p(data, base=2)  
    return data



def read_samples(dir_sample, sample_name):
    adata = sc.read_10x_mtx(dir_sample, var_names='gene_symbols', cache=True)
    adata.obs['batch'] = sample_name
    adata = data_filtering(adata, min_genes=200, min_cells=3, min_counts=400,  mt_pct=0.3, mt_pct2=0.3, ncount2=1000, max_fraction=0.1, plot=True)
    data_array = adata.X.T
    if isinstance(data_array, csr_matrix):
        data_array = data_array.toarray()  
    return adata, data_array


#read and filter samples in parallel
results = Parallel(n_jobs=25)( 
    delayed(read_samples)(dir_sample, sample_name) 
    for dir_sample, sample_name in tqdm(zip(dir_samples, sample_names)))

adata_list, data_array_list = zip(*results)
adata_list = list(adata_list)
data_array_list = list(data_array_list)




"""
#######################################################################################################
-------------------------------5. Doublets Detection via scDblFinder in R------------------------------
#######################################################################################################
"""


import rpy2.robjects as ro
from rpy2.robjects import pandas2ri
pandas2ri.activate()

ro.r('library(Seurat)')
ro.r('library(scater)')
ro.r('library(scDblFinder)')
ro.r('library(BiocParallel)')
ro.r('library(future)')
ro.r('library(future.apply)')


j = len(data_array_list)
ro.globalenv['data_array_list'] = ro.r('list()')  #R list
for i, data_array in tqdm(enumerate(data_array_list)):
    i = i+1
    print(f"Converting pandas Dataframe to R Dataframe {i}/{j}")
    ro.globalenv['data_array'] = data_array  
    ro.globalenv['i'] = i
    ro.r("""data_array_list[[{i}]] <- data_array""") 



#R code scDblFinder
ro.r('''
doublet_score_list <- list()
doublet_class_list <- list()
set.seed(123)

plan(multisession, workers = 25)
results <- future_lapply(data_array_list, function(obj) {
    sce <- scDblFinder(SingleCellExperiment(list(counts=obj)))
    doublet_score <- sce$scDblFinder.score
    doublet_class <- sce$scDblFinder.class
    list(doublet_score = doublet_score, doublet_class = doublet_class)
})

doublet_score_list <- lapply(results, function(x) x$doublet_score)
doublet_class_list <- lapply(results, function(x) x$doublet_class)
plan(sequential)
''')

# Convert the R lists to NumPy array
doublet_score_list = [np.array(r_doublet_score) for r_doublet_score in ro.globalenv['doublet_score_list']]
doublet_class_list = [np.array(r_doublet_class) for r_doublet_class in ro.globalenv['doublet_class_list']]

adata_list_QC_doublets = []
for i in range(len(adata_list)):
    doublet_score = doublet_score_list[i]
    doublet_class = doublet_class_list[i]
    adata = adata_list[i]
    adata.obs["scDblFinder_score"] = doublet_score
    adata.obs["scDblFinder_class"] = doublet_class
    adata_list_QC_doublets.append(adata)
 
    
adata_combined = adata_list_QC_doublets[0].concatenate(*adata_list_QC_doublets[1:], batch_key="batch", batch_categories=sample_names, join='outer') 



"""
###########################################################################################################
----------------------------------------------6-11----------Level_1_Clusters-------------------------------
###########################################################################################################
"""

sc.pp.highly_variable_genes(adata_combined, n_top_genes = 6000, flavor ='seurat', batch_key='batch')

adata_combined.layers["counts_normalized"] = adata_combined.X.copy()  
adata_combined = sc.pp.scale(adata_combined, max_value=10, copy=True, zero_center=True )  

sc.tl.pca(adata_combined, svd_solver='arpack', use_highly_variable = True)
sc.pp.neighbors(adata_combined, n_pcs=10, n_neighbors=500)
sc.tl.umap(adata_combined, init_pos=adata_combined.obsm['X_pca'][:,:2] *0.0001)
sc.tl.leiden(adata_combined, resolution=1, key_added = 'cluster_lvl1' )

sc.pl.umap(adata_combined, color=['batch'],size=1)
sc.pl.umap(adata_combined, color=['cluster_lvl1'],size=1, legend_loc='on data', legend_fontsize=10)


"""
###########################################################################################################
----------------------------------------------12---------------------------------------------------------
###########################################################################################################
"""

sc.pl.umap(adata_combined, color=["cluster_lvl1", "scDblFinder_score", "scDblFinder_class"], wspace=0.5,size=3)
sc.pl.umap(adata_combined, color=["cluster_lvl1", "percent_mito", "n_counts", "n_genes"], wspace=0.5,size=3)

sc.pl.umap(adata_combined, color=cell_type_marker_lvl1["Epithelial Cells"], wspace=0.5,size=3)



"""
###########################################################################################################
----------------------------------------------13---------------------------------------------------------
###########################################################################################################
"""

path_save_adata_combined = "/mnt/cluster/environments/willmstid/Projekte/RNA_data/scRNA_GC_David_Tim3_Gal9/AnnData_objects/adata_combined.h5ad"
adata_combined.write(filename=path_save_adata_combined, compression=None, compression_opts=None, as_dense=())



