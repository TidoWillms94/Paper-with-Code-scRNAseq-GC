
# PYTHONHASHSEED=0 CUDA_VISIBLE_DEVICES=3 pixi run -e gpu ipython
SEED = 0
import os
os.environ["PYTHONHASHSEED"] = str(SEED)
os.environ["CUBLAS_WORKSPACE_CONFIG"] = ":4096:8"
import random
import numpy as np
import cupy as cp
random.seed(SEED)
np.random.seed(SEED)
cp.random.seed(SEED)
import rpy2.robjects as ro
ro.r(f"set.seed({SEED})")

from pathlib import Path
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt 
from tqdm import tqdm
from joblib import Parallel, delayed
import anndata as ad
import rapids_singlecell as rsc
import infercnvpy as cnv
import scanpy.external as sce

# CUDA_VISIBLE_DEVICES=0 pixi run -e gpu ipython


sc.set_figure_params( frameon=True,color_map='Spectral_r') 
sc.verbosity = 3    
plt.rcParams["axes.grid"] = False


def to_cpu(x):
    return cp.asnumpy(x) if isinstance(x, cp.ndarray) else x

def cluster_function(adata, num_genes, num_pca, num_neighbor, min_dist, spread, resolution, path_img):

    sc.pp.highly_variable_genes(adata, n_top_genes=num_genes, flavor="seurat", batch_key="sample")

    hvg = adata[:, adata.var["highly_variable"]].copy()
    rsc.get.anndata_to_GPU(hvg)
    rsc.pp.scale(hvg, max_value=10, zero_center=True)
    rsc.pp.pca(hvg, n_comps=50, random_state=SEED)

    adata.obsm["X_pca"] = to_cpu(hvg.obsm["X_pca"])
    adata.uns["pca"] = {k: to_cpu(v) for k, v in hvg.uns["pca"].items()}
    del hvg
    cp.get_default_memory_pool().free_all_blocks()

    sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
    plt.savefig(os.path.join(path_img, "elbow_plot_log.png"), dpi=300, bbox_inches="tight")
    plt.close()

    rsc.pp.neighbors(adata, n_pcs=num_pca, n_neighbors=num_neighbor, random_state=SEED)
    rsc.tl.umap(adata, min_dist=min_dist, spread=spread, random_state=SEED)
    rsc.tl.leiden(adata, resolution=resolution, key_added="cluster_lvl1", random_state=SEED)

    return adata

marker_genes_lvl1 = {
    "Bcells":            ["MS4A1", "CD19", "BANK1", "VPREB3", "CD79A"],
    "PlasmaCells":       ["JCHAIN", "MZB1", "DERL3", "TNFRSF17"],
    "EndocrineCells":    ["CHGA", "CHGB", "SCG3", "SCG5", "TTR"],
    "EpithelialCells":   ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1", "ELF3"],
    "T&NKcells":         ["CD2", "CD3D", "CD3E", "CD3G", "CD7", "NKG7", "KLRD1"],
    "Erythrocytes":      ["ALAS2", "GYPA", "CA1", "SLC4A1", "AHSP", "HBB", "HBA1"],
    "MastCells":         ["TPSAB1", "TPSB2", "CPA3", "MS4A2", "HPGDS"],
    "MyeloidCells":      ["CD14", "CD68", "AIF1", "C1QA", "CSF1R", "FCN1"],
    "EndothelialCells":  ["PLVAP", "VWF", "CLDN5", "CDH5", "ACKR1"],
    "Fibroblasts":       ["PDGFRA", "DCN", "LUM", "DPT", "TRPA1"],
    "SmoothMuscleCells": ["ACTA2", "ACTG2", "MYH11", "DES", "CNN1"],}


marker_genes_epithilial_cells = {
    "Chief cell":    ["PGA4", "PGA3", "LIPF"],
    "GMC":           ["MUC6", "FUT9"],
    "PMC":           ["MUC5AC", "TFF1", "TFF2", "GKN1"],
    "Parietal cell": ["ATP4A", "ATP4B", "CBLIF"],             
    "Goblet cell":   ["MUC2", "ATOH1", "TFF3", "SPINK4", "CLCA1", "FCGBP"],
    "Enterocyte":    ["FABP1", "VIL1", "CDX1", "CDX2", "REG4", "KRT20"],
    "Tuft cell":     ["POU2F3", "TRPM5", "AVIL"],}


"""
marker_genes_epithilial_cells = {
    # normal gastric lineages
    "PMC":             ["MUC5AC", "TFF1", "GKN1", "GKN2"],
    "GMC":             ["MUC6", "TFF2", "FUT9"],
    "Chief":           ["PGA3", "PGA4", "PGA5", "LIPF", "CHIA"],
    "Parietal":        ["ATP4A", "ATP4B", "CBLIF"],
    "Tuft":            ["POU2F3", "TRPM5", "AVIL"],
    # intestinal metaplasia
    "IM_enterocyte":   ["CDX1", "CDX2", "CDH17", "GUCY2C", "GPA33", "FABP1", "VIL1"],
    "IM_goblet":       ["MUC2", "ATOH1", "SPINK4", "CLCA1"],
    # tumor-associated states (descriptive, not decisive)
    "Normal_lost":     ["CLDN18", "MT1G", "MT1X"],
    "Tumor_up":        ["CLDN3", "CLDN4", "CLDN7", "CEACAM5", "CEACAM6"],
    "Esophageal-like": ["KRT7", "KRT17", "LYPD2"],
    "Stem/prolif":     ["OLFM4", "LGR5", "CD44", "MKI67"],}
"""

gene_groups = {
    "IEG_stress": ["FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "ATF3", "IER2", "IER3", "DUSP1", "NR4A1", "ZFP36"],
    "HSP_stress": ["HSPA1A", "HSPA1B", "HSPH1", "DNAJB1"],
    "sex_female": ["XIST", "TSIX"],
    "sex_male":   ["RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY"],
    "s_genes" :   ["MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4", "RRM1", "UNG", "GINS2", "MCM6",
                    "CDCA7", "DTL", "PRIM1", "UHRF1", "CENPU", "HELLS", "RFC2", "RPA2", "NASP",
                    "RAD51AP1", "GMNN", "WDR76", "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2", "ATAD2",
                    "RAD51", "RRM2", "CDC45", "CDC6", "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
                    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1", "E2F8"],
    "g2m_genes": ["HMGB2", "CDK1", "NUSAP1", "UBE2C", "BIRC5", "TPX2", "TOP2A", "NDC80", "CKS2",
                    "NUF2", "CKS1B", "MKI67", "TMPO", "CENPF", "TACC3", "PIMREG", "SMC4", "CCNB2",
                    "CKAP2L", "CKAP2", "AURKB", "BUB1", "KIF11", "ANP32E", "TUBB4B", "GTSE1",
                    "KIF20B", "HJURP", "CDCA3", "JPT1", "CDC20", "TTK", "CDC25C", "KIF2C",
                    "RANGAP1", "NCAPD2", "DLGAP5", "CDCA2", "CDCA8", "ECT2", "KIF23", "HMMR",
                    "AURKA", "PSRC1", "ANLN", "LBR", "CKAP5", "CENPE", "CTCF", "NEK2", "G2E3",
                    "GAS2L3", "CBX5", "CENPA"],
    "cycling":   ["MKI67", "TOP2A"],
    "nuclear":   ["MALAT1", "NEAT1"]}


path_save_adata = Path("/mnt/cluster/datasets/TIDO_Data/scRNAseq_GC/AnnData/")
adata_filtered = sc.read_h5ad(path_save_adata/"adata_filtered_doublets_lvl1.h5ad")

"""
#########################################################################################################
---------------------------3----Tumor-Score-Calculation--ALL-TISSUE -------------------------------------
#########################################################################################################
"""

output_dir_tumorscore = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/tumorscore_tissue/res_1"
os.makedirs(output_dir_tumorscore, exist_ok=True) 

samples_tissue = ["GC01T","GC02T","GC03T1","GC03T2","GC04T", "GC05T","GC06T","GC07T","GC08T1","GC08T2","GC09T","GC10T",
                    "GC03P","GC04P","GC05P","GC06P","GC07P","GC08P","GC09P","GC10P1","GC10P2"]

adata_tissue = adata_filtered[adata_filtered.obs["sample"].isin(samples_tissue)]
adata_tissue = cluster_function(adata_tissue, num_genes=4000, num_pca=26, num_neighbor=30, min_dist=0.5, spread=2, resolution=1.0, path_img=output_dir_tumorscore)

qc_plots = [
    # QC Group 1: General & QC metrics
    ("cluster_lvl1",    dict(legend_loc="on data", legend_fontsize=8, legend_fontoutline=2)),
    ("sample",          dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("sample_type",     dict(legend_loc="right margin", legend_fontsize=10)),
    ("pct_counts_mt",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_ribo", dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_hb",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 2: Sex & Stress scores
    ("SexScore",        dict(cmap="RdBu_r", vcenter=0)),
    ("IEG_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("HSP_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Nuclear",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 3: Cell cycle metrics
    ("S_phase",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("G2M_phase",       dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Cycling",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("total_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("n_genes_by_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),]


for col, kw in qc_plots:
    sc.pl.umap(adata_tissue, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(output_dir_tumorscore, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  



cols = ["n_genes_by_counts", "total_counts", "pct_counts_mt", "Nuclear", "HSP_stress", "scDblFinder_score"]
qc = adata_tissue.obs.groupby("cluster_lvl1", observed=True)[cols].median().round(2)
print(qc.loc[["6", "8", "14", "17", "3"]])         


cell_types = ["EpithelialCells", "EndocrineCells", "EndothelialCells", "Fibroblasts", "SmoothMuscleCells"]
for ct in cell_types:
    sc.pl.umap(adata_tissue, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_tumorscore, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


cell_types = ["T&NKcells", "Erythrocytes", "MastCells", "MyeloidCells", "PlasmaCells", "Bcells"]
for ct in cell_types:
    sc.pl.umap(adata_tissue, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_tumorscore, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


#annotale clusters 
cluster_celltypes = {
    "T&NKcells":         ["3", "5", "8", "10", "11", "14", "17", "18", "21", "23", "25"],
    "Bcells":            ["16", "28"],
    "PlasmaCells":       ["4"],
    "MyeloidCells":      ["7", "9"],
    "MastCells":         ["15"],
    "EndothelialCells":  ["13"],
    "Fibroblasts":       ["1", "31"],
    "SmoothMuscleCells": ["19"],
    "EndocrineCells":    ["20", "22"],
    "EpithelialCells":   ["0", "2", "12", "24", "27", "29", "30", "26"],
    "Remove":            ["6"],}

mapping = {
    cluster: cell_type
    for cell_type, clusters in cluster_celltypes.items()
    for cluster in clusters}

adata_tissue.obs["cell_type"] = adata_tissue.obs["cluster_lvl1"].astype(str).map(mapping)
print("unmapped:", sorted(set(adata_tissue.obs["cluster_lvl1"].astype(str)) - set(mapping), key=int))
adata_tissue.obs["cell_type"] = adata_tissue.obs["cell_type"].astype("category")


adata_tissue_epi = adata_tissue[adata_tissue.obs["cell_type"] == "EpithelialCells"].copy()
adata_tissue_tnk = adata_tissue[adata_tissue.obs["cell_type"] == "T&NKcells"].copy()



##################################################################################################
#######################################################
# LVL 2 SUBCLUSTER MALIGNANT AND NON-MALIGNANT  TUMOR SCORE

output_dir_tumorscore_lvl2 = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/tumorscore_tissue/lvl2"
os.makedirs(output_dir_tumorscore_lvl2, exist_ok=True) 


malignant_genes = [
  "REG4", "TFF3", "PHGR1", "FABP1", "PRAP1", "APOA1", "CLDN4", "LEFTY1", 
  "SPINK4", "ANPEP", "CLDN7", "ALDOB", "S100A10", "IL32", "MALAT1", "LGALS4", 
  "CLDN3", "MDK", "OLFM4", "PRSS3", "LGALS3", "DMBT1", "CES2", "HLA-A", 
  "CDHR5", "FABP2", "CCT2", "PPP1R1B", "AOC1", "TSPAN7", "KRT7", "NFKBIA", 
  "ISG15", "CCL20", "TCIM", "FTL", "AGR3", "CXCL1", "RBP2", "TMC5", "PCK1", 
  "TRIM31", "CDH17", "PI3", "HES1", "LCN2", "TMEM17", "TM4SF4", "SERINC2", "MISP"]

nonmalignant_genes = [
  "REP15", "TM9SF3", "SLC12A2", "IL33", "RPS4Y1", "BRI3", "GHRL", "S100P", 
  "CXCL17", "PAPSS1", "FOXQ1", "FABP5", "CD63", "CLDN18", "MUC1", "KLF2", 
  "CTSE", "SOSTDC1", "EEF1A1", "ANXA10", "KCNE2", "TAGLN2", "TSPO", "VSIG1", 
  "AGR2", "METTL7B", "ALDH1A1", "VSIG2", "PSCA", "ALDH3A1", "PGA4", "TMSB4X", 
  "MUCL3", "SPINK1", "IGFBP2", "CA9", "GLUL", "PGA3", "TFF1", "MSMB", "CA2", 
  "CYSTM1", "RNASE1", "MUC6", "MUC5AC", "TFF2", "GKN2", "GKN1", "PGC", "LIPF"]

sc.tl.score_genes(adata_tissue_epi, malignant_genes, score_name="MalignantScore", random_state=0)
sc.tl.score_genes(adata_tissue_epi, nonmalignant_genes, score_name="NonMalignantScore", random_state=0)
adata_tissue_epi.obs["TumorScore"] = adata_tissue_epi.obs["MalignantScore"] - adata_tissue_epi.obs["NonMalignantScore"]


adata_tissue_epi = cluster_function(adata_tissue_epi, num_genes=3000, num_pca=32, num_neighbor=20, min_dist=0.5, spread=1.8, resolution=1.0, path_img=output_dir_tumorscore_lvl2)



qc_plots = [
    # QC Group 1: General & QC metrics
    ("cluster_lvl1",    dict(legend_loc="on data", legend_fontsize=8, legend_fontoutline=2)),
    ("sample",          dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("sample_type",     dict(legend_loc="right margin", legend_fontsize=10)),
    ("pct_counts_mt",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_ribo", dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_hb",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 2: Sex & Stress scores
    ("SexScore",        dict(cmap="RdBu_r", vcenter=0)),
    ("TumorScore",      dict(cmap="RdBu_r", vcenter=0)),
    ("IEG_stress",      dict(cmap="RdBu_r", vmin=0)),
    ("HSP_stress",      dict(cmap="RdBu_r", vmin=0)),
    ("Nuclear",         dict(cmap="RdBu_r", vmin=0)),
    # QC Group 3: Cell cycle metrics
    ("S_phase",         dict(cmap="RdBu_r", vmin=0)),
    ("G2M_phase",       dict(cmap="RdBu_r", vmin=0)),
    ("Cycling",         dict(cmap="RdBu_r", vmin=0)),]

for col, kw in qc_plots:
    sc.pl.umap(adata_tissue_epi, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(output_dir_tumorscore_lvl2, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory


################################################################################
# COPY NUMBER VARIATION 

output_dir_cnv = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/cnv_tissue"
os.makedirs(output_dir_cnv, exist_ok=True) 

adata_cnv_tissue = adata_tissue[~adata_tissue.obs["cell_type"].isin(["EndocrineCells", "Remove"])].copy()
adata_cnv_tissue.var = adata_cnv_tissue.var.rename(columns={"chr_start": "start", "chr_end": "end"})
adata_cnv_tissue = adata_cnv_tissue[:, adata_cnv_tissue.var["chromosome"].isin([f"chr{i}" for i in range(1, 23)])].copy()
names = adata_cnv_tissue.var_names
keep = (
    (adata_cnv_tissue.var["gene_type"] == "protein_coding")                              # also removes IG/TR genes
    & ~names.str.startswith("HLA-")                                                      # HLA region (chr6)
    & ~names.str.match(r"^(RP[SL]\d+[AXYL]?\d*|RPLP\d|RPSA)$",  case=False)              # ribosomal proteins
    & ~names.str.match(r"^HB[ABDEGMQZ]\d?$",  case=False)                                # hemoglobins
)

print("removing", (~keep).sum(), "genes; keeping", keep.sum())
adata_cnv_tissue = adata_cnv_tissue[:, keep].copy()

pd.set_option("future.infer_string", False)
adata_cnv_tissue.var_names = pd.Index(adata_cnv_tissue.var_names.to_numpy(dtype=object))
print(adata_cnv_tissue.var_names.dtype) 
reference = ['T&NKcells', 'MyeloidCells', 'Bcells', 'PlasmaCells', 'Fibroblasts',
             'SmoothMuscleCells', 'MastCells', 'EndothelialCells']
cnv.tl.infercnv(adata_cnv_tissue, reference_key="cell_type", reference_cat=reference, n_jobs=32, chunksize=625,
                layer="lognorm", window_size=250, calculate_gene_values=False)




barcodes_epi = adata_tissue_epi.obs.index.values
adata_cnv_tissue_epi = adata_cnv_tissue[adata_cnv_tissue.obs.index.isin(barcodes_epi)]
assert (adata_tissue_epi.obs_names == adata_cnv_tissue_epi.obs_names).all()

adata_tissue_epi.obsm["X_cnv"] = adata_cnv_tissue_epi.obsm["X_cnv"] 
adata_tissue_epi.uns["cnv"] = adata_cnv_tissue_epi.uns["cnv"] 

cnv.pl.chromosome_heatmap(adata_tissue_epi, groupby="cluster_lvl1")
output_file = os.path.join(output_dir_cnv, "heatmap_chromosomes_window_250.png") 
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()



adata_tissue_epi_test = adata_tissue_epi.copy()
X = adata_tissue_epi_test.obsm["X_cnv"]
X2 = X.multiply(X) if hasattr(X, "multiply") else X ** 2
adata_tissue_epi_test.obs["cnv_score_sq"] = np.asarray(X2.mean(axis=1)).ravel()
adata_tissue_epi_test.obs["cnv_score_sq_cluster"] = (adata_tissue_epi_test.obs.groupby("cluster_lvl1")["cnv_score_sq"].transform("mean"))

vmax = np.percentile(adata_tissue_epi_test.obs["cnv_score_sq_cluster"], 99) 
fig, axes = plt.subplots(1, 5, figsize=(30, 5))
plt.subplots_adjust(wspace=0.5)
sc.pl.umap(adata_tissue_epi_test, color="sample", size=3,
           legend_fontsize=8, ax=axes[0], show=False)
sc.pl.umap(adata_tissue_epi_test, color="sample_type", size=3,
           legend_fontsize=8, ax=axes[1], show=False)

sc.pl.umap(adata_tissue_epi_test, color="cluster_lvl1", size=3,
           legend_loc="on data", legend_fontsize=10, ax=axes[2], show=False)
sc.pl.umap(adata_tissue_epi_test, color="TumorScore", size=3,
           cmap="RdBu_r", vcenter=0, ax=axes[3], show=False)
sc.pl.umap(adata_tissue_epi_test, color="cnv_score_sq_cluster", size=3,
           cmap="RdBu_r", vmin=0, vmax=vmax, ax=axes[4], show=False)

output_file = os.path.join(output_dir_cnv, "umap_cluster_tumorscore_CNV_mean.png")
fig.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close(fig)



marker_genes_epithilial_cells_new =  {
    "PMC":        ["MUC5AC", "TFF1", "GKN1", "GKN2"],
    "GMC":        ["MUC6", "TFF2", "FUT9"],
    "Chief":      ["PGA3", "PGA4", "PGA5", "LIPF", "CHIA"],
    "Parietal":   ["ATP4A", "ATP4B", "CBLIF"],
    "Enterocyte": ["CDX1", "CDX2", "CDH17", "FABP1", "VIL1", "KRT20"],
    "Goblet":     ["MUC2", "ATOH1", "SPINK4", "CLCA1"],
    "Tuft":       ["POU2F3", "TRPM5", "AVIL"],}


cell_types = ["EpithelialCells", "EndocrineCells", "EndothelialCells", "Fibroblasts", "SmoothMuscleCells"]
for ct in cell_types:
    sc.pl.umap(adata_tissue_epi, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_tumorscore_lvl2, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


cell_types = ["Chief", "GMC", "PMC", "Parietal", "Goblet", "Tuft", "Enterocyte"]
for ct in cell_types:
    sc.pl.umap(adata_tissue_epi, color=marker_genes_epithilial_cells_new[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_tumorscore_lvl2, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


sc.pl.umap(adata_tissue_epi, color=["TumorScore"], size=1, legend_fontsize=10)
output_file = os.path.join(output_dir_tumorscore_lvl2, "umap_lvl1_cluster_TumorScore.png") 
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()



#annotale clusters 
epi_celltypes = {
    "Tumor_GC10":      ["9", "23"],
    "Tumor_GC08_1":    ["17"],
    "Tumor_GC08_2":    ["20", "24"],
    "Tumor_GC07":      ["16"],
    "Uncertain_GC02":  ["11"],
    "IM_enterocyte":   ["15", "6", "2"],
    "IM_goblet":       ["19"],
    "Normal_PMC":      ["21", "7", "3", "4"],
    "Normal_GMC":      ["10", "13", "14", "5"],
    "Normal_chief":    ["18", "8", "25", "1"],
    "Normal_parietal": ["22"],
    "Tuft":            ["12"],
    "Uncertain_check": ["0"],
}


mapping = {
    cluster: cell_type
    for cell_type, clusters in epi_celltypes.items()
    for cluster in clusters}

adata_tissue_epi.obs["cell_type"] = adata_tissue_epi.obs["cluster_lvl1"].astype(str).map(mapping)
print("unmapped:", sorted(set(adata_tissue_epi.obs["cluster_lvl1"].astype(str)) - set(mapping), key=int))
adata_tissue_epi.obs["cell_type"] = adata_tissue_epi.obs["cell_type"].astype("category")

sc.pl.umap(adata_tissue_epi, color=["cell_type"], size=3, legend_fontsize=7, legend_loc="on data")
output_file = os.path.join(output_dir_tumorscore_lvl2, "umap_tissue_epi_cell_type.png") 
plt.savefig(output_file, dpi=300, bbox_inches="tight")
plt.close()



genes = adata_tissue_epi.var.index.to_list()
barcodes = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs.index.to_list()
matrix = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].X.toarray()

df_tissue_epi_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_tissue_epi_analyse["sample"] = [x.split("-")[-1] for x in df_tissue_epi_analyse.index]
df_tissue_epi_analyse["cell_type"] = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["cell_type"].to_list()
df_tissue_epi_analyse["sample_type"] = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["sample_type"].to_list()

df_POI_tissue_epi = df_tissue_epi_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]




########################################################################################################################################################
# CLUSTER T and NK CELLS 



output_dir_TNK_blood = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/cluster_TNK/blood"
os.makedirs(output_dir_TNK_blood, exist_ok=True) 



adata_blood_tnk = adata_filtered[adata_filtered.obs["sample_type"] == "blood"]
adata_blood_tnk = cluster_function(adata_blood_tnk, num_genes=2000, num_pca=20, num_neighbor=30, min_dist=0.5, spread=2, resolution=1.0, path_img=output_dir_TNK_blood)


qc_plots = [
    # QC Group 1: General & QC metrics
    ("cluster_lvl1",    dict(legend_loc="on data", legend_fontsize=8, legend_fontoutline=2)),
    ("sample",          dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("sample_type",     dict(legend_loc="right margin", legend_fontsize=10)),
    ("pct_counts_mt",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_ribo", dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_hb",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 2: Sex & Stress scores
    ("SexScore",        dict(cmap="RdBu_r", vcenter=0, vmax="p99")),
    ("IEG_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("HSP_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Nuclear",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 3: Cell cycle metrics
    ("S_phase",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("G2M_phase",       dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Cycling",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("total_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("n_genes_by_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),]


for col, kw in qc_plots:
    sc.pl.umap(adata_blood_tnk, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(output_dir_TNK_blood, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory



cell_types = ["T&NKcells", "Erythrocytes", "MastCells", "MyeloidCells", "PlasmaCells", "Bcells"]
for ct in cell_types:
    sc.pl.umap(adata_blood_tnk, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_TNK_blood, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)







cluster_tnk = ["0", "1", "14","5","4","12","6","7","13","8","10"]
adata_blood_tnk = adata_blood_tnk[adata_blood_tnk.obs["cluster_lvl1"].isin(cluster_tnk)]


output_dir_TNK_blood_lvl2 = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/cluster_TNK/blood_lvl2"
os.makedirs(output_dir_TNK_blood_lvl2, exist_ok=True) 

adata_blood_tnk = cluster_function(adata_blood_tnk, num_genes=1000, num_pca=18, num_neighbor=30, min_dist=0.5, spread=2, resolution=1.0, path_img=output_dir_TNK_blood_lvl2)


tnk_markers = {
    "T":      ["CD3D", "CD3E", "CD2"],
    "CD4":    ["CD4", "CD40LG", "IL7R"],
    "CD8":    ["CD8A", "CD8B"],
    "Treg":   ["FOXP3", "IL2RA", "IKZF2", "CTLA4"],
    "NK":     ["NCR1", "KLRF1", "FCGR3A", "NCAM1", "TYROBP"],
    "gdT":    ["TRDC", "TRGC1", "TRGC2"],
    "MAIT":   ["SLC4A10", "TRAV1-2", "KLRB1"],
    "Naive":  ["CCR7", "SELL", "LEF1", "TCF7"],
    "Effector": ["GZMB", "PRF1", "GNLY", "GZMK", "CX3CR1"],
    "Cycling":  ["MKI67", "TOP2A"],}

tnk_markers = {k: [g for g in v if g in adata_blood_tnk.var_names] for k, v in tnk_markers.items()}

sc.pl.umap(adata_blood_tnk, color=["cluster_lvl1"], size=2, legend_loc="on data",
           legend_fontsize=8, wspace=0.4, show=False)
plt.savefig(os.path.join(output_dir_TNK_blood_lvl2, "umap_tnk_clusters.png"), dpi=300, bbox_inches="tight")
plt.close()

sc.tl.dendrogram(adata_blood_tnk, groupby="cluster_lvl1")

sc.pl.dotplot(adata_blood_tnk, tnk_markers, groupby="cluster_lvl1", layer="lognorm",
              standard_scale="var", dendrogram=True, show=False)
plt.savefig(os.path.join(output_dir_TNK_blood_lvl2, "dotplot_tnk_clusters.png"),
            dpi=300, bbox_inches="tight")
plt.close()




cluster_celltypes = {
    "CD4":   ["7", "11", "8", "3", "0"],
    "Treg":  ["10"],
    "CD8":   ["1", "15"],
    "NK":    ["2", "6", "13", "9"],
    "Other": ["14", "4", "5", "12"],     # MAIT (14), gamma-delta T (4, 5, 12)
}


mapping = {
    cluster: cell_type
    for cell_type, clusters in cluster_celltypes.items()
    for cluster in clusters}

adata_blood_tnk.obs["cell_type"] = adata_blood_tnk.obs["cluster_lvl1"].astype(str).map(mapping)
print("unmapped:", sorted(set(adata_blood_tnk.obs["cluster_lvl1"].astype(str)) - set(mapping), key=int))
adata_blood_tnk.obs["cell_type"] = adata_blood_tnk.obs["cell_type"].astype("category")



genes = adata_blood_tnk.var.index.to_list()
barcodes = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs.index.to_list()
matrix = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].X.toarray()

df_blood_tnk_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_blood_tnk_analyse["sample"] = [x.split("-")[-1] for x in df_blood_tnk_analyse.index]
df_blood_tnk_analyse["cell_type"] = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs["cell_type"].to_list()
df_blood_tnk_analyse["sample_type"] = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs["sample_type"].to_list()

df_POI_blood_tkn = df_blood_tnk_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]
df_POI_blood_tkn = df_POI_blood_tkn[df_POI_blood_tkn["cell_type"]!="Other"]



#############################################################################################################################################################

output_dir_TNK_tissue = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/cluster_TNK/tissue"
os.makedirs(output_dir_TNK_tissue, exist_ok=True) 
adata_tissue_tnk = adata_tissue[adata_tissue.obs["cell_type"] == "T&NKcells"].copy()
adata_tissue_tnk = cluster_function(adata_tissue_tnk, num_genes=2500, num_pca=30, num_neighbor=30, min_dist=0.4, spread=2.5, resolution=1.2 ,path_img=output_dir_TNK_tissue)


qc_plots = [
    # QC Group 1: General & QC metrics
    ("cluster_lvl1",    dict(legend_loc="on data", legend_fontsize=8, legend_fontoutline=2)),
    ("sample",          dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("sample_type",     dict(legend_loc="right margin", legend_fontsize=10)),
    ("pct_counts_mt",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_ribo", dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_hb",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 2: Sex & Stress scores
    ("SexScore",        dict(cmap="RdBu_r", vcenter=0, vmax="p99")),
    ("IEG_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("HSP_stress",      dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Nuclear",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 3: Cell cycle metrics
    ("S_phase",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("G2M_phase",       dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("Cycling",         dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("total_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("n_genes_by_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),]


for col, kw in qc_plots:
    sc.pl.umap(adata_tissue_tnk, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(output_dir_TNK_tissue, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory


cell_types = ["T&NKcells", "Erythrocytes", "MastCells", "MyeloidCells", "PlasmaCells", "Bcells"]
for ct in cell_types:
    sc.pl.umap(adata_tissue_tnk, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(output_dir_TNK_tissue, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


adata_tissue_tnk = adata_tissue_tnk[adata_tissue_tnk.obs["cluster_lvl1"] != "0"]



output_dir_TNK_tissue_lvl2 = "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/cluster_TNK/tissue_lvl2"
os.makedirs(output_dir_TNK_tissue_lvl2, exist_ok=True) 

sc.pp.highly_variable_genes(adata_tissue_tnk, n_top_genes=1500, flavor="seurat_v3")


stress_genes = [# heat shock
        "HSPA1A", "HSPA1B", "HSPA6", "HSPA8", "HSPH1", "HSP90AA1", "HSP90AB1", "HSPB1", "HSPD1", "HSPE1",
        "DNAJB1", "DNAJA1", "BAG3",
        # immediate early / dissociation
        "FOS", "FOSB", "JUN", "JUNB", "JUND", "EGR1", "EGR2", "EGR3", "ATF3", "IER2", "IER3", "IER5",
        "DUSP1", "NR4A1", "NR4A2", "NR4A3", "ZFP36", "KLF2", "KLF4", "KLF6", "PPP1R15A", "GADD45B",]

sex_genes = ["XIST", "TSIX","RPS4Y1", "DDX3Y", "KDM5D", "UTY", "EIF1AY", "USP9Y", "ZFY", "NLGN4Y", "TXLNGY", "PRKY"]



v = adata_tissue_tnk.var_names
exclude = (
    v.str.match(r"^TR[ABDG][VJ]")
    | v.str.match(r"^IG[HKL][VJ]")
    | v.str.match(r"^MT-", case=False)
    | v.str.match(r"^HB[ABDEGMQZ]\d?$", case=False)
    | v.str.match(r"^(RP[SL]\d+[AXYL]?\d*|RPLP\d|RPSA)$", case=False)
    | v.isin(sex_genes)
    | v.isin(stress_genes))
print("excluded genes:", int(exclude.sum()))

adata_sub = adata_tissue_tnk[:, ~exclude]

sc.pp.highly_variable_genes(adata_sub, n_top_genes=1500, flavor="seurat",batch_key="sample",)

adata_tissue_tnk.var["highly_variable"] = False
adata_tissue_tnk.var.loc[adata_sub.var_names, "highly_variable"] = ( adata_sub.var["highly_variable"])

print(f"HVGs: {adata_tissue_tnk.var['highly_variable'].sum()}")
print(f"Excluded genes among HVGs: {(adata_tissue_tnk.var['highly_variable'] & exclude).sum()}")

hvg = adata_tissue_tnk[:, adata_tissue_tnk.var["highly_variable"]].copy()

rsc.get.anndata_to_GPU(hvg)
rsc.pp.scale(hvg, max_value=10)
rsc.pp.pca(hvg, n_comps=50, random_state=SEED)

adata_tissue_tnk.obsm["X_pca"] = to_cpu(hvg.obsm["X_pca"])
adata_tissue_tnk.uns["pca"] = {k: to_cpu(v) for k, v in hvg.uns["pca"].items()}
del hvg
cp.get_default_memory_pool().free_all_blocks()


sc.pl.pca_variance_ratio(adata_tissue_tnk, n_pcs=50, log=True, show=False)
plt.savefig(os.path.join(output_dir_TNK_tissue_lvl2, "elbow_plot_log.png"), dpi=300, bbox_inches="tight")
plt.close()

rsc.pp.harmony_integrate(adata_tissue_tnk, key="patient", basis="X_pca",adjusted_basis="X_pca_harmony", theta=0.5)

rsc.pp.neighbors(adata_tissue_tnk, use_rep="X_pca_harmony", n_pcs=21, n_neighbors=30, random_state=SEED)
rsc.tl.umap(adata_tissue_tnk, random_state=SEED)
rsc.tl.leiden(adata_tissue_tnk, resolution=1.0, key_added="cluster_lvl1", random_state=SEED)
#sc.tl.dendrogram(adata_tissue_tnk, groupby="cluster_lvl1")


qc_plots = [
    # QC Group 1: General & QC metrics
    ("cluster_lvl1",    dict(legend_loc="on data", legend_fontsize=8, legend_fontoutline=2)),
    ("sample",          dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("sample_type",     dict(legend_loc="right margin", legend_fontsize=10)),
    ("pct_counts_mt",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_ribo", dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("pct_counts_hb",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    # QC Group 2: Sex & Stress scores
    ("SexScore",        dict(cmap="RdBu_r", vcenter=0)),
    ("IEG_stress",      dict(cmap="RdBu_r", vmin=0)),
    ("HSP_stress",      dict(cmap="RdBu_r", vmin=0)),
    ("Nuclear",         dict(cmap="RdBu_r", vmin=0)),
    # QC Group 3: Cell cycle metrics
    ("S_phase",         dict(cmap="RdBu_r", vmin=0)),
    ("G2M_phase",       dict(cmap="RdBu_r", vmin=0)),
    ("Cycling",         dict(cmap="RdBu_r", vmin=0)),
    ("total_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("n_genes_by_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),]

for col, kw in qc_plots:
    sc.pl.umap(adata_tissue_tnk, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(output_dir_TNK_tissue_lvl2, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory




tnk_markers = {
    "T":      ["CD3D", "CD3E", "CD2"],
    "CD4":    ["CD4", "CD40LG", "IL7R"],
    "CD8":    ["CD8A", "CD8B"],
    "Treg":   ["FOXP3", "IL2RA", "IKZF2", "CTLA4"],
    "NK":     ["NCR1", "KLRF1", "FCGR3A", "NCAM1", "TYROBP"],
    "gdT":    ["TRDC", "TRGC1", "TRGC2"],
    "MAIT":   ["SLC4A10", "TRAV1-2", "KLRB1"],
    "Naive":  ["CCR7", "SELL", "LEF1", "TCF7"],
    "Effector": ["GZMB", "PRF1", "GNLY", "GZMK", "CX3CR1"],
    "Cycling":  ["MKI67", "TOP2A"],}

tnk_markers = {k: [g for g in v if g in adata_tissue_tnk.var_names] for k, v in tnk_markers.items()}

sc.pl.umap(adata_tissue_tnk, color=["cluster_lvl1"], size=2, legend_loc="on data",
           legend_fontsize=8, wspace=0.4, show=False)
plt.savefig(os.path.join(output_dir_TNK_tissue_lvl2, "umap_tnk_clusters.png"), dpi=300, bbox_inches="tight")
plt.close()

sc.tl.dendrogram(adata_tissue_tnk, groupby="cluster_lvl1")
sc.pl.dotplot(adata_tissue_tnk, tnk_markers, groupby="cluster_lvl1", layer="lognorm",
              standard_scale="var", dendrogram=True, show=False)
plt.savefig(os.path.join(output_dir_TNK_tissue_lvl2, "dotplot_tnk_clusters.png"),
            dpi=300, bbox_inches="tight")
plt.close()



cluster_celltypes = {
    "CD4":    ["1", "2", "9", "11"],
    "Treg":   ["12"],
    "CD8":    ["0", "3", "4", "6", "7"],
    "NK":     ["13", "15"],
    "Other":  ["5", "8", "10"],   # 5 MAIT, 8 γδ T, 10 cycling T
    "Remove": ["14"],
}

mapping = {
    cluster: cell_type
    for cell_type, clusters in cluster_celltypes.items()
    for cluster in clusters}

adata_tissue_tnk.obs["cell_type"] = adata_tissue_tnk.obs["cluster_lvl1"].astype(str).map(mapping)
print("unmapped:", sorted(set(adata_tissue_tnk.obs["cluster_lvl1"].astype(str)) - set(mapping), key=int))
adata_tissue_tnk.obs["cell_type"] = adata_tissue_tnk.obs["cell_type"].astype("category")




genes = adata_tissue_tnk.var.index.to_list()
barcodes = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs.index.to_list()
matrix = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].X.toarray()

df_tissue_tnk_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_tissue_tnk_analyse["sample"] = [x.split("-")[-1] for x in df_tissue_tnk_analyse.index]
df_tissue_tnk_analyse["cell_type"] = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["cell_type"].to_list()
df_tissue_tnk_analyse["sample_type"] = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["sample_type"].to_list()

df_POI_tissue_tkn = df_tissue_tnk_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]
df_POI_tissue_tkn = df_POI_tissue_tkn[~df_POI_tissue_tkn["cell_type"].isin(["Other", "Remove"])]
df_POI_tissue_tkn["cell_type"].unique()









"""
#########################################################################################################
---------------------------SAVE ANNOTATED OBJECTS--------------------------------------------------------
#########################################################################################################
"""

import anndata as ad
ad.settings.allow_write_nullable_strings = False

def prepare_for_saving(adata):
    """Make an AnnData object safe for writing to h5ad."""
    adata = adata.copy()                                   # resolve views
    for col in adata.obs.columns:                          # pandas string/object columns -> category
        if adata.obs[col].dtype == object or str(adata.obs[col].dtype).startswith("str"):
            adata.obs[col] = adata.obs[col].astype(str).astype("category")
    adata.var_names = pd.Index(adata.var_names.astype(str), dtype=object)
    adata.obs_names = pd.Index(adata.obs_names.astype(str), dtype=object)
    for key in list(adata.uns.keys()):                     # drop GPU arrays left in .uns
        if isinstance(adata.uns[key], cp.ndarray):
            adata.uns[key] = to_cpu(adata.uns[key])
    return adata

objects = {
    "adata_tissue_annotated_lvl1.h5ad":   adata_tissue,       # level-1 cell types (all tissue)
    "adata_tissue_epi_annotated.h5ad":    adata_tissue_epi,   # epithelial / tumor + X_cnv
    "adata_tissue_tnk_annotated.h5ad":    adata_tissue_tnk,   # tissue T/NK subtypes
    "adata_blood_tnk_annotated.h5ad":     adata_blood_tnk,    # blood T/NK subtypes
}
for fname, obj in objects.items():
    out = prepare_for_saving(obj)
    out.write_h5ad(path_save_adata / fname, compression="gzip")
    print(f"saved {fname}: {out.n_obs} cells, cell types: {out.obs['cell_type'].nunique()}")


# ---- combined annotation table: level-1 label + final fine label per cell ----
ann = pd.DataFrame(index=adata_filtered.obs_names)
ann["sample"] = adata_filtered.obs["sample"].astype(str)
ann["sample_type"] = adata_filtered.obs["sample_type"].astype(str)
ann["cell_type_lvl1"] = adata_tissue.obs["cell_type"].astype(str).reindex(ann.index)
ann.loc[adata_blood_tnk.obs_names, "cell_type_lvl1"] = "T&NKcells"

ann["cell_type_final"] = ann["cell_type_lvl1"]
ann.loc[adata_tissue_epi.obs_names, "cell_type_final"] = adata_tissue_epi.obs["cell_type"].astype(str)
ann.loc[adata_tissue_tnk.obs_names, "cell_type_final"] = adata_tissue_tnk.obs["cell_type"].astype(str)
ann.loc[adata_blood_tnk.obs_names,  "cell_type_final"] = adata_blood_tnk.obs["cell_type"].astype(str)
ann["cell_type_final"] = ann["cell_type_final"].fillna("Unannotated")   # e.g. blood non-T/NK, removed clusters

ann.to_csv(path_save_adata / "cell_annotation_final.csv")
print(ann["cell_type_final"].value_counts())

# ---- full dataset with final labels ----
adata_filtered.obs["cell_type_lvl1"] = ann["cell_type_lvl1"].fillna("Unannotated").astype("category")
adata_filtered.obs["cell_type_final"] = ann["cell_type_final"].astype("category")
prepare_for_saving(adata_filtered).write_h5ad(path_save_adata / "adata_annotated_final.h5ad",
                                              compression="gzip")



path_load = Path("/mnt/cluster/datasets/TIDO_Data/scRNAseq_GC/AnnData")

files = {
    "adata_tissue":     "adata_tissue_annotated_lvl1.h5ad",  # level-1 cell types (all tissue)
    "adata_tissue_epi": "adata_tissue_epi_annotated.h5ad",   # epithelial / tumor + X_cnv
    "adata_tissue_tnk": "adata_tissue_tnk_annotated.h5ad",   # tissue T/NK subtypes
    "adata_blood_tnk":  "adata_blood_tnk_annotated.h5ad",    # blood T/NK subtypes
}

adatas = {name: sc.read_h5ad(path_load / fname) for name, fname in files.items()}

adata_tissue     = adatas["adata_tissue"]
adata_tissue_epi = adatas["adata_tissue_epi"]
adata_tissue_tnk = adatas["adata_tissue_tnk"]
adata_blood_tnk  = adatas["adata_blood_tnk"]




###########################################################################
genes = adata_tissue_epi.var.index.to_list()
barcodes = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs.index.to_list()
matrix = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].X.toarray()

df_tissue_epi_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_tissue_epi_analyse["sample"] = [x.split("-")[-1] for x in df_tissue_epi_analyse.index]
df_tissue_epi_analyse["cell_type"] = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["cell_type"].to_list()
df_tissue_epi_analyse["sample_type"] = adata_tissue_epi[adata_tissue_epi.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["sample_type"].to_list()

df_POI_tissue_epi = df_tissue_epi_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]

###########################################################################

genes = adata_tissue_tnk.var.index.to_list()
barcodes = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs.index.to_list()
matrix = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].X.toarray()

df_tissue_tnk_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_tissue_tnk_analyse["sample"] = [x.split("-")[-1] for x in df_tissue_tnk_analyse.index]
df_tissue_tnk_analyse["cell_type"] = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["cell_type"].to_list()
df_tissue_tnk_analyse["sample_type"] = adata_tissue_tnk[adata_tissue_tnk.obs["sample"].isin(['GC07P','GC07T','GC08P','GC08T1','GC07P', "GC10P1", 'GC10T',])].obs["sample_type"].to_list()

df_POI_tissue_tkn = df_tissue_tnk_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]
df_POI_tissue_tkn = df_POI_tissue_tkn[~df_POI_tissue_tkn["cell_type"].isin(["Other", "Remove"])]

#####################################################################

genes = adata_blood_tnk.var.index.to_list()
barcodes = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs.index.to_list()
matrix = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].X.toarray()

df_blood_tnk_analyse = pd.DataFrame(data=matrix, index=barcodes, columns=genes)
df_blood_tnk_analyse["sample"] = [x.split("-")[-1] for x in df_blood_tnk_analyse.index]
df_blood_tnk_analyse["cell_type"] = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs["cell_type"].to_list()
df_blood_tnk_analyse["sample_type"] = adata_blood_tnk[adata_blood_tnk.obs["sample"].isin(['GC06B','GC07B','GC08B','GC10B'])].obs["sample_type"].to_list()

df_POI_blood_tkn = df_blood_tnk_analyse[["sample_type", "sample", "cell_type", "TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]]
df_POI_blood_tkn = df_POI_blood_tkn[df_POI_blood_tkn["cell_type"]!="Other"]








########################################
#for df_POI_tissue_epi


coarse_map = {
    "Normal_PMC":      "Normal_epithelium",
    "Normal_GMC":      "Normal_epithelium",
    "Normal_chief":    "Normal_epithelium",
    "Normal_parietal": "Normal_epithelium",
    "Tuft":            "Normal_epithelium",
    "IM_enterocyte":   "Intestinal_metaplasia",
    "IM_goblet":       "Intestinal_metaplasia",
    "Tumor_GC07":      "Tumor",
    "Tumor_GC08_1":    "Tumor",
    "Tumor_GC08_2":    "Tumor",
    "Tumor_GC10":      "Tumor",
    "Uncertain_check": "Uncertain",}

medium_map = {**coarse_map,
              "Tumor_GC07": "Tumor_GC07",
              "Tumor_GC08_1": "Tumor_GC08",          
              "Tumor_GC08_2": "Tumor_GC08",
              "Tumor_GC10": "Tumor_GC10"}

df_POI_tissue_epi = df_POI_tissue_epi
df_POI_tissue_epi["cell_group"] = df_POI_tissue_epi["cell_type"].map(coarse_map)

df_POI_tissue_epi = df_POI_tissue_epi[~(df_POI_tissue_epi["cell_group"].isin(["Other","Intestinal_metaplasia", "Uncertain"]))]

df_POI_tissue_epi["cell_group"].unique()
df_POI_tissue_tkn["cell_type"].unique()
df_POI_blood_tkn["cell_type"].unique()




output_dir_publication= "/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/plot_publication"
os.makedirs(output_dir_publication, exist_ok=True) 



import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

PATIENTS = ["GC07", "GC08", "GC10"]
GENES = ["TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]
CMAPS = ["Oranges", "Greens", "Blues", "Purples", "Reds", "Greys"]
MIN_CELLS = 10
SIZE = 300

# ---------- 1. combine the three tables ----------
blood = df_POI_blood_tkn.assign(compartment="Blood", cell=df_POI_blood_tkn["cell_type"])
tnk   = df_POI_tissue_tkn.assign(compartment=df_POI_tissue_tkn["sample_type"].str.capitalize(),cell=df_POI_tissue_tkn["cell_type"])
epi   = df_POI_tissue_epi.assign(compartment=df_POI_tissue_epi["sample_type"].str.capitalize(),cell=df_POI_tissue_epi["cell_group"].map({"Normal_epithelium": "Epithelial cells","Tumor": "Tumor cells"}))
df = pd.concat([blood, tnk, epi], ignore_index=True)

df["cell"] = df["cell"].replace({"NK": "NK cells", "CD4": "CD4⁺ T cells","CD8": "CD8⁺ T cells", "Treg": "Treg cells"})
df["patient"] = df["sample"].astype(str).str.extract(r"(GC\d+)")[0]
df = df[df["patient"].isin(PATIENTS)]


rows = ([("Blood", c) for c in ["NK cells", "CD4⁺ T cells", "CD8⁺ T cells", "Treg cells"]] +
        [("Paratumor", c) for c in ["Epithelial cells", "NK cells", "CD4⁺ T cells",
                                    "CD8⁺ T cells", "Treg cells"]] +
        [("Tumor", c) for c in ["Epithelial cells", "NK cells", "CD4⁺ T cells", "CD8⁺ T cells",
                                "Treg cells", "Tumor cells"]])


# ---------- 2. statistics per group ----------
records = []
for comp, cell in rows:
    base = df[(df["compartment"] == comp) & (df["cell"] == cell)]
    for gene in GENES:
        for p in PATIENTS + ["Total"]:
            sub = base if p == "Total" else base[base["patient"] == p]
            if len(sub) < MIN_CELLS:
                m = f = np.nan
            else:
                v = sub[gene].values
                m, f = v.mean(), (v > 0).mean()
            records.append(dict(comp=comp, cell=cell, gene=gene, col=p, mean=m, frac=f))
st = pd.DataFrame(records)
st["mean_scaled"] = st.groupby("gene")["mean"].transform(lambda x: x / x.max())

# ---------- 3. x positions ----------
xpos, x, panel_ends, panel_centers = {}, 0, [], []
for gene in GENES:
    xs = []
    for p in PATIENTS:
        xpos[(gene, p)] = x; xs.append(x); x += 1
    panel_centers.append(np.mean(xs)); panel_ends.append(x - 0.5 + 0.3)
    x += 0.6
total_start = x
for gene in GENES:
    xpos[(gene, "Total")] = x; x += 1
ypos = {r: i for i, r in enumerate(rows)}

# ---------- 4. plot ----------
fig, ax = plt.subplots(figsize=(22, 7))
for gene, cmap in zip(GENES, CMAPS):
    s = st[st["gene"] == gene].dropna()
    ax.scatter([xpos[(gene, c)] for c in s["col"]],
               [ypos[(a, b)] for a, b in zip(s["comp"], s["cell"])],
               s=s["frac"] * SIZE, c=s["mean_scaled"], cmap=cmap,
               vmin=0, vmax=1, edgecolors="grey", linewidths=0.4)

for xe in panel_ends:
    ax.axvline(xe, color="black", lw=0.8)
for i in [3.5, 8.5]:
    ax.axhline(i, color="lightgrey", lw=0.6)

ax.set_xticks([xpos[(g, p)] for g in GENES for p in PATIENTS] + [xpos[(g, "Total")] for g in GENES])
ax.set_xticklabels(PATIENTS * len(GENES) + GENES, rotation=90)
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([c for _, c in rows])
ax.invert_yaxis()
ax.set_xlim(-0.6, x - 0.4)

# gene names above the per-patient panels, "Total" above the last panel
for gene, xc in zip(GENES, panel_centers):
    ax.text(xc, -0.9, gene, ha="center", va="bottom", fontsize=11, style="italic")
ax.text(total_start + (len(GENES) - 1) / 2, -0.9, "Total", ha="center", va="bottom", fontsize=11)

# compartment brackets
for comp in ["Blood", "Paratumor", "Tumor"]:
    idx = [i for i, (c, _) in enumerate(rows) if c == comp]
    ax.annotate("", xy=(-0.11, idx[0] - 0.4), xycoords=("axes fraction", "data"),
                xytext=(-0.11, idx[-1] + 0.4), arrowprops=dict(arrowstyle="-", lw=1.2))
    ax.text(-0.125, np.mean(idx), comp, transform=ax.get_yaxis_transform(),
            rotation=90, va="center", ha="center", fontsize=11)

# legends: dot size + one colorbar per gene
for f in [0.15, 0.30, 0.45, 0.60]:
    ax.scatter([], [], s=f * SIZE, c="grey", label=f"{int(f*100)}")
ax.legend(title="Fraction of cells\nin group (%)", loc="upper left",
          bbox_to_anchor=(1.01, 1.0), frameon=False, labelspacing=1.2)
for gi, (gene, cmap) in enumerate(zip(GENES, CMAPS)):
    cax = fig.add_axes([0.92, 0.50 - gi * 0.075, 0.07, 0.018])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(f"$\\it{{{gene}}}$", fontsize=8, labelpad=1)
    cb.set_ticks([0, 0.5, 1]); cb.ax.tick_params(labelsize=7)

fig.savefig(os.path.join(output_dir_publication, "dotplot_TIGIT_ligands_per_patient.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)








import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

PATIENTS = ["GC07", "GC08", "GC10"]
GENES = ["TIGIT", "PVR", "NECTIN1", "NECTIN2", "NECTIN3", "NECTIN4"]
CMAPS = ["Oranges", "Greens", "Blues", "Purples", "Reds", "Greys"]
MIN_CELLS = 50           
MIN_PATIENTS_TOTAL = 2  
SIZE = 300


def truncate_cmap(name, lo=0.2):
    """Drop the near-white bottom of a sequential colormap so low values stay visible."""
    base = plt.get_cmap(name)
    return LinearSegmentedColormap.from_list(f"{name}_trunc", base(np.linspace(lo, 1, 256)))


CMAPS = [truncate_cmap(c) for c in CMAPS]

# ---------- 1. combine the three tables ----------
blood = df_POI_blood_tkn.assign(compartment="Blood", cell=df_POI_blood_tkn["cell_type"])
tnk = df_POI_tissue_tkn.assign(compartment=df_POI_tissue_tkn["sample_type"].str.capitalize(),
                               cell=df_POI_tissue_tkn["cell_type"])
epi = df_POI_tissue_epi.assign(compartment=df_POI_tissue_epi["sample_type"].str.capitalize(),
                               cell=df_POI_tissue_epi["cell_group"].map(
                                   {"Normal_epithelium": "Epithelial cells", "Tumor": "Tumor cells"}))
df = pd.concat([blood, tnk, epi], ignore_index=True)

df["cell"] = df["cell"].replace({"NK": "NK cells", "CD4": "CD4⁺ T cells",
                                 "CD8": "CD8⁺ T cells", "Treg": "Treg cells"})
df["patient"] = df["sample"].astype(str).str.extract(r"(GC\d+)")[0]
df = df[df["patient"].isin(PATIENTS)]

rows = ([("Blood", c) for c in ["NK cells", "CD4⁺ T cells", "CD8⁺ T cells", "Treg cells"]] +
        [("Paratumor", c) for c in ["Epithelial cells", "NK cells", "CD4⁺ T cells",
                                    "CD8⁺ T cells", "Treg cells"]] +
        [("Tumor", c) for c in ["Epithelial cells", "NK cells", "CD4⁺ T cells", "CD8⁺ T cells",
                                "Treg cells", "Tumor cells"]])

# report anything that is not plotted (e.g. CNV-high cells in paratumor)
counts = df.groupby(["compartment", "cell", "patient"]).size().unstack(fill_value=0)
excluded = counts.loc[~counts.index.isin(pd.MultiIndex.from_tuples(rows))]
if len(excluded):
    print("Cells not shown in the plot:\n", excluded)

# ---------- 2. statistics per group ----------
# Per patient: mean expression and fraction of expressing cells (log-normalised data).
# Total: mean of the per-patient values, so every patient is weighted equally
# (pooling cells would let the patient with the most cells dominate).
records = []
for comp, cell in rows:
    base = df[(df["compartment"] == comp) & (df["cell"] == cell)]
    for gene in GENES:
        pm, pf = [], []
        for p in PATIENTS:
            sub = base[base["patient"] == p]
            if len(sub) < MIN_CELLS:
                m = f = np.nan
            else:
                v = sub[gene].to_numpy()
                m, f = v.mean(), (v > 0).mean()
                pm.append(m); pf.append(f)
            records.append(dict(comp=comp, cell=cell, gene=gene, col=p, mean=m, frac=f, n=len(sub)))
        if len(pm) >= MIN_PATIENTS_TOTAL:
            m, f = np.mean(pm), np.mean(pf)
        else:
            m = f = np.nan
        records.append(dict(comp=comp, cell=cell, gene=gene, col="Total", mean=m, frac=f, n=len(pm)))
st = pd.DataFrame(records)

# scale mean per gene to its maximum; guard against genes never detected (max = 0)
st["mean_scaled"] = st.groupby("gene")["mean"].transform(
    lambda x: x / x.max() if x.max() > 0 else x * 0)

# ---------- 3. x positions ----------
xpos, x, panel_ends, panel_centers = {}, 0, [], []
for gene in GENES:
    xs = []
    for p in PATIENTS:
        xpos[(gene, p)] = x; xs.append(x); x += 1
    panel_centers.append(np.mean(xs)); panel_ends.append(x - 0.5 + 0.3)
    x += 0.6
total_start = x
for gene in GENES:
    xpos[(gene, "Total")] = x; x += 1
ypos = {r: i for i, r in enumerate(rows)}

# ---------- 4. plot ----------
fig, ax = plt.subplots(figsize=(22, 7))
for gene, cmap in zip(GENES, CMAPS):
    s = st[(st["gene"] == gene) & st["mean"].notna()]
    ax.scatter([xpos[(gene, c)] for c in s["col"]],
               [ypos[(a, b)] for a, b in zip(s["comp"], s["cell"])],
               s=s["frac"] * SIZE, c=s["mean_scaled"], cmap=cmap,
               vmin=0, vmax=1, edgecolors="dimgrey", linewidths=0.4)

# groups with too few cells: small grey cross, so they are not confused with zero expression
miss = st[st["mean"].isna()]
ax.scatter([xpos[(g, c)] for g, c in zip(miss["gene"], miss["col"])],
           [ypos[(a, b)] for a, b in zip(miss["comp"], miss["cell"])],
           marker="x", s=18, c="lightgrey", linewidths=0.8)

for xe in panel_ends:
    ax.axvline(xe, color="black", lw=0.8)
for i in [3.5, 8.5]:
    ax.axhline(i, color="lightgrey", lw=0.6)

ax.set_xticks([xpos[(g, p)] for g in GENES for p in PATIENTS] + [xpos[(g, "Total")] for g in GENES])
ax.set_xticklabels(PATIENTS * len(GENES) + GENES, rotation=90)
for lab in ax.get_xticklabels()[-len(GENES):]:
    lab.set_fontstyle("italic")
ax.set_yticks(range(len(rows)))
ax.set_yticklabels([c for _, c in rows])
ax.set_ylim(len(rows) - 0.5, -0.5)   # inverted y-axis
ax.set_xlim(-0.6, x - 0.4)

# gene names above the per-patient panels, "Total" above the last panel
for gene, xc in zip(GENES, panel_centers):
    ax.text(xc, -0.7, gene, ha="center", va="bottom", fontsize=11, style="italic")
ax.text(total_start + (len(GENES) - 1) / 2, -0.7, "Total (mean of patients)",
        ha="center", va="bottom", fontsize=11)

# compartment brackets
for comp in ["Blood", "Paratumor", "Tumor"]:
    idx = [i for i, (c, _) in enumerate(rows) if c == comp]
    ax.annotate("", xy=(-0.11, idx[0] - 0.4), xycoords=("axes fraction", "data"),
                xytext=(-0.11, idx[-1] + 0.4), arrowprops=dict(arrowstyle="-", lw=1.2))
    ax.text(-0.125, np.mean(idx), comp, transform=ax.get_yaxis_transform(),
            rotation=90, va="center", ha="center", fontsize=11)

# legends: colorbars top right, dot size + missing marker bottom right
pos = ax.get_position()
cb_left, cb_w, cb_h, cb_step = pos.x1 + 0.015, 0.06, 0.018, 0.07
for gi, (gene, cmap) in enumerate(zip(GENES, CMAPS)):
    cax = fig.add_axes([cb_left, pos.y1 - cb_h - gi * cb_step, cb_w, cb_h])
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(0, 1))
    cb = fig.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label(f"$\\it{{{gene}}}$", fontsize=8, labelpad=2)
    cb.ax.xaxis.set_label_position("top")   # gene name above its own bar
    cb.set_ticks([0, 0.5, 1]); cb.ax.tick_params(labelsize=7)
fig.text(cb_left, pos.y1 + 0.045, "Mean expression\n(scaled to max per gene)", fontsize=8, va="bottom")

handles = [ax.scatter([], [], s=f * SIZE, facecolor="white", edgecolor="black", lw=0.6)
           for f in [0.2, 0.4, 0.6, 0.8]]
labels = ["20", "40", "60", "80"]
handles.append(ax.scatter([], [], marker="x", s=18, c="lightgrey", linewidths=0.8))
labels.append(f"< {MIN_CELLS} cells")
ax.legend(handles, labels, title="Fraction of cells\nin group (%)", loc="lower left",
          bbox_to_anchor=(1.01, 0.0), frameon=False, labelspacing=1.2, fontsize=8, title_fontsize=8)

fig.savefig(os.path.join(output_dir_publication, "dotplot_TIGIT_ligands_per_patient.png"),
            dpi=300, bbox_inches="tight")
plt.close(fig)
