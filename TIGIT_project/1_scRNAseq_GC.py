
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

from rpy2.robjects.conversion import localconverter
import anndata2ri

ro.r("""
suppressPackageStartupMessages({
    library(Seurat)
    library(scater)
    library(scDblFinder)
    library(BiocParallel)
    library(future)
    library(future.apply)})
  """)


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

    #sc.pl.umap(adata, color=["cluster_lvl1"], size=3, legend_loc="on data", legend_fontsize=10)
    #output_file = os.path.join(path_img, "umap_lvl1_cluster.png") 
    #plt.savefig(output_file, dpi=300, bbox_inches="tight")
    #plt.close()

    #sc.pl.umap(adata, color=["sample"], size=1, legend_fontsize=10)
    #output_file = os.path.join(path_img, "umap_lvl1_cluster_batch.png") 
    #plt.savefig(output_file, dpi=300, bbox_inches="tight")
    #plt.close()

    #sc.pl.umap(adata, color=["sample_type"], size=1, legend_fontsize=10)
    #output_file = os.path.join(path_img, "umap_lvl1_cluster_sample_type.png") 
    #plt.savefig(output_file, dpi=300, bbox_inches="tight")
    #plt.close()

    return adata


# GET REFERENCE TABLE FOR GENE NAMES AND LOCATIONS
gtf_path = "/mnt/cluster/datasets/TIDO_Data/scRNAseq_GC/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz"
cols = ["chromosome", "source", "feature", "start", "end", "score", "strand", "frame", "attributes"]
gtf = pd.read_csv(gtf_path, sep="\t", comment="#", header=None, names=cols,
                  usecols=["chromosome", "feature", "start", "end", "strand", "attributes"],
                  dtype={"chromosome": str})

gtf["gene_id"] = gtf["attributes"].str.extract(r'gene_id "([^"]+)"', expand=False)
gtf["gene_type"] = gtf["attributes"].str.extract(r'gene_type "([^"]+)"', expand=False)
gtf["gene_name"] = gtf["attributes"].str.extract(r'gene_name "([^"]+)"', expand=False)

# Strand is used for RNA velocity 
gtf_gene = gtf[gtf["feature"] == "gene"][["gene_name", "gene_type", 'chromosome', 'start', 'end','gene_id']]
gtf_gene = gtf_gene.set_index("gene_id")



dir_samples = [
#Tumortissue
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207224/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207225/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207228/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207229/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207231/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207233/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207236/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207239/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207242/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207243/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207245/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207249/outs/filtered_feature_bc_matrix",
#Paratumor Tissue
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207227/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207230/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207232/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207235/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207238/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207241/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207244/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207247/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207248/outs/filtered_feature_bc_matrix",
#Blood
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207234/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207237/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207240/outs/filtered_feature_bc_matrix",
"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/data/HRR207246/outs/filtered_feature_bc_matrix",]

sample_names = [
    "GC01T","GC02T","GC03T1","GC03T2","GC04T", "GC05T","GC06T","GC07T","GC08T1","GC08T2","GC09T","GC10T",
    "GC03P","GC04P","GC05P","GC06P","GC07P","GC08P","GC09P","GC10P1","GC10P2",
    "GC06B","GC07B","GC08B","GC10B"             
        ]
   

"""
###########################################################################################################
----------------------------------------------- ---------------------------------------------------------
###########################################################################################################
"""


qc_dir = Path("/mnt/cluster/environments/willmstid/GitHub/PDAC_TME/figures/QC_lvl1")
os.makedirs(qc_dir, exist_ok=True) 

def data_filtering(adata, sample_name, gtf_gene, min_genes, min_counts, mt_pct, qc_dir):

    qc_dir_samples = qc_dir/sample_name
    os.makedirs(qc_dir_samples, exist_ok=True)

    adata.var_names_make_unique()   # adds  -1 -2 to gene names is duplets                      
    adata.var["gene_type"] = adata.var["gene_ids"].map(gtf_gene["gene_type"])
    adata.var["chromosome"] = adata.var["gene_ids"].map(gtf_gene["chromosome"])
    adata.var["chr_start"] = adata.var["gene_ids"].map(gtf_gene["start"])
    adata.var["chr_end"] = adata.var["gene_ids"].map(gtf_gene["end"])
    keep = adata.var["chromosome"].astype(str).str.startswith("chr")
    adata = adata[:, keep].copy()
    
    adata.var["mt"] = adata.var_names.str.match(r"^MT-", case=False)
    adata.var["ribo"] = adata.var_names.str.match(r"^(RP[SL]\d+[AXYL]?\d*|RPLP\d|RPSA)$",  case=False)  #r (^...$) forces an exact, full match for regex
    adata.var["hb"] = adata.var_names.str.match(r"^HB[ABDEGMQZ]\d?$",  case=False)
    adata.var["nuclear"] = adata.var_names.str.match("^(MALAT1|NEAT1)$")

    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt", "ribo", "hb", "nuclear"], percent_top=None, log1p=False, inplace=True)

    adata = adata[
        (adata.obs.n_genes_by_counts > min_genes) &
        (adata.obs.total_counts > min_counts) &
        (adata.obs.pct_counts_mt < mt_pct)].copy()

    g = sc.pl.violin(adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"],jitter=0.4, multi_panel=True, show=False)
    g.figure.savefig(qc_dir_samples / f"{sample_name}_violin_QC.png", dpi=150, bbox_inches="tight")
    g = sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", show=False)
    g.figure.savefig(qc_dir_samples /f"_{sample_name}_counts_vs_mito.png", dpi=150, bbox_inches="tight")
    g = sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", show=False)
    g.figure.savefig(qc_dir_samples /f"_{sample_name}_counts_vs_genes.png", dpi=150, bbox_inches="tight")
    return adata


def parse_sample_type(batch_name):
    if "T" in batch_name:
        return "tumor"
    elif "B" in batch_name:
        return "blood"
    elif "P" in batch_name:
        return "paratumor"
    else:
        return "unknown"


def read_samples(dir_sample, sample_name):
    #sample_name = sample_names[0]
    #dir_sample = dir_samples[0]
    adata = sc.read_10x_mtx(dir_sample, var_names='gene_symbols', cache=True)
    adata.layers["raw_counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4, exclude_highly_expressed=True, max_fraction=0.1)
    sc.pp.log1p(adata, base=2)
    adata.layers["lognorm"] = adata.X.copy()
    adata.X = adata.layers["raw_counts"].copy()
    adata.obs['sample'] = sample_name
    adata.obs["sample_type"] = adata.obs["sample"].astype(str).apply(parse_sample_type).astype("category")
    adata.obs["patient"] = adata.obs["sample"].astype(str).str.extract(r"(GC\d+)")[0]

    adata = data_filtering(adata, sample_name, gtf_gene, min_genes=200, min_counts=400, mt_pct=30, qc_dir=qc_dir)

    return adata


adata_list = Parallel(n_jobs=25)(delayed(read_samples)(dir_sample, sample_name) 
    for dir_sample, sample_name in tqdm(zip(dir_samples, sample_names)))


#DOUPLET DETECTION
# copy raw counts (genes x cells) of each sample into an R list
j = len(adata_list)
ro.r("data_array_list <- list()")
for i, adata in enumerate(tqdm(adata_list), start=1):
    counts = adata.layers["raw_counts"].T.tocsc()
    with localconverter(ro.default_converter + anndata2ri.converter):
        ro.globalenv["data_array"] = counts
    ro.globalenv["cell_names"] = ro.StrVector(adata.obs_names.tolist())
    ro.globalenv["i"] = i
    ro.r("colnames(data_array) <- cell_names; data_array_list[[i]] <- data_array")

print(list(ro.r("ls()")))                           # all objects in R's workspace
print(ro.r("length(data_array_list)")[0])          # number of samples
print(list(ro.r("dim(data_array_list[[1]])")))     # [genes, cells] of sample 1
print(ro.r("class(data_array_list[[1]])")[0])      # should be dgCMatrix


ro.r(f"""
plan(multisession, workers = 25)
results <- future_lapply(data_array_list, function(obj) {{
    sce <- scDblFinder::scDblFinder(SingleCellExperiment::SingleCellExperiment(list(counts = obj)))
    list(cells = colnames(sce),
         score = sce$scDblFinder.score,
         class = as.character(sce$scDblFinder.class))
}}, future.seed = {SEED}, future.packages = c("scDblFinder", "SingleCellExperiment"))
plan(sequential)
""")


# bring results back and add them to each AnnData object
results = ro.globalenv["results"]
for adata, res in zip(adata_list, results):
    cells = list(res.rx2("cells"))
    assert cells == adata.obs_names.tolist(), "cell order mismatch"
    adata.obs["scDblFinder_score"] = np.array(res.rx2("score"))
    adata.obs["scDblFinder_class"] = list(res.rx2("class"))


adata_all = ad.concat(adata_list, join="outer", label="sample", keys=sample_names, index_unique="-", fill_value=0, merge="same",)

path_save_adata = Path("/mnt/cluster/datasets/TIDO_Data/scRNAseq_GC/AnnData/")
adata_all.write_h5ad(path_save_adata/"adata_all.h5ad", compression="gzip")


#######################################################################################
# QC CLUSTERS lvl 1 all 

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




qc_dir_cluster_lvl1 = qc_dir/"clusters_lvl1"
os.makedirs(qc_dir_cluster_lvl1, exist_ok=True) 

adata_all.X = adata_all.layers["lognorm"]
adata_all = cluster_function(adata_all, num_genes=6000, num_pca=30, num_neighbor=20, min_dist=0.5, spread=2, resolution=1, path_img=qc_dir_cluster_lvl1)



qc_dir_cluster_lvl1_cell_types = qc_dir_cluster_lvl1/"cell_types"
os.makedirs(qc_dir_cluster_lvl1_cell_types, exist_ok=True) 

cell_types = ["T&NKcells", "Erythrocytes", "MastCells", "MyeloidCells", "PlasmaCells", "Bcells"]
for ct in cell_types:
    sc.pl.umap(adata_all, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(qc_dir_cluster_lvl1_cell_types, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)

cell_types = ["EpithelialCells", "EndocrineCells", "EndothelialCells", "Fibroblasts", "SmoothMuscleCells"]
for ct in cell_types:
    sc.pl.umap(adata_all, color=marker_genes_lvl1[ct], wspace=0.5, size=3, show=False)
    fig = plt.gcf()
    fig.suptitle(f"{ct} markers", fontsize=16, y=1.02)
    fig.savefig(os.path.join(qc_dir_cluster_lvl1_cell_types, f"umap_{ct}.png"),
                dpi=300, bbox_inches="tight")
    plt.close(fig)


sc.tl.score_genes(adata_all, [g for g in gene_groups["IEG_stress"] if g in adata_all.var_names], score_name="IEG_stress")
sc.tl.score_genes(adata_all, [g for g in gene_groups["HSP_stress"] if g in adata_all.var_names], score_name="HSP_stress")
sc.tl.score_genes(adata_all, [g for g in gene_groups["s_genes"] if g in adata_all.var_names], score_name="S_phase")
sc.tl.score_genes(adata_all, [g for g in gene_groups["g2m_genes"] if g in adata_all.var_names], score_name="G2M_phase")
sc.tl.score_genes(adata_all, [g for g in gene_groups["cycling"] if g in adata_all.var_names], score_name="Cycling")
sc.tl.score_genes(adata_all, [g for g in gene_groups["nuclear"] if g in adata_all.var_names], score_name="Nuclear")
sc.tl.score_genes(adata_all, [g for g in gene_groups["sex_male"] if g in adata_all.var_names], score_name="MaleScore")
sc.tl.score_genes(adata_all, [g for g in gene_groups["sex_female"] if g in adata_all.var_names], score_name="FemaleScore")
adata_all.obs["SexScore"] = adata_all.obs["MaleScore"] - adata_all.obs["FemaleScore"]



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
    ("Cycling",         dict(cmap="RdBu_r", vmin=0)),]

for col, kw in qc_plots:
    sc.pl.umap(adata_all, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(qc_dir_cluster_lvl1, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory

##################################################################################################
#Removing douples via TRA and TRB 


def get_Tcell_doublets_and_subtypes(VDJ_samples):
    list_barcodes_singles = []
    list_barcodes_remove = []

    for sample, path in VDJ_samples.items():
        print(f"Processing sample: {sample}")
        df = pd.read_csv(path)

        # Keep productive AND non-productive rearrangements (both occupy one allele)
        valid_contigs = df[
            (df["is_cell"].astype(str).str.lower() == "true")
            & (df["high_confidence"].astype(str).str.lower() == "true")
            & (df["chain"].isin(["TRA", "TRB", "TRG", "TRD"]))].copy()

        # UMI filtering (removes ambient background contigs)
        valid_contigs["umis"] = valid_contigs["umis"].astype(int)
        valid_contigs = valid_contigs[valid_contigs["umis"] >= 2]
        max_umi_per_chain = valid_contigs.groupby(["barcode", "chain"])["umis"].transform("max")
        valid_contigs = valid_contigs[valid_contigs["umis"] >= 0.10 * max_umi_per_chain].copy()

        # One entry per distinct rearrangement (V + J + CDR3), not per contig
        valid_contigs["rearrangement"] = (valid_contigs["v_gene"].astype(str) + "|"
                                          + valid_contigs["j_gene"].astype(str) + "|"
                                          + valid_contigs["cdr3_nt"].astype(str))

        chain_counts = (valid_contigs.groupby(["barcode", "chain"])["rearrangement"]
                        .nunique().unstack(fill_value=0))
        for c in ["TRA", "TRB", "TRG", "TRD"]:
            if c not in chain_counts.columns:
                chain_counts[c] = 0

        chain_counts["has_ab"] = (chain_counts["TRA"] + chain_counts["TRB"]) > 0
        chain_counts["has_gd"] = (chain_counts["TRG"] + chain_counts["TRD"]) > 0

        def classify_tcr(row):
            tra = row["TRA"]
            trb = row["TRB"]
            trg = row["TRG"]
            trd = row["TRD"]

            # Rule 1: 3+ rearrangements at any locus exceeds two alleles -> doublet
            if tra >= 3 or trb >= 3 or trg >= 3 or trd >= 3:
                return "doublet"

            # Rule 2: Mixed ab + gd chains (kept as single cell)
            if row["has_ab"] and row["has_gd"]:
                return "single"

            # Rule 3: Alpha/Beta
            if row["has_ab"]:
                # Two TRB rearrangements (e.g. 1 productive + 1 non-productive allele)
                if trb == 2 and tra <= 2:
                    return "single"
                # Standard paired
                if trb == 1 and tra in [1, 2]:
                    return "single"
                # Incomplete capture / single chain
                if (trb == 1 and tra == 0) or (trb == 0 and tra in [1, 2]):
                    return "single"

            # Rule 4: Gamma/Delta
            if row["has_gd"]:
                # Standard paired
                if trg in [1, 2] and trd in [1, 2]:
                    return "single"
                # Incomplete capture
                if (trg in [1, 2] and trd == 0) or (trg == 0 and trd in [1, 2]):
                    return "single"

            return "other"

        chain_counts["call"] = chain_counts.apply(classify_tcr, axis=1)
        print(chain_counts["call"].value_counts().to_dict())

        sample_singles = (chain_counts[chain_counts["call"] == "single"].index + f"_{sample}").tolist()
        sample_remove = (chain_counts[chain_counts["call"].isin(["doublet", "other"])].index + f"-{sample}").tolist()

        list_barcodes_singles.extend(sample_singles)
        list_barcodes_remove.extend(sample_remove)

    return list_barcodes_remove, list_barcodes_singles


VDJ_samples = {
    "GC01T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207251/outs/all_contig_annotations.csv",
    "GC02T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207252/outs/all_contig_annotations.csv",
    "GC04P" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207254/outs/all_contig_annotations.csv",
    "GC04T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207255/outs/all_contig_annotations.csv",
    "GC05P" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207256/outs/all_contig_annotations.csv",
    "GC05T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207257/outs/all_contig_annotations.csv",
    "GC06B" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207258/outs/all_contig_annotations.csv",
    "GC06P" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207259/outs/all_contig_annotations.csv",
    "GC06T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207260/outs/all_contig_annotations.csv",
    "GC07B" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207261/outs/all_contig_annotations.csv",
    #"GC07P" :"/mnt/cluster/environments/willmstid/Projekte/RNA_data/scGC/VDJ/output/HRR207262/outs/all_contig_annotations.csv",  #No alignment posssible
    "GC07T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207263/outs/all_contig_annotations.csv",
    "GC08B" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207264/outs/all_contig_annotations.csv",
    "GC08P" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207265/outs/all_contig_annotations.csv",
    "GC08T1" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207266/outs/all_contig_annotations.csv",
    "GC08T2" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207267/outs/all_contig_annotations.csv",
    "GC09P" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207268/outs/all_contig_annotations.csv",
    "GC09T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207269/outs/all_contig_annotations.csv",
    "GC10B" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207270/outs/all_contig_annotations.csv",
    "GC10P1" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207271/outs/all_contig_annotations.csv",
    "GC10P2" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207272/outs/all_contig_annotations.csv",
    "GC10T" :"/mnt/nct-zfs/TCO-Test/willmstid/data_2/Projekte/RNA_data/scGC/VDJ/output/HRR207273/outs/all_contig_annotations.csv",}


list_barcodes_remove, list_barcodes_Tcells = get_Tcell_doublets_and_subtypes(VDJ_samples)


adata_all.obs["doublets_VDJ"] = adata_all.obs.index.isin(list_barcodes_remove)
sum(adata_all.obs["doublets_VDJ"])


qc_plots = [
    # QC Group 1: General & QC metrics
    ("scDblFinder_class",dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("doublets_VDJ",     dict(legend_loc="right margin", legend_fontsize=8, legend_fontoutline=2)),
    ("total_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),
    ("n_genes_by_counts",   dict(cmap="RdBu_r", vmin=0, vmax="p99")),]

for col, kw in qc_plots:
    sc.pl.umap(adata_all, color=col,  size=1, title=col,show=False, **kw)
    output_file = os.path.join(qc_dir_cluster_lvl1, f"umap_{col}.png")
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()  # Close the plot to free memory



adata_filtered = adata_all[~adata_all.obs["cluster_lvl1"].isin(["17", "31", "36"])]
adata_filtered = adata_filtered[~adata_filtered.obs["doublets_VDJ"]==True]
adata_filtered = adata_filtered[~((adata_filtered.obs["scDblFinder_class"] == "doublet") & (adata_filtered.obs["cluster_lvl1"].astype(str) != "23"))]

path_save_adata = Path("/mnt/cluster/datasets/TIDO_Data/scRNAseq_GC/AnnData/")
adata_filtered.write_h5ad(path_save_adata/"adata_filtered_doublets_lvl1.h5ad", compression="gzip")
