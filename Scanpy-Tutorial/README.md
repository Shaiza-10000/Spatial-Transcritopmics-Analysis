# Spatial Transcriptomics with Scanpy

This notebook walks through a complete spatial transcriptomics analysis pipeline — from raw count data to biologically meaningful maps of gene expression across tissue. It covers two technologies: **10x Genomics Visium** and **MERFISH**, each with its own strengths and quirks.

If you are new to spatial transcriptomics, think of it this way: regular single-cell RNA sequencing tells you *what* genes a cell expresses, but throws away *where* that cell was sitting in the tissue. Spatial transcriptomics keeps that location information. The result is something closer to a map than a list — you can see which cell types cluster near blood vessels, which genes are active at the tissue edge, and how cell identity changes across anatomy.

---

## What You Need

```bash
pip install scanpy igraph
```

The notebook also relies on `matplotlib`, `pandas`, `seaborn`, and `openpyxl`. The MERFISH section requires two supplementary files from the original paper — details in that section.

---

## The Dataset

The main dataset is a **human lymph node** captured with 10x Visium. The lymph node is a good test case for spatial analysis because it has well-defined anatomical zones — germinal centers packed with B cells, T cell zones, sinuses — so clustering results can be sanity-checked against known biology.

The second dataset is a **mouse brain** section profiled with MERFISH, which operates at true single-cell resolution (as opposed to Visium's ~55 µm spots that may contain a few cells each).

---

## Step-by-Step Breakdown

### Environment Setup

```python
sc.logging.print_versions()
sc.set_figure_params(facecolor="white", figsize=(8, 8))
sc.settings.verbosity = 3
```

Nothing special here, just making sure figure backgrounds are clean for export and that Scanpy logs its progress. Verbosity 3 is useful during long steps like neighbor graph construction; it tells you what's happening instead of going silent for two minutes.

---

### Loading the Data

```python
adata = sc.datasets.visium_sge(sample_id="V1_Human_Lymph_Node")
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True)
```

`sc.datasets.visium_sge()` downloads the dataset directly, no manual file handling needed. After loading, two things happen immediately:

- Gene names are deduplicated (`var_names_make_unique`) because some genes appear under multiple probe IDs in Visium data.
- Mitochondrial genes are flagged. Genes that start with `MT-` encode proteins in the mitochondria. A spot with a very high proportion of mitochondrial reads is usually a sign of a dying or damaged cell, the nuclear RNA has leaked out but the mitochondria are still intact enough to be sequenced.

`calculate_qc_metrics` computes per-spot statistics (total counts, number of genes, % mitochondrial) that you will use to filter in the next step.

---

### Quality Control Plots

```python
fig, axs = plt.subplots(1, 4, figsize=(15, 4))
sns.histplot(adata.obs["total_counts"], ...)
sns.histplot(adata.obs["total_counts"][adata.obs["total_counts"] < 10000], ...)
sns.histplot(adata.obs["n_genes_by_counts"], ...)
sns.histplot(adata.obs["n_genes_by_counts"][adata.obs["n_genes_by_counts"] < 4000], ...)
```

Before filtering anything, you look at the distributions. Each panel answers a specific question:

| Panel | Question it answers |
|---|---|
| Total counts (full range) | Is there a long tail of very high-count spots? |
| Total counts (< 10,000) | Where does the low-count cliff actually fall? |
| Genes per spot (full range) | Do most spots detect a reasonable number of genes? |
| Genes per spot (< 4,000) | What does the bottom of the distribution look like? |

The zoomed-in panels are the important ones. A cliff or gap in the histogram is what you are looking for — it tells you where "real tissue" ends and "background" begins.

---

### Filtering

```python
sc.pp.filter_cells(adata, min_counts=5000)
sc.pp.filter_cells(adata, max_counts=35000)
adata = adata[adata.obs["pct_counts_mt"] < 20].copy()
sc.pp.filter_genes(adata, min_cells=10)
```

Four filters, each with a biological rationale:

- **min_counts = 5,000** — spots below this are likely off-tissue or have poor capture
- **max_counts = 35,000** — spots above this may be ink spots, folded tissue, or other artifacts
- **pct_counts_mt < 20%** — removes spots where too much of the signal is coming from mitochondria (stressed or dying tissue)
- **min_cells = 10 per gene** — drops genes seen in fewer than 10 spots; these are too sparse to contribute anything meaningful to clustering

These thresholds are not universal. They should always be chosen by looking at your QC plots, not copied blindly from a tutorial.

---

### Normalization and Feature Selection

```python
sc.pp.normalize_total(adata, inplace=True)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, flavor="seurat", n_top_genes=2000)
```

Raw counts are noisy and incomparable across spots — one spot might have 8,000 total counts and another 30,000, which has nothing to do with biology. Normalization fixes that by scaling each spot to the same total.

The log transformation then compresses the dynamic range. A gene expressed at 10 counts vs. 1,000 counts is not actually 100× more biologically interesting, the log scale reflects that.

Finally, 2,000 highly variable genes (HVGs) are selected. These are genes that vary meaningfully across spots, not just noisy low-expressed genes, but genes that genuinely differ between cell types or regions. Everything downstream uses only these 2,000 genes.

---

### Dimensionality Reduction and Clustering

```python
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata, key_added="clusters", flavor="igraph", directed=False, n_iterations=2)
```

This is the analytical core of the pipeline. Four steps in sequence:

**PCA** reduces the 2,000 HVGs down to a manageable number of principal components — linear combinations of genes that capture the most variance. Think of it as compression that preserves the important differences between spots.

**Neighbor graph** connects each spot to its most similar spots in PCA space. This graph is the data structure that both UMAP and Leiden operate on.

**UMAP** takes the neighbor graph and creates a 2D layout where similar spots appear close together. It is a visualization tool, the distances in UMAP space are not exact, but the groupings are meaningful.

**Leiden clustering** finds communities in the neighbor graph. Unlike k-means, you do not need to specify the number of clusters in advance, the algorithm finds natural groupings in the data. The results are stored as `adata.obs["clusters"]`.

---

### UMAP Visualization

```python
sc.pl.umap(adata, color=["total_counts", "n_genes_by_counts", "clusters"], wspace=0.4)
```

Three UMAP panels side by side. The first two are QC checks, if `total_counts` or `n_genes_by_counts` strongly separates the UMAP, that is a warning sign that technical variation is driving your clusters rather than biology. The third shows the Leiden clusters, which ideally should look like coherent islands.

---

### Spatial Plots — Projecting Results Back to Tissue

```python
sc.pl.spatial(adata, img_key="hires", color=["total_counts", "n_genes_by_counts"])
sc.pl.spatial(adata, img_key="hires", color="clusters", size=1.5)
sc.pl.spatial(
    adata, img_key="hires", color="clusters",
    groups=["0", "6"], crop_coord=[7000, 10000, 0, 6000],
    alpha=0.5, size=1.3,
)
```

This is where spatial transcriptomics earns its name. Instead of plotting in UMAP space, you plot on top of the actual H&E tissue image. Each Visium spot is placed at its physical location.

The cropped view isolates clusters 0 and 6 in a specific tissue region, useful when you want to examine a particular anatomical structure without the full tissue overwhelming the view.

---

### Differential Expression and Marker Genes

```python
sc.tl.rank_genes_groups(adata, "clusters", method="t-test")
sc.pl.rank_genes_groups_heatmap(adata, groups="9", n_genes=10, groupby="clusters")
```

For each cluster, Scanpy runs a one-vs-rest t-test to find genes that are significantly more expressed in that cluster than in all others combined. These marker genes are what you use to identify *what cell type* a cluster represents.

The heatmap shows the top 10 marker genes for cluster 9 and their expression across every cluster. If cluster 9 lights up and all others stay dark, you have a clean marker. If the expression is scattered, that cluster may need re-examination.

---

### Gene Expression on Tissue

```python
sc.pl.spatial(adata, img_key="hires", color=["clusters", "CR2"])
sc.pl.spatial(adata, img_key="hires", color=["COL1A2", "SYPL1"], alpha=0.7)
```

Individual genes plotted on tissue. `CR2` (Complement Receptor 2) is a well-established B cell marker, in the lymph node it should concentrate in follicles and germinal centers. If your B cell cluster spatially overlaps with high `CR2` expression, your clustering is capturing real biology.

`COL1A2` is a collagen gene expressed by fibroblasts and stromal cells. `SYPL1` marks specific compartments in lymphoid tissue. These serve as additional validation landmarks.

---

### MERFISH Dataset

```python
coordinates = pd.read_excel("./pnas.1912459116.sd15.xlsx", index_col=0)
counts = sc.read_csv("./pnas.1912459116.sd12.csv").transpose()

adata_merfish = counts[coordinates.index, :].copy()
adata_merfish.obsm["spatial"] = coordinates.to_numpy()
```

This section loads data from [Moffitt et al., Science 2018](https://doi.org/10.1126/science.aau5324) — a landmark MERFISH study of the mouse hypothalamus. The two supplementary files need to be downloaded manually from the paper.

Unlike Visium, MERFISH captures individual cells rather than spots that average over a few cells. The coordinates go directly into `adata.obsm["spatial"]`, which is the standard slot Scanpy uses for spatial coordinates.

The preprocessing pipeline is the same in structure but with two differences:

- **15 PCA components** instead of more — MERFISH panels are smaller (hundreds of genes vs. the full transcriptome), so fewer components are appropriate
- **Leiden resolution 0.5** — a lower value produces fewer, broader clusters, which suits the coarser gene panel

```python
sc.pp.normalize_per_cell(adata_merfish, counts_per_cell_after=1e6)
sc.pp.log1p(adata_merfish)
sc.pp.pca(adata_merfish, n_comps=15)
sc.pp.neighbors(adata_merfish)
sc.tl.umap(adata_merfish)
sc.tl.leiden(adata_merfish, key_added="clusters", resolution=0.5, n_iterations=2, flavor="igraph", directed=False)
```

---

### MERFISH Visualization

```python
sc.pl.umap(adata_merfish, color="clusters")
sc.pl.embedding(adata_merfish, basis="spatial", color="clusters")
```

Two views of the same clustering. The UMAP shows which cell types are transcriptionally similar. The spatial embedding shows where those cells physically sit in the mouse brain section. In a well-structured dataset like this, you expect to see anatomically coherent regions — neurons of the same type grouped together, not scattered randomly.

---

## Common Pitfalls

**Copying QC thresholds from tutorials.** The filters in this notebook (5,000–35,000 counts, < 20% MT) are reasonable starting points for this dataset. For a different tissue or sequencing depth, they may be completely wrong. Always look at your own QC histograms first.

**Treating UMAP distances as meaningful.** UMAP preserves local structure (what is near what) but not global distances (how far cluster A is from cluster B). Two clusters that look far apart in UMAP may be more similar than two that look close.

**Ignoring the spatial plots.** A clustering that looks clean in UMAP may reveal a spatial artifact when projected onto tissue — for example, a cluster that perfectly follows a tissue fold or an ink mark. Always cross-check.

---

## References

- Wolf F.A. et al. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology* 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0
- Moffitt J.R. et al. **Molecular, spatial, and functional single-cell profiling of the hypothalamic preoptic region.** *Science* 362 (2018). https://doi.org/10.1126/science.aau5324
- Scanpy documentation: https://scanpy.readthedocs.io
- 10x Genomics Visium: https://www.10xgenomics.com/spatial-transcriptomics
