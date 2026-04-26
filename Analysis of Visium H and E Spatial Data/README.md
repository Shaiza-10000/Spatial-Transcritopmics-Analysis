# Analysis of Visium H&E Spatial Transcriptomics Data

## Overview

This notebook performs **spatial transcriptomics analysis** on a 10x Genomics Visium dataset paired with a Hematoxylin & Eosin (H&E) stained tissue image. The workflow integrates gene expression data with tissue morphology to uncover spatially organized biological patterns in mouse brain tissue, particularly the **hippocampus** and surrounding cortical regions.

Spatial transcriptomics allows researchers to measure gene expression while preserving the physical location of cells within a tissue, bridging the gap between histology (what the tissue looks like) and genomics (what genes are active).

---

## Biological Context

### What is Visium?
[10x Genomics Visium](https://www.10xgenomics.com/spatial-transcriptomics) captures mRNA from tissue sections placed on a slide containing spatially barcoded spots (~55 µm diameter each). Each spot captures transcripts from the cells above it, giving both gene expression and spatial coordinates simultaneously.

### What is H&E staining?
Hematoxylin & Eosin (H&E) is the standard histological stain used in pathology. **Hematoxylin** stains nuclei blue/purple; **Eosin** stains cytoplasm and extracellular matrix pink. The resulting image reveals tissue architecture (cell density, layer structure, morphology) that complements molecular data.

### Why Mouse Brain?
The mouse hippocampus and cortex are well-characterized brain regions with distinct laminar organization (layers of cell types). This makes them ideal for validating spatial methods, expected biology (e.g., pyramidal neuron layers, dentate gyrus) should be recoverable from the data.

---

## Dependencies

| Package | Role |
|---|---|
| `scanpy` | Single-cell / transcriptomics analysis (clustering, PCA, preprocessing) |
| `squidpy` | Spatial transcriptomics analysis (image features, spatial graphs, ligand-receptor) |
| `anndata` | Data structure holding expression matrix + metadata |
| `pandas` | Tabular data manipulation |
| `numpy` | Numerical computation |
| `leidenalg` | Community detection algorithm for clustering |

Install all dependencies:
```bash
pip install anndata scanpy squidpy leidenalg
```

---

## Analysis Workflow

### 1. Data Loading

```python
img = sq.datasets.visium_hne_image()
adata = sq.datasets.visium_hne_adata()
```

**What happens:** Loads a pre-processed `AnnData` object (`adata`) containing the gene expression matrix and spatial coordinates, alongside the high-resolution H&E tissue image (`img`).

**Biology:** The `adata` object already has cells assigned to clusters (brain regions) from prior gene-expression-based clustering (e.g., Hippocampus, Pyramidal_layer, Pyramidal_layer_dentate_gyrus, Cortex, etc.).

---

### 2. Spatial Cluster Visualization

```python
sq.pl.spatial_scatter(adata, color="cluster")
```

**What happens:** Overlays gene-expression-derived cluster labels onto the H&E image. Each spot on the tissue is colored by its assigned cell/region type.

**Biology:** Validates that transcriptomic clusters match expected tissue anatomy. For example, hippocampal spots should spatially co-localize with the hippocampal structure visible in the H&E image.

---

### 3. Image Feature Extraction (Morphology Analysis)

```python
for scale in [1.0, 2.0]:
    sq.im.calculate_image_features(
        adata, img.compute(),
        features="summary",
        key_added=f"features_summary_scale{scale}",
        n_jobs=4, scale=scale,
    )
```

**What happens:** For each Visium spot, extracts **summary statistics** (mean, standard deviation, percentiles) of pixel intensities from the H&E image at two zoom levels (scales). This converts raw image patches into numerical feature vectors — one per spot.

**Biology:** H&E pixel statistics encode morphological information: dense nuclei (dark patches) appear in neuron-rich layers, while white matter or sparse regions look different. Scale 1.0 captures local texture; scale 2.0 captures broader tissue context. Together, they allow the tissue image itself to "explain" spatial organization, independently of gene expression.

```python
adata.obsm["features"] = pd.concat([...])
```
Combines multi-scale features into a single feature matrix stored per spot.

---

### 4. Morphology-Based Clustering

```python
def cluster_features(features, like=None):
    adata = ad.AnnData(features)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    return adata.obs["leiden"]

adata.obs["features_cluster"] = cluster_features(adata.obsm["features"], like="summary")
sq.pl.spatial_scatter(adata, color=["features_cluster", "cluster"])
```

**What happens:**
1. **Scaling** - normalizes image features to zero mean and unit variance (critical since pixel statistics have very different ranges).
2. **PCA** - reduces feature dimensionality, capturing the major axes of morphological variation.
3. **Neighborhood graph** - connects spots with similar morphological profiles.
4. **Leiden clustering** - community detection on the graph, grouping spots with similar H&E appearance.
5. **Side-by-side visualization** - compares morphology-based clusters (`features_cluster`) with gene-expression-based clusters (`cluster`).

**Biology:** If H&E-derived morphological clusters agree with gene-expression clusters, it confirms that tissue architecture and molecular identity co-vary - a fundamental principle of tissue biology. Discordances may reveal cases where cells look similar but are molecularly distinct (or vice versa), which is scientifically interesting.

---

### 5. Spatial Neighborhood Enrichment

```python
sq.gr.spatial_neighbors(adata)
sq.gr.nhood_enrichment(adata, cluster_key="cluster")
sq.pl.nhood_enrichment(adata, cluster_key="cluster")
```

**What happens:**
1. **`spatial_neighbors`** - builds a spatial graph where spots that are physically adjacent on the tissue are connected.
2. **`nhood_enrichment`** - for each pair of cluster types, tests whether they are found neighboring each other more (enriched) or less (depleted) than expected by chance, using permutation testing.
3. **Heatmap** — displays the enrichment Z-scores between all cluster pairs.

**Biology:** Reveals which brain regions are spatially adjacent in the tissue. For example, the Pyramidal layer and Hippocampus should show strong neighborhood enrichment because they are anatomically contiguous. This is essentially a data-driven reconstruction of brain region borders from spatial proximity.

---

### 6. Co-occurrence Analysis

```python
sq.gr.co_occurrence(adata, cluster_key="cluster")
sq.pl.co_occurrence(
    adata, cluster_key="cluster",
    clusters="Hippocampus", figsize=(8, 4),
)
```

**What happens:** Calculates how the probability of observing a given cluster changes as a function of **distance** from the Hippocampus cluster. Unlike neighborhood enrichment (which is binary - adjacent or not), co-occurrence gives a continuous, distance-resolved view of spatial relationships.

**Biology:** Shows the spatial "sphere of influence" of hippocampal tissue. Clusters that peak at short distances are immediate neighbors; those peaking further away reflect more distal organizational relationships. This captures the layered organization of the hippocampus and dentate gyrus across spatial scales.

---

### 7. Ligand-Receptor Interaction Analysis

```python
sq.gr.ligrec(adata, n_perms=100, cluster_key="cluster")
sq.pl.ligrec(
    adata, cluster_key="cluster",
    source_groups="Hippocampus",
    target_groups=["Pyramidal_layer", "Pyramidal_layer_dentate_gyrus"],
    means_range=(3, np.inf),
    alpha=1e-4,
    swap_axes=True,
)
```

**What happens:** Tests for statistically significant **ligand-receptor (L-R) interactions** between pairs of clusters using permutation testing (100 permutations). Filters results to show only high-expression interactions (`means_range=(3, np.inf)`) that are highly significant (`alpha=1e-4`). Plots a dot plot of significant L-R pairs between the Hippocampus and its neighboring pyramidal layers.

**Biology:** Cell-to-cell communication relies on ligands (signaling molecules) secreted by one cell binding to receptors on another. By linking L-R analysis to spatial data, we can identify **which signaling pathways operate at specific tissue interfaces** - e.g., between hippocampal interneurons and pyramidal neurons. This goes beyond simple co-expression and tests whether the molecular machinery for communication is present in physically adjacent regions.

---

### 8. Spatial Autocorrelation (Moran's I)

```python
genes = adata[:, adata.var.highly_variable].var_names.values[:1000]
sq.gr.spatial_autocorr(
    adata, mode="moran",
    genes=genes, n_perms=100, n_jobs=1,
)
```

**What happens:** Computes **Moran's I** statistic for up to 1,000 highly variable genes. Moran's I measures whether a gene's expression is spatially autocorrelated — i.e., whether nearby spots tend to have similar expression levels. Uses permutation testing (100 permutations) to assess significance.

**Biology:**
- **High Moran's I (close to +1):** Gene expression is spatially clustered — the gene is a marker of a specific anatomical region or layer. These are the most biologically interpretable spatially variable genes.
- **Moran's I near 0:** Gene expression is randomly distributed across the tissue — likely housekeeping genes with no spatial pattern.
- **Negative Moran's I:** Expression is spatially dispersed (rare in practice).

Identifying spatially variable genes is foundational for understanding which molecular programs are spatially organized in the tissue, and for discovering region-specific markers beyond what is known from bulk RNA-seq.

---

## Key Data Structures

| Object | Description |
|---|---|
| `adata` | `AnnData` object: gene expression matrix (spots × genes) + spatial coordinates + cluster labels + computed features |
| `img` | `ImageContainer`: high-resolution H&E tissue image used for morphology feature extraction |
| `adata.obsm["features"]` | Per-spot morphological feature matrix extracted from H&E image |
| `adata.obs["features_cluster"]` | H&E morphology-based cluster assignments |
| `adata.obs["cluster"]` | Gene-expression-based cluster assignments (pre-computed) |
| `adata.uns["leiden_nhood_enrichment"]` | Neighborhood enrichment results |
| `adata.uns["leiden_co_occurrence"]` | Co-occurrence scores across distance bands |
| `adata.uns["leiden_ligrec"]` | Ligand-receptor interaction test results |
| `adata.uns["moranI"]` | Moran's I spatial autocorrelation results per gene |

---

## Analysis Summary

```
Load Data
    │
    ▼
Visualize Gene-Expression Clusters on Tissue
    │
    ▼
Extract H&E Image Features (multi-scale)
    │
    ▼
Cluster by Morphology → Compare to Gene Clusters
    │
    ▼
Build Spatial Graph → Neighborhood Enrichment
    │
    ▼
Distance-Resolved Co-occurrence (Hippocampus focus)
    │
    ▼
Ligand-Receptor Communication (Hippocampus → Pyramidal layers)
    │
    ▼
Moran's I Spatial Autocorrelation (top 1000 HVGs)
```

---

## Key Biological Findings This Workflow Enables

1. **Tissue architecture validation** - Confirms that transcriptomic clusters map to anatomically expected regions.
2. **Morphology–expression concordance** - Tests whether H&E appearance and gene expression agree in their tissue segmentation.
3. **Spatial organization of brain regions** = Identifies which regions border each other and at what distances.
4. **Cell-cell communication at region interfaces** = Pinpoints active signaling between hippocampal and pyramidal layer neurons.
5. **Spatially variable genes** - Discovers which genes are most spatially patterned, providing leads for region-specific biology.

---

## References

- Palla et al. (2022). *Squidpy: a scalable framework for spatial omics analysis.* **Nature Methods**. https://doi.org/10.1038/s41592-021-01358-2
- Wolf et al. (2018). *SCANPY: large-scale single-cell gene expression data analysis.* **Genome Biology**. https://doi.org/10.1186/s13059-017-1382-0
- 10x Genomics Visium Documentation: https://www.10xgenomics.com/spatial-transcriptomics
- Squidpy Tutorials: https://squidpy.readthedocs.io

