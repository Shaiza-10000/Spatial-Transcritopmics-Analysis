# Visium Fluorescence Image Analysis with Squidpy

This notebook extends the standard spatial transcriptomics workflow by bringing **fluorescence image data** into the analysis alongside gene expression. Instead of treating the tissue image as a passive background, it actively extracts quantitative features from it — cell counts, intensity patterns, texture — and asks whether those image-derived features tell the same biological story as the gene expression clusters.

The key question this notebook answers: **can you recover biologically meaningful tissue structure from the image alone, without looking at any gene expression data?**

---

## What You Need

```bash
pip install anndata scanpy==1.11.0 squidpy python-igraph leidenalg
```

This notebook uses [Squidpy](https://squidpy.readthedocs.io/) — a library built on top of Scanpy specifically for spatial data, with dedicated tools for image processing, segmentation, and feature extraction.

---

## The Dataset

The dataset is a **cropped Visium fluorescence** section, loaded directly from Squidpy's built-in datasets. Unlike the standard Visium workflow that uses H&E staining, this section uses **fluorescence imaging** — where specific molecules are tagged with fluorescent dyes and imaged across separate channels. Each channel captures a different biological signal (e.g., DAPI for nuclei, another channel for a protein of interest).

The dataset comes pre-processed with cluster labels already assigned from gene expression, which serves as the ground truth for comparison later.

---

## Step-by-Step Breakdown

### Setup and Data Loading

```python
import anndata as ad
import scanpy as sc
import squidpy as sq

img = sq.datasets.visium_fluo_image_crop()
adata = sq.datasets.visium_fluo_adata_crop()
```

Two objects are loaded here and it is important to understand what each one is:

- `adata` — the standard AnnData object containing gene expression counts, spot metadata, and pre-computed cluster labels. This is what you would have after running a normal Scanpy pipeline.
- `img` — an `ImageContainer` object, which is Squidpy's way of storing the multi-channel fluorescence image alongside spatial coordinates. It knows how the image aligns to the spots in `adata`.

Having both objects linked together is what makes the image feature analysis possible.

---

### Visualizing Gene Expression Clusters on Tissue

```python
sq.pl.spatial_scatter(adata, color="cluster")
```

The first plot maps the gene-expression-derived clusters onto the tissue. This is the **reference** — the biologically grounded picture of cell type organization that the rest of the notebook will try to reproduce using only image features. Keep this in mind as you go through the analysis. At the end, you compare back to it.

---

### Exploring the Fluorescence Image

```python
img.show(channelwise=True)
```

Displays each fluorescence channel separately. This is an important sanity check before doing any analysis — you want to confirm that the channels look as expected, that there is no obvious imaging artifact, and that the signal is distributed across the tissue in a way that makes biological sense. A channel that is uniformly bright or uniformly dark is not going to contribute useful features downstream.

---

### Image Smoothing and Cell Segmentation

```python
sq.im.process(img=img, layer="image", method="smooth")

sq.im.segment(img=img, layer="image_smooth", method="watershed", channel=0, chunks=1000)
```

Before segmenting cells, the image is smoothed. Smoothing reduces high-frequency noise — small pixel-level variation that has nothing to do with biology — so that the segmentation algorithm can find real cell boundaries rather than chasing noise.

**Watershed segmentation** is then applied to the smoothed image. Watershed is a classic image segmentation algorithm that treats pixel intensities like a topographic map — it "floods" the image from local intensity minima and places boundaries where the flooding from different minima would meet. In fluorescence microscopy, where nuclei (stained with DAPI, for example) appear as bright islands on a dark background, watershed reliably separates individual cells.

The crop-and-compare plot that follows shows a 500×500 pixel region of the tissue before and after segmentation — you can visually confirm that the algorithm is correctly outlining individual cells.

---

### Extracting Segmentation Features per Spot

```python
features_kwargs = {"segmentation": {"label_layer": "segmented_watershed"}}
sq.im.calculate_image_features(
    adata, img,
    features="segmentation",
    layer="image",
    key_added="features_segmentation",
    n_jobs=1,
    features_kwargs=features_kwargs,
)
```

Now that individual cells are segmented, Squidpy counts how many cells fall within each Visium spot and measures the mean fluorescence intensity per channel within those segmented cells. These become per-spot features stored in `adata.obsm["features_segmentation"]`.

This is biologically meaningful because cell density varies across tissue regions — germinal centers are densely packed, while stromal areas are sparse. If the segmentation is working correctly, cell count per spot should correlate with known cell-dense regions in your gene expression clusters.

The four-panel plot compares:

| Panel | What it shows |
|---|---|
| `segmentation_label` | Number of segmented cells per spot |
| `cluster` | Gene expression cluster (ground truth) |
| `segmentation_ch-0_mean_intensity_mean` | Mean fluorescence of channel 0 per spot |
| `segmentation_ch-1_mean_intensity_mean` | Mean fluorescence of channel 1 per spot |

If the intensity maps and cell density map visually resemble the cluster map, image and transcriptomics are telling the same story.

---

### Extracting Richer Image Features

```python
params = {
    "features_orig": {
        "features": ["summary", "texture", "histogram"],
        "scale": 1.0,
        "mask_circle": True,
    },
    "features_context": {"features": ["summary", "histogram"], "scale": 1.0},
    "features_lowres": {"features": ["summary", "histogram"], "scale": 0.25},
}

for feature_name, cur_params in params.items():
    sq.im.calculate_image_features(adata, img, layer="image", key_added=feature_name, n_jobs=1, **cur_params)
```

This is the most feature-rich part of the notebook. Three different feature extraction strategies are run, each capturing a different aspect of the image:

**features_orig** — extracts summary statistics, texture features, and pixel intensity histograms from the exact circular area underneath each Visium spot (`mask_circle=True`). This is the most precise — only pixels that the spot actually captures are used.

**features_context** — same features but without the circular mask, so it includes some surrounding tissue context. Useful when a biological signal (like a diffuse extracellular matrix pattern) extends beyond the spot boundary.

**features_lowres** — uses a downscaled version of the image (`scale=0.25`), which captures broader spatial patterns at the cost of fine detail. A cluster of cells looks like a blob at low resolution, which can sometimes be more informative than pixel-level noise.

The three feature types are:

- **Summary** — mean, standard deviation, percentiles of pixel intensities within the spot. Captures overall brightness and spread.
- **Histogram** — distribution of pixel intensity values binned into buckets. More detailed than summary stats — captures whether a spot has a bimodal intensity distribution (two cell populations?) or a uniform one.
- **Texture** — features derived from the Grey-Level Co-occurrence Matrix (GLCM), which describe how pixel intensities relate to their neighbors. Captures whether a region looks grainy, smooth, striated, or patchy — structural properties invisible to intensity-only features.

All three sets are concatenated into a single feature matrix stored in `adata.obsm["features"]`.

---

### Image-Based Clustering

```python
def cluster_features(features: pd.DataFrame, like=None):
    if like is not None:
        features = features.filter(like=like)
    adata = ad.AnnData(features)
    sc.pp.scale(adata)
    sc.pp.pca(adata, n_comps=min(10, features.shape[1] - 1))
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata)
    return adata.obs["leiden"]

adata.obs["features_summary_cluster"] = cluster_features(adata.obsm["features"], like="summary")
adata.obs["features_histogram_cluster"] = cluster_features(adata.obsm["features"], like="histogram")
adata.obs["features_texture_cluster"] = cluster_features(adata.obsm["features"], like="texture")
```

The `cluster_features` function runs the standard Scanpy clustering pipeline — scale, PCA, neighbors, Leiden — but on image features instead of gene expression. Each feature type is clustered independently so you can see which one best captures tissue organization.

One detail worth noting: the features are **scaled before PCA**. This is important because image features have very different numerical ranges — a mean intensity might be 0.4 while a texture contrast value might be 120. Without scaling, PCA would be dominated by whichever feature has the largest absolute values, regardless of biological relevance.

---

### Comparing Image Clusters to Gene Expression Clusters

```python
sq.pl.spatial_scatter(
    adata,
    color=["features_summary_cluster", "features_histogram_cluster", "features_texture_cluster", "cluster"],
    ncols=3,
)
```

The final plot puts everything side by side — three image-derived clusterings and the gene expression clustering. This is the payoff of the whole notebook.

If the image-based clusters spatially resemble the gene expression clusters, it means the tissue's visual appearance carries real biological information. A region that looks different under fluorescence microscopy genuinely has different cell types or states — and vice versa. This kind of cross-modality validation is valuable because it builds confidence that the clusters represent real biology rather than computational artifacts.

In practice, no image clustering will perfectly match gene expression — the modalities capture different things. But partial agreement, especially in anatomically distinct regions, is a meaningful result.

---

## Key Concepts

| Concept | Description |
|---|---|
| `ImageContainer` | Squidpy's object for storing multi-channel tissue images aligned to spatial spots |
| Watershed segmentation | Algorithm that detects cell boundaries by treating image intensity as a topographic map |
| GLCM texture features | Describe structural patterns in an image — graininess, smoothness, directionality |
| Histogram features | Per-spot distribution of pixel intensities across intensity bins |
| Summary features | Per-spot statistics: mean, std, percentiles of pixel intensities |
| `mask_circle` | Restricts feature extraction to the circular capture area of each Visium spot |
| Multi-scale analysis | Extracting features at different image resolutions to capture both fine and coarse structure |

---

## How This Differs from the Standard Scanpy Workflow

| Standard Scanpy (Visium) | This Notebook |
|---|---|
| Gene expression matrix is the input | Fluorescence image features are the input |
| Clusters reflect transcriptional similarity | Clusters reflect visual/morphological similarity |
| H&E image used as background only | Fluorescence image actively mined for features |
| One modality | Two modalities compared and cross-validated |

---

## Common Pitfalls

**Skipping the per-channel visualization.** Always look at `img.show(channelwise=True)` before extracting features. A saturated or empty channel will produce useless features and misleading clusters.

**Not scaling before PCA.** Image features span very different numerical ranges. Forgetting `sc.pp.scale()` before running PCA on image features is one of the most common mistakes — it will let one feature type dominate the entire embedding.

**Expecting a perfect match between image and gene clusters.** The two modalities are complementary, not identical. Gene expression captures cell identity; image features capture morphology and protein levels. Partial overlap is the norm and is still informative.

---

## References

- Palla G. et al. **Squidpy: a scalable framework for spatial omics analysis.** *Nature Methods* 19, 171–178 (2022). https://doi.org/10.1038/s41592-021-01358-2
- Wolf F.A. et al. **SCANPY: large-scale single-cell gene expression data analysis.** *Genome Biology* 19, 15 (2018). https://doi.org/10.1186/s13059-017-1382-0
- Squidpy documentation: https://squidpy.readthedocs.io
- Visium Fluorescence dataset: https://squidpy.readthedocs.io/en/stable/api.html#squidpy.datasets.visium_fluo_adata_crop
