# Scanpy scRNA-seq Pipeline

> Bulk RNA-seq gives you the average expression across millions of cells. That's like measuring the average temperature of a hospital and concluding nobody has a fever. Single-cell measures each cell individually, and that changes everything.

## Why Single-Cell Matters

A tumour isn't one thing. It's a heterogeneous mix of cancer cells, immune cells, fibroblasts, and more. Bulk RNA-seq blends all of these. scRNA-seq lets you ask: Which immune cells are infiltrating? Are there rare cancer stem cells? How do cells communicate?

## What This Pipeline Does

1. **QC**: Remove dead cells (high mito), doublets, low-quality cells
2. **Normalisation**: Account for sequencing depth differences per cell
3. **HVG selection**: Keep the ~2,000 most variable genes
4. **Dimensionality reduction**: PCA then UMAP for visualisation
5. **Leiden clustering**: Group similar cells
6. **Marker genes**: Wilcoxon tests per cluster
7. **Trajectory analysis**: PAGA for developmental paths

## How to Read a UMAP Plot

Each dot is a cell. Cells close together have similar expression. Distinct clusters = distinct cell types.

But be careful: **distances between clusters on UMAP are not meaningful.** Two far-apart clusters might be very similar. UMAP preserves local structure, not global distances. Use it to identify groups, not to measure relationships between groups.

## The Most Common Mistakes

1. **Not filtering aggressively**: Dead cells form their own clusters and confuse your biology.
2. **Over-clustering**: Just because Leiden finds 30 clusters doesn't mean there are 30 cell types. Start with resolution 0.5.
3. **Calling cell types without markers**: If your "T cells" don't express CD3, something is wrong.

## Usage
```bash
python scrna_pipeline.py --input data/filtered_feature_bc_matrix.h5 --output results/
```
