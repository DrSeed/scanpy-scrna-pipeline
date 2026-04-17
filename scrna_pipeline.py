#!/usr/bin/env python3
import scanpy as sc, argparse
from pathlib import Path
import warnings; warnings.filterwarnings('ignore')
sc.settings.verbosity = 2

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', required=True)
    parser.add_argument('--output', default='results')
    parser.add_argument('--resolution', type=float, default=0.8)
    args = parser.parse_args()
    output_dir = Path(args.output); output_dir.mkdir(parents=True, exist_ok=True)
    adata = sc.read_10x_h5(args.input) if args.input.endswith('.h5') else sc.read_10x_mtx(args.input)
    adata.var_names_make_unique()
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], inplace=True)
    adata = adata[adata.obs.pct_counts_mt < 20, :]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
    adata.raw = adata
    adata = adata[:, adata.var.highly_variable]
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata)
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)
    sc.tl.leiden(adata, resolution=args.resolution)
    sc.tl.rank_genes_groups(adata, groupby='leiden', method='wilcoxon', use_raw=True)
    adata.write(output_dir / 'processed.h5ad')
    print('scRNA-seq pipeline complete.')

if __name__ == '__main__':
    main()
