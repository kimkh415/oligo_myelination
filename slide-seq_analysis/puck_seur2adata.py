import os
import gzip
import scipy.io
import pandas as pd
import anndata as ad
import numpy as np
import tacco as tc


coords_file = "data/rotated_coordinates.csv"
cells_file = "data/barcodes.tsv"
genes_file = "data/genes.tsv"
counts_file = "data/matrix.mtx"

#with gzip.open(counts_file, 'r') as f:
X = scipy.io.mmread(counts_file).T


cells = pd.read_csv(cells_file, sep='\t', header=None).to_numpy().flatten()
genes = pd.read_csv(genes_file, sep='\t', header=None).iloc[:,0].to_numpy().flatten()

adata = ad.AnnData(X.tocsr(),obs=pd.DataFrame(index=cells),var=pd.DataFrame(index=genes),dtype=np.float32)

df = pd.read_csv(coords_file)
adata.obs['x'] = df.x.values
adata.obs['y'] = df.y.values

# make adata conform to the expected conventions
#adata.obs[['x','y']] /= 0.65 # go back from µm to pixel units
adata.X = adata.X.tocsc()
adata = adata[:,tc.sum(adata.X,axis=0) != 0].copy()

adata.write("data/slideseq.h5ad", compression='gzip', compression_opts=9)

print("DONE")

