"""Reproducible example: scanpy.pp.combat with a categorical covariate leaves the batch effect in.

Input: data/combat_mre_plate6.npz, Tahoe-100M plate 6 as pseudobulk log expression, one row per
(well, cell line), 2000 genes (written by section 18 of notebooks/aim2_baselines.ipynb).
Batch = well, covariate = cell line.

scanpy.pp.combat encodes a categorical covariate with a full set of dummy columns. Together with
one column per batch, the design matrix is one rank short, and scanpy inverts X'X with
numpy.linalg.inv. The same model is fitted twice:
  (1) covariates=["line"]            categorical, full dummy set   -> singular design
  (2) the same lines as k - 1 numeric 0/1 columns (one dropped)    -> full-rank design
The two should give the same corrected data. They do not, and only (2) removes the well effect.

Score: share of the variance (after removing each line's mean) that the well explains.
Close to 0 when the well effect is removed.

On synthetic data with the same structure the singular version usually still works, so the
failure depends on the data; this file uses the real input for that reason.

Run from the repository root:  python scripts/combat_singular_design_mre.py
"""
import sys
import warnings

import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
from scanpy.preprocessing._combat import _design_matrix

warnings.filterwarnings("ignore")
path = sys.argv[1] if len(sys.argv) > 1 else "data/combat_mre_plate6.npz"
z = np.load(path, allow_pickle=True)
X, well, line = z["X"].astype(np.float64), z["well"].astype(str), z["line"].astype(str)
obs = pd.DataFrame({"well": well, "line": line})
obs.index = obs.index.astype(str)
dum = pd.get_dummies(obs.line, drop_first=True).astype(float)
dum.columns = [f"L{i}" for i in range(dum.shape[1])]
ENC = {"categorical covariate": obs, "k-1 numeric columns": pd.concat([obs[["well"]], dum], axis=1)}


def well_share(Z):
    R = Z - pd.DataFrame(Z).groupby(line).transform("mean").to_numpy()   # remove line means
    W = pd.DataFrame(R).groupby(well).transform("mean").to_numpy()        # well means of the rest
    return (W ** 2).sum() / (R ** 2).sum()


print(f"scanpy {sc.__version__}, numpy {np.__version__}; {X.shape[0]} rows, {len(set(well))} wells, "
      f"{len(set(line))} lines, {X.shape[1]} genes\n")
print(f"{'':24s}{'columns':>8s}{'rank':>6s}{'cond(XtX)':>11s}{'well share':>12s}")
print(f"{'uncorrected':24s}{'':25s}{well_share(X):12.3f}")
OUT = {}
for name, o in ENC.items():
    D = _design_matrix(o, "well", tuple(o.groupby("well", observed=True).groups.keys())).to_numpy(float)
    A = ad.AnnData(X.copy(), obs=o)
    sc.pp.combat(A, key="well", covariates=[c for c in o.columns if c != "well"])
    OUT[name] = np.asarray(A.X)
    print(f"{name:24s}{D.shape[1]:8d}{np.linalg.matrix_rank(D):6d}{np.linalg.cond(D.T @ D):11.1e}"
          f"{well_share(OUT[name]):12.3f}")
d = np.abs(OUT["categorical covariate"] - OUT["k-1 numeric columns"]).max()
print(f"\nmax |difference| between the two encodings of the same model: {d:.3g}")
