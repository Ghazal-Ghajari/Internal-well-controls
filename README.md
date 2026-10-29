# Internal well controls

Code for the paper

> Ghajari G., Amsaad F., Thirunarayan K. *Plate layout limits which perturbation responses can be separated from well effects in multiplexed single-cell screens.* (preprint link to be added)

In multiplexed single-cell screens, several cell types or cell lines share a well. The paper shows that a technical signal shared by everything in a well is present in Parse PBMC and Tahoe-100M. It then asks what within-well centering removes, when it removes biology, and which drug responses the plate layout can separate from well effects at all.

## Contents

```
notebooks/   analysis notebooks (one per part of the paper, see below)
scripts/     stand-alone scripts
data/        small derived inputs (large data are downloaded, see below)
```

Each notebook writes its figures and tables to `results/cache/<notebook name>/`.

| Notebook | What it does |
|---|---|
| `aim2_download.ipynb`, `aim2_extend_download.ipynb` | Download of Tahoe-100M plates and metadata (run first) |
| `aim2_parse_shared_well_effect.ipynb` | Shared well signal in Parse PBMC (Tests A to D), ambient RNA, doublets, noise ceiling |
| `aim2_tahoe_shared_well_effect.ipynb` | Shared well signal in Tahoe-100M (plates 6, 13, 14), ambient RNA |
| `aim2_model_simulation.ipynb` | First simulator (low-rank prototype) |
| `aim2_model_simulation_v2.ipynb` | Main simulator, calibration of calls, misspecified simulator |
| `aim2_model_real_data.ipynb` | Estimators on Parse and Tahoe-100M, two-part signature, independent tests |
| `aim2_baselines.ipynb` | Comparison with ComBat, Harmony, RUV and the well mean; leakage theory; ComBat diagnosis |
| `aim2_extend_analysis.ipynb` | OP3; all 14 Tahoe plates; identifiability and variance test; confirmatory test |

`notebooks/fig1b_inert_size.csv` lists the inert cytokines and their response sizes from the companion paper (its Figure 1B); `aim2_parse_shared_well_effect.ipynb` and `aim2_model_real_data.ipynb` read it.

| Figure | Notebook | Output file |
|---|---|---|
| 1 | `aim2_tahoe_shared_well_effect` | `aim2_tahoe_shared_well_effect.png` |
| 2 | `aim2_model_simulation_v2` | `aim2_model_simulation_v2.png` |
| 3 | `aim2_model_real_data` | `independent_tests.png` |
| 4 | `aim2_extend_analysis` | `aim2_tahoe_all_plates.png` |
| S1 | `aim2_parse_shared_well_effect` | `aim2_shared_well_effect.png`, `aim2_testD_robustness.png` |
| S2 | `aim2_parse_shared_well_effect`, `aim2_tahoe_shared_well_effect` | `aim2_parse_ambient_check.png`, `aim2_tahoe_ambient_check.png` |
| S3 | `aim2_model_simulation` | `aim2_model_simulation.png` |
| S4 | `aim2_model_simulation_v2`, `aim2_model_real_data` | `aim2_calibration.png`, `aim2_two_part_calibration_variants.png` |
| S5 | `aim2_extend_analysis` | `aim2_op3_tests.png` |
| S6 | `aim2_extend_analysis` (nine-plate option, see the first cells) | `aim2_tahoe_nine_plates.png` |
| S7 | `aim2_baselines` | `aim2_baselines.png` |
| S8 | `aim2_baselines` | `leakage_theory.png` |
| S9 | `aim2_extend_analysis` | `variance_test.png` |
| S10 | `aim2_model_simulation_v2` | `simulation_v3_misspecified.png` |
| S11 | `aim2_model_real_data` | `real_data.png` |



### ComBat example

`scripts/combat_singular_design_mre.py` shows that `scanpy.pp.combat` with a categorical covariate leaves the well effect in Tahoe-100M plate 6. With the same covariate encoded as k − 1 numeric columns (a full-rank design for the same model) it is removed. Input: `data/combat_mre_plate6.npz`.

```
python scripts/combat_singular_design_mre.py
```

## Data

| Data set | Source |
|---|---|
| Parse PBMC cytokine screen (10M cells) | Parse Biosciences, https://www.parsebiosciences.com/datasets/10-million-human-pbmcs-in-a-single-experiment/ |
| Tahoe-100M | HuggingFace `tahoebio/Tahoe-100M`; h5ad files from the Arc Institute Virtual Cell Atlas bucket |
| OP3 compound screen | GEO GSE279945 |

## Requirements

Python 3.12. Exact versions used for the paper are in `requirements.txt`.

```
pip install -r requirements.txt
```

Related package: [donorvar](https://github.com/Ghazal-Ghajari/plate-perturbation-confounders) (companion paper).

## Pre-specification

The confirmatory test on Tahoe plates 4, 5 and 7 to 9 was fixed before it was run: commit `02bcb19` (26 October 2026, 21:05), tagged [`prespecification`](../../tree/prespecification). The same files were sent by e-mail to a co-author before the new plates were processed.

## Citation

Please cite the paper above (details to be added on publication). An archived copy with a DOI will be on Zenodo.

## License

MIT (see `LICENSE`).
