# armd-common

A harmonization layer for the four ARMD antimicrobial resistance datasets: Stanford, MGB, UTSW and
ECUH. The package ships code and split manifests only. It never ships data.

Harmonization decisions and data quirks are documented in `docs/HARMONIZATION.md`.

## Setup

```bash
uv venv --python 3.11 .venv
VIRTUAL_ENV=.venv uv pip install -e ".[model,dev]"
cp .env.example .env    # set ARMD_DATA_ROOT, and Dryad API credentials for UTSW and ECUH
```

## Data

| site | source | folder under `ARMD_DATA_ROOT` |
|---|---|---|
| Stanford | Dryad, doi:10.5061/dryad.jq2bvq8kp | `armd-stanford` |
| MGB | PhysioNet credentialed access, doi:10.13026/2r5k-b955 | `armd-mgb/physionet.org/files/armd-mgb/1.0.0` |
| UTSW | `python -m armd_common.download --site utsw` | `armd-utsw` |
| ECUH | `python -m armd_common.download --site ecuh` | `armd-ecuh` |

Every raw file is checked against the sizes and SHA-256 digests in
`src/armd_common/data_sources.yaml` before it is staged.

## Pipeline

```bash
python -m armd_common.convert --site stanford    # stage the CSVs as Parquet, once per site
python analysis/00_gate.py                       # Phase 0 feasibility gate across staged sites
python -m armd_common.qc --site stanford         # per-site QC report
python -m armd_common.sex_inference              # evidence behind mappings/sex_codes.yaml
python -m armd_common.splits --site stanford     # generate and freeze the split manifest
python -m armd_common.features --site stanford   # harmonized 48 features, written under the site data folder
python analysis/01_build_cohorts.py              # E1, E2 and sensitivity cohorts plus Analysis 1, derivation only
python analysis/02_train_e1.py                   # models, cross-fitted predictions and recalibration maps, derivation only
python analysis/03_evaluate_ladder.py            # the only script that reads the test partition
python analysis/04_cascade_sensitivity.py        # sentinel reporting and the evaluable subset, derivation only
python analysis/05_e3_bounds.py                  # population bounds from the 03 results
.venv/bin/pytest
```

Scripts write provenance to `results/<script>/<timestamp>/manifest.json`.

## Canonical tables

```python
from armd_common import labels
from armd_common.sites import get_loader

loader = get_loader("stanford")
tables = loader.build()   # TEMP tables inside loader.con
e1 = labels.esbl_phenotype_composite(loader.con, tables["culture_episodes"], tables["isolate_ast"])
e2 = labels.esbl_decision_time(loader.con, tables["culture_episodes"], tables["isolate_ast"])
```

| table | grain |
|---|---|
| `culture_episodes` | one row per culture |
| `isolate_ast` | one row per distinct culture, organism, antibiotic and result |
| `enzyme_tests` | one row per culture, organism and enzyme test (MGB only) |
| `culture_keys` | raw IDs and order date for joining staged feature tables |

Every label constructor returns a frame with `culture_id`, `label`, `evaluable` and `reason`, plus a
provenance dictionary.

## Data use

All four releases are under data use terms, and MGB is credentialed. Never print or commit row-level
data. Reports carry aggregate counts only, and counts from 1 to 10 are suppressed.
