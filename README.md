# CRC Mutational Signature Pipeline

Replication and expansion of Battuello et al. (*Brief Bioinform* 2024, [bbae249](https://doi.org/10.1093/bib/bbae249)):

| Phase | Question | Data |
|-------|----------|------|
| **1 — Algorithm + Reference** | How do fitting tool and COSMIC catalogue affect cosine, ΔMMR, ΔPOLE? | Preclinical WES (n=230) + Clinical WES (n=152) |
| **2 — NGS workflow** | How does sequencing scope affect the same readouts? | Preclinical WES / WGS (n=63) / TSO-500 (n=230) |

**Status: phases 1 and 2 complete.** Master cosine table: **33,750** rows (10 tools × 5 COSMIC × all available datasets).

## Layout

| Path | Role |
|------|------|
| `input_files/` | SBS96 catalogues + Table S1–S3 |
| `cosmic_files/` | COSMIC SBS references (GRCh38) |
| `config/` | Tool registry, COSMIC list, signature-group maps |
| `metadata/` | MSI / MSS / POLE sample lists |
| `scripts/r/` | R fitters + merge + downstream |
| `scripts/python/` | SPA / MuSiCal runners |
| `scripts/_archive/` | Previous exploratory scripts (do not run) |
| `results/` | Local run outputs (not stored in this repository) |

## Prerequisites

- Run all commands from **`mywork/`** as the working directory.
- R packages: MutationalPatterns, deconstructSigs, signature.tools.lib, siglasso, sigfit, sigminer, SigsPack, YAPSA, readxl, jsonlite, ggplot2, dplyr
- Python (conda):
  - MuSiCal in env `musical` (`/opt/anaconda3/envs/musical`)
  - SigProfilerAssignment in env `pyenv` (`/opt/anaconda3/envs/pyenv`)
  - SignatureAnalyzer is **disabled** in config (not equivalent catalogue fitting)

## Phase 1 run order (WES)

```bash
cd mywork

Rscript scripts/r/00_build_metadata.R
Rscript scripts/r/run_r_fitting.R
bash scripts/python/run_py_fitting.sh
Rscript scripts/r/00_combine_cosine.R
Rscript scripts/r/downstream/run_downstream.R
```

Optional smoke (MutationalPatterns × v2 + v3.2, preclinical WES):

```bash
Rscript scripts/r/run_mp_smoke.R
```

## Phase 2 run order (NGS)

```bash
cd mywork

Rscript scripts/r/run_ngs_smoke.R          # MP × v3.2 × WGS+TSO
Rscript scripts/r/run_r_fitting.R --phase 2
bash scripts/python/run_py_fitting.sh --phase 2
Rscript scripts/r/00_combine_cosine.R
Rscript scripts/r/downstream/run_downstream.R
Rscript scripts/r/downstream/run_ngs_workflow.R
```

## Phase 3 — optional expansion

```bash
cd mywork

# 1) CRC tissue-specific reference (paper GitHub CRC_Specific ⊂ COSMIC v3.2)
Rscript scripts/r/run_tissue_specific.R
Rscript scripts/r/00_combine_cosine.R
Rscript scripts/r/downstream/run_downstream.R
Rscript scripts/r/downstream/run_reference_expansion.R

# 2) Mutation-count threshold Fig. 6A–F across ALL COSMIC versions in config
Rscript scripts/r/run_mutation_downsampling.R
# optional subset: Rscript scripts/r/run_mutation_downsampling.R --cosmic v2,v3.2,v3.4,v3.5,v3.6,TS

# Problematic-sample discard (cosine / mut-threshold / flat ≥ 0.25)
Rscript scripts/r/downstream/run_problematic_discard.R
```

| Analysis | Status | Output |
|----------|--------|--------|
| Tissue-specific (TS) + Fig. 5-style panels | **Done** | `results/downstream/reference_expansion/` |
| Mutation downsampling Fig. **6A–F** (all COSMIC) | **Runnable** | `results/downstream/mutation_threshold/plots/Fig6{A–F}_*` |
| Problematic-sample discard + TCGA OS | **Done** | `results/downstream/problematic_discard/` |
| Paper metanormal rematch (original 6C–D) | **Blocked** | needs 21 PBMC WES (not in workspace) |
| Exon / intron / extragenic WGS pies | **Blocked** | needs regional SBS96; use WES vs WGS flat inflation as proxy |

TS catalogue: `cosmic_files/COSMIC_CRC_tissue_specific_SBS_GRCh38.txt` (24 SBS; paper text says 26 — we use the published GitHub list).

**Fig. 6 panel map (this pipeline):** A/B cosine curves by COSMIC; C/D thresholds by COSMIC; E/F assigned-fraction curves. Paper’s metanormal C/D are not reproducible here.
## Deliverables

**Cosine / QC**

- `results/cosine/CosineResults_AllTools_AllCOSMIC.tsv`
- `results/cosine/cosine_qc_summary.json`

**Phase 1 (Algorithm + Reference)**

- Tables: `results/downstream/tables/` — `delta_MMR_POLE.tsv`, `assigned_fraction.tsv`, `MMR_composition_MSI.tsv`, Wilcoxon pairs
- Plots: `results/downstream/plots/` — cosine, ΔMMR, ΔPOLE, MMR composition

**Phase 2 (NGS)**

- Tables: `results/downstream/ngs_workflow/tables/`
- Plots: `results/downstream/ngs_workflow/plots/`

**Phase 3 (expansion)**

- `results/downstream/reference_expansion/` — v2 / v3.2 / TS comparison
- `results/downstream/mutation_threshold/` — Fig. 6A–F downsampling
- `results/downstream/problematic_discard/` — Maura-style QC discard experiment (ΔMMR/ΔPOLE + TCGA OS)
- `metadata/clinical_survival.tsv` — TCGA COAD/READ overall survival (152 samples, 36 events)
## Output contracts

- Contributions: `signatures × samples` under `results/contributions/{Tool}/{COSMIC_Version}/`
- Reconstructed: `contexts × samples` under `results/reconstructed/...`
- Cosine long TSV: `Sample, Cosine_Similarity, Dataset, Data_Type, Tool, COSMIC_Version, Mutation_Count`
- Cosine threshold for “assigned”: **0.9**

## Notes / caveats

- **sigLASSO** needs `nnls`, `glmnet`, and `RColorBrewer` attached (handled in the fitter).
- **sigfit** expects catalogues as samples×contexts and signatures as signatures×contexts (handled in the fitter).
- **SigsPack** on overcomplete catalogues (e.g. COSMIC v3.6, K>96) falls back from QP to NNLS when needed.
- **SignatureAnalyzer** remains disabled (not catalogue fitting).
- A few samples have `NA` cosine (notably **sigLASSO × v2**); they are dropped from boxplots.
- Larger COSMIC catalogues generally raise cosine and shrink ΔMMR for several tools.
