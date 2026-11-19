# Thesis draft — Results & gap analysis

Replication and expansion of Battuello et al., *Briefings in Bioinformatics* (2024), [bbae249](https://doi.org/10.1093/bib/bbae249).  
Figures refer to files in `results/thesis/figures/` (see `FIGURES.md`). Numeric checks: `paper_vs_replicate.tsv` and `results/downstream/reference_expansion/tables/ref_paper_vs_replicate.tsv`.

---

## 1. Aims and scope

This work asked how **fitting algorithm**, **COSMIC reference**, and **NGS scope** change mutational-signature readouts in colorectal cancer (CRC), using two complementary criteria from the paper:

1. **Mathematical:** cosine similarity between observed SBS96 profiles and reconstructions (assigned if cosine ≥ 0.9).
2. **Biological:** ΔMMR (median MSI − MSS MMR-signature contribution) and ΔPOLE (median POLE-MUT − POLE-WT among MSS).

**Cohorts:** preclinical CRC cell lines (WES/TSO-500 n=230; WGS n=63) and TCGA clinical WES (n=152).  
**Grid:** 10 fitting tools × COSMIC v2, v3.2, v3.4, v3.5, v3.6 (+ CRC tissue-specific subset, TS) where applicable. SignatureAnalyzer was disabled (not equivalent catalogue fitting). Master cosine table: **34,132** rows.

---

## 2. Results

### 2.1 Algorithm choice (paper Fig. 4 analogue)

On preclinical WES with **COSMIC v2**, tools differed in both cosine distributions and biological separation (**Fig01–Fig05**).

- Most tools reached median cosine ≥ 0.9; SPA left a substantial fraction unassigned (**182/230 = 79.1%** assigned), matching the paper’s ~79% claim.
- All tools separated MSI from MSS on MMR contributions (Wilcoxon), but **ΔMMR magnitude** varied strongly. Among the paper’s tool set, **SPA had the highest ΔMMR (0.675)**; MutationalPatterns ΔMMR was **0.344** (paper ~0.34).
- With the expanded tool set, **sigLASSO ΔMMR (0.734)** exceeded SPA — an expansion finding, not a contradiction of the paper’s five-tool ranking.
- POLE separation was significant across tools; SPA again among the strongest (ΔPOLE ≈ 0.70; MP ≈ 0.59, paper ~0.59).

**Claim for thesis:** Algorithm choice is a first-order driver of both assignment rate and subtype Δ scores; SPA maximises biological contrast among classical tools at the cost of more unassigned samples.

### 2.2 Reference catalogue (paper Fig. 5 analogue)

MutationalPatterns on preclinical WES across **v2 / v3.2 / TS** (**Fig06, Fig11–Fig13**):

| Reference | Paper ΔMMR | This run | Match |
|-----------|------------|----------|-------|
| COSMIC v2 | ~0.34 | **0.344** | yes |
| COSMIC v3.2 | ~0.26 | **0.257** | yes |
| Tissue-specific (TS) | ~0.27 | **0.272** | yes |

Larger catalogues raised cosine but **shrunk ΔMMR**, consistent with signal dilution / overfitting discussed in the paper. Per-SBS MMR composition in MSI shifted with reference (e.g. SBS6 high on v2, redistributed on v3.2/TS) — **Fig13**.

**Expansion:** COSMIC v3.4–v3.6 were included in the full grid; trends continue toward higher cosine and attenuated ΔMMR for several tools (supplement **FigS1–FigS2**).

### 2.3 NGS workflow (paper Fig. 3 analogue)

MutationalPatterns · COSMIC v3.2 on preclinical WES / WGS / TSO-500 (**Fig07–Fig10**):

| Data type | Median cosine (smoke) | ΔMMR | Flat in MSI (median) |
|-----------|----------------------|------|----------------------|
| TSO-500 | ~0.960 | **0.323** | ~0 |
| WES | ~0.988 | **0.257** | ~0 |
| WGS | ~0.993 | **0.115** | **0.269** |

Ordering **TSO > WES > WGS** for ΔMMR matches the paper (0.32 / 0.26 / 0.11). High technical cosine on WGS coexists with **weaker biological MMR separation**, explained by elevated **flat** signatures in MSI on WGS (**Fig10**).

**Claim:** Large panels and WES remain effective for CRC MMR/POLE stratification; genome-wide catalogues do not automatically improve subtype readout.

### 2.4 Minimum mutation count (paper Fig. 6 expansion)

Subsampling 5–95% of mutations (5 replicates), refitting MutationalPatterns for **every COSMIC version** in the pipeline (**Fig6A–F**):

| Dataset | COSMIC | Inferred threshold (cos ≥ 0.9) | Paper (v2 context) |
|---------|--------|--------------------------------|--------------------|
| Preclinical | v2 | **~351** | 323 |
| Preclinical | TS | ~284 | — |
| Preclinical | v3.2–v3.6 | ~146 | — |
| Clinical | v2 | **~49** | 64 |
| Clinical | v3.2–v3.6 | ~17–21 | — |

**Interpretation:** Threshold depends on **both** mutation quality/quantity **and** reference size. Richer catalogues reach cosine ≥ 0.9 at lower mutation counts. Clinical WES needs far fewer mutations than the unmatched preclinical cell-line cohort — same direction as the paper.

Paper **metanormal** panels (original Fig. 6C–D: 323→145 after PBMC-matched calling; artefact drop) were **not** reproduced (no local metanormal). Our Fig6C–F report the **multi-COSMIC** expansion instead.

---

## 3. Gap analysis (replication vs expansion)

### 3.1 Successfully replicated (within tolerance)

| Item | Paper | Ours | Notes |
|------|-------|------|-------|
| MP ΔMMR preclinical WES v2 / v3.2 / TS | 0.34 / 0.26 / 0.27 | 0.344 / 0.257 / 0.272 | Abs diff &lt; 0.01 |
| SPA assigned ≥0.9 (preclinical, v2) | ~79% | 79.1% | Exact count 182/230 |
| SPA highest ΔMMR among paper tools | yes | yes | Among MP, DS, STL, SPA, MuSiCal |
| NGS ΔMMR order (MP·v3.2) | TSO&gt;WES&gt;WGS | 0.323&gt;0.257&gt;0.115 | Pattern match |
| Flat signatures ↑ on WGS MSI | yes | median 0.269 | WES/TSO ~0 |
| Mutation threshold order of magnitude (v2) | 323 / 64 | ~351 / ~49 | Same cohort contrast |

### 3.2 Intentional expansions (own beyond the paper)

1. **Tools:** +sigLASSO, sigfit, sigminer, SigsPack, YAPSA, MuSiCal (paper five-tool core retained).
2. **References:** +COSMIC v3.4, v3.5, v3.6 and CRC **TS** subset; full Algorithm×Reference grids.
3. **Fig. 6 multi-COSMIC:** thresholds and assigned-fraction curves for all configured catalogues (**Fig6A–F**), not only v2.
4. **Unified dual pipeline:** R + Python fitters → shared cosine schema → R downstream and paper-styled thesis figures.
5. **Problematic-sample discard experiment:** cosine / mutation-threshold / flat-load flags; ΔMMR–ΔPOLE and TCGA OS before vs after discard (`results/downstream/problematic_discard/`).

### 3.3 Known gaps / limitations (state explicitly in thesis)

| Gap | Impact | Mitigation / wording |
|-----|--------|----------------------|
| **No metanormal** (paper Fig. 6C–D) | Cannot claim artefact halving or 323→145 rematch | State blocked; cite WES vs WGS flat inflation as related biology |
| **No exon/intron/extragenic SBS96** | Cannot draw paper Fig. 3D *right* pies | Use class stacked bars by NGS type (**Fig10**) as proxy |
| **TS = 24 SBS** (GitHub list) vs paper text “26” | Minor catalogue difference | Disclose; ΔMMR still matches ~0.27 |
| **SignatureAnalyzer off** | Not in algorithm ranking | Justify: not catalogue fitting equivalent |
| **Few NA cosines** (e.g. sigLASSO×v2) | Dropped from boxplots | Report in Methods/QC |
| **SigsPack×v3.6 NNLS fallback** | Solver, not biology | Methods note |
| **Downsampling method** | Multinomial SBS96 subsampling + MP refit; paper used `shuf` on variant lists | Same scientific question; thresholds are approximate |
| **Clinical ΔMMR** | Paper quotes higher clinical ΔMMR (e.g. MP v2 ~0.61); we get ~0.61 | Aligns; keep table in supplement |
| **Expanded-tool ranking** | sigLASSO &gt; SPA on ΔMMR | Always qualify “among paper tools” vs “all tools” |
| **Survival power** | Only 36 OS events in clinical WES | Report as underpowered exploratory Cox/log-rank |
| **Flat cut = 0.25** | Aggressive for SPA on large catalogues | Sensitivity analyses; disclose discard rates |

### 3.4 What not to over-claim

- Do **not** claim exact reproduction of paper Fig. 6C–D or regional WGS pies.
- Do **not** say SPA is universally “best” — it is best **biological separator among the paper tool set**, with lower assignment rate.
- Do **not** imply WGS is “worse technically” — cosine is high; the deficit is **subtype contrast**, not reconstruction fidelity.
- Do **not** claim that discarding problematic samples broadly “improves” clinical associations — median ΔMMR change was zero; SPA gains cost most samples.

---

## 4. Suggested Results subsection map (copy into thesis)

1. **Study design** — three axes; cohorts; cosine / ΔMMR / ΔPOLE definitions.  
2. **Algorithms** — Fig01–05; SPA vs MP; expanded tools.  
3. **References** — Fig11–13; table v2/v3.2/TS.  
4. **NGS scope** — Fig07–10; TSO&gt;WES&gt;WGS; flat on WGS.  
5. **Mutation burden** — Fig6A–F; multi-COSMIC thresholds; metanormal gap.  
6. **Problematic-sample discard** — flags; ΔMMR change; TCGA OS.  
7. **Synthesis** — workflow choices dominate; QC discard is not a free robustness win.
---

## 5. One-paragraph abstract-style summary

We replicated Battuello et al. on CRC mutational-signature fitting and extended the analysis across additional tools and COSMIC releases (v3.4–v3.6) plus a CRC tissue-specific catalogue. Key paper benchmarks reproduced within tight tolerance (MutationalPatterns ΔMMR 0.344/0.257/0.272 for v2/v3.2/TS; SPA assignment 79.1%; NGS ΔMMR ordering TSO &gt; WES &gt; WGS with flat-signature inflation on WGS). Mutation-downsampling across all configured references showed that the mutation count needed for cosine ≥ 0.9 falls as catalogues grow. Metanormal rematching and true regional WGS matrices were unavailable and are reported as limitations. Overall, algorithm and reference choice remain dominant sources of variability, and sequencing breadth does not substitute for biologically appropriate genomic scope when the goal is MMR/POLE stratification.

---

## 6. File pointers

| Need | Path |
|------|------|
| Figures | `results/thesis/figures/` |
| Figure captions | `results/thesis/FIGURES.md` |
| Replication table | `results/thesis/paper_vs_replicate.tsv` |
| Reference Δ table | `results/downstream/reference_expansion/tables/ref_paper_vs_replicate.tsv` |
| Mutation thresholds | `results/downstream/mutation_threshold/tables/mutation_count_threshold.tsv` |
| Pipeline README | `mywork/README.md` |
