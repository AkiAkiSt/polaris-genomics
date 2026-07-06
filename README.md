# POLARIS
### Posterior-calibrated, Latent-truth Attribution of Regulatory-variant Influence and target-gene Selection

### 🔭 **Live interactive site → https://akiakist.github.io/polaris-genomics/**
An interactive, self-teaching explorer of the real results (829 variants, two diseases) with the honesty
layer front-and-centre — transparent confidence, visible "uncertain" flags, and disagreements shown, not
hidden. Screenshots in [`docs/screenshots/`](docs/screenshots/).

**Broad question (Candidate A from the proposal):** *How can we understand the function of non-coding DNA
variants to predict whether they may cause disease?* → For a single complex disease, convert its non-coding
GWAS loci into **ranked, mechanistically-annotated causal-variant → target-gene hypotheses.**

POLARIS is a transparent, calibrated framework that turns the field's two guardrails into mathematics:
**functional ≠ pathogenic** (a factorized posterior) and **distrust any single oracle** (a reliability-weighted
latent-truth model). It is built entirely from a *white-box* oracle stack and is fully reproducible from
free public APIs — providing a scaffold into which black-box predictors (AlphaGenome, Enformer) can be
*calibrated* rather than trusted blindly.

---

## Headline results
| | Type 2 diabetes | Coronary artery disease |
|---|---|---|
| Loci / fine-mapped variants | 30 / 569 | 14 / 260 |
| Target-gene top-1 (controls) | **0.85** (vs 0.69 nearest) | **0.83** (vs 0.67 nearest) |
| Conservation enrichment (causal vs passenger) | AUC 0.68, p=0.011, CI [0.55,0.80] | AUC 0.87, p=5e-4 |
| Rarity **inverted** (functional≠pathogenic) | AUC 0.28, p=2e-3 | AUC 0.32 |
| EM-inferred top channel | phyloP (wt 0.62, perm p=0.013) | phyloP (wt 0.66) |
| Calibration | ECE 0.005, Brier 0.034 | replicated |
| Evidence-channel independence (Marchenko–Pastur) | 2/7 shared dims | — |

- **Mechanism recovered from first principles:** FTO rs1421085 → **ARID5B** motif (ΔLLR −8.6 bits) → IRX3
  (Claussnitzer 2015); SORT1 (ΔLLR −15.4); distal PHACTR1 (−13.7).
- **Honest uncertainty:** confidence-stratified abstention raises gene precision 0.85→0.89; the two misses
  are the canonical hard cases (distal FTO→IRX3, flagged; imprinted DUSP8, lowest confidence).

---

## Deliverables
| Artifact | Path |
|---|---|
| **Manuscript** (Nature-style, 13 pp, 7 figs) | [`manuscript/polaris.pdf`](manuscript/polaris.pdf) |
| **Foundational explainer** (no-prior-knowledge) | [`explainer/explainer.pdf`](explainer/explainer.pdf) |
| **Presentation** (10-min, figure-driven) | [`presentation/POLARIS_deck.pptx`](presentation/POLARIS_deck.pptx) |
| **Website** (interactive) | [`website/index.html`](website/index.html) |
| **Figures** (PNG + SVG) | [`results/figures/`](results/figures/) |
| **Ranked hypothesis table** | [`results/tables/ranked_hypotheses.csv`](results/tables/ranked_hypotheses.csv) |
| Literature review · harsh self-review · plan | [`docs/`](docs/), [`PROJECT_PLAN.md`](PROJECT_PLAN.md) |

## Code (`polaris/`) — every engine from scratch, validated on simulation
- `apis.py`, `cache.py` — cached clients for Open Targets, Ensembl, GTEx, gnomAD, UCSC, JASPAR, ClinVar.
- `finemap.py` — **SuSiE-RSS** fine-mapping (sim: PIP 0.96, FDR 0.00).
- `motifs.py` — **information-theoretic TF-binding disruption** (ΔLLR ≈ ΔΔG, exact analytic null p-value).
- `molecular.py` — multimodal molecular-effect orchestrator (motif + conservation + cCRE + splice + AF).
- `coloc.py` — **Giambartolomei colocalization** (sim: shared/distinct PP4 0.999/0.000; real-data via dyneqtl).
- `integrate.py` — **latent-truth EM** (continuous Dawid–Skene) + calibration (Platt/isotonic, Brier, ECE).
- `geometry.py` — diffusion maps, entropic optimal transport (Sinkhorn), Marchenko–Pastur null.
- `dynamics.py` — gene-regulatory **ODE** (ΔLLR → directional expression prediction).
- `viz.py` — publication figure identity + TextPath sequence logos.

## Reproduce
```bash
python scripts/01_fetch_loci.py        # GWAS loci + fine-mapped credible sets (Open Targets)
python scripts/04_genes.py             # gene TSS distances per locus (Ensembl)
python scripts/03_molecular_effect.py  # transparent molecular feature matrix (~15 min, cached)
python scripts/05_integrate.py         # POLARIS posterior, gene model, validation, ranked table
python scripts/06_real_coloc.py        # real-data colocalization (GTEx dyneqtl)
python scripts/07_second_disease.py    # CAD generalization (full pipeline)
python scripts/08_robustness.py        # bootstrap CIs, permutation null, abstention
python scripts/fig_*.py                # all figures
python scripts/pptx_build.py           # presentation
# manuscript / explainer: data/external/tec_msvc/tectonic.exe via manuscript/compile.ps1
```
All external calls are cached under `data/cache/`, so reruns are offline-fast and every number is reproducible.

**Data sources (all free, no institutional approval):** GWAS Catalog · Open Targets Platform · Ensembl REST ·
GTEx Portal v2 · gnomAD · NCBI/ClinVar · UCSC · JASPAR 2024.
