# POLARIS — Project Plan
### Posterior-calibrated, Latent-truth Attribution of Regulatory-variant Influence and target-gene Selection

**Candidate A** from the proposal: *"GWAS hit → mechanistic hypothesis."* For a single complex disease,
convert non-coding GWAS loci into **ranked, mechanistically-annotated causal-variant → target-gene hypotheses.**

**Disease scope:** Type 2 Diabetes (T2D) — extremely well-powered GWAS (DIAMANTE/Mahajan 2022, ~1.3M),
rich islet/adipose/liver eQTLs, and *textbook non-coding mechanistic positive controls* for validation:
`TCF7L2` (rs7903146, islet enhancer), `FTO/IRX3-IRX5` (rs1421085, ARID5B motif disruption — Claussnitzer 2015),
`MTNR1B` (rs10830963), `ADCY5`, `CDKN2A/B` (ANRIL), `KCNQ1`, `SLC30A8`.

---

## Guiding philosophy (the two PDF guardrails, made mathematical)
1. **Functional ≠ Pathogenic.** POLARIS *factorizes* the posterior into a molecular-function factor and a
   disease-linkage factor, keeping them conceptually + numerically separate.
2. **Triangulate; distrust any single oracle.** POLARIS treats every predictor as a *noisy annotator* and
   infers its reliability via a Dawid–Skene-style latent-truth EM. Credibility = convergence, not any one score.

## Novel contributions (the "interesting part" = the integration logic in Step 5)
- **C1 — Transparent, white-box molecular-effect module.** Instead of a black-box deep-learning oracle, an
  information-theoretic, allele-resolved TF-binding disruption model (JASPAR PWMs + order-k Markov background +
  calibrated nulls), giving *auditable mechanism* annotations ("disrupts a GATA motif").
- **C2 — Calibrated Bayesian evidence fusion** with explicit functional/pathogenic factorization and full
  uncertainty quantification; calibrated against known mechanisms (Brier/ECE/reliability curves).
- **C3 — Reliability-weighted oracle adjudication** (latent-truth EM) — learns which oracle to trust where.
- **C4 (stretch / NeurIPS-level)** — geometry & dynamics: diffusion-map manifold of the evidence space
  (graph-Laplacian heat kernel), entropic optimal-transport tissue attribution (Sinkhorn), a Marchenko–Pastur
  random-matrix null for oracle redundancy, and an ODE regulatory-dynamics model for mechanistic vignettes.

## Pipeline (extends Candidate A's 6 steps)
0. **Scope & loci** — T2D, 15–25 loci incl. positive controls.
1. **GWAS signal** — Open Targets Platform + GWAS Catalog → lead SNPs, ±500 kb windows.
2. **Fine-mapping → credible sets** — (a) Open Targets credible sets+PIPs; (b) **SuSiE reimplemented in Python**
   on ≥1 locus using Ensembl LD (methodological depth).
3. **Molecular effect (multimodal, transparent)** — motif disruption (C1) + conservation (phyloP/phastCons/GERP,
   UCSC) + regulatory context (Ensembl regulatory build / ENCODE cCREs) + CADD v1.7 + splicing (SpliceAI lookup or
   MaxEntScan-style from scratch). → variant × modality feature matrix.
4. **Variant→gene→tissue** — (a) Open Targets colocalisation + L2G; (b) **Giambartolomei coloc from scratch** on
   GWAS×GTEx eQTL/sQTL; triangulate with nearest gene + enhancer-target.
5. **Integrate & rank** — the **POLARIS** model (C2+C3); baselines: naive weighted sum, L2G, distance, best oracle.
6. **Output + sanity check** — ranked hypothesis table; recovery of known mechanisms + ClinVar; vignettes.

## Data sources (all free, no institutional approval, API/installable)
GWAS Catalog REST · Open Targets Platform GraphQL · Ensembl REST (VEP/LD/sequence/regulatory) · GTEx Portal API v2
(eQTL/sQTL) · gnomAD GraphQL (AF/constraint) · NCBI E-utilities (ClinVar) · UCSC API (sequence/conservation/ENCODE)
· JASPAR API (TF PWMs) · Broad SpliceAI-lookup API.

## Deliverables
Reproducible code + notebook · ranked hypothesis table + figures · breathtaking EDA · LaTeX manuscript (→PDF via
tectonic) · harsh self-review + revision · 10-min PPTX (figure-driven) · static website · foundational explainer PDF.

## Risk register & pivots
- **AlphaGenome key unavailable** → C1 transparent module is the *headline*, not a fallback.
- **No R** → SuSiE + coloc reimplemented in Python (depth, not deficit).
- **No LaTeX preinstalled** → fetch tectonic; fallback to matplotlib/reportlab PDF if needed.
- **Python 3.14, no torch** → classical ML (sklearn) + custom numpy/scipy Bayesian models (aligns with white-box thesis).
- **Any API failing** → cache aggressively; degrade gracefully; always pivot and continue.
