# Literature Review & Novelty Positioning — POLARIS

## 1. The problem
~90% of GWAS associations for complex disease fall in **non-coding** DNA. The interpretive chain is
`variant → molecular effect → gene → tissue → trait`. Each link has dedicated tools, but they disagree and
are individually unreliable — so the field's open problem is **principled integration with calibrated confidence.**

## 2. The state of the art, by pipeline link

### Sequence → molecular effect (the "oracles")
- **AlphaGenome** (Avsec/DeepMind, *Nature* 2025; *Nat Struct Mol Biol* 2026). Largest multimodal DNA model:
  1 Mb context, base-pair resolution, predicts accessibility, TF binding, histone marks, splicing, expression,
  3D contacts; SOTA on 24/26 human variant-effect benchmarks. **Non-commercial API + Python SDK, key-gated.**
  *DeepMind's own stated limits:* distal elements >100 kb remain hard; tissue specificity needs work;
  "not designed/validated for individual-level prediction"; "cannot fully explain how genetic variation leads
  to complex disease." → motivates triangulation + functional≠pathogenic.
- **Enformer** (Avsec et al., *Nature Methods* 2021) — 196 kb transformer, expression/epigenome tracks.
- **Borzoi** (Linder et al., *Nature Genetics* 2025) — RNA-seq coverage prediction; strong expression baseline.
- **CADD v1.7** (Rentzsch et al., *Genome Medicine* 2021/2019) — genome-wide ensemble deleteriousness.
- **SpliceAI** (Jaganathan et al., *Cell* 2019) & **Pangolin** (Zeng & Li, *Genome Biology* 2022) — splicing.
- Reviews (Briefings in Bioinformatics 2024, bbae446; arXiv 2411.11158; Frontiers Mol Biosci 2026) converge on
  one message: deep models excel at **local, high-effect-size** variants but **degrade on distal regulation,
  direction-of-effect, and tissue/individual specificity** — confirmed by CRISPR-based evaluations.

### GWAS → causal variants (fine-mapping)
- **SuSiE** (Wang, Sarkar, Carbonetto, Stephens, *JRSS-B* 2020) — Sum of Single Effects regression → credible
  sets + PIPs. **SuSiE-RSS** (Zou et al., *PLoS Genet* 2022) for summary statistics + LD.
- **PolyFun+SuSiE** (Weissbrod et al., *Nat Genet* 2020) — functionally-informed priors.

### Variant → gene → tissue (linking)
- **coloc** (Giambartolomei et al., *PLoS Genet* 2014) and **coloc-SuSiE** (Wallace, *PLoS Genet* 2021).
- **Open Targets L2G** (Mountjoy et al., *Nat Genet* 2021) — ML over genetic+functional features.
- **cS2G** (Gazal et al., *Nat Genet* 2022) — combined SNP-to-gene strategy.
- **ABC model** (Fulco et al., *Nat Genet* 2019; Nasser et al., *Nature* 2021) — activity-by-contact enhancer→gene.
- **FLAMES** (*Nat Genet*, Oct 2025) — closest competitor: fuses SNP-based + "convergence" evidence for the
  effector gene; benchmarked on expert-curated locus-gene pairs. **Reports no calibrated posterior, no
  transparent per-variant mechanism, no per-oracle reliability inference** — the POLARIS gaps.

### Functional ground truth (validation)
- **MPRA** (Tewhey et al., *Cell* 2016; Abell et al., *Nature* 2022) — measured regulatory activity.
- **ClinVar** — clinical pathogenic/benign labels (non-coding subset).
- Textbook mechanistic positive controls used here: **FTO/IRX3-IRX5** rs1421085 ARID5B-motif disruption
  (Claussnitzer et al., *NEJM* 2015), **SORT1** rs12740374 C/EBP site (Musunuru et al., *Nature* 2010),
  **TCF7L2** islet enhancer (Gaulton; Grant et al. 2006), **MTNR1B** rs10830963.

### Data resources
GWAS Catalog (Sollis 2023); Open Targets Platform (Ochoa 2023); DIAMANTE T2D GWAS (Mahajan et al., *Nat Genet*
2022); GTEx v8 (GTEx Consortium, *Science* 2020); gnomAD v4 (Chen et al., *Nature* 2024); Ensembl (2024);
JASPAR 2024 (Rauluseviciute et al., *NAR* 2024); UCSC; ENCODE cCREs (Moore et al., *Nature* 2020).

## 3. Methodological lineage for POLARIS's novel core
- **Latent-truth from noisy annotators:** Dawid & Skene (*JRSS-C* 1979, EM for observer error rates);
  Raykar et al. (*JMLR* 2010). POLARIS adapts this to **variant-effect oracles**.
- **Calibration:** Platt (1999); isotonic (Zadrozny & Elkan 2002); ECE/Brier; Guo et al. (*ICML* 2017).
- **Geometry/dynamics (stretch):** diffusion maps (Coifman & Lafon, *ACHA* 2006 — graph-Laplacian heat kernel);
  entropic optimal transport / Sinkhorn (Cuturi, *NeurIPS* 2013; Peyré & Cuturi 2019); Marchenko–Pastur (1967)
  random-matrix null; gene-regulatory ODE steady-state analysis (Alon, *Systems Biology* 2006).

## 4. The gap → POLARIS's contribution
| Axis | L2G | cS2G | ABC | FLAMES | AlphaGenome | **POLARIS** |
|---|---|---|---|---|---|---|
| Calibrated posterior probability | partial | no | no | no | no | **yes** |
| Transparent, auditable mechanism | no | no | partial | no | **no (black box)** | **yes (biophysical)** |
| Per-oracle reliability *inferred* | no | no | no | no | n/a | **yes (latent-truth EM)** |
| Functional vs pathogenic *factorized* | no | no | no | no | no | **yes** |
| Uncertainty quantification | no | no | no | no | no | **yes** |
| Joint variant **and** gene posterior | no | no | no | gene-only | variant-only | **yes** |

**POLARIS thesis:** the credibility of a non-coding-variant hypothesis comes from *calibrated convergence across
independent, individually-fallible evidence channels* — not from any single oracle. We operationalize the
proposal's two guardrails as mathematics: a factorized posterior (functional ≠ pathogenic) and a
reliability-weighted latent-truth model (distrust any single oracle), validated by recovery of textbook
mechanisms and ClinVar, with full calibration diagnostics.
