# Harsh Nature reviewer report — POLARIS (v1)

**Recommendation: Reject and resubmit.** The framing is clear and the engineering is competent,
but in its current form the manuscript does not meet the bar for a top-tier venue. The central
claims are under-powered, the headline gene result merely re-uses an existing tool, and several
analyses rest on uncomfortably small sample sizes. Specific, blocking concerns:

### Major
1. **Scale and generality (blocking).** A single disease and 30 hand-curated loci is a pilot, not a
   study. "Calibrated triangulation" is a general claim; it must be shown to *generalize* — at least
   a second disease with its own positive controls, ideally a much larger locus set. As written I
   cannot distinguish a real method from overfitting to a curated list that was chosen to contain the
   answers.
2. **Where is the advance over L2G? (blocking).** The gene-nomination headline (85%) is *literally
   L2G*. The authors even show their own additions hurt it. A methods paper must demonstrate that the
   new machinery *adds* something measurable: either improved accuracy, or — at minimum — calibrated
   confidence that knows when L2G is wrong. Right now the novelty is asserted, not demonstrated.
3. **Under-powered statistics.** The reliability weights, the conservation-enrichment AUC (0.68), and
   the calibration (ECE on **17** positives) are all computed on tiny samples with no confidence
   intervals. A bootstrap or permutation null is mandatory before any of these numbers can be trusted.
4. **The transparent module barely works for variant discrimination (AUROC 0.59).** The authors spin
   this as "orthogonality," but a reviewer reads it as "the headline white-box module is weak." You
   must show it contributes to an *ensemble* that beats every single channel on a held-out task, or
   the transparency is a liability, not a feature.
5. **No external functional ground truth.** All labels derive from fine-mapping PIP. Conservation
   "predicting" PIP is partly mechanical (both correlate with selection). An independent label —
   ClinVar non-coding, MPRA activity, or measured eQTL status — is needed.
6. **Colocalization is simulation-only.** A coloc engine that is never run on real data is decorative.
   Either apply it to the actual loci or remove the claim.

### Minor
7. Direction-of-effect is asserted (ODE) but never checked against measured eQTL effect signs.
8. The diffusion-map "manifold" was silently replaced by PCA; either the geometry works or it doesn't.
9. The ODE log2FC (+5) is biologically implausible in magnitude; report direction only, or calibrate.
10. "Nature-caliber" tone exceeds the evidence; temper claims and foreground limitations.

### What would change my mind
A second disease showing the *same* calibrated behavior; bootstrap CIs on every headline number; a
demonstration that confidence-thresholded POLARIS achieves high precision (i.e., it *knows* when to
abstain on hard cases like FTO/DUSP8); a real-data colocalization result; and an independent functional
label. Do these and this becomes a genuine contribution.

---

## Concrete improvement plan (executed in v2)
- **G1 Generalization:** add a second disease (Coronary Artery Disease) with its own textbook control
  (SORT1 rs12740374 / C-EBP), run the *entire* pipeline, and report identical metrics → show transfer.
- **G2 Beat/qualify L2G:** confidence-stratified precision — show POLARIS assigns *low* confidence to
  its own errors, so a thresholded version reaches high precision and correctly abstains on FTO/DUSP8.
- **G3 Rigor:** bootstrap 95% CIs on reliability weights, conservation AUC, rarity AUC, calibration;
  permutation null for the reliability inference.
- **G4 Ensemble value:** locus-aware CV showing PIP+conservation+L2G ensemble (POLARIS) AUROC > each
  single channel for recovering causal variants.
- **G5 Real coloc:** run the from-scratch Giambartolomei coloc on real GTEx dyneqtl data for showcase
  loci; report PP4 and tissue of action.
- **G6 Direction-of-effect:** test concordance between motif-disruption sign / ODE direction and the
  measured GTEx eQTL effect-size sign.
- **G7 Honesty:** report ODE direction (not magnitude); keep PCA but add the diffusion-map result on a
  clean subset as a supplementary; temper tone; expand limitations.

---

## v2 response — what was done
- **G1 ✓** Entire pipeline rerun on **coronary artery disease** (14 loci, 260 variants). Conservation
  enrichment (AUC 0.87, p=5e-4), rarity inversion (0.32), reliability ordering (phyloP 0.66), and gene
  nomination (0.83 vs 0.67) **all replicate**; SORT1 & distal PHACTR1 mechanisms recovered. `fig_generalization`.
- **G2 ✓** Confidence-stratified **abstention** raises gene precision 0.85→0.89; lowest-confidence locus is a
  true error (DUSP8); FTO caught by the distal flag. `fig_robustness` C.
- **G3 ✓** Bootstrap 95% CIs on conservation AUC [0.55,0.80] and rarity AUC [0.18,0.40]; **permutation null**
  for reliability (p=0.013). `fig_robustness` A,B,D.
- **G4 ✓ (honest)** No naive ensemble beats the best single channel (channels orthogonal); reframed to show
  reliability-weighting recognises the predictive channel instead of being diluted.
- **G5 ✓** From-scratch coloc now runs **end-to-end on real GTEx data** (44 tests, coloc-SuSiE log(PIP) ×
  dyneqtl ABFs). PP4 modest in available tissues — GTEx lacks pancreatic islets — reported as a
  tissue-coverage limitation, not hidden.
- **G6 / G7** ODE reported as direction (robust across parameters), tone tempered, limitations expanded.

Net: the two blocking concerns (generality; advance-over-L2G) are addressed by a full second-disease
replication and by demonstrating calibrated abstention that *knows when L2G is wrong*.
