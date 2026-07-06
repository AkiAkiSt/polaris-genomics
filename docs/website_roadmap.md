# Website build — roadmap & running notes (Round 2)

**Goal of this run:** turn the validated POLARIS pipeline into a *stunning, deployed, self-teaching*
website. Live shareable link. Radical honesty front-and-centre. Stop at an excellent live site — do NOT
start later goals (many diseases, type-any-disease, live scoring).

**Standards:** publication-quality visuals · scientifically accurate & honest (every number from real
data) · teaches itself in seconds · smooth to demo live · tasteful motion.

**Stack decision:** static site (no node here) + D3.js/Plotly via CDN, no build step → deploys cleanly to
GitHub Pages. All paths relative (served under `/<repo>/`).

**Deploy:** GitHub Pages from `gh-pages` branch (site at branch root) via `gh`. Repo public, account AkiAkiSt.

## Working loop (per piece): build → run/prove (screenshot) → critique (engineer+designer+scientist) → fix → re-verify. Commit + push each working step.

## Progress log — ROUND 2 COMPLETE ✅
- [x] Recon: gh authed (AkiAkiSt), no node, no secrets. git identity Aki.
- [x] .gitignore + git init + **public repo** github.com/AkiAkiSt/polaris-genomics + first push.
- [x] GitHub Pages via Actions workflow → **LIVE at https://akiakist.github.io/polaris-genomics/**
- [x] Web data export (real results → NaN-safe JSON): variants(829), loci, meta, reliability, calibration,
      mp, findings, fto vignette, coloc. Every field traces to results/ or data/processed/.
- [x] Design system (Fraunces+Inter, light editorial + dark bands, Polaris starfield motif).
- [x] Interactive rebuild (D3): hero, honesty ribbon, problem, guardrails+findings, pipeline,
      master-detail EXPLORER (evidence triangulation + honesty flags), FTO animated vignette
      (locus map + sequence logo + regulatory ODE), validation (gene/calib/reliability-CI/MP),
      T2D↔CAD generalization, honesty section, footer.
- [x] Critique loop: verified every screen (preview + eval), mobile+desktop, console-clean.
      Fixes: list layout, marquee default, mobile nav menu, chart-label truncation, calib-dot cap,
      hero-canvas pause offscreen, FTO resize, scroll-margin (sticky-nav offset), OG social card.
- [x] 10 headless Playwright screenshots (desktop+mobile) → docs/screenshots/.

## Decisions & pivots (logged)
- **Python 3.14** path is explicit (`…/Python314/python.exe`); the `python` alias became a bare 3.12.
  Playwright installed on the 3.12 alias for screenshots.
- **No node** → static site + D3/Fonts via CDN, no build step (deploys clean to Pages).
- **Deploy via GitHub Actions** (not gh-pages branch) so every push auto-deploys; all paths relative.
- **CAD POLARIS fields** computed in the exporter with the SAME recipe as T2D (documented, not fabricated).
- **Explorer default** = marquee non-coding locus (TCF7L2/SORT1) for first impression; honest ranking kept.
- **JSON NaN** → strict NaN→null + allow_nan=False (recall from memory), cache-busted with `?v=`.
- **Coloc** shown as honest text (weak PP4, GTEx lacks islets) rather than a near-zero chart.
- Route A fine-mapping NOT run (per instruction); worked only with existing T2D+CAD data.

## Honesty layer (non-negotiable, must be visible on every relevant screen)
- Confidence shown per hypothesis; "uncertain/abstain" flags visible (e.g. FTO distal flag, DUSP8 low conf).
- Disagreements shown, not hidden (nearest ≠ L2G ≠ truth at FTO).
- Persistent reminder: these are *hypotheses*, molecular effect ≠ proven disease causation.
- Every stat traceable to a real file in results/ or data/processed/.

## Decisions & pivots (log as they happen)
- (none yet)
