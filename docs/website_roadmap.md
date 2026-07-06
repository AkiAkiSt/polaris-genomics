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

## Progress log
- [x] Recon: gh authed (AkiAkiSt), no node, no secrets. git identity Aki.
- [x] .gitignore (exclude .env, data/cache 34M, data/external 119M, build artifacts).
- [ ] Git init + public repo + first push.
- [ ] Web data export (real results → NaN-safe JSON).
- [ ] Design system + interactive rebuild (D3): hero, honesty, pipeline, locus/variant explorer,
      FTO animated vignette, calibration/reliability/MP, T2D↔CAD.
- [ ] Deploy to Pages, verify live.
- [ ] Cold-eyes review + screenshots into repo + report.

## Honesty layer (non-negotiable, must be visible on every relevant screen)
- Confidence shown per hypothesis; "uncertain/abstain" flags visible (e.g. FTO distal flag, DUSP8 low conf).
- Disagreements shown, not hidden (nearest ≠ L2G ≠ truth at FTO).
- Persistent reminder: these are *hypotheses*, molecular effect ≠ proven disease causation.
- Every stat traceable to a real file in results/ or data/processed/.

## Decisions & pivots (log as they happen)
- (none yet)
