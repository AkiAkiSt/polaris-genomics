"""Export ALL real pipeline results to clean, NaN-safe JSON for the interactive website.
Every field traces to a real file in results/ or data/processed/. Nothing is invented.

CAD's POLARIS fields are computed here with the SAME recipe as T2D (05_integrate) so the
two diseases are directly comparable; this is documented, not fabricated.
"""
from __future__ import annotations
import sys, json, math, pathlib
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
import numpy as np, pandas as pd
from scipy.special import expit
from polaris import config, integrate as IN, motifs as M, dynamics as DY, apis, utils

OUT = config.ROOT / "website" / "data"
OUT.mkdir(parents=True, exist_ok=True)
PROC, TAB = config.PROC, config.TAB
FCOLS = ["motif_motif_disruption", "phylop_max", "phastcons470", "ccre_enhancer", "splice_delta"]


def clean(o):
    if isinstance(o, dict):
        return {k: clean(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [clean(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating, float)):
        return None if (o is None or math.isnan(o) or math.isinf(o)) else round(float(o), 5)
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if o is None or (isinstance(o, float) and math.isnan(o)):
        return None
    return o


def dump(name, obj):
    p = OUT / name
    json.dump(clean(obj), open(p, "w"), allow_nan=False)
    return p.stat().st_size


def add_polaris_fields(df):
    """Compute P_functional (relative molecular EM score) + POLARIS for a disease frame."""
    for c in FCOLS:
        df[c] = df.get(c, np.nan)
    X = df[FCOLS].astype(float).values
    Xz = (X - np.nanmean(X, 0)) / (np.nanstd(X, 0) + 1e-9)
    y = np.full(len(df), np.nan)
    y[df.pip.values > 0.5] = 1.0
    y[df.pip.values < 0.05] = 0.0
    em = IN.latent_class_em(Xz, y)
    zl = (em["logit"] - np.nanmean(em["logit"])) / (np.nanstd(em["logit"]) + 1e-9)
    df["P_functional"] = np.clip(expit(zl), 1e-3, 1 - 1e-3)
    df["S_causal_functional"] = df["P_functional"] * df["pip"].fillna(0)
    df["polaris_gene"] = df["best_l2g_gene"]
    df["P_gene_best"] = df["best_l2g_score"].fillna(0)
    df["POLARIS"] = df["S_causal_functional"] * df["P_gene_best"]
    df["distal_flag"] = 0
    return df


def mech(r):
    tf, d, cc, g = (r.get("motif_top_disrupt_tf"), r.get("motif_top_disrupt_dLLR"),
                    r.get("ccre_class"), r.get("polaris_gene"))
    bits = []
    if isinstance(tf, str) and pd.notna(d):
        bits.append(f"{'abolishes' if d < -2 else 'weakens'} a {tf} motif ({d:+.1f} bits)")
    if isinstance(cc, str):
        bits.append(f"in {cc}")
    if isinstance(g, str):
        bits.append(f"→ {g}")
    return "; ".join(bits) if bits else "no transparent motif mechanism"


def variant_records(df, disease):
    recs = []
    for _, r in df.iterrows():
        vid = r["vid"]
        try:
            ch, pos, ref, alt = utils.parse_ot_vid(vid)
        except Exception:
            ch, pos, ref, alt = r.get("chrom"), r.get("pos"), r.get("ref"), r.get("alt")
        known = r.get("known_effector")
        gene = r.get("polaris_gene")
        hit = None
        if isinstance(known, str) and isinstance(gene, str):
            hit = (known == gene)
        recs.append(dict(
            disease=disease, locus=r.get("locus"), rsid=r.get("rsid"), vid=vid,
            chrom=str(ch), pos=int(pos) if pd.notna(pos) else None, ref=ref, alt=alt,
            pip=r.get("pip"), is95=bool(r.get("is95")) if pd.notna(r.get("is95")) else None,
            is_lead=bool(r.get("is_lead")) if pd.notna(r.get("is_lead")) else None,
            is_control=bool(r.get("is_control")) if pd.notna(r.get("is_control")) else None,
            coding=bool(r.get("coding_lead")) if pd.notna(r.get("coding_lead")) else None,
            P_functional=r.get("P_functional"), POLARIS=r.get("POLARIS"),
            polaris_gene=gene, known_effector=known if isinstance(known, str) else None,
            nearest_gene=r.get("nearest_gene") if isinstance(r.get("nearest_gene"), str) else None,
            l2g_gene=r.get("best_l2g_gene") if isinstance(r.get("best_l2g_gene"), str) else None,
            l2g_score=r.get("best_l2g_score"), gene_hit=hit,
            motif_tf=r.get("motif_top_disrupt_tf") if isinstance(r.get("motif_top_disrupt_tf"), str) else None,
            motif_dLLR=r.get("motif_top_disrupt_dLLR"), motif_pref=r.get("motif_top_disrupt_pref"),
            create_tf=r.get("motif_top_create_tf") if isinstance(r.get("motif_top_create_tf"), str) else None,
            phylop=r.get("phylop"), phylop_max=r.get("phylop_max"), phastcons470=r.get("phastcons470"),
            ccre=r.get("ccre_class") if isinstance(r.get("ccre_class"), str) else None,
            in_ccre=bool(r.get("in_ccre")) if pd.notna(r.get("in_ccre")) else None,
            gnomad_af=r.get("gnomad_af"), rarity=r.get("rarity"),
            distal_flag=int(r.get("distal_flag")) if pd.notna(r.get("distal_flag")) else 0,
            mechanism=r.get("mechanism") if isinstance(r.get("mechanism"), str) else mech(r),
        ))
    return recs


def export_fto():
    """FTO rs1421085 vignette: gene positions, ARID5 motif PWM + binding-site alleles, ODE dynamics."""
    POS, REF, ALT, F = 53767042, "T", "C", 25
    fg = pd.read_csv(PROC / "genes.csv")
    fg = fg[fg.locus == "FTO"].sort_values("dist_to_tss")
    genes = [dict(gene=r.gene, tss=int(r.tss), dist=int(r.dist_to_tss), nearest=bool(r.is_nearest),
                  role=("truth (distal)" if r.gene == "IRX3" else ("L2G call" if r.gene == "FTO"
                        else ("nearest" if r.is_nearest else "")))) for _, r in fg.iterrows()]
    jpath = str(config.EXT / "JASPAR2024_CORE_vertebrates_nr.txt")
    pwms = M.load_pwms(jpath); mo = M.load_jaspar(jpath)
    arid = next(mid for mid, p in pwms.items() if p["name"].lower().startswith("arid5"))
    pwm = pwms[arid]; w = pwm["w"]
    seq = apis.ucsc_sequence("chr16", POS - 1 - F, POS - 1 + F + 1)
    si = M.seq_to_idx(seq)
    best = max(range(max(0, F - w + 1), min(F, len(seq) - w) + 1),
               key=lambda s: pwm["llr"][si[s:s + w], np.arange(w)].sum())
    llr_ref = float(pwm["llr"][si[best:best + w], np.arange(w)].sum())
    si_alt = si.copy(); si_alt[F] = M.B2I[ALT]
    llr_alt = float(pwm["llr"][si_alt[best:best + w], np.arange(w)].sum())
    counts = mo[arid][1]; prob = (counts + 0.01) / (counts + 0.01).sum(0)
    ic = np.sum(prob * np.log2(prob / M.HG38_BG[:, None]), 0)
    voff = F - best
    site_ref = seq[best:best + w]; site_alt = site_ref[:voff] + ALT + site_ref[voff + 1:]
    res = DY.allelic_response(llr_ref, llr_alt, mode="repressor", tf_conc=1.0, Kd0=12.0)
    Kd0s = np.logspace(0, 1.7, 30)
    l2fc = [DY.allelic_response(llr_ref, llr_alt, mode="repressor", tf_conc=1.0, Kd0=float(k))["log2fc"] for k in Kd0s]
    sub = slice(None, None, 3)  # subsample ODE for lighter payload
    dump("fto.json", dict(
        variant=dict(rsid="rs1421085", chrom="16", pos=POS, ref=REF, alt=ALT),
        genes=genes,
        motif=dict(name=pwm["name"], id=arid, w=int(w),
                   prob=[[float(prob[b, i]) for i in range(w)] for b in range(4)], ic=[float(x) for x in ic],
                   var_offset=int(voff), site_ref=site_ref, site_alt=site_alt,
                   llr_ref=llr_ref, llr_alt=llr_alt, dLLR=llr_alt - llr_ref, p_ref=M.match_pvalue(pwm, llr_ref)),
        ode=dict(t=[float(x) for x in res["ref"]["t"][sub]],
                 x_ref=[float(x) for x in res["ref"]["x"][sub]], x_alt=[float(x) for x in res["alt"]["x"][sub]],
                 occ_ref=float(res["O_ref"]), occ_alt=float(res["O_alt"]),
                 log2fc=float(res["log2fc"]), direction=res["direction"],
                 sens_Kd0=[float(x) for x in Kd0s], sens_log2fc=[float(x) for x in l2fc])))


def main():
    t2d = pd.read_csv(PROC / "integrated.csv")
    cad = add_polaris_fields(pd.read_csv(PROC / "cad_integrated.csv"))
    cad["mechanism"] = cad.apply(mech, axis=1)
    m = json.load(open(TAB / "polaris_metrics.json"))
    cadm = json.load(open(TAB / "cad_metrics.json"))
    rob = json.load(open(TAB / "robustness.json"))

    # ---- variants ----
    variants = variant_records(t2d, "T2D") + variant_records(cad, "CAD")
    dump("variants.json", variants)

    # ---- per-locus summary ----
    loci = []
    for disease, df in [("T2D", t2d), ("CAD", cad)]:
        for L, g in df.groupby("locus"):
            g = g.sort_values("POLARIS", ascending=False)
            top = g.iloc[0]
            known = top.get("known_effector")
            gene = top.get("polaris_gene")
            loci.append(dict(
                disease=disease, locus=L, lead_rsid=g[g.is_lead == True]["rsid"].iloc[0] if (g.is_lead == True).any() else top["rsid"],
                n_variants=int(len(g)), top_rsid=top["rsid"], top_pip=top.get("pip"),
                top_POLARIS=top.get("POLARIS"), polaris_gene=gene,
                known_effector=known if isinstance(known, str) else None,
                nearest_gene=top.get("nearest_gene") if isinstance(top.get("nearest_gene"), str) else None,
                l2g_score=top.get("best_l2g_score"),
                gene_hit=(known == gene) if (isinstance(known, str) and isinstance(gene, str)) else None,
                is_control=bool(top.get("is_control")) if pd.notna(top.get("is_control")) else None,
                distal_flag=int(top.get("distal_flag")) if pd.notna(top.get("distal_flag")) else 0,
                mechanism=top.get("mechanism") if isinstance(top.get("mechanism"), str) else mech(top),
            ))
    dump("loci.json", loci)

    # ---- meta / headline ----
    def dstats(df, met, cadmet=None):
        POS, NEG = df.pip > 0.5, df.pip < 0.05
        return dict(
            n_loci=int(df.locus.nunique()), n_variants=int(len(df)),
            n_in95=int((df.is95 == True).sum()),
            n_controls=int(df[df.is_control == True].locus.nunique()),
        )
    meta = dict(
        generated_note="All numbers computed by the POLARIS pipeline from public data; see README.",
        diseases=dict(
            T2D=dict(name="Type 2 diabetes", **dstats(t2d, m),
                     gene_top1=m["gene_top1_POLARIS(=L2G)"], gene_nearest=m["gene_top1_nearest"],
                     ece=m["calib_ece"], brier=m["calib_brier"],
                     cons_auc=m["finding_conservation_enrichment"]["auc"],
                     cons_p=m["finding_conservation_enrichment"]["p"], cons_ci=rob["conservation_AUC_CI"],
                     rarity_auc=m["finding_rarity_inversion"]["auc"],
                     rarity_p=m["finding_rarity_inversion"]["p"], rarity_ci=rob["rarity_AUC_CI"],
                     mp_signal=m["oracle_spectrum"]["n_signal_dims"], mp_total=m["oracle_spectrum"]["n_channels"],
                     phylop_perm_p=rob["phyloP_weight_perm_p"]),
            CAD=dict(name="Coronary artery disease", n_loci=int(cad.locus.nunique()), n_variants=int(len(cad)),
                     n_controls=int(cad[cad.is_control == True].locus.nunique()),
                     gene_top1=cadm["gene_top1_L2G"], gene_nearest=cadm["gene_top1_nearest"],
                     cons_auc=cadm["conservation_auc"], cons_p=cadm["conservation_p"],
                     rarity_auc=cadm["rarity_auc"], rarity_p=cadm["rarity_p"]),
        ),
        sources=dict(gwas="Open Targets Platform", eqtl="GTEx v8", conservation="UCSC phyloP/phastCons",
                     ccre="ENCODE SCREEN", motifs="JASPAR 2024", af="gnomAD v4"),
    )
    dump("meta.json", meta)

    # ---- reliability (weights + CIs + perm p), both diseases ----
    labels = {"motif_motif_disruption": "TF motif disruption", "phylop_max": "conservation (phyloP)",
              "phastcons470": "conservation (phastCons470)", "ccre_enhancer": "enhancer (cCRE)",
              "splice_delta": "cryptic splice"}
    ci = rob["reliability_weight_CI"]
    rel = dict(perm_p=rob["phyloP_weight_perm_p"], channels=[
        dict(key=c, label=labels[c], weight=m["em_reliability"][c]["weight"],
             dprime=m["em_reliability"][c]["dprime"],
             ci_lo=ci[c][0], ci_med=ci[c][1], ci_hi=ci[c][2],
             cad_weight=cadm["em_reliability"].get(c)) for c in FCOLS])
    dump("reliability.json", rel)

    # ---- calibration reliability curve ----
    oof = np.load(PROC / "calib_oof.npy")
    xs, ys, ns = IN.reliability_curve(oof[:, 0], oof[:, 1], bins=6)
    dump("calibration.json", dict(ece=m["calib_ece"], brier=m["calib_brier"],
         points=[dict(pred=float(a), obs=float(b), n=int(c)) for a, b, c in zip(xs, ys, ns)]))

    # ---- Marchenko-Pastur oracle spectrum ----
    sp = m["oracle_spectrum"]
    dump("mp.json", dict(evals=sp["evals"], edge=sp["mp_hi"], gamma=sp["gamma"],
         n_signal=sp["n_signal_dims"], n_channels=sp["n_channels"]))

    # ---- findings: conservation & rarity distributions (both diseases) ----
    def dist(df, col):
        POS, NEG = df.pip > 0.5, df.pip < 0.05
        return dict(causal=[x for x in df.loc[POS, col].dropna().round(3).tolist()],
                    passenger=[x for x in df.loc[NEG, col].dropna().round(3).tolist()])
    findings = dict(
        T2D=dict(phylop=dist(t2d, "phylop_max"), rarity=dist(t2d, "rarity"),
                 cons_auc=m["finding_conservation_enrichment"]["auc"], rarity_auc=m["finding_rarity_inversion"]["auc"]),
        CAD=dict(phylop=dist(cad, "phylop_max"), rarity=dist(cad, "rarity"),
                 cons_auc=cadm["conservation_auc"], rarity_auc=cadm["rarity_auc"]))
    dump("findings.json", findings)

    # ---- real colocalization (honest: modest PP4, islet caveat) ----
    try:
        rc = pd.read_csv(TAB / "real_coloc.csv")
        best = rc.sort_values("PP4", ascending=False).groupby("locus").first().reset_index()
        dump("coloc.json", dict(note="GWAS log(PIP) x GTEx eQTL ABFs (coloc-SuSiE style). PP4 modest in "
             "available tissues — GTEx lacks pancreatic islets, the disease-relevant T2D tissue.",
             tests=int(len(rc)), loci=int(rc.locus.nunique()),
             rows=[dict(locus=r.locus, gene=r.gene, tissue=r.tissue, n=int(r.n_common), PP4=r.PP4)
                   for _, r in best.iterrows()]))
    except Exception as e:
        dump("coloc.json", dict(note="unavailable", err=str(e)[:80]))

    export_fto()
    print("exported:", sorted(p.name for p in OUT.glob("*.json")))
    print("variants:", len(variants), "| loci:", len(loci), "| T2D+CAD POLARIS fields computed")


if __name__ == "__main__":
    main()
