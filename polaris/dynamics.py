"""Gene-regulatory dynamics: turn a motif dLLR into a predicted expression change.

A single-base variant changes a TF's binding log-odds LLR (bits). By the Berg-von Hippel
relation, the dissociation constant scales as Kd ~ 2^{-LLR}, so the variant changes the
TF's fractional occupancy at the enhancer. We propagate that through a small ODE cascade
(enhancer activity -> target mRNA) to a new steady state, giving a mechanistic, dynamical
prediction of the allelic expression fold-change.

For ARID5B (a repressor) at FTO rs1421085: the risk (alt) allele weakens ARID5B binding
(dLLR<0) -> derepression -> higher IRX3/IRX5 -> reduced adipocyte browning (Claussnitzer 2015).
"""
from __future__ import annotations
import numpy as np
from scipy.integrate import solve_ivp


def occupancy(tf_conc, llr, Kd0=1.0):
    """Fractional TF occupancy given binding log-odds (bits): Kd = Kd0 * 2^{-LLR}."""
    Kd = Kd0 * 2.0 ** (-llr)
    return tf_conc / (tf_conc + Kd)


def _rhs(t, y, O, mode, alpha, delta, beta, gamma):
    E, x = y
    drive = (1 - O) if mode == "repressor" else O          # repressor derepresses
    dE = alpha * drive - delta * E
    dx = beta * E - gamma * x
    return [dE, dx]


def simulate(llr, tf_conc=1.0, Kd0=1.0, mode="repressor",
             alpha=1.0, delta=1.0, beta=1.0, gamma=1.0, t_max=20.0, n=400, y0=None):
    O = occupancy(tf_conc, llr, Kd0)
    if y0 is None:
        y0 = [0.0, 0.0]
    sol = solve_ivp(_rhs, (0, t_max), y0, args=(O, mode, alpha, delta, beta, gamma),
                    t_eval=np.linspace(0, t_max, n), rtol=1e-7, atol=1e-9)
    drive = (1 - O) if mode == "repressor" else O
    x_ss = (beta / gamma) * (alpha / delta) * drive
    return dict(t=sol.t, E=sol.y[0], x=sol.y[1], occupancy=O, x_steady=x_ss)


def allelic_response(llr_ref, llr_alt, mode="repressor", tf_conc=1.0, Kd0=1.0, **kw):
    r = simulate(llr_ref, tf_conc, Kd0, mode, **kw)
    a = simulate(llr_alt, tf_conc, Kd0, mode, **kw)
    fc = (a["x_steady"] + 1e-9) / (r["x_steady"] + 1e-9)
    return dict(ref=r, alt=a, log2fc=float(np.log2(fc)), fold_change=float(fc),
                O_ref=r["occupancy"], O_alt=a["occupancy"],
                direction=("up" if fc > 1 else "down"))


if __name__ == "__main__":
    # FTO rs1421085: measured dLLR for ARID5 ~ -8.63 (llr_ref ~ 8.69, llr_alt ~ 0.06)
    # calibrate Kd0 so the reference repressor sits at intermediate-high occupancy
    res = allelic_response(8.69, 0.06, mode="repressor", tf_conc=1.0, Kd0=12.0)
    print("[regulatory ODE | FTO/ARID5B repressor]")
    print(f"  occupancy ref={res['O_ref']:.3f}  alt={res['O_alt']:.3f}")
    print(f"  predicted IRX3 steady state log2FC (alt/ref) = {res['log2fc']:+.2f} "
          f"({res['direction']}-regulated on risk allele)")
    print(f"  -> risk allele derepresses enhancer => IRX3/IRX5 UP (matches Claussnitzer 2015): "
          f"{'CONSISTENT' if res['direction']=='up' else 'INCONSISTENT'}")
