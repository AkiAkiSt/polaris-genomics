"""Small shared helpers: variant-id parsing & cross-resource conversion."""
from __future__ import annotations
import math


def parse_ot_vid(vid):
    """'10_112998590_C_T' -> ('10', 112998590, 'C', 'T')."""
    c, p, r, a = vid.split("_")
    return c, int(p), r, a


def to_gnomad(vid):
    c, p, r, a = parse_ot_vid(vid)
    return f"{c}-{p}-{r}-{a}"


def to_gtex(vid):
    c, p, r, a = parse_ot_vid(vid)
    return f"chr{c}_{p}_{r}_{a}_b38"


def ucsc_chrom(chrom):
    return chrom if str(chrom).startswith("chr") else f"chr{chrom}"


def neglog10p(mantissa, exponent):
    if mantissa is None or exponent is None:
        return None
    try:
        return -(math.log10(mantissa) + exponent)
    except Exception:
        return None


FINEMAP_PRIORITY = {"SuSiE": 4, "SuSiE-inf": 4, "SuSie": 4, "PICS": 2, "conditional": 1}


def finemap_rank(method):
    if not method:
        return 0
    for k, v in FINEMAP_PRIORITY.items():
        if k.lower() in method.lower():
            return v
    return 1


def safe_log10(x, floor=1e-300):
    return math.log10(max(x, floor))
