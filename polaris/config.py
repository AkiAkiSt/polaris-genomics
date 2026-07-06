"""Project constants: paths, T2D scope, curated loci, tissue panel, motif set."""
from __future__ import annotations
import pathlib

ROOT = pathlib.Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
RAW = DATA / "raw"
PROC = DATA / "processed"
EXT = DATA / "external"
FIG = ROOT / "results" / "figures"
TAB = ROOT / "results" / "tables"
for p in (RAW, PROC, EXT, FIG, TAB):
    p.mkdir(parents=True, exist_ok=True)

GENOME = "GRCh38"
T2D_EFO = "MONDO_0005148"   # type 2 diabetes mellitus (Open Targets)

# trait keywords used to keep only T2D / glycemic GWAS credible sets
T2D_TRAIT_KEYWORDS = [
    "type 2 diabetes", "type ii diabetes", "t2d", "diabetes mellitus",
    "fasting glucose", "fasting insulin", "glycated", "hba1c", "hemoglobin a1c",
    "2-hour glucose", "glucose", "insulin", "homa", "beta-cell", "proinsulin",
]

# Disease-relevant tissues/cell contexts for T2D (GTEx ids where available)
T2D_TISSUES = [
    "Pancreas", "Adipose_Subcutaneous", "Adipose_Visceral_Omentum", "Liver",
    "Muscle_Skeletal", "Whole_Blood", "Artery_Tibial", "Nerve_Tibial",
]

# ---------------------------------------------------------------------------
# Curated T2D lead variants (established loci). These guarantee inclusion of
# textbook NON-CODING mechanistic positive controls used for validation.
# 'control' = literature-known regulatory mechanism + effector gene.
# coding=True flags missense leads kept only as contrast (not the focus).
# ---------------------------------------------------------------------------
CURATED_LOCI = [
    # rsid          locus            effector(known/likely)  control  coding  note
    ("rs7903146",  "TCF7L2",        "TCF7L2",  True,  False, "islet enhancer; strongest T2D signal"),
    ("rs1421085",  "FTO",           "IRX3",    True,  False, "ARID5B motif disruption -> IRX3/IRX5 (Claussnitzer 2015)"),
    ("rs10830963", "MTNR1B",        "MTNR1B",  True,  False, "FOXA2-bound islet enhancer; melatonin receptor"),
    ("rs972283",   "KLF14",         "KLF14",   True,  False, "imprinted master trans-eQTL regulator (Small 2011)"),
    ("rs11708067", "ADCY5",         "ADCY5",   True,  False, "islet enhancer; insulin secretion"),
    ("rs10811661", "CDKN2A_CDKN2B", "CDKN2B",  True,  False, "9p21 ANRIL/CDKN2B-AS1 regulatory"),
    ("rs13266634", "SLC30A8",       "SLC30A8", True,  True,  "R325W (coding) zinc transporter; islet"),
    ("rs4402960",  "IGF2BP2",       "IGF2BP2", False, False, "intronic; well-replicated"),
    ("rs7754840",  "CDKAL1",        "CDKAL1",  False, False, "intronic; islet"),
    ("rs7578597",  "THADA",         "THADA",   False, False, "missense/regulatory"),
    ("rs1111875",  "HHEX_IDE",      "HHEX",    False, False, "intergenic HHEX/IDE"),
    ("rs864745",   "JAZF1",         "JAZF1",   False, False, "intronic"),
    ("rs1801214",  "WFS1",          "WFS1",    False, False, "Wolfram syndrome gene; regulatory"),
    ("rs4430796",  "HNF1B",         "HNF1B",   True,  False, "TCF2; MODY gene; regulatory"),
    ("rs340874",   "PROX1",         "PROX1",   False, False, "fasting glucose; islet"),
    ("rs560887",   "G6PC2",         "G6PC2",   True,  False, "splice/regulatory; fasting glucose"),
    ("rs10923931", "NOTCH2",        "NOTCH2",  False, False, "intronic"),
    ("rs12779790", "CDC123_CAMK1D", "CAMK1D",  True,  False, "enhancer -> CAMK1D (Fogarty 2014)"),
    ("rs896854",   "TP53INP1",      "TP53INP1",False, False, "intergenic"),
    ("rs11063069", "CCND2",         "CCND2",   False, False, "cell cycle; islet"),
    ("rs1552224",  "ARAP1_CENTD2",  "ARAP1",   True,  False, "proinsulin; ARAP1/STARD10 enhancer"),
    ("rs8042680",  "PRC1",          "PRC1",    False, False, "intergenic"),
    ("rs17791513", "TLE1",          "TLE1",    False, False, "intergenic"),
    ("rs7177055",  "HMG20A",        "HMG20A",  False, False, "islet TF"),
    ("rs2334499",  "DUSP8_HCCA2",   "DUSP8",   True,  False, "imprinted; parent-of-origin"),
    ("rs231362",   "KCNQ1",         "KCNQ1",   True,  False, "imprinted KCNQ1OT1 regulatory"),
    ("rs5215",     "KCNJ11_ABCC8",  "KCNJ11",  False, True,  "E23K (coding) K-ATP channel"),
    ("rs459193",   "ANKRD55",       "ANKRD55", False, False, "intergenic insulin resistance"),
    ("rs10401969", "CILP2_TM6SF2",  "TM6SF2",  False, False, "lipid/T2D pleiotropy"),
    ("rs6815464",  "MAEA",          "MAEA",    False, False, "East-Asian T2D"),
]

# ---- Second disease: Coronary Artery Disease (for generalization / G1) ----
CAD_EFO = "EFO_0001645"
CAD_TRAIT_KEYWORDS = ["coronary artery disease", "coronary heart disease", "myocardial infarction",
                      "coronary", "cad", "atherosclerosis", "ischemic heart"]
CAD_TISSUES = ["Artery_Aorta", "Artery_Coronary", "Artery_Tibial", "Liver",
               "Whole_Blood", "Adipose_Subcutaneous"]
CAD_LOCI = [
    ("rs12740374", "SORT1",     "SORT1",   True,  False, "1p13 C/EBP site creation -> hepatic SORT1 (Musunuru 2010)"),
    ("rs9349379",  "PHACTR1",   "EDN1",    True,  False, "distal: regulates EDN1, not nearest PHACTR1 (Gupta 2017)"),
    ("rs1333049",  "CDKN2B_AS1","CDKN2B",  True,  False, "9p21 ANRIL locus"),
    ("rs11206510", "PCSK9",     "PCSK9",   False, False, "LDL-C; near PCSK9"),
    ("rs6511720",  "LDLR",      "LDLR",    True,  False, "LDLR regulatory; LDL-C"),
    ("rs10455872", "LPA",       "LPA",     True,  False, "Lp(a); strong CAD"),
    ("rs2107595",  "HDAC9",     "HDAC9",   False, False, "HDAC9/TWIST1"),
    ("rs17114036", "PLPP3",     "PLPP3",   False, False, "PPAP2B endothelial enhancer"),
    ("rs1746048",  "CXCL12",    "CXCL12",  False, False, "intergenic CXCL12"),
    ("rs1412444",  "LIPA",      "LIPA",    False, False, "lysosomal acid lipase"),
    ("rs7692387",  "GUCY1A1",   "GUCY1A1", False, False, "sGC; vascular"),
    ("rs4845625",  "IL6R",      "IL6R",    False, False, "IL6R inflammation"),
    ("rs2891168",  "CDKN2B_AS1b","CDKN2B", True,  False, "9p21 second signal"),
    ("rs12936587", "UBE2Z",     "UBE2Z",   False, False, "RAI1/UBE2Z region"),
]

# JASPAR motif families especially relevant to islet/adipose/metabolic enhancers.
# (used to prioritise which PWMs to scan first; full CORE set also scanned)
PRIORITY_TF_FAMILIES = [
    "FOXA1", "FOXA2", "FOXO1", "PDX1", "NEUROD1", "NKX6-1", "MAFA", "MAFB",
    "PAX6", "RFX6", "HNF1A", "HNF1B", "HNF4A", "GATA4", "GATA6", "TCF7L2",
    "ISL1", "PPARG", "CEBPA", "CEBPB", "ARID5B", "NFKB1", "STAT3", "CTCF",
    "TBX", "ETS1", "SP1", "KLF", "EGR1", "USF1", "MEF2A",
]

VIS = dict(
    navy="#0b1f3a", ink="#10263f", teal="#1f9e89", amber="#f0a202",
    crimson="#b3123a", slate="#5b6b7b", mist="#dce6ef", paper="#fbfcfe",
    violet="#6a4c93", grass="#3a8d5b",
)
