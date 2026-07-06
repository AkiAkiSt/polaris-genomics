"""Validated clients for all public genomics APIs used by POLARIS.

Sources (all free, no institutional approval):
  Open Targets Platform GraphQL  - credible sets/PIPs, L2G, colocalisation, variant index
  Ensembl REST                   - sequence, VEP, LD (1000G), regulatory overlap
  GTEx Portal API v2             - cis-eQTL / sQTL effect sizes by gene+tissue / variant
  gnomAD GraphQL                 - allele frequency & constraint
  UCSC API                       - reference sequence, phyloP/phastCons, ENCODE cCREs
  JASPAR API                     - TF position frequency matrices (CORE 2024)
  NCBI E-utilities               - ClinVar clinical labels
"""
from __future__ import annotations
from . import cache

OT = "https://api.platform.opentargets.org/api/v4/graphql"
ENS = "https://rest.ensembl.org"
GTEX = "https://gtexportal.org/api/v2"
GNOMAD = "https://gnomad.broadinstitute.org/api"
UCSC = "https://api.genome.ucsc.edu"
JASPAR = "https://jaspar.elixir.no/api/v1"
EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# --------------------------------------------------------------------------- #
# Open Targets
# --------------------------------------------------------------------------- #
def ot_resolve_variant(rsid):
    """rsID -> Open Targets variant id (chrom_pos_ref_alt, GRCh38)."""
    q = '{search(queryString:"%s",entityNames:["variant"]){hits{id name entity}}}' % rsid
    hits = cache.graphql(OT, q)["data"]["search"]["hits"]
    for h in hits:
        if h["entity"] == "variant":
            return h["id"]
    return None


VARIANT_FULL_Q = """
query V($vid:String!){
 variant(variantId:$vid){
  id rsIds chromosome position referenceAllele alternateAllele
  mostSevereConsequence{id label}
  alleleFrequencies{populationName alleleFrequency}
  variantEffect{method assessment score normalisedScore}
  credibleSets(page:{index:0,size:50}){
    count
    rows{
      studyType finemappingMethod confidence credibleSetIndex
      study{id traitFromSource publicationFirstAuthor publicationDate studyType nCases nControls nSamples pubmedId}
      locus(page:{index:0,size:200}){count rows{posteriorProbability is95CredibleSet logBF beta standardError variant{id rsIds position}}}
      l2GPredictions(page:{index:0,size:8}){rows{score target{id approvedSymbol}}}
      colocalisation(page:{index:0,size:25}){rows{colocalisationMethod h3 h4 clpp rightStudyType
        otherStudyLocus{studyType study{id traitFromSource projectId} variant{id rsIds}}}}
    }
  }
 }
}"""


def ot_variant_full(vid):
    return cache.graphql(OT, VARIANT_FULL_Q, {"vid": vid})["data"]["variant"]


CORE_Q = """query V($vid:String!){variant(variantId:$vid){
  id rsIds chromosome position referenceAllele alternateAllele
  mostSevereConsequence{id label}
  alleleFrequencies{populationName alleleFrequency}
  variantEffect{method assessment score normalisedScore}}}"""


def ot_variant_core(vid):
    return cache.graphql(OT, CORE_Q, {"vid": vid})["data"]["variant"]


# Two-tier credible-set retrieval (OT enforces a query-cost cap):
#   light pass  = cheap scalar page -> locate the right T2D credible set's index
#   detail pass = targeted size:1 page at that index -> full members/L2G/coloc
CS_LIGHT_Q = """query L($vid:String!,$i:Int!,$sz:Int!){variant(variantId:$vid){
 credibleSets(page:{index:$i,size:$sz}){count rows{
   studyLocusId studyType finemappingMethod confidence credibleSetIndex
   study{id traitFromSource publicationFirstAuthor studyType nSamples nCases nControls pubmedId}}}}}"""

CS_DETAIL_Q = """query D($vid:String!,$k:Int!){variant(variantId:$vid){credibleSets(page:{index:$k,size:1}){rows{
   studyLocusId studyType finemappingMethod confidence
   study{id traitFromSource publicationFirstAuthor publicationDate studyType nSamples nCases nControls pubmedId}
   locus(page:{index:0,size:150}){count rows{posteriorProbability is95CredibleSet logBF beta standardError variant{id rsIds position}}}
   l2GPredictions(page:{index:0,size:10}){rows{score target{id approvedSymbol}}}
   colocalisation(page:{index:0,size:40}){rows{colocalisationMethod h3 h4 clpp rightStudyType
     otherStudyLocus{studyType study{id traitFromSource projectId} variant{id rsIds}}}}}}}}"""


def ot_cs_light(vid, i=0, sz=100):
    return cache.graphql(OT, CS_LIGHT_Q, {"vid": vid, "i": i, "sz": sz})["data"]["variant"]["credibleSets"]


def ot_cs_light_all(vid, sz=100, cap=2000):
    """All credible sets (scalar fields) with a stable global index '_idx'."""
    first = ot_cs_light(vid, 0, sz)
    total = min(first.get("count", 0), cap)
    rows = []
    for off, r in enumerate(first["rows"]):
        r["_idx"] = off; rows.append(r)
    i = sz
    while i < total:
        page = ot_cs_light(vid, i // sz, sz)
        for off, r in enumerate(page["rows"]):
            r["_idx"] = i + off; rows.append(r)
        i += sz
    return rows


def ot_cs_detail_at(vid, k):
    rows = cache.graphql(OT, CS_DETAIL_Q, {"vid": vid, "k": int(k)})["data"]["variant"]["credibleSets"]["rows"]
    return rows[0] if rows else None


def ot_cs_detail_verified(vid, k, slid):
    """Fetch detail at global index k; verify studyLocusId, scanning +-2 on drift."""
    for kk in (k, k - 1, k + 1, k - 2, k + 2):
        if kk is None or kk < 0:
            continue
        d = ot_cs_detail_at(vid, kk)
        if d and d.get("studyLocusId") == slid:
            return d
    return ot_cs_detail_at(vid, k)


def ot_disease_targets(efo, size=50):
    q = """query D($efo:String!,$n:Int!){disease(efoId:$efo){name
      associatedTargets(page:{index:0,size:$n}){count rows{score target{id approvedSymbol}}}}}"""
    return cache.graphql(OT, q, {"efo": efo, "n": size})["data"]["disease"]


# --------------------------------------------------------------------------- #
# Ensembl REST
# --------------------------------------------------------------------------- #
def ens_sequence(chrom, start, end, strand=1):
    """1-based inclusive reference sequence (GRCh38)."""
    url = f"{ENS}/sequence/region/human/{chrom}:{start}..{end}:{strand}"
    return cache.get_text(url, headers={"Content-Type": "text/plain"})


def ens_vep_rsid(rsid):
    url = f"{ENS}/vep/human/id/{rsid}"
    return cache.get_json(url, headers={"Content-Type": "application/json"})


def ens_vep_region(chrom, pos, ref, alt):
    # region form: chr:start-end:strand/allele
    url = f"{ENS}/vep/human/region/{chrom}:{pos}-{pos}:1/{alt}"
    return cache.get_json(url, headers={"Content-Type": "application/json"})


def ens_ld(rsid, pop="1000GENOMES:phase_3:EUR", r2=0.2):
    url = f"{ENS}/ld/human/{rsid}/{pop}"
    return cache.get_json(url, params={"r2": r2}, headers={"Content-Type": "application/json"})


def ens_ld_region(chrom, start, end, pop="1000GENOMES:phase_3:EUR", r2=0.1):
    url = f"{ENS}/ld/human/region/{chrom}:{start}..{end}/{pop}"
    return cache.get_json(url, params={"r2": r2}, headers={"Content-Type": "application/json"})


def ens_regulatory(chrom, start, end):
    url = f"{ENS}/overlap/region/human/{chrom}:{start}-{end}"
    return cache.get_json(url, params={"feature": "regulatory"},
                          headers={"Content-Type": "application/json"})


def ens_overlap_genes(chrom, start, end):
    url = f"{ENS}/overlap/region/human/{chrom}:{start}-{end}"
    return cache.get_json(url, params={"feature": "gene"},
                          headers={"Content-Type": "application/json"})


# --------------------------------------------------------------------------- #
# GTEx
# --------------------------------------------------------------------------- #
def gtex_gene(symbol):
    d = cache.get_json(f"{GTEX}/reference/gene", params={"geneId": symbol})
    return d["data"][0] if d.get("data") else None


def gtex_eqtl_gene_tissue(gencode_id, tissue, max_items=2000):
    url = f"{GTEX}/association/singleTissueEqtl"
    return cache.get_json(url, params={"gencodeId": gencode_id, "tissueSiteDetailId": tissue,
                                       "datasetId": "gtex_v8", "page": 0, "itemsPerPage": max_items})


def gtex_eqtl_variant(variant_id, max_items=2000):
    """variant_id like chr10_112998590_C_T_b38"""
    url = f"{GTEX}/association/singleTissueEqtl"
    return cache.get_json(url, params={"variantId": variant_id, "datasetId": "gtex_v8",
                                       "page": 0, "itemsPerPage": max_items})


def gtex_dyneqtl(variant_id, gencode_id, tissue):
    """On-demand cis-eQTL summary for ANY variant-gene-tissue: nes, error(SE), pValue."""
    d = cache.get_json(f"{GTEX}/association/dyneqtl",
                       params={"variantId": variant_id, "gencodeId": gencode_id,
                               "tissueSiteDetailId": tissue})
    return {k: d.get(k) for k in ("nes", "error", "pValue", "tStatistic", "maf")}


def gtex_sqtl_variant(variant_id, max_items=2000):
    url = f"{GTEX}/association/singleTissueSqtl"
    return cache.get_json(url, params={"variantId": variant_id, "datasetId": "gtex_v8",
                                       "page": 0, "itemsPerPage": max_items})


# --------------------------------------------------------------------------- #
# gnomAD
# --------------------------------------------------------------------------- #
def gnomad_variant(variant_id, dataset="gnomad_r4"):
    """variant_id like 10-112998590-C-T"""
    q = """query V($v:String!,$d:DatasetId!){variant(variantId:$v,dataset:$d){
      variantId genome{ac an af} exome{ac an af}}}"""
    try:
        d = cache.graphql(GNOMAD, q, {"v": variant_id, "d": dataset})
        return d.get("data", {}).get("variant")
    except Exception:
        return None


# --------------------------------------------------------------------------- #
# UCSC
# --------------------------------------------------------------------------- #
def ucsc_sequence(chrom, start, end, genome="hg38"):
    """0-based half-open. Returns uppercase DNA string."""
    d = cache.get_json(f"{UCSC}/getData/sequence",
                       params={"genome": genome, "chrom": chrom, "start": start, "end": end})
    return d["dna"].upper()


def ucsc_track(track, chrom, start, end, genome="hg38"):
    d = cache.get_json(f"{UCSC}/getData/track",
                       params={"genome": genome, "track": track, "chrom": chrom,
                               "start": start, "end": end})
    return d.get(track, d)


def ucsc_phylop(chrom, start, end):
    return ucsc_track("phyloP100way", chrom, start, end)


def ucsc_phastcons(chrom, start, end):
    return ucsc_track("phastCons100way", chrom, start, end)


def ucsc_ccre(chrom, start, end):
    return ucsc_track("encodeCcreCombined", chrom, start, end)


# --------------------------------------------------------------------------- #
# JASPAR
# --------------------------------------------------------------------------- #
def jaspar_list(tax_group="vertebrates", collection="CORE", page_size=1000, max_pages=6):
    out = []
    url = f"{JASPAR}/matrix/?page_size={page_size}&tax_group={tax_group}&collection={collection}&version=latest"
    for _ in range(max_pages):
        d = cache.get_json(url)
        out.extend(d.get("results", []))
        url = d.get("next")
        if not url:
            break
    return out


def jaspar_pfm(matrix_id):
    return cache.get_json(f"{JASPAR}/matrix/{matrix_id}/")


# --------------------------------------------------------------------------- #
# ClinVar via NCBI E-utilities
# --------------------------------------------------------------------------- #
def clinvar_region(chrom, start, end, retmax=200):
    term = f"{chrom}[chr] AND {start}:{end}[chrpos38] AND (\"clinsig pathogenic\"[Properties] OR \"clinsig likely pathogenic\"[Properties])"
    s = cache.get_json(f"{EUTILS}/esearch.fcgi",
                       params={"db": "clinvar", "term": term, "retmode": "json", "retmax": retmax})
    return s.get("esearchresult", {})


def clinvar_summary(ids):
    if not ids:
        return {}
    d = cache.get_json(f"{EUTILS}/esummary.fcgi",
                       params={"db": "clinvar", "id": ",".join(ids), "retmode": "json"})
    return d.get("result", {})
