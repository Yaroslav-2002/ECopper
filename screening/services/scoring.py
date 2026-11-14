"""Rule-based ranking of BLAST hits against the copper-resistance reference set (deterministic, no LLM)."""
import re
from dataclasses import dataclass

from ..models import BacMetGene
from .blast import run_blast_many

CU_RE = re.compile(r"(?i:\bcopper\b)|\bCu\b")

# Step 1: significance gate on the best hit (highest bitscore).
MAX_EVALUE = 1e-10
MIN_BITSCORE = 50  # independent of database size, unlike the E-value
MIN_COVERAGE = 80  # % of the reference gene covered by the alignment
MIN_IDENTITY = 40  # below this, function cannot be inferred from homology

NO_EVIDENCE = "No significant evidence"
NOT_TRANSFERABLE = "Function not transferable"
CATEGORIES = ("High potential", "Moderate", "Low", NOT_TRANSFERABLE, NO_EVIDENCE)  # best first (Step 5)
UNCERTAINTIES = ("Low", "Moderate", "High")

FLAG_INDIRECT = "indirect role"
FLAG_DISTANT = "evidence from distant homolog"  # criterion D; never raised while D is skipped
FLAG_SPECIFICITY = "possible different metal specificity"

CRITERIA = {  # key: (label, max points); maximum 11
    "A": ("Identity", 3), "B": ("Copper specificity", 2), "C": ("Mechanism", 2),
    "D": ("Experimental evidence", 2), "E": ("Phylogenetic consistency", 1), "F": ("Taxonomic proximity", 1),
}

# Protein family and mechanism points, matched against the protein DESCRIPTION (never the gene name).
# First matching rule wins, so regulators are recognised before the Cu or transport words they mention.
FAMILIES = [
    ("regulator/sensor", r"regulat|repressor|activator|sensor|sensing|kinase|transcription", 0),
    ("multicopper oxidase", r"copper oxidase", 2),
    ("metallothionein", r"metallothionein", 2),
    ("Cu chaperone", r"chaperone", 2),
    ("Cu-binding protein", r"(copper|\bcu\S*)[- ]binding|cuproprotein|periplasmic binding", 2),
    ("P-type ATPase", r"p-type|atpase", 1),
    ("RND/CBA efflux", r"\brnd\b|efflux system|membrane fusion", 1),
    ("transporter", r"transport|efflux|pump|permease|export(?!ed)", 1),
    ("homeostasis protein", r"homeostasis", 0),
]
SPECIFICITY_FAMILIES = {"P-type ATPase": 60, "multicopper oxidase": 60}  # min identity for criterion E

LAB_CHECK = "Recommended lab check: growth at elevated Cu (MIC) and Cu removal measured by ICP-MS."
DISCLAIMER = (
    "Annotation uses the standard best-hit method with literature thresholds (E <= 1e-10, coverage >= 80%, "
    "identity 40/50/80%). The score is an expert heuristic for prioritising candidates for lab validation, "
    "not an established standard. The best-hit method can miss novel genes with low similarity. "
    "Results must be confirmed experimentally."
)


def has_copper(text):
    return bool(CU_RE.search(text or ""))


@dataclass
class Ref:
    """What the database knows about one reference gene."""
    gene_name: str
    organism: str = ""
    compounds: str = ""
    descriptions: tuple = ()  # UniProt protein name, NCBI description
    ncbi_gene_url: str = ""


GENERIC_FAMILIES = {"transporter", "homeostasis protein"}  # used only if no description names a specific family


def classify(descriptions):
    """(family, mechanism points, description used): specific families first, then generic ones."""
    for generic in (False, True):
        for text in descriptions:
            for family, pattern, points in FAMILIES:
                if text and (family in GENERIC_FAMILIES) == generic and re.search(pattern, text, re.I):
                    return family, points, text
    return None, None, next((t for t in descriptions if t), "")


def metals(compounds):
    """Element symbols in a BacMet compound list, e.g. 'Copper (Cu), Zinc (Zn)' -> {'Cu', 'Zn'}."""
    return set(re.findall(r"\(([A-Z][a-z]?)\)", compounds or ""))


def format_evalue(evalue, digits=1):
    """BLAST prints E-values below ~1e-180 as 0.0 (they are too small to print), so show 0 as '< 1e-180'."""
    if evalue is None:
        return ""
    return "< 1e-180" if evalue == 0 else f"{evalue:.{digits}e}"


def gate_checks(evalue, bitscore, ref_coverage):
    """Step 1 checks as (label, value, requirement, passed)."""
    return [
        ("E-value", format_evalue(evalue),f"≤ {MAX_EVALUE:.0e}", evalue <= MAX_EVALUE),
        ("Bitscore", f"{bitscore:.1f}", f"≥ {MIN_BITSCORE}", bitscore >= MIN_BITSCORE),
        ("Reference coverage", f"{ref_coverage:.0f}%", f"≥ {MIN_COVERAGE}%", ref_coverage >= MIN_COVERAGE),
    ]


def ranking_key(category, score_pct, uncertainty, identity):
    """Step 5: category, then score %, then lower uncertainty, then higher identity."""
    return (
        CATEGORIES.index(category) if category in CATEGORIES else len(CATEGORIES),
        -(score_pct or 0),
        UNCERTAINTIES.index(uncertainty) if uncertainty in UNCERTAINTIES else len(UNCERTAINTIES),
        -(identity or 0),
    )


def screen_many(sequences, seq_type="protein", query_ids=None):
    """Screen sequences with one BLAST run. Returns one result dict per sequence, in order."""
    query_ids = query_ids or [f"Query {i + 1}" for i in range(len(sequences))]
    raw_lists = run_blast_many(sequences, seq_type)
    refs = _references({h["bacmet_id"] for raw in raw_lists for h in raw})
    return [evaluate(qid, raw, refs) for qid, raw in zip(query_ids, raw_lists)]


def _references(ids):
    return {
        bacmet_id: Ref(
            gene.gene_name, gene.organism, gene.compounds,
            (gene.protein_name, gene.ncbi_description) if gene.protein_name or gene.ncbi_description else (),
            gene.ncbi_gene_url,
        )
        for bacmet_id, gene in BacMetGene.objects.in_bulk(ids, field_name="bacmet_id").items()
    }


def _annotate(h, refs):
    ref = refs.get(h["bacmet_id"]) or Ref(h["bacmet_id"])
    family, _, description = classify(ref.descriptions)
    ref_cov = 100 * (abs(h["subject_end"] - h["subject_start"]) + 1) / h["subject_length"]
    return {
        "bacmet_id": h["bacmet_id"], "gene_name": ref.gene_name, "organism": ref.organism,
        "description": description[:255], "family": family or "", "ncbi_gene_url": ref.ncbi_gene_url,
        "identity": h["identity"], "align_length": h["align_length"], "evalue": h["evalue"],
        "bitscore": h["bitscore"], "query_start": h["query_start"], "query_end": h["query_end"],
        "ref_coverage": min(100.0, ref_cov),
        "passes_gate": all(ok for *_, ok in gate_checks(h["evalue"], h["bitscore"], ref_cov)),
    }


def evaluate(query_id, raw_hits, refs):
    """Apply the ranking specification to one query's BLAST hits."""
    hits = sorted((_annotate(h, refs) for h in raw_hits), key=lambda h: -h["bitscore"])
    for rank, h in enumerate(hits, 1):
        h["rank"] = rank
    result = {"hits": hits, "score": None, "max_score": None, "score_pct": None,
              "uncertainty": "", "flags": [], "breakdown": []}
    best = hits[0] if hits else None
    name = _name(best) if best else ""

    if not best or not best["passes_gate"]:
        return {**result, "category": NO_EVIDENCE, "conclusion": (
            f"{query_id} has no significant match to a Cu-resistance reference (gate: E ≤ 1e-10, bitscore ≥ 50, "
            f"≥ 80% coverage of the reference), so no Cu bioremediation function can be inferred.")}
    if best["identity"] < MIN_IDENTITY:
        return {**result, "category": NOT_TRANSFERABLE, "conclusion": (
            f"{query_id} matches {name} at only {best['identity']:.0f}% identity: below 40%, "
            f"its function cannot be reliably inferred from homology.")}

    ref = refs.get(best["bacmet_id"]) or Ref(best["bacmet_id"])
    flags = []
    points = {
        "A": _identity(best),
        "B": _specificity(ref),
        "C": _mechanism(ref, flags),
        "D": (None, "Skipped: the database has no experimental bioremediation evidence reference."),
        "E": _consistency(hits, best, flags),
        "F": (None, "Skipped: no user organism is collected and the database has no taxonomic class "
                    "for reference organisms."),
    }
    breakdown = [{"key": k, "label": CRITERIA[k][0], "points": p, "max": CRITERIA[k][1], "note": note}
                 for k, (p, note) in points.items()]
    score = sum(r["points"] for r in breakdown if r["points"] is not None)
    max_score = sum(r["max"] for r in breakdown if r["points"] is not None)
    pct = 100 * score / max_score
    category = "High potential" if pct >= 75 else "Moderate" if pct >= 50 else "Low"
    uncertainty = _uncertainty(best["identity"], flags)
    level = {"High potential": "high", "Moderate": "moderate", "Low": "low"}[category]
    return {
        **result, "score": score, "max_score": max_score, "score_pct": pct, "category": category,
        "uncertainty": uncertainty, "flags": flags, "breakdown": breakdown,
        "conclusion": (f"{query_id} is {best['identity']:.0f}% identical to {name}: {level} Cu bioremediation "
                       f"potential, {uncertainty.lower()} uncertainty. {LAB_CHECK}"),
    }


def _name(hit):
    desc, gene = hit["description"], hit["gene_name"]
    if not desc:
        return gene
    return desc if gene.split("/")[0].lower() in desc.lower() else f"{desc} ({gene})"


def _identity(best):
    i = best["identity"]
    pts = 3 if i >= 80 else 2 if i >= 50 else 1
    return pts, f"{i:.1f}% identity to the best hit (≥ 80% = 3, 50–80% = 2, 40–50% = 1)."


def _specificity(ref):
    found = metals(ref.compounds)
    if not found:
        return None, "Skipped: no metal list for the best hit in the database."
    listed = ", ".join(sorted(found))
    if found == {"Cu"}:
        return 2, "Metal list of the best hit: Cu only."
    return (1 if "Cu" in found else 0), f"Metal list of the best hit: {listed}."


def _mechanism(ref, flags):
    family, pts, desc = classify(ref.descriptions)
    if family is None:
        return None, (f"Skipped: no mechanism-describing protein description for {ref.gene_name} in the database"
                      + (f" ('{desc}')" if desc else "") + "; gene names are not used.")
    if pts == 0:
        flags.append(FLAG_INDIRECT)
    return pts, f"Description '{desc}' → {family}."


def _consistency(hits, best, flags):
    top = [h for h in hits[:5] if h["passes_gate"]]
    families = {h["family"] for h in top}
    if "" in families:
        return None, "Skipped: protein family unknown for at least one significant top-5 hit (no usable description)."
    if len(families) > 1:
        flags.append(FLAG_SPECIFICITY)
        return 0, f"Significant top-5 hits span several families: {', '.join(sorted(families))}."
    family = families.pop()
    min_id = SPECIFICITY_FAMILIES.get(family)
    if min_id and best["identity"] < min_id:
        flags.append(FLAG_SPECIFICITY)
        return 0, f"All {len(top)} significant top-5 hits are {family}, but best-hit identity is below {min_id}%."
    return 1, f"All {len(top)} significant top-5 hits are {family}."


def _uncertainty(identity, flags):
    if identity < 50 or FLAG_SPECIFICITY in flags or len(flags) >= 2:
        return "High"
    if identity < 80 or flags:
        return "Moderate"
    return "Low"
