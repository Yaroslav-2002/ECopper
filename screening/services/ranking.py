"""The 10-point ranking of copper-resistance candidates.

Plain Python with no Django or pandas imports, so the web app (services/scoring.py) and the standalone
script (scripts/rank_copper_genes.py) share exactly the same rules. Edit SCORING to change the weights.
"""
import math
import re

# ============================================================================
# SCORING RULES  (edit the weights here; maximum total = 5 + 3 + 2 = 10)
# ============================================================================
SCORING = {
    # 1. BLAST similarity = identity x coverage / 100  -> 0-5 points
    #    (minimum similarity, points), checked top-down
    "blast_bins": [(80, 5), (60, 4), (40, 3), (20, 2), (0, 0)],
    "evalue_cutoff": 1e-5,  # if blast_evalue > cutoff (or missing) ...
    "evalue_cap": 1,        # ... the BLAST points are capped at this value
    # 2. BacMet evidence -> 0-3 points
    "bacmet": {"experimental": 3, "predicted": 1},
    # 3. Other metals -> 0-2 points
    "other_metal": 2,       # BacMet lists at least one metal besides copper
    "copper_only": 1,
    # Tiers: (minimum total score, label), checked top-down
    "tiers": [(8, "High"), (5, "Medium"), (float("-inf"), "Low")],
}
MAX_POINTS = {
    "blast": max(p for _, p in SCORING["blast_bins"]),
    "bacmet": max(SCORING["bacmet"].values()),
    "metals": max(SCORING["other_metal"], SCORING["copper_only"]),
}
MAX_SCORE = sum(MAX_POINTS.values())
TIERS = tuple(label for _, label in SCORING["tiers"])  # best first

# Metals recognised in BacMet compound lists (BacMet's metal list, including the metalloids
# As, Sb, Se, Te). Name -> chemical symbol.
METALS = {
    "copper": "Cu", "zinc": "Zn", "cadmium": "Cd", "cobalt": "Co", "nickel": "Ni",
    "lead": "Pb", "mercury": "Hg", "silver": "Ag", "gold": "Au", "iron": "Fe",
    "manganese": "Mn", "chromium": "Cr", "aluminium": "Al", "aluminum": "Al",
    "magnesium": "Mg", "molybdenum": "Mo", "tungsten": "W", "vanadium": "V",
    "bismuth": "Bi", "gallium": "Ga", "palladium": "Pd", "platinum": "Pt",
    "ruthenium": "Ru", "titanium": "Ti", "uranium": "U", "caesium": "Cs",
    "cesium": "Cs", "arsenic": "As", "antimony": "Sb", "selenium": "Se",
    "tellurium": "Te",
}

EMPTY_REF = {"", "-", "na", "n/a", "none"}
PMID_RE = re.compile(r"^(?:PMID:?\s*|https?://pubmed\.ncbi\.nlm\.nih\.gov/)?(\d{1,9})/?$", re.I)
DOI_RE = re.compile(r"10\.\d{4,9}/[^\s;]+")

LIMITATIONS = (
    "High sequence similarity to a BacMet gene suggests a similar function but does not prove it, and copper "
    "resistance is not the same as bioremediation ability (uptake, accumulation, growth in real effluents). "
    "The weights, bins and tier cut-offs are arbitrary choices for prioritising candidates, not a validated "
    "model. Confirm results experimentally."
)


def is_missing(value):
    return value is None or (isinstance(value, float) and math.isnan(value)) or not str(value).strip()


def metals_in(compounds):
    """Chemical symbols of the metals named in a BacMet compounds string."""
    text = "" if is_missing(compounds) else str(compounds)
    return {sym for name, sym in METALS.items()
            if re.search(rf"\b{name}\b", text, re.I) or re.search(rf"\b{sym}\b", text)}


def parse_references(refs):
    """Return (number of articles, list of links). Titles without a PMID/DOI get no link."""
    if is_missing(refs):
        return 0, []
    items = [r.strip() for r in str(refs).split(";") if r.strip().lower() not in EMPTY_REF]
    links = []
    for item in items:
        if m := PMID_RE.match(item):
            links.append(f"https://pubmed.ncbi.nlm.nih.gov/{m[1]}/")
        elif m := DOI_RE.search(item):
            links.append(f"https://doi.org/{m[0].rstrip('.,)')}")
    return len(items), links


def tier(total):
    return next(label for minimum, label in SCORING["tiers"] if total >= minimum)


def score(*, dataset, compounds, identity, coverage, evalue, references):
    """Score one gene. Numbers may be None (missing); missing data scores 0 and is listed in notes."""
    notes = []

    # 1. BLAST similarity
    if is_missing(identity) or is_missing(coverage):
        similarity, s_blast = None, 0
        notes.append("BLAST missing")
    else:
        similarity = identity * coverage / 100
        s_blast = next((p for cutoff, p in SCORING["blast_bins"] if similarity >= cutoff), 0)
        if is_missing(evalue):
            s_blast = min(s_blast, SCORING["evalue_cap"])
            notes.append("e-value missing (BLAST capped)")
        elif evalue > SCORING["evalue_cutoff"]:
            s_blast = min(s_blast, SCORING["evalue_cap"])
            notes.append(f"e-value > {SCORING['evalue_cutoff']:g} (BLAST capped)")

    # 2. BacMet evidence
    dataset = "" if is_missing(dataset) else str(dataset).strip().lower()
    s_bacmet = SCORING["bacmet"].get(dataset, 0)
    if dataset not in SCORING["bacmet"]:
        notes.append(f"unknown dataset '{dataset}'" if dataset else "BacMet dataset missing")

    # 3. Other metals
    if is_missing(compounds):
        s_metals = 0
        notes.append("compounds missing")
    else:
        metals = metals_in(compounds)
        s_metals = (SCORING["other_metal"] if metals - {"Cu"} else SCORING["copper_only"]) if "Cu" in metals else 0
        if "Cu" not in metals:
            notes.append("copper not listed in compounds")

    # Article links (information only, not scored)
    n_refs, links = parse_references(references)
    if n_refs == 0:
        notes.append("no references")
    elif len(links) < n_refs:
        notes.append(f"{n_refs - len(links)} reference(s) without PMID/DOI")

    total = s_blast + s_bacmet + s_metals
    return {
        "similarity": similarity, "score_blast": s_blast, "score_bacmet": s_bacmet, "score_metals": s_metals,
        "total_score": total, "tier": tier(total), "n_references": n_refs, "article_links": links,
        "notes": notes,
    }


def ranking_key(total_score, similarity, evalue):
    """Higher total first; ties broken by higher BLAST similarity, then lower e-value. Missing values last."""
    return (
        -(total_score if total_score is not None else -math.inf),
        -(similarity if similarity is not None else -math.inf),
        evalue if evalue is not None else math.inf,
    )
