import re

from ..models import BacMetGene, NucleotideGene
from . import europepmc
from .blast import run_blast_many

CU_RE = re.compile(r"(?i:\bcopper\b)|\bCu\b")
SYSTEM_PREFIXES = ("cop", "cus", "pco", "sil", "cue", "cut")
NO_HIT_VERDICT = "No Cu-resistance homolog found"


def has_copper(text):
    return bool(CU_RE.search(text or ""))


def primary_gene(gene_name):
    return gene_name.split("/")[0].strip()


# C1 identity thresholds (max 3 points). DNA identity runs higher than protein identity for the
# same degree of relatedness, so blastn hits need more identity for the same points.
IDENTITY_STEPS = {
    "protein": ((80, 3), (65, 2), (50, 1)),
    "nucleotide": ((90, 3), (80, 2), (70, 1)),
}
# blastn hits also need the alignment to cover at least this much of the reference gene to earn C1,
# so a short conserved stretch is not mistaken for the whole gene.
MIN_NT_COVERAGE = 50


def score_identity(identity, seq_type="protein", coverage=None):
    if seq_type == "nucleotide" and (coverage or 0) < MIN_NT_COVERAGE:
        return 0
    return next((pts for threshold, pts in IDENTITY_STEPS[seq_type] if identity >= threshold), 0)


def verdict_for(score, has_hits=True):
    if not has_hits:
        return NO_HIT_VERDICT
    if score >= 6:
        return "Strong Cu bioremediation marker"
    if score >= 3:
        return "Possible marker"
    return "Weak / unlikely"


def screen(sequence, seq_type="protein"):
    return screen_many([sequence], seq_type)[0]


def screen_many(sequences, seq_type="protein"):
    """Screen sequences with one BLAST run. Returns (best, verdict, hits) per sequence, in order."""
    raw_lists = run_blast_many(sequences, seq_type)
    model = NucleotideGene if seq_type == "nucleotide" else BacMetGene
    genes = model.objects.in_bulk({h["bacmet_id"] for raw in raw_lists for h in raw}, field_name="bacmet_id")
    pmc_cache = {}  # one Europe PMC lookup per gene, even if it fails
    return [_score(raw, genes, pmc_cache, seq_type) for raw in raw_lists]


def _score(raw, genes, pmc_cache, seq_type):
    hits = []
    for h in raw:
        h = dict(h)
        subject_length = h.pop("subject_length")
        gene = genes.get(h["bacmet_id"])
        gene_name = gene.gene_name if gene else h["bacmet_id"]
        key = primary_gene(gene_name)
        if key not in pmc_cache:
            pmc_cache[key] = europepmc.hit_count(key)
        pmc = pmc_cache[key]
        ref_length = getattr(gene, "length_nt", None) or subject_length
        coverage = min(100.0, 100.0 * h["align_length"] / ref_length) if ref_length else None
        c1 = score_identity(h["identity"], seq_type, coverage)
        c2 = 3 if gene_name.lower().startswith(SYSTEM_PREFIXES) else 0
        c3 = 1 if gene and has_copper(gene.compounds) else 0
        c4 = 1 if pmc is not None and pmc >= 5 else 0
        hits.append({
            **h,
            "gene_name": gene_name,
            "organism": gene.organism if gene else "",
            "coverage": coverage,
            "description": (getattr(gene, "protein_name", "") or getattr(gene, "ncbi_description", ""))[:255],
            "ncbi_gene_url": getattr(gene, "ncbi_gene_url", ""),
            "c1": c1, "c2": c2, "c3": c3, "c4": c4,
            "total": c1 + c2 + c3 + c4,
            "pmc_hits": pmc,
        })
    hits.sort(key=lambda x: (-x["total"], -x["identity"], -x["bitscore"]))
    for rank, h in enumerate(hits, 1):
        h["rank"] = rank
    best = hits[0]["total"] if hits else 0
    return best, verdict_for(best, bool(hits)), hits


def explain(hit, seq_type="protein"):
    pmc = "Europe PMC request failed" if hit.pmc_hits is None else f"{hit.pmc_hits} Europe PMC articles"
    (t3, _), (t2, _), (t1, _) = IDENTITY_STEPS[seq_type]
    unit = "nt" if seq_type == "nucleotide" else "aa"
    c1_text = (f"{hit.identity:.1f}% identity over {hit.align_length} {unit} "
               f"(≥{t3}% = 3, {t2}–{t3 - 0.1}% = 2, {t1}–{t2 - 0.1}% = 1).")
    if seq_type == "nucleotide":
        cov = "unknown" if hit.coverage is None else f"{hit.coverage:.0f}%"
        c1_text += f" Alignment covers {cov} of the reference gene (≥{MIN_NT_COVERAGE}% needed for any points)."
    return [
        ("C1 Identity", hit.c1, 3, c1_text),
        ("C2 Multi-gene resistance system", hit.c2, 3,
         f"Gene '{hit.gene_name}' {'belongs' if hit.c2 else 'does not belong'} to the cop/cus/pco/sil/cue/cut systems."),
        ("C3 Activity against Cu", hit.c3, 1,
         "Copper is listed in the BacMet compound field." if hit.c3 else "Copper is not listed in the BacMet compound field."),
        ("C4 Literature support", hit.c4, 1,
         f'{pmc} for "{primary_gene(hit.gene_name)}" AND bioremediation (≥5 needed).'),
    ]
