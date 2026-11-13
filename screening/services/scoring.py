import re

from ..models import BacMetGene
from . import europepmc
from .blast import run_blastp

CU_RE = re.compile(r"(?i:\bcopper\b)|\bCu\b")
SYSTEM_PREFIXES = ("cop", "cus", "pco", "sil", "cue", "cut")
NO_HIT_VERDICT = "No Cu-resistance homolog found"


def has_copper(text):
    return bool(CU_RE.search(text or ""))


def primary_gene(gene_name):
    return gene_name.split("/")[0].strip()


def score_identity(identity):
    if identity >= 80:
        return 3
    if identity >= 65:
        return 2
    if identity >= 50:
        return 1
    return 0


def verdict_for(score, has_hits=True):
    if not has_hits:
        return NO_HIT_VERDICT
    if score >= 6:
        return "Strong Cu bioremediation marker"
    if score >= 3:
        return "Possible marker"
    return "Weak / unlikely"


def screen(sequence):
    raw = run_blastp(sequence)
    genes = BacMetGene.objects.in_bulk([h["bacmet_id"] for h in raw], field_name="bacmet_id")
    hits = []
    for h in raw:
        gene = genes.get(h["bacmet_id"])
        gene_name = gene.gene_name if gene else h["bacmet_id"]
        pmc = europepmc.hit_count(primary_gene(gene_name))
        c1 = score_identity(h["identity"])
        c2 = 3 if gene_name.lower().startswith(SYSTEM_PREFIXES) else 0
        c3 = 1 if gene and has_copper(gene.compounds) else 0
        c4 = 1 if pmc is not None and pmc >= 5 else 0
        hits.append({
            **h,
            "gene_name": gene_name,
            "organism": gene.organism if gene else "",
            "c1": c1, "c2": c2, "c3": c3, "c4": c4,
            "total": c1 + c2 + c3 + c4,
            "pmc_hits": pmc,
        })
    hits.sort(key=lambda x: (-x["total"], -x["identity"], -x["bitscore"]))
    for rank, h in enumerate(hits, 1):
        h["rank"] = rank
    best = hits[0]["total"] if hits else 0
    return best, verdict_for(best, bool(hits)), hits


def explain(hit):
    pmc = "Europe PMC request failed" if hit.pmc_hits is None else f"{hit.pmc_hits} Europe PMC articles"
    return [
        ("C1 Identity", hit.c1, 3,
         f"{hit.identity:.1f}% identity over {hit.align_length} aa (≥80% = 3, 65–79.9% = 2, 50–64.9% = 1)."),
        ("C2 Multi-gene resistance system", hit.c2, 3,
         f"Gene '{hit.gene_name}' {'belongs' if hit.c2 else 'does not belong'} to the cop/cus/pco/sil/cue/cut systems."),
        ("C3 Activity against Cu", hit.c3, 1,
         "Copper is listed in the BacMet compound field." if hit.c3 else "Copper is not listed in the BacMet compound field."),
        ("C4 Literature support", hit.c4, 1,
         f'{pmc} for "{primary_gene(hit.gene_name)}" AND bioremediation (≥5 needed).'),
    ]
