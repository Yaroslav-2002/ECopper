"""Score each query's BLAST hits with the 10-point ranking (services/ranking.py); the top-ranked hit is the result."""
import re
from dataclasses import dataclass

from ..models import BacMetGene
from . import ranking
from .blast import run_blast_many

CU_RE = re.compile(r"(?i:\bcopper\b)|\bCu\b")
NO_HIT = "No hit"  # tier of a query without BLAST hits ("" marks a result from an older scoring scheme)

DISCLAIMER = (
    "Score out of 10: BLAST similarity = identity × reference coverage / 100 (0–5), BacMet evidence (0–3) and "
    "other metals besides copper (0–2). " + ranking.LIMITATIONS
)


def has_copper(text):
    return bool(CU_RE.search(text or ""))


@dataclass
class Ref:
    """What the database knows about one reference gene."""
    gene_name: str
    organism: str = ""
    compounds: str = ""
    dataset: str = "experimental"
    references: str = ""
    description: str = ""
    ncbi_gene_url: str = ""
    bacmet_url: str = ""
    annotation_note: str = ""


def format_evalue(evalue, digits=1):
    """BLAST prints E-values below ~1e-180 as 0.0 (they are too small to print), so show 0 as '< 1e-180'."""
    if evalue is None:
        return ""
    return "< 1e-180" if evalue == 0 else f"{evalue:.{digits}e}"


def screen_many(sequences, seq_type="protein", query_ids=None):
    """Screen sequences with one BLAST run. Returns one result dict per sequence, in order."""
    query_ids = query_ids or [f"Query {i + 1}" for i in range(len(sequences))]
    raw_lists = run_blast_many(sequences, seq_type)
    refs = _references({h["bacmet_id"] for raw in raw_lists for h in raw})
    return [evaluate(qid, raw, refs) for qid, raw in zip(query_ids, raw_lists)]


def _references(ids):
    return {
        bacmet_id: Ref(
            gene.gene_name, gene.organism, gene.compounds, gene.dataset, gene.references,
            gene.protein_name or gene.ncbi_description, gene.ncbi_gene_url, gene.source_link, gene.annotation_note,
        )
        for bacmet_id, gene in BacMetGene.objects.in_bulk(ids, field_name="bacmet_id").items()
    }


def score_hit(h, refs):
    """One raw BLAST hit -> the fields of a Hit, including its score."""
    ref = refs.get(h["bacmet_id"]) or Ref(h["bacmet_id"], dataset="", compounds="")
    ref_cov = min(100.0, 100 * (abs(h["subject_end"] - h["subject_start"]) + 1) / h["subject_length"])
    s = ranking.score(dataset=ref.dataset, compounds=ref.compounds, references=ref.references,
                      identity=h["identity"], coverage=ref_cov, evalue=h["evalue"])
    if ref.annotation_note:
        s["notes"].append(ref.annotation_note)
    return {
        "dataset": ref.dataset,
        "bacmet_id": h["bacmet_id"], "gene_name": ref.gene_name, "organism": ref.organism,
        "compounds": ref.compounds, "description": ref.description[:255], "ncbi_gene_url": ref.ncbi_gene_url,
        "bacmet_url": ref.bacmet_url,
        "identity": h["identity"], "align_length": h["align_length"], "evalue": h["evalue"],
        "bitscore": h["bitscore"], "query_start": h["query_start"], "query_end": h["query_end"],
        "ref_coverage": ref_cov, **s, "notes": "; ".join(s["notes"]),
    }


def evaluate(query_id, raw_hits, refs):
    """Score and rank one query's hits; the query takes the score and tier of its top-ranked hit."""
    # Experimental matches take priority; search both databases independently so predicted hits cannot hide them.
    eligible = [h for h in raw_hits if h["bacmet_id"] in refs
                and refs[h["bacmet_id"]].dataset in ("experimental", "predicted")
                and h["evalue"] <= ranking.SCORING["evalue_cutoff"]]
    hits = sorted((score_hit(h, refs) for h in eligible),
                  key=lambda h: (h["dataset"] != "experimental",
                                 *ranking.ranking_key(h["total_score"], h["similarity"], h["evalue"])))
    for rank, h in enumerate(hits, 1):
        h["rank"] = rank
    if not hits:
        return {"hits": [], "tier": NO_HIT, "total_score": 0, "conclusion": (
            f"{query_id} has no qualifying match in the experimental or predicted BacMet copper references "
            "(E-value cutoff 1e-5). No candidate result; score 0.")}
    top = hits[0]
    return {"hits": hits, "tier": top["tier"], "total_score": top["total_score"], "conclusion": (
        f"{query_id}: top candidate {top['gene_name']} ({top['bacmet_id']}), {top['identity']:.0f}% identity over "
        f"{top['ref_coverage']:.0f}% of the reference, scores {top['total_score']:g}/{ranking.MAX_SCORE} "
        f"({top['tier']} tier).")}
