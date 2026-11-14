#!/usr/bin/env python3
"""Rank candidate copper-resistance genes from BacMet with a simple 0-10 score.

Usage:
    python scripts/rank_copper_genes.py genes.tsv            # writes genes_ranked.csv
    python scripts/rank_copper_genes.py genes.csv -o out.csv

Input columns: bacmet_id, gene_name, organism, bacmet_dataset, compounds,
blast_identity, blast_coverage, blast_evalue, references (";"-separated).

The scoring rules (weights, bins, tiers, metal list) live in ONE dictionary, SCORING in
screening/services/ranking.py, which the CuRemediate web app uses too.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from screening.services import ranking  # noqa: E402  (plain Python, no Django needed)

REQUIRED = ["bacmet_id", "gene_name", "organism", "bacmet_dataset", "compounds",
            "blast_identity", "blast_coverage", "blast_evalue", "references"]


def to_number(value):
    """'85%', '85', 85 -> 85.0; empty or unparsable -> None."""
    number = pd.to_numeric(str(value).replace("%", "").strip(), errors="coerce")
    return None if pd.isna(number) else float(number)


def score_row(row):
    evalue = to_number(row["blast_evalue"])
    r = ranking.score(
        dataset=row["bacmet_dataset"], compounds=row["compounds"], references=row["references"],
        identity=to_number(row["blast_identity"]), coverage=to_number(row["blast_coverage"]), evalue=evalue,
    )
    return pd.Series({
        "score_blast": r["score_blast"], "score_bacmet": r["score_bacmet"], "score_metals": r["score_metals"],
        "total_score": r["total_score"], "tier": r["tier"], "n_references": r["n_references"],
        "article_links": " ; ".join(r["article_links"]), "notes": "; ".join(r["notes"]),
        "_key": ranking.ranking_key(r["total_score"], r["similarity"], evalue),
    })


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("table", type=Path, help="input CSV or TSV, one row per gene")
    parser.add_argument("-o", "--output", type=Path, help="output CSV (default: <input>_ranked.csv)")
    args = parser.parse_args()

    with open(args.table, encoding="utf-8-sig") as fh:
        sep = "\t" if "\t" in fh.readline() else ","
    genes = pd.read_csv(args.table, sep=sep, dtype=str, encoding="utf-8-sig")
    missing = [c for c in REQUIRED if c not in genes.columns]
    if missing:
        raise SystemExit(f"Missing column(s) in {args.table}: {', '.join(missing)}")

    ranked = genes.join(genes.apply(score_row, axis=1)).sort_values("_key", kind="stable")
    ranked.insert(len(genes.columns) + 4, "rank", range(1, len(ranked) + 1))
    ranked = ranked.drop(columns="_key")

    out = args.output or args.table.with_name(f"{args.table.stem}_ranked.csv")
    ranked.to_csv(out, index=False, encoding="utf-8-sig")  # utf-8-sig opens cleanly in Excel

    top = ranked.head(10)[["rank", "bacmet_id", "gene_name", "organism", "total_score", "tier", "n_references"]]
    top = top.assign(organism=top["organism"].fillna("").str.slice(0, 35))
    print(f"Top {len(top)} of {len(ranked)} genes (full table: {out})\n")
    print(top.to_string(index=False))


if __name__ == "__main__":
    main()

# LIMITATIONS
# - High sequence similarity to a BacMet gene suggests a similar function but does not prove it;
#   only lab assays can confirm copper resistance in the candidate organism.
# - Copper resistance (efflux, sequestration) is not the same as bioremediation ability, which also
#   depends on uptake/accumulation, growth in real effluents, and gene expression.
# - The weights, bins and tier cut-offs are arbitrary expert choices, not a validated model;
#   the ranking should be used to prioritise candidates, not to draw conclusions.
