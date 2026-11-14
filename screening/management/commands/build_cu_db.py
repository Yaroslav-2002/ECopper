import csv
import gzip
import hashlib
import re
from collections import defaultdict
from functools import lru_cache

from Bio import SeqIO
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from screening.models import BacMetGene
from screening.services.blast import BlastError, make_db
from screening.services.scoring import has_copper
from screening.services.ranking import metals_in

COLUMN_PREFIXES = {
    "bacmet_id": "bacmet",
    "gene_name": "gene",
    "accession": "accession",
    "organism": "organism",
    "compounds": "compound",
}
REQUIRED = ("bacmet_id", "gene_name", "compounds")


def norm(name):
    return re.sub(r"[^a-z0-9]+", "_", name.strip().lstrip("#").lower()).strip("_")


def read_mapping(path):
    lines = path.read_text(encoding="utf-8-sig", errors="replace").splitlines()
    for start, line in enumerate(lines):
        header = [norm(c) for c in line.split("\t")]
        if any(h.startswith("bacmet") for h in header) and any(h.startswith("compound") for h in header):
            break
    else:
        raise CommandError(f"No header with BacMet_ID and Compound columns found in {path.name}.")
    idx = {}
    for key, prefix in COLUMN_PREFIXES.items():
        idx[key] = next((i for i, h in enumerate(header) if h.startswith(prefix)), None)
    missing = [k for k in REQUIRED if idx[k] is None]
    if missing:
        raise CommandError(f"Missing columns in mapping file: {', '.join(missing)}")
    rows = []
    for line in lines[start + 1:]:
        cols = [c.strip() for c in line.split("\t")]
        get = lambda k: cols[idx[k]] if idx[k] is not None and idx[k] < len(cols) else ""
        if get("bacmet_id"):
            rows.append({k: get(k) for k in COLUMN_PREFIXES})
    return rows


def read_annotation(path):
    """Annotation per BacMet ID: UniProt / NCBI Gene descriptions, the NCBI Gene link and, from optional columns,
    the gene's literature references (";"-separated PubMed IDs, DOIs or titles) and its BacMet EXP page."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "bacmet_id" not in (reader.fieldnames or []):
            raise CommandError(f"No bacmet_id column in {path.name}.")
        return {
            r["bacmet_id"].strip(): {
                "protein_name": (r.get("uniprot_protein_name") or "").strip()[:255],
                "ncbi_description": (r.get("ncbi_description") or "").strip()[:255],
                "ncbi_gene_url": (r.get("ncbi_gene_url") or "").strip()[:255],
                "references": (r.get("references") or "").strip(),
                "bacmet_url": (r.get("bacmet_url") or "").strip()[:255],
            }
            for r in reader if r["bacmet_id"].strip()
        }


def predicted_genes(path, rows):
    """Link predicted header families to EXP mapping, without calling their evidence experimental.

    The supplied predicted FASTA has no per-protein compound table. Use only metals common to all
    mapping entries for the named family. Never copy an experimental organism, PMID or gene link.
    """
    families, aliases = defaultdict(list), defaultdict(list)
    for row in rows:
        families[row["gene_name"].lower()].append(row)
        for alias in row["gene_name"].lower().split("/"):
            aliases[alias.strip()].append(row)

    @lru_cache(maxsize=None)
    def family_annotation(family):
        matches = families.get(family) or aliases.get(family, [])
        if not any(has_copper(row["compounds"]) for row in matches):
            return None
        metal_sets = [metals_in(row["compounds"]) for row in matches]
        common_metals = set.intersection(*metal_sets)
        compounds = ", ".join(sorted(common_metals)) if "Cu" in common_metals else ""
        note = "Predicted family annotation inferred from EXP mapping: " + ", ".join(row["bacmet_id"] for row in matches)
        if len({frozenset(metals) for metals in metal_sets}) > 1:
            note += "; differing family metal annotations: only shared metals used"
        if not compounds:
            note += "; copper association ambiguous, compounds left unscored"
        return compounds, note

    seen = set()
    with gzip.open(path, "rt", encoding="utf-8-sig") as handle:
        for rec in SeqIO.parse(handle, "fasta"):
            labels = re.findall(r"\[([^\]]+)\]", rec.description)
            if not labels:
                raise CommandError(f"Predicted FASTA header lacks a gene-family label: {rec.description}")
            family = labels[-1].strip()
            annotation = family_annotation(family.lower())
            if annotation is None:
                continue
            bid = "PRE" + hashlib.sha256(rec.description.encode()).hexdigest()[:24]
            if bid in seen:
                continue
            seen.add(bid)
            compounds, note = annotation
            fields = rec.id.split("|")
            accession = fields[3] if len(fields) > 3 else rec.id
            description = rec.description[len(rec.id):].split("[")[0].strip()
            yield BacMetGene(
                bacmet_id=bid, gene_name=family[:128], accession=accession[:64],
                organism=labels[-2][:255] if len(labels) > 1 else "",
                dataset="predicted", sequence=str(rec.seq).upper().rstrip("*"), compounds=compounds,
                protein_name=description[:255], annotation_note=note,
            )


class Command(BaseCommand):
    help = "Build experimental and predicted copper reference tables and separate BLAST protein databases."

    def handle(self, *args, **options):
        for path in (settings.BACMET_FASTA, settings.BACMET_MAPPING, settings.CU_ANNOTATION, settings.BACMET_PREDICTED_FASTA):
            if not path.exists():
                raise CommandError(f"Missing file: {path}")

        rows = read_mapping(settings.BACMET_MAPPING)
        cu_rows = {r["bacmet_id"]: r for r in rows if has_copper(r["compounds"])}
        self.stdout.write(f"Mapping entries: {len(rows)}, copper entries: {len(cu_rows)}")

        sequences = {}
        for rec in SeqIO.parse(str(settings.BACMET_FASTA), "fasta"):
            bacmet_id = rec.id.split("|")[0].strip()
            if bacmet_id in cu_rows and bacmet_id not in sequences:
                sequences[bacmet_id] = str(rec.seq).upper().rstrip("*")

        missing = sorted(set(cu_rows) - set(sequences))
        if missing:
            self.stdout.write(self.style.WARNING(f"No sequence for {len(missing)} entries: {', '.join(missing[:10])}"))
        if not sequences:
            raise CommandError("No copper entries matched FASTA headers.")

        annotation = read_annotation(settings.CU_ANNOTATION)
        unannotated = sorted(set(sequences) - set(annotation))
        if unannotated:
            self.stdout.write(self.style.WARNING(
                f"No annotation for {len(unannotated)} entries: {', '.join(unannotated[:10])}"
            ))
        # Every entry comes from BacMet EXP, the experimentally confirmed dataset.
        genes = [
            BacMetGene(sequence=sequences[bid], dataset="experimental", **cu_rows[bid], **annotation.get(bid, {}))
            for bid in sorted(sequences)
        ]
        predicted = list(predicted_genes(settings.BACMET_PREDICTED_FASTA, rows))
        if not predicted:
            raise CommandError("No predicted copper families matched the mapping file.")
        for entries, filename in ((genes, settings.CU_FASTA_NAME), (predicted, settings.CU_PRED_FASTA_NAME)):
            out_fasta = settings.DATA_DIR / filename
            with open(out_fasta, "w", newline="\n") as fh:
                for g in entries:
                    seq = "\n".join(g.sequence[i:i + 60] for i in range(0, len(g.sequence), 60))
                    fh.write(f">{g.bacmet_id} {g.gene_name}\n{seq}\n")
            self.stdout.write(f"Wrote {len(entries)} sequences to {out_fasta}")

        with transaction.atomic():
            BacMetGene.objects.all().delete()
            BacMetGene.objects.bulk_create(genes + predicted, batch_size=500)

        try:
            make_db()
            make_db(settings.CU_PRED_FASTA_NAME, settings.CU_PRED_DB_NAME)
        except BlastError as exc:
            raise CommandError(f"makeblastdb failed: {exc}")
        self.stdout.write(self.style.SUCCESS(f"BLAST databases ready: {len(genes)} experimental, {len(predicted)} predicted"))
