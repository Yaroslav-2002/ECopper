import csv
import re

from Bio import SeqIO
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from screening.models import BacMetGene
from screening.services.blast import BlastError, make_db
from screening.services.scoring import has_copper

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
    """UniProt / NCBI Gene annotation per BacMet ID: protein descriptions and the NCBI Gene link."""
    with open(path, encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if "bacmet_id" not in (reader.fieldnames or []):
            raise CommandError(f"No bacmet_id column in {path.name}.")
        return {
            r["bacmet_id"].strip(): {
                "protein_name": (r.get("uniprot_protein_name") or "").strip()[:255],
                "ncbi_description": (r.get("ncbi_description") or "").strip()[:255],
                "ncbi_gene_url": (r.get("ncbi_gene_url") or "").strip()[:255],
            }
            for r in reader if r["bacmet_id"].strip()
        }


class Command(BaseCommand):
    help = "Build the copper subset of BacMet EXP with its annotation, and its BLAST protein database."

    def handle(self, *args, **options):
        for path in (settings.BACMET_FASTA, settings.BACMET_MAPPING, settings.CU_ANNOTATION):
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
                f"No annotation for {len(unannotated)} entries (mechanism criterion skipped): {', '.join(unannotated[:10])}"
            ))
        genes = [
            BacMetGene(sequence=sequences[bid], **cu_rows[bid], **annotation.get(bid, {}))
            for bid in sorted(sequences)
        ]
        out_fasta = settings.DATA_DIR / settings.CU_FASTA_NAME
        with open(out_fasta, "w", newline="\n") as fh:
            for g in genes:
                seq = "\n".join(g.sequence[i:i + 60] for i in range(0, len(g.sequence), 60))
                fh.write(f">{g.bacmet_id} {g.gene_name}\n{seq}\n")
        self.stdout.write(f"Wrote {len(genes)} sequences to {out_fasta}")

        with transaction.atomic():
            BacMetGene.objects.all().delete()
            BacMetGene.objects.bulk_create(genes)

        try:
            make_db()
        except BlastError as exc:
            raise CommandError(f"makeblastdb failed: {exc}")
        self.stdout.write(self.style.SUCCESS(f"BLAST database ready: {settings.DATA_DIR / settings.CU_DB_NAME}"))
