import csv

from Bio import SeqIO
from django.conf import settings
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from screening.models import NucleotideGene
from screening.services.blast import BAC_RE, BlastError, make_nt_db

REQUIRED_COLUMNS = ("bacmet_id", "bacmet_gene", "bacmet_compounds", "status")


def _int(value):
    try:
        return int(value)
    except (TypeError, ValueError):
        return None


class Command(BaseCommand):
    help = "Load the copper genes' nucleotide annotation and build the blastn database."

    def handle(self, *args, **options):
        fasta = settings.DATA_DIR / settings.CU_NT_FASTA_NAME
        for path in (fasta, settings.CU_NT_ANNOTATION):
            if not path.exists():
                raise CommandError(f"Missing file: {path}")

        with open(settings.CU_NT_ANNOTATION, encoding="utf-8-sig", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            missing = [c for c in REQUIRED_COLUMNS if c not in (reader.fieldnames or [])]
            if missing:
                raise CommandError(f"Missing columns in {settings.CU_NT_ANNOTATION.name}: {', '.join(missing)}")
            rows = {r["bacmet_id"].strip(): r for r in reader if r["bacmet_id"].strip()}

        sequences = {}
        for rec in SeqIO.parse(str(fasta), "fasta"):
            m = BAC_RE.match(rec.id)
            if not m:
                raise CommandError(f"FASTA header without a BacMet ID: {rec.id}")
            if m.group(0) in sequences:
                raise CommandError(f"Duplicate BacMet ID in FASTA: {m.group(0)}")
            sequences[m.group(0)] = str(rec.seq).upper()

        no_annotation = sorted(set(sequences) - set(rows))
        if no_annotation:
            raise CommandError(f"FASTA sequences missing from the annotation: {', '.join(no_annotation)}")
        unresolved = sorted(set(rows) - set(sequences))
        self.stdout.write(
            f"Annotation rows: {len(rows)}, sequences: {len(sequences)}"
            + (f", annotated without sequence (not searchable): {', '.join(unresolved)}" if unresolved else "")
        )

        genes = []
        for bid in sorted(sequences):
            r, seq = rows[bid], sequences[bid]
            if _int(r.get("length_nt")) not in (None, len(seq)):
                self.stdout.write(self.style.WARNING(
                    f"{bid}: annotation length {r['length_nt']} nt, FASTA sequence {len(seq)} nt"
                ))
            genes.append(NucleotideGene(
                bacmet_id=bid,
                gene_name=r["bacmet_gene"].strip(),
                uniprot_accession=r.get("bacmet_accession", "").strip(),
                organism=(r.get("ncbi_organism") or r.get("bacmet_organism") or "").strip()[:255],
                location=r.get("bacmet_location", "").strip()[:255],
                compounds=r["bacmet_compounds"].strip(),
                protein_name=r.get("uniprot_protein_name", "").strip()[:255],
                ncbi_description=r.get("ncbi_description", "").strip()[:255],
                ncbi_gene_url=r.get("ncbi_gene_url", "").strip()[:255],
                genomic_accession=r.get("genomic_accession", "").strip()[:64],
                start=_int(r.get("start")),
                end=_int(r.get("end")),
                strand=r.get("strand", "").strip()[:1],
                length_nt=len(seq),
                sequence=seq,
            ))

        with transaction.atomic():
            NucleotideGene.objects.all().delete()
            NucleotideGene.objects.bulk_create(genes)
        self.stdout.write(f"Loaded {len(genes)} nucleotide genes")

        try:
            make_nt_db()
        except BlastError as exc:
            raise CommandError(f"makeblastdb failed: {exc}")
        self.stdout.write(self.style.SUCCESS(f"blastn database ready: {settings.DATA_DIR / settings.CU_NT_DB_NAME}"))
