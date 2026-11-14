from urllib.parse import quote

from django.conf import settings
from django.db import models


class BacMetGene(models.Model):
    bacmet_id = models.CharField(max_length=32, unique=True)
    gene_name = models.CharField(max_length=128)
    accession = models.CharField(max_length=64, blank=True)
    organism = models.CharField(max_length=255, blank=True)
    compounds = models.TextField(blank=True)
    dataset = models.CharField(max_length=16, default="experimental")  # BacMet "experimental" or "predicted"
    annotation_note = models.TextField(blank=True)
    sequence = models.TextField()
    # From data/copper_genes_annotation.tsv: protein description, literature references and the NCBI Gene link.
    references = models.TextField(blank=True)  # ";"-separated PubMed IDs, DOIs or titles
    bacmet_url = models.CharField(max_length=255, blank=True)  # entry in the BacMet EXP database
    protein_name = models.CharField(max_length=255, blank=True)  # UniProt
    ncbi_description = models.CharField(max_length=255, blank=True)  # NCBI Gene
    ncbi_gene_url = models.CharField(max_length=255, blank=True)

    def __str__(self):
        return f"{self.bacmet_id} {self.gene_name}"

    @property
    def source_link(self):
        if self.dataset == "predicted":
            return f"https://www.ncbi.nlm.nih.gov/protein/{quote(self.accession, safe='')}" if self.accession else ""
        return self.bacmet_url


SEQ_TYPES = [("protein", "Protein (blastp)"), ("nucleotide", "Nucleotide (blastx)")]


class Batch(models.Model):
    created_at = models.DateTimeField(auto_now_add=True)
    filename = models.CharField(max_length=255)
    seq_type = models.CharField(max_length=16, choices=SEQ_TYPES, default="protein")
    sequence_count = models.PositiveIntegerField()
    screened_count = models.PositiveIntegerField()
    skipped = models.JSONField(default=list, blank=True)  # [[query_id, reason], ...]

    class Meta:
        ordering = ["-created_at"]

    @property
    def unit(self):
        return "nt" if self.seq_type == "nucleotide" else "aa"


class Analysis(models.Model):
    batch = models.ForeignKey(Batch, null=True, blank=True, on_delete=models.CASCADE, related_name="analyses")
    created_at = models.DateTimeField(auto_now_add=True)
    filename = models.CharField(max_length=255)
    seq_type = models.CharField(max_length=16, choices=SEQ_TYPES, default="protein")
    query_id = models.CharField(max_length=255)
    sequence = models.TextField()
    length = models.PositiveIntegerField()
    # Score and tier of the top-ranked hit (services/ranking.py); tier "" marks a result from an older scheme.
    tier = models.CharField(max_length=16, blank=True)
    total_score = models.FloatField(null=True, blank=True)
    conclusion = models.TextField(blank=True)

    class Meta:
        ordering = ["-created_at"]

    @property
    def unit(self):
        return "nt" if self.seq_type == "nucleotide" else "aa"


class Hit(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE, related_name="hits")
    rank = models.PositiveSmallIntegerField()
    dataset = models.CharField(max_length=16, blank=True)
    bacmet_id = models.CharField(max_length=32)
    gene_name = models.CharField(max_length=128)
    organism = models.CharField(max_length=255, blank=True)
    compounds = models.TextField(blank=True)  # BacMet metal/compound list of the reference gene
    identity = models.FloatField()
    align_length = models.PositiveIntegerField()
    evalue = models.FloatField()
    bitscore = models.FloatField()
    query_start = models.PositiveIntegerField(null=True, blank=True)
    query_end = models.PositiveIntegerField(null=True, blank=True)
    ref_coverage = models.FloatField(null=True, blank=True)  # % of the reference gene covered by the alignment
    description = models.CharField(max_length=255, blank=True)
    ncbi_gene_url = models.CharField(max_length=255, blank=True)
    bacmet_url = models.CharField(max_length=255, blank=True)  # entry in the BacMet EXP database
    # 10-point ranking (services/ranking.py); total_score None marks a hit from an older scoring scheme.
    similarity = models.FloatField(null=True, blank=True)  # identity x ref_coverage / 100
    score_blast = models.FloatField(null=True, blank=True)
    score_bacmet = models.FloatField(null=True, blank=True)
    score_metals = models.FloatField(null=True, blank=True)
    total_score = models.FloatField(null=True, blank=True)
    tier = models.CharField(max_length=16, blank=True)
    n_references = models.PositiveSmallIntegerField(null=True, blank=True)
    article_links = models.JSONField(default=list, blank=True)
    notes = models.TextField(blank=True)  # missing data and caps, "; "-separated

    class Meta:
        ordering = ["rank"]

    @property
    def bacmet_link(self):
        """BacMet EXP entry or the source protein record for a predicted reference."""
        if self.dataset == "predicted":
            if self.bacmet_url:
                return self.bacmet_url
            # Older saved predicted hits did not store a link; resolve their reference accession.
            reference = BacMetGene.objects.filter(bacmet_id=self.bacmet_id, dataset="predicted").first()
            return reference.source_link if reference else ""
        return self.bacmet_url or settings.BACMET_EXP_URL.format(query=quote(self.gene_name.split("/")[0].strip()))
