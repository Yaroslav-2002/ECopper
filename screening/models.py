from urllib.parse import quote, urlencode

from django.conf import settings
from django.db import models


class BacMetGene(models.Model):
    bacmet_id = models.CharField(max_length=32, unique=True)
    gene_name = models.CharField(max_length=128)
    accession = models.CharField(max_length=64, blank=True)
    organism = models.CharField(max_length=255, blank=True)
    compounds = models.TextField(blank=True)
    sequence = models.TextField()

    def __str__(self):
        return f"{self.bacmet_id} {self.gene_name}"


class LiteratureCache(models.Model):
    gene_name = models.CharField(max_length=128, unique=True)
    hit_count = models.IntegerField()
    fetched_at = models.DateTimeField(auto_now=True)


class Analysis(models.Model):
    created_at = models.DateTimeField(auto_now_add=True)
    filename = models.CharField(max_length=255)
    query_id = models.CharField(max_length=255)
    sequence = models.TextField()
    length = models.PositiveIntegerField()
    best_score = models.PositiveSmallIntegerField(default=0)
    verdict = models.CharField(max_length=64)

    class Meta:
        ordering = ["-created_at"]

    @property
    def level(self):
        if self.best_score >= 6:
            return "strong"
        if self.best_score >= 3:
            return "possible"
        return "weak"


class Hit(models.Model):
    analysis = models.ForeignKey(Analysis, on_delete=models.CASCADE, related_name="hits")
    rank = models.PositiveSmallIntegerField()
    bacmet_id = models.CharField(max_length=32)
    gene_name = models.CharField(max_length=128)
    organism = models.CharField(max_length=255, blank=True)
    identity = models.FloatField()
    align_length = models.PositiveIntegerField()
    evalue = models.FloatField()
    bitscore = models.FloatField()
    c1 = models.PositiveSmallIntegerField()
    c2 = models.PositiveSmallIntegerField()
    c3 = models.PositiveSmallIntegerField()
    c4 = models.PositiveSmallIntegerField()
    total = models.PositiveSmallIntegerField()
    pmc_hits = models.IntegerField(null=True, blank=True)

    class Meta:
        ordering = ["rank"]

    @property
    def query_gene(self):
        return self.gene_name.split("/")[0].strip()

    @property
    def bacmet_url(self):
        return settings.BACMET_GENE_URL.format(gene=quote(self.query_gene))

    @property
    def pmc_url(self):
        return "https://europepmc.org/search?" + urlencode({"query": f'"{self.query_gene}" AND bioremediation'})
