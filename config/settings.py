import os
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent.parent

SECRET_KEY = os.environ.get("DJANGO_SECRET_KEY", "dev-only-cu-remediate-key")
DEBUG = True
ALLOWED_HOSTS = ["127.0.0.1", "localhost"]

INSTALLED_APPS = [
    "django.contrib.contenttypes",
    "django.contrib.staticfiles",
    "screening",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

ROOT_URLCONF = "config.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [],
        "APP_DIRS": True,
        "OPTIONS": {"context_processors": ["django.template.context_processors.request"]},
    }
]

WSGI_APPLICATION = "config.wsgi.application"

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": BASE_DIR / "db.sqlite3",
    }
}

LANGUAGE_CODE = "en-us"
TIME_ZONE = "UTC"
USE_I18N = False
USE_TZ = True

STATIC_URL = "static/"
DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

DATA_DIR = Path(os.environ.get("CUREMEDIATE_DATA_DIR", BASE_DIR / "data"))
BACMET_FASTA = DATA_DIR / "BacMet2_EXP_database.fasta"
BACMET_MAPPING = DATA_DIR / "BacMet2_EXP.753.mapping.txt"
CU_FASTA_NAME = "bacmet_cu.fasta"
CU_DB_NAME = "bacmet_cu"
# Nucleotide (blastn) database: copper genes' DNA sequences plus their annotation table.
CU_NT_FASTA_NAME = "copper_genes_nucleotide.fasta"
CU_NT_ANNOTATION = DATA_DIR / "copper_genes_annotation_nucleotide.tsv"
CU_NT_DB_NAME = "bacmet_cu_nt"

BLAST_BIN_DIR = os.environ.get("BLAST_BIN_DIR") or next(
    (str(p) for p in sorted(BASE_DIR.glob("tools/ncbi-blast-*/bin"), reverse=True)), ""
)
BLASTP_BIN = os.path.join(BLAST_BIN_DIR, "blastp") if BLAST_BIN_DIR else "blastp"
BLASTN_BIN = os.path.join(BLAST_BIN_DIR, "blastn") if BLAST_BIN_DIR else "blastn"
MAKEBLASTDB_BIN = os.path.join(BLAST_BIN_DIR, "makeblastdb") if BLAST_BIN_DIR else "makeblastdb"

MAX_UPLOAD_BYTES = 5 * 1024 * 1024  # per FASTA file or pasted sequence text
# Django rejects form bodies over 2.5 MB by default; leave room for a full-size pasted sequence.
DATA_UPLOAD_MAX_MEMORY_SIZE = MAX_UPLOAD_BYTES + 1024 * 1024
EUROPEPMC_URL = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
EUROPEPMC_TIMEOUT = 10
BACMET_GENE_URL = "http://bacmet.biomedicine.gu.se/advanced_BacMet_EXP_connect_to_ALL.pl?gene={gene}"
