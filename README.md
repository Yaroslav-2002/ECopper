# CuRemediate

A small web app that screens protein or DNA sequences for **copper bioremediation potential**.
Each sequence is searched with NCBI BLAST+ against experimental and predicted copper-resistance proteins from
[BacMet](http://bacmet.biomedicine.gu.se/), and the hits are ranked with a transparent 10-point score
so that candidates can be prioritised for lab validation.

- One sequence → a detailed report (score breakdown, article links, ranked hits, CSV export).
- Upload or paste exactly one sequence. Multiple FASTA entries are rejected before screening.
- Protein input is searched with `blastp`. Nucleotide input (genes, contigs) is translated in all six reading frames
  with the bacterial genetic code and searched with `blastx`, against the same copper-resistance proteins, so both
  modes are scored on protein identity.

> The score is an expert heuristic for prioritising candidates, not an established standard.
> The best-hit method can miss novel genes with low similarity. Confirm results experimentally.

## Quick start with Docker (recommended)

Requires [Docker Desktop](https://www.docker.com/products/docker-desktop/) (Windows, macOS) or Docker Engine (Linux).

```bash
git clone https://github.com/Yaroslav-2002/ECopper.git
cd ECopper
docker compose up --build
```

Open <http://localhost:8000>. The first build downloads BLAST+ 2.17.0 from NCBI (about 300 MB) and takes a few minutes;
later starts take seconds. On every start the container builds the copper reference table and BLAST protein database
from `data/`, then serves the app.

- Stop: `Ctrl+C`, or `docker compose down`.
- Your analyses are kept in the `curemediate-data` Docker volume. `docker compose down -v` deletes them.
- Use another port: `CUREMEDIATE_PORT=8080 docker compose up` (PowerShell: `$env:CUREMEDIATE_PORT=8080; docker compose up`).
- Run the tests: `docker compose run --rm web python manage.py test screening`.

Images build for both x86-64 and ARM (Apple Silicon) machines.

## Running without Docker

Requires Python 3.11+ and [NCBI BLAST+](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) (2.17.0 is the tested version).

1. Install BLAST+ so that `blastp`, `blastx` and `makeblastdb` are on your PATH. Alternatively:
   - unpack a BLAST+ release into `tools/` (the app finds `tools/ncbi-blast-*/bin` automatically), or
   - set `BLAST_BIN_DIR` to the folder that contains the programs.

   With conda you can install it with `conda install -c bioconda -c conda-forge blast`.

2. Create a virtual environment and install the Python packages:

   ```bash
   python -m venv .venv
   source .venv/bin/activate          # Windows: .venv\Scripts\activate
   pip install -r requirements.txt
   ```

3. Create the database and build the copper reference database:

   ```bash
   python manage.py migrate
   python manage.py build_cu_db       # copper reference table + BLAST protein database
   ```

4. Start the development server and open <http://127.0.0.1:8000>:

   ```bash
   python manage.py runserver
   ```

Tests: `python manage.py test screening` (BLAST integration tests are skipped if the databases are not built).

## Configuration

All settings have working defaults for local use. Environment variables:

| Variable | Default | Purpose |
|---|---|---|
| `DJANGO_DEBUG` | `1` (`0` in Docker) | Django debug mode. Turn off whenever the app is reachable by others. |
| `DJANGO_SECRET_KEY` | development key (Docker: generated once and kept in the volume) | Set your own for any shared deployment. |
| `DJANGO_ALLOWED_HOSTS` | `127.0.0.1,localhost` | Comma-separated host names the app answers to. |
| `DJANGO_CSRF_TRUSTED_ORIGINS` | empty | Needed behind a proxy or on another domain, e.g. `https://cu.example.org`. |
| `CUREMEDIATE_DB_PATH` | `db.sqlite3` in the project folder | SQLite database file. |
| `CUREMEDIATE_DATA_DIR` | `data/` | Folder with the BacMet input files and the generated BLAST databases. |
| `BLAST_BIN_DIR` | auto-detected | Folder containing `blastp`, `blastx` and `makeblastdb`. |

The Docker setup is meant for running on your own machine. To expose it on a network, set your own
`DJANGO_SECRET_KEY` and `DJANGO_ALLOWED_HOSTS`, and put it behind an HTTPS reverse proxy.

## How results are ranked

The experimental and predicted copper databases are searched separately (up to 5 hits per protein or 10 per DNA
sequence from each database; E-value cutoff 1e-5). Each qualifying hit gets a score out of 10:

1. **BLAST similarity (0–5)**: identity × reference coverage / 100. ≥ 80 → 5, 60–79 → 4, 40–59 → 3, 20–39 → 2,
   below 20 → 0. If the E-value is above 1e-5 (or missing), this criterion is capped at 1.
2. **BacMet evidence (0–3)**: experimentally confirmed → 3, predicted → 1. An experimental match takes priority
   over a predicted one. Neither database matches → no candidate result, score 0.
3. **Other metals (0–2)**: copper plus at least one other metal → 2, copper only → 1;
   absent, unrecognized or non-copper annotations → 0.

The supplied predicted FASTA contains gene-family names but no per-protein compounds. Its families are matched
to the EXP mapping (full family name first, then an exact alias). Only metals shared by **all** matching mapping
entries are used. Conflicting copper associations leave compounds empty and unscored. Notes identify the mapping
IDs and mark the annotation as inferred; EXP organisms, descriptions, article references and NCBI links are not
copied to predicted proteins. Predicted IDs starting with `PRE` are stable local identifiers, not official BacMet IDs.

Tiers: **High** ≥ 8, **Medium** 5–7.9, **Low** below 5. Missing data counts as 0 points and is listed in the notes.
Literature references are parsed into PubMed / DOI links and counted (`n_references`), but are not scored.

Hits are ranked with experimental matches first, then total score, higher similarity and lower E-value. A sequence takes the score of its
top-ranked hit, and batch results are ordered the same way.

All weights, bins, tier cut-offs and the metal list are in one dictionary, `SCORING`, in
[`screening/services/ranking.py`](screening/services/ranking.py). The same rules rank a ready-made table of BacMet genes
from the command line (needs pandas):

```bash
python scripts/rank_copper_genes.py genes.tsv -o genes_ranked.csv
```

The input needs the columns `bacmet_id, gene_name, organism, bacmet_dataset, compounds, blast_identity, blast_coverage,
blast_evalue, references`. The script prints the top 10 and writes the full ranked table.

**Limitations.** Sequence similarity to a BacMet gene and copper resistance do not prove bioremediation ability,
and the weights are arbitrary choices for prioritising candidates, not a validated model.

## Project layout

```
config/                  Django settings and URLs
screening/               the app: models, views, forms, templates, static files
  services/blast.py      runs blastp / blastx
  services/ranking.py    the 10-point scoring rules (SCORING dictionary), shared with scripts/
  services/scoring.py    scores and ranks each sequence's BLAST hits
  management/commands/   build_cu_db
data/                    reference data (the BLAST database is generated here)
scripts/                 rank_copper_genes.py: rank a table of BacMet genes with the same rules
examples/                example FASTA files, also offered as one-click examples on the upload page
docker/entrypoint.sh     container start-up: migrate, build databases, start gunicorn
```

## Data sources

- **BacMet** – Antibacterial Biocide and Metal Resistance Genes Database, experimental and predicted sets (BacMet 2).
  Pal C., Bengtsson-Palme J., Rensing C., Kristiansson E., Larsson D.G.J. (2014).
  *BacMet: antibacterial biocide and metal resistance genes database.* Nucleic Acids Research 42(D1): D737–D743.
- **UniProt** and **NCBI Gene** – protein names, descriptions and NCBI Gene links for the copper genes
  (`data/copper_genes_annotation.tsv`), shown with each hit. UniProt data is CC BY 4.0.
- **NCBI BLAST+** – Camacho C. et al. (2009). *BLAST+: architecture and applications.* BMC Bioinformatics 10: 421.

Files in `data/`:

| File | Source | Used for |
|---|---|---|
| `BacMet2_EXP_database.fasta` | BacMet 2 (EXP) | experimental copper proteins for the `bacmet_cu` BLAST database |
| `BacMet2_predicted_database.fasta.gz` | BacMet 2 (predicted), supplied locally | copper-associated gene families for the separate `bacmet_cu_predicted` BLAST database; read directly from gzip |
| `BacMet2_EXP.753.mapping.txt` | BacMet 2 (EXP) | gene, organism and metal list per entry (copper selection, other-metals criterion) |
| `copper_genes_annotation.tsv` | derived from UniProt and NCBI Gene | protein descriptions and NCBI Gene links; `references` (PubMed IDs) and `bacmet_url` (the entry in the BacMet experimentally confirmed database) were taken from the BacMet website, which lacks the 37 entries added in BacMet 2 (those are left empty and link to a BacMet EXP search by gene name) |
