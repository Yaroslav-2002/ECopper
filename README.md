# CuRemediate

A small web app that screens protein or DNA sequences for **copper bioremediation potential**.
Each sequence is searched with NCBI BLAST+ against experimentally confirmed copper-resistance genes from
[BacMet](http://bacmet.biomedicine.gu.se/), and the best hit is ranked with a transparent, rule-based score
so that candidates can be prioritised for lab validation.

- One sequence → a detailed report (significance gate, score breakdown, ranked hits, CSV export).
- A multi-sequence FASTA (a proteome, a gene set, contigs) → a batch screen with a filterable, ranked table.
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

1. **Significance gate** on the best hit (highest bitscore): E-value ≤ 1e-10, bitscore ≥ 50 and
   ≥ 80% of the reference gene covered by the alignment. Failing it → *No significant evidence*.
2. **Function transfer**: below 40% identity → *Function not transferable*.
3. **Score** from criteria A–F: Identity (3), Copper specificity (2), Mechanism (2), Experimental evidence (2),
   Phylogenetic consistency (1), Taxonomic proximity (1). Criteria without data in the database are skipped and
   lower the maximum instead of counting as zero.
4. **Category** from the percentage of the available maximum: ≥ 75% *High potential*, ≥ 50% *Moderate*, otherwise *Low*,
   with an uncertainty level (Low / Moderate / High) and flags such as *possible different metal specificity*.

Batch results are ordered by category, then score %, then lower uncertainty, then higher identity.
The rules live in [`screening/services/scoring.py`](screening/services/scoring.py).

## Project layout

```
config/                  Django settings and URLs
screening/               the app: models, views, forms, templates, static files
  services/blast.py      runs blastp / blastx
  services/scoring.py    significance gate, criteria and ranking
  management/commands/   build_cu_db
data/                    reference data (the BLAST database is generated here)
examples/                example FASTA files, also offered as one-click examples on the upload page
docker/entrypoint.sh     container start-up: migrate, build databases, start gunicorn
```

## Data sources

- **BacMet** – Antibacterial Biocide and Metal Resistance Genes Database, experimentally confirmed set (BacMet2 EXP).
  Pal C., Bengtsson-Palme J., Rensing C., Kristiansson E., Larsson D.G.J. (2014).
  *BacMet: antibacterial biocide and metal resistance genes database.* Nucleic Acids Research 42(D1): D737–D743.
- **UniProt** and **NCBI Gene** – protein names, descriptions and NCBI Gene links for the copper genes
  (`data/copper_genes_annotation.tsv`), used to classify each protein's mechanism. UniProt data is CC BY 4.0.
- **NCBI BLAST+** – Camacho C. et al. (2009). *BLAST+: architecture and applications.* BMC Bioinformatics 10: 421.

Files in `data/`:

| File | Source | Used for |
|---|---|---|
| `BacMet2_EXP_database.fasta` | BacMet 2 (EXP) | protein sequences; the 133 copper entries form the BLAST database |
| `BacMet2_EXP.753.mapping.txt` | BacMet 2 (EXP) | gene, organism and metal list per entry (copper selection, criterion B) |
| `copper_genes_annotation.tsv` | derived from UniProt and NCBI Gene | protein descriptions (criterion C) and NCBI Gene links |
