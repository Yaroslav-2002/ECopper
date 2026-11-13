import re
import subprocess

from django.conf import settings

BAC_RE = re.compile(r"BAC\d+")
OUT_FIELDS = "sseqid pident length evalue bitscore stitle"


class BlastError(Exception):
    pass


def _run(cmd, input_text=None):
    try:
        proc = subprocess.run(
            cmd, cwd=settings.DATA_DIR, input=input_text, capture_output=True, text=True, timeout=300
        )
    except FileNotFoundError:
        raise BlastError(f"Executable not found: {cmd[0]}. Install BLAST+ or set BLAST_BIN_DIR.")
    except subprocess.TimeoutExpired:
        raise BlastError(f"{cmd[0]} timed out.")
    if proc.returncode != 0:
        raise BlastError(proc.stderr.strip() or f"{cmd[0]} exited with code {proc.returncode}.")
    return proc.stdout


def db_exists():
    return any(settings.DATA_DIR.glob(settings.CU_DB_NAME + ".p*"))


def make_db():
    return _run([
        settings.MAKEBLASTDB_BIN, "-in", settings.CU_FASTA_NAME,
        "-dbtype", "prot", "-out", settings.CU_DB_NAME,
    ])


def run_blastp(sequence, max_hits=5):
    if not db_exists():
        raise BlastError("Cu database not found. Run: python manage.py build_cu_db")
    out = _run([
        settings.BLASTP_BIN, "-query", "-", "-db", settings.CU_DB_NAME,
        "-outfmt", f"6 {OUT_FIELDS}", "-evalue", "1e-5",
        "-max_target_seqs", str(max_hits), "-max_hsps", "1",
    ], f">query\n{sequence}\n")
    hits = {}
    for line in out.splitlines():
        cols = line.split("\t")
        if len(cols) < 5:
            continue
        title = cols[5] if len(cols) > 5 else ""
        m = BAC_RE.search(cols[0]) or BAC_RE.search(title)
        bacmet_id = m.group(0) if m else cols[0]
        if bacmet_id in hits:
            continue
        hits[bacmet_id] = {
            "bacmet_id": bacmet_id,
            "identity": float(cols[1]),
            "align_length": int(cols[2]),
            "evalue": float(cols[3]),
            "bitscore": float(cols[4]),
        }
        if len(hits) == max_hits:
            break
    return list(hits.values())
