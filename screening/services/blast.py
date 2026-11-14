import os
import re
import subprocess

from django.conf import settings

BAC_RE = re.compile(r"BAC\d+")
OUT_FIELDS = "qseqid sseqid pident length evalue bitscore qstart qend qlen sstart send slen stitle"
EVALUE = "1e-5"


class BlastError(Exception):
    pass


# Per sequence type: BLAST program, extra options and hits kept per query. Both search the same copper
# protein database, so identity and coverage always refer to the reference PROTEIN.
SEARCHES = {
    "protein": {"bin": lambda: settings.BLASTP_BIN, "options": [], "max_hits": 5},
    "nucleotide": {
        # DNA is translated in all six reading frames with the bacterial genetic code (table 11). A contig can
        # carry several copper genes, so more hits are kept. Query coordinates are in nucleotides.
        "bin": lambda: settings.BLASTX_BIN, "options": ["-query_gencode", "11"], "max_hits": 10,
    },
}
BUILD_CMD = "python manage.py build_cu_db"


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


def run_blast_many(sequences, seq_type="protein"):
    """BLAST all sequences in one run (blastp or blastx). Returns one hit list per input sequence, in order."""
    s = SEARCHES[seq_type]
    if not db_exists():
        raise BlastError(f"Cu protein database not found. Run: {BUILD_CMD}")
    max_hits = s["max_hits"]
    query = "".join(f">q{i}\n{seq}\n" for i, seq in enumerate(sequences))
    out = _run([
        s["bin"](), "-query", "-", "-db", settings.CU_DB_NAME, *s["options"],
        "-outfmt", f"6 {OUT_FIELDS}", "-evalue", EVALUE,
        "-max_target_seqs", str(max_hits), "-max_hsps", "1",
        "-num_threads", str(max(1, min(8, os.cpu_count() or 1))),
    ], query)
    hits = [{} for _ in sequences]
    for line in out.splitlines():
        cols = line.split("\t")
        if len(cols) < 12 or not cols[0].startswith("q"):
            continue
        qseqid, sseqid, pident, length, evalue, bitscore, *coords = cols[:12]
        title = cols[12] if len(cols) > 12 else ""
        per_query = hits[int(qseqid[1:])]
        if len(per_query) == max_hits:
            continue
        m = BAC_RE.search(sseqid) or BAC_RE.search(title)
        bacmet_id = m.group(0) if m else sseqid
        if bacmet_id in per_query:
            continue
        qstart, qend, qlen, sstart, send, slen = map(int, coords)
        per_query[bacmet_id] = {
            "bacmet_id": bacmet_id,
            "identity": float(pident),
            "align_length": int(length),
            "evalue": float(evalue),
            "bitscore": float(bitscore),
            "query_start": qstart, "query_end": qend, "query_length": qlen,
            "subject_start": sstart, "subject_end": send, "subject_length": slen,
        }
    return [list(h.values()) for h in hits]
