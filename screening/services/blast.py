import os
import re
import subprocess

from django.conf import settings

BAC_RE = re.compile(r"BAC\d+")
OUT_FIELDS = "qseqid sseqid pident length evalue bitscore qstart qend slen stitle"
EVALUE = "1e-5"


class BlastError(Exception):
    pass


# Per sequence type: BLAST program, database, database file suffix, how to (re)build it, extra options.
SEARCHES = {
    "protein": {
        "bin": lambda: settings.BLASTP_BIN, "db": lambda: settings.CU_DB_NAME, "suffix": ".p*",
        "build_cmd": "python manage.py build_cu_db", "options": [], "max_hits": 5,
    },
    "nucleotide": {
        # Full blastn (not megablast) so diverged homologs from other species are still found.
        "bin": lambda: settings.BLASTN_BIN, "db": lambda: settings.CU_NT_DB_NAME, "suffix": ".n*",
        "build_cmd": "python manage.py build_cu_nt_db", "options": ["-task", "blastn"], "max_hits": 10,
    },
}


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


def db_exists(seq_type="protein"):
    s = SEARCHES[seq_type]
    return any(settings.DATA_DIR.glob(s["db"]() + s["suffix"]))


def make_db():
    return _run([
        settings.MAKEBLASTDB_BIN, "-in", settings.CU_FASTA_NAME,
        "-dbtype", "prot", "-out", settings.CU_DB_NAME,
    ])


def make_nt_db():
    return _run([
        settings.MAKEBLASTDB_BIN, "-in", settings.CU_NT_FASTA_NAME,
        "-dbtype", "nucl", "-out", settings.CU_NT_DB_NAME,
    ])


def run_blast_many(sequences, seq_type="protein"):
    """BLAST all sequences in one run (blastp or blastn). Returns one hit list per input sequence, in order."""
    s = SEARCHES[seq_type]
    if not db_exists(seq_type):
        raise BlastError(f"Cu {seq_type} database not found. Run: {s['build_cmd']}")
    max_hits = s["max_hits"]
    query = "".join(f">q{i}\n{seq}\n" for i, seq in enumerate(sequences))
    out = _run([
        s["bin"](), "-query", "-", "-db", s["db"](), *s["options"],
        "-outfmt", f"6 {OUT_FIELDS}", "-evalue", EVALUE,
        "-max_target_seqs", str(max_hits), "-max_hsps", "1",
        "-num_threads", str(max(1, min(8, os.cpu_count() or 1))),
    ], query)
    hits = [{} for _ in sequences]
    for line in out.splitlines():
        cols = line.split("\t")
        if len(cols) < 9 or not cols[0].startswith("q"):
            continue
        qseqid, sseqid, pident, length, evalue, bitscore, qstart, qend, slen = cols[:9]
        title = cols[9] if len(cols) > 9 else ""
        per_query = hits[int(qseqid[1:])]
        if len(per_query) == max_hits:
            continue
        m = BAC_RE.search(sseqid) or BAC_RE.search(title)
        bacmet_id = m.group(0) if m else sseqid
        if bacmet_id in per_query:
            continue
        per_query[bacmet_id] = {
            "bacmet_id": bacmet_id,
            "identity": float(pident),
            "align_length": int(length),
            "evalue": float(evalue),
            "bitscore": float(bitscore),
            "query_start": int(qstart),
            "query_end": int(qend),
            "subject_length": int(slen),
        }
    return [list(h.values()) for h in hits]
