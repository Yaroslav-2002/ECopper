import io
import os

from Bio import SeqIO
from django import forms
from django.conf import settings

from .models import SEQ_TYPES

ALLOWED_EXT = {".fasta", ".fa", ".faa", ".fna", ".ffn", ".fas"}
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
NUCLEOTIDE_LETTERS = set("ACGTUN")
IUPAC_NUCLEOTIDE_LETTERS = set("ACGTUNRYKMSWBDHV")
PROTEIN_LENGTH = (30, 5000)  # aa
MIN_NUCLEOTIDE_LENGTH = 30  # nt; the upper bound is the file size limit
MAX_UPLOAD_MB = f"{settings.MAX_UPLOAD_BYTES / (1024 * 1024):g} MB"


class UploadForm(forms.Form):
    # Declared before `file` so clean_file() can read it.
    seq_type = forms.ChoiceField(
        label="Sequence type", choices=SEQ_TYPES, initial="protein", widget=forms.RadioSelect
    )
    file = forms.FileField(
        label=f"FASTA file (.fasta, .fa, .faa, .fna; max {MAX_UPLOAD_MB})", required=False
    )
    sequence_text = forms.CharField(
        label="…or paste sequence(s)",
        required=False,
        strip=True,
        widget=forms.Textarea(attrs={
            "rows": 8, "spellcheck": "false",
            "placeholder": ">my_sequence\nMKVLAAGIVG...\n\nFASTA with one or more '>' headers, or just the bare sequence.",
        }),
    )

    def clean_file(self):
        f = self.cleaned_data.get("file")
        if not f:
            return None
        if os.path.splitext(f.name)[1].lower() not in ALLOWED_EXT:
            raise forms.ValidationError(f"Only {', '.join(sorted(ALLOWED_EXT))} files are accepted.")
        if f.size > settings.MAX_UPLOAD_BYTES:
            raise forms.ValidationError(f"File is larger than {MAX_UPLOAD_MB}.")
        try:
            text = f.read().decode("utf-8-sig")
        except UnicodeDecodeError:
            raise forms.ValidationError("File is not a readable text file.")
        if not text.lstrip().startswith(">"):
            raise forms.ValidationError("Not a FASTA file: the first line must start with '>'.")
        return {"filename": f.name[:255], "text": text}

    def clean_sequence_text(self):
        text = self.cleaned_data.get("sequence_text", "")
        if not text:
            return None
        if len(text.encode("utf-8")) > settings.MAX_UPLOAD_BYTES:
            raise forms.ValidationError(f"Pasted text is larger than {MAX_UPLOAD_MB}.")
        if not text.startswith(">"):
            # A bare sequence: give it a header so it parses as FASTA.
            text = ">pasted_sequence\n" + text
        return {"filename": "Pasted sequence", "text": text}

    def clean(self):
        cleaned = super().clean()
        if self.errors:
            return cleaned
        source = cleaned.get("file") or cleaned.get("sequence_text")
        if cleaned.get("file") and cleaned.get("sequence_text"):
            raise forms.ValidationError("Upload a file or paste a sequence, not both.")
        if not source:
            raise forms.ValidationError("Upload a FASTA file or paste a sequence.")
        cleaned["input"] = parse_input(source["text"], source["filename"], cleaned.get("seq_type", "protein"))
        return cleaned


def parse_input(text, filename, seq_type):
    """Parse FASTA text into validated records (single sequence) or records + skipped (batch)."""
    try:
        records = list(SeqIO.parse(io.StringIO(text), "fasta"))
    except Exception:
        raise forms.ValidationError("The FASTA input could not be parsed.")
    if not records:
        raise forms.ValidationError("The input contains no sequences.")

    if len(records) == 1:
        # Single sequence: any problem rejects the input.
        record = _check_record(records[0], 1, seq_type)
        return {"filename": filename, "seq_type": seq_type, "batch": False, "records": [record], "skipped": []}

    # Several sequences: screen the valid ones and report the rest as skipped.
    valid, skipped = [], []
    for i, rec in enumerate(records, 1):
        try:
            valid.append(_check_record(rec, i, seq_type))
        except forms.ValidationError as exc:
            skipped.append([_query_id(rec, i), exc.messages[0]])
    if not valid:
        raise forms.ValidationError(
            f"None of the {len(records)} sequences could be screened. First problem: "
            f"{skipped[0][0]}: {skipped[0][1]}"
        )
    return {"filename": filename, "seq_type": seq_type, "batch": True, "records": valid, "skipped": skipped}


def _query_id(rec, index):
    return (rec.id or f"sequence_{index}")[:255]


def _check_record(rec, index, seq_type):
    seq = str(rec.seq).upper().replace(" ", "")
    if seq_type == "nucleotide":
        seq = _check_nucleotide(seq.replace("-", ""))
    else:
        seq = _check_protein(seq.rstrip("*"))
    return {"query_id": _query_id(rec, index), "sequence": seq}


def _check_protein(seq):
    if not seq:
        raise forms.ValidationError("The sequence is empty.")
    bad = sorted(set(seq) - PROTEIN_LETTERS)
    if bad:
        raise forms.ValidationError(f"Invalid amino-acid characters: {', '.join(bad)}.")
    if set(seq) <= NUCLEOTIDE_LETTERS:
        raise forms.ValidationError(
            "This looks like a nucleotide sequence; choose 'Nucleotide (blastx)' as the sequence type."
        )
    low, high = PROTEIN_LENGTH
    if not low <= len(seq) <= high:
        raise forms.ValidationError(f"Sequence length must be {low}–{high} aa (got {len(seq)}).")
    return seq


def _check_nucleotide(seq):
    if not seq:
        raise forms.ValidationError("The sequence is empty.")
    bad = sorted(set(seq) - IUPAC_NUCLEOTIDE_LETTERS)
    if bad:
        if set(seq) <= PROTEIN_LETTERS:
            raise forms.ValidationError(
                "This looks like a protein sequence; choose 'Protein (blastp)' as the sequence type."
            )
        raise forms.ValidationError(f"Invalid nucleotide characters: {', '.join(bad)}.")
    if len(seq) < MIN_NUCLEOTIDE_LENGTH:
        raise forms.ValidationError(f"Sequence must be at least {MIN_NUCLEOTIDE_LENGTH} nt (got {len(seq)}).")
    if sum(seq.count(c) for c in "ACGTU") < 0.5 * len(seq):
        raise forms.ValidationError("Sequence is mostly ambiguous bases (N or IUPAC codes).")
    return seq.replace("U", "T")
