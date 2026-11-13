import io
import os

from Bio import SeqIO
from django import forms
from django.conf import settings

ALLOWED_EXT = {".fasta", ".fa", ".faa"}
PROTEIN_LETTERS = set("ACDEFGHIKLMNPQRSTVWYBXZJUO")
NUCLEOTIDE_LETTERS = set("ACGTUN")


class UploadForm(forms.Form):
    file = forms.FileField(label="Protein FASTA (.fasta, .fa, .faa; max 1 MB)")

    def clean_file(self):
        f = self.cleaned_data["file"]
        if os.path.splitext(f.name)[1].lower() not in ALLOWED_EXT:
            raise forms.ValidationError("Only .fasta, .fa or .faa files are accepted.")
        if f.size > settings.MAX_UPLOAD_BYTES:
            raise forms.ValidationError("File is larger than 1 MB.")
        try:
            text = f.read().decode("utf-8-sig")
        except UnicodeDecodeError:
            raise forms.ValidationError("File is not a readable text file.")
        if not text.lstrip().startswith(">"):
            raise forms.ValidationError("Not a FASTA file: the first line must start with '>'.")
        try:
            records = list(SeqIO.parse(io.StringIO(text), "fasta"))
        except Exception:
            raise forms.ValidationError("The FASTA file could not be parsed.")
        if len(records) != 1:
            raise forms.ValidationError(f"The file must contain exactly one sequence (found {len(records)}).")
        rec = records[0]
        seq = str(rec.seq).upper().replace(" ", "").rstrip("*")
        if not seq:
            raise forms.ValidationError("The sequence is empty.")
        bad = sorted(set(seq) - PROTEIN_LETTERS)
        if bad:
            raise forms.ValidationError(f"Invalid amino-acid characters: {', '.join(bad)}.")
        if set(seq) <= NUCLEOTIDE_LETTERS:
            raise forms.ValidationError("This looks like a nucleotide sequence; upload a protein sequence.")
        if not 30 <= len(seq) <= 5000:
            raise forms.ValidationError(f"Sequence length must be 30–5000 aa (got {len(seq)}).")
        return {"filename": f.name, "query_id": rec.id or "query", "sequence": seq}
