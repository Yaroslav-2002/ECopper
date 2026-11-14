import csv

from django.conf import settings
from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render

from .forms import UploadForm
from .models import Analysis, Batch, Hit
from .services.blast import BlastError
from .services.ranking import MAX_POINTS, MAX_SCORE, TIERS, ranking_key
from .services.scoring import DISCLAIMER, screen_many

OLD_SCHEME = "Old scoring – re-run"


def upload(request):
    form = UploadForm(request.POST, request.FILES) if request.method == "POST" else UploadForm()
    if form.is_bound and form.is_valid():
        data = form.cleaned_data["input"]
        records = data["records"]
        try:
            results = screen_many([r["sequence"] for r in records], data["seq_type"], [r["query_id"] for r in records])
        except BlastError as exc:
            form.add_error(None, f"BLAST error: {exc}")
        else:
            with transaction.atomic():
                if not data["batch"]:
                    analysis = _save_analysis(None, data, records[0], results[0])
                    return redirect("result", pk=analysis.pk)
                batch = Batch.objects.create(
                    filename=data["filename"],
                    seq_type=data["seq_type"],
                    sequence_count=len(records) + len(data["skipped"]),
                    screened_count=len(records),
                    skipped=data["skipped"],
                )
                # Only sequences with at least one copper hit are stored.
                for record, res in zip(records, results):
                    if res["hits"]:
                        _save_analysis(batch, data, record, res)
            return redirect("batch", pk=batch.pk)
    return render(request, "screening/upload.html", {"form": form, "examples": _examples()})


# One-click example inputs shown on the upload page: (file in examples/, button label, sequence type).
EXAMPLES = [
    ("copA_ecoli.fasta", "CopA · E. coli", "protein"),
    ("cueO_ecoli.fasta", "CueO · E. coli", "protein"),
    ("copZ_bsubtilis.fasta", "CopZ · B. subtilis", "protein"),
    ("zntA_ecoli.fasta", "ZntA · E. coli", "protein"),
    ("hba_human.fasta", "Hemoglobin α · human", "protein"),
    ("copA_ecoli_gene.fna", "copA gene (DNA) · E. coli", "nucleotide"),
]


def _examples():
    examples = {}
    for filename, label, seq_type in EXAMPLES:
        path = settings.BASE_DIR / "examples" / filename
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        examples[filename] = {"label": label, "seq_type": seq_type, "text": text}
    return examples


def _save_analysis(batch, data, record, res):
    res = dict(res)
    hits = res.pop("hits")
    analysis = Analysis.objects.create(
        batch=batch, filename=data["filename"], seq_type=data["seq_type"], query_id=record["query_id"],
        sequence=record["sequence"], length=len(record["sequence"]), **res,
    )
    Hit.objects.bulk_create([Hit(analysis=analysis, **h) for h in hits])
    return analysis


def _batch_rows(batch):
    """Analyses of a batch with their top-ranked hit, best first (same ranking as the hits of one sequence)."""
    rows = [(a, a.hits.all()[0]) for a in batch.analyses.prefetch_related("hits")]
    rows.sort(key=lambda r: ranking_key(r[1].total_score, r[1].similarity, r[1].evalue))
    return rows


def batch_detail(request, pk):
    batch = get_object_or_404(Batch, pk=pk)
    rows = _batch_rows(batch)
    tiers = [a.tier or OLD_SCHEME for a, _ in rows]
    return render(request, "screening/batch.html", {
        "b": batch,
        "rows": rows,
        "counts": [(t, tiers.count(t)) for t in (*TIERS, OLD_SCHEME) if t in tiers],
        "no_hit_count": batch.screened_count - len(rows),
        "max_points": MAX_POINTS,
        "max_score": MAX_SCORE,
        "disclaimer": DISCLAIMER,
    })


# One CSV row per hit: the input columns of the ranking specification, then its score columns.
HIT_HEADER = [
    "Rank", "BacMet ID", "Gene", "Description", "Organism", "Compounds",
    "Identity %", "Reference coverage %", "E-value", "Bitscore", "Similarity (identity x coverage / 100)",
    f"score_blast (/{MAX_POINTS['blast']})", f"score_bacmet (/{MAX_POINTS['bacmet']})",
    f"score_metals (/{MAX_POINTS['metals']})", f"total_score (/{MAX_SCORE})", "Tier",
    "n_references", "article_links", "Notes",
]


def _hit_row(h):
    return [
        h.rank, h.bacmet_id, h.gene_name, h.description, h.organism, h.compounds,
        f"{h.identity:.1f}", _num(h.ref_coverage, ".0f"), f"{h.evalue:.2e}", f"{h.bitscore:.1f}",
        _num(h.similarity, ".1f"), _num(h.score_blast), _num(h.score_bacmet), _num(h.score_metals),
        _num(h.total_score), h.tier or OLD_SCHEME, _num(h.n_references), " ; ".join(h.article_links), h.notes,
    ]


def _num(value, fmt="g"):
    return "" if value is None else format(value, fmt)


def _csv_response(name):
    resp = HttpResponse(content_type="text/csv")
    resp["Content-Disposition"] = f'attachment; filename="{name}"'
    return resp, csv.writer(resp)


def batch_csv(request, pk):
    b = get_object_or_404(Batch, pk=pk)
    rows = _batch_rows(b)
    resp, w = _csv_response(f"curemediate_batch_{b.pk}.csv")
    w.writerow(["File", b.filename])
    w.writerow(["Sequence type", b.get_seq_type_display()])
    w.writerow(["Date", b.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow(["Sequences in file", b.sequence_count])
    w.writerow(["Screened", b.screened_count])
    w.writerow(["Skipped", len(b.skipped)])
    w.writerow(["With Cu hits", len(rows)])
    w.writerow(["Disclaimer", DISCLAIMER])
    w.writerow([])
    w.writerow(["Rank", "Query ID", "Length", *HIT_HEADER[1:]])  # hit columns describe the top-ranked hit
    for rank, (a, h) in enumerate(rows, 1):
        w.writerow([rank, a.query_id, a.length, *_hit_row(h)[1:]])
    return resp


def result(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    hits = list(analysis.hits.all())
    return render(request, "screening/result.html", {
        "a": analysis, "hits": hits, "top": hits[0] if hits else None,
        "max_points": MAX_POINTS, "max_score": MAX_SCORE, "disclaimer": DISCLAIMER,
    })


def result_csv(request, pk):
    a = get_object_or_404(Analysis, pk=pk)
    resp, w = _csv_response(f"curemediate_analysis_{a.pk}.csv")
    w.writerow(["File", a.filename])
    w.writerow(["Query ID", a.query_id])
    w.writerow(["Sequence type", a.get_seq_type_display()])
    w.writerow(["Date", a.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow(["Tier", a.tier or OLD_SCHEME])
    w.writerow(["Conclusion", a.conclusion])
    w.writerow(["Disclaimer", DISCLAIMER])
    w.writerow([])
    w.writerow(HIT_HEADER)
    for h in a.hits.all():
        w.writerow(_hit_row(h))
    return resp


def history(request):
    return render(request, "screening/history.html", {
        "batches": Batch.objects.all(),
        "analyses": Analysis.objects.filter(batch__isnull=True),
        "max_score": MAX_SCORE,
    })
