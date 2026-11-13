import csv

from django.conf import settings
from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render

from .forms import UploadForm
from .models import Analysis, Batch, Hit
from .services.blast import BlastError
from .services.scoring import CATEGORIES, CRITERIA, DISCLAIMER, gate_checks, ranking_key, screen_many


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


# One-click example inputs shown on the upload page: (file in examples/, button label).
EXAMPLES = [
    ("copA_ecoli.fasta", "CopA · E. coli"),
    ("cueO_ecoli.fasta", "CueO · E. coli"),
    ("copZ_bsubtilis.fasta", "CopZ · B. subtilis"),
    ("zntA_ecoli.fasta", "ZntA · E. coli"),
    ("hba_human.fasta", "Hemoglobin α · human"),
]


def _examples():
    examples = {}
    for filename, label in EXAMPLES:
        path = settings.BASE_DIR / "examples" / filename
        try:
            text = path.read_text(encoding="utf-8")
        except OSError:
            continue
        examples[filename] = {"label": label, "seq_type": "protein", "text": text}
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
    """Analyses of a batch with their best hit, ordered by the ranking specification (Step 5)."""
    rows = [(a, a.hits.all()[0]) for a in batch.analyses.prefetch_related("hits")]
    rows.sort(key=lambda r: ranking_key(r[0].category, r[0].score_pct, r[0].uncertainty, r[1].identity))
    return rows


def batch_detail(request, pk):
    batch = get_object_or_404(Batch, pk=pk)
    rows = _batch_rows(batch)
    categories = [a.category for a, _ in rows]
    return render(request, "screening/batch.html", {
        "b": batch,
        "rows": rows,
        "counts": [(c, categories.count(c)) for c in CATEGORIES if c in categories],
        "no_hit_count": batch.screened_count - len(rows),
        "disclaimer": DISCLAIMER,
    })


# Report columns of the ranking specification, for the best hit of each sequence.
REPORT_HEADER = [
    "Query ID", "Length", "Best hit gene", "BacMet ID", "Description", "Organism",
    "Identity %", "Reference coverage %", "E-value", "Bitscore",
    *(f"{k} {label} (/{mx})" for k, (label, mx) in CRITERIA.items()),
    "Total", "Max available", "%", "Category", "Uncertainty", "Flags", "Evidence reference", "Conclusion",
]


def _report_row(a, h):
    points = a.points
    return [
        a.query_id, a.length, h.gene_name, h.bacmet_id, h.description, h.organism,
        f"{h.identity:.1f}", _pct(h.ref_coverage), f"{h.evalue:.2e}", f"{h.bitscore:.1f}",
        *(("skipped" if points[k] is None else points[k]) if k in points else "" for k in CRITERIA),
        _blank(a.score), _blank(a.max_score), _pct(a.score_pct), a.category or "old scoring - re-run",
        a.uncertainty, "; ".join(a.flags), "not in database", a.conclusion,
    ]


def _pct(value):
    return "" if value is None else f"{value:.0f}"


def _blank(value):
    return "" if value is None else value


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
    w.writerow(REPORT_HEADER)
    for a, h in rows:
        w.writerow(_report_row(a, h))
    return resp


def result(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    hits = list(analysis.hits.all())
    top = hits[0] if hits else None
    gate = (gate_checks(top.evalue, top.bitscore, top.ref_coverage)
            if top and top.ref_coverage is not None else [])
    return render(request, "screening/result.html", {
        "a": analysis, "hits": hits, "top": top, "gate": gate, "disclaimer": DISCLAIMER,
    })


def result_csv(request, pk):
    a = get_object_or_404(Analysis, pk=pk)
    hits = list(a.hits.all())
    resp, w = _csv_response(f"curemediate_analysis_{a.pk}.csv")
    w.writerow(["File", a.filename])
    w.writerow(["Sequence type", a.get_seq_type_display()])
    w.writerow(["Date", a.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow(["Category", a.category or "old scoring - re-run"])
    w.writerow(["Conclusion", a.conclusion])
    w.writerow(["Disclaimer", DISCLAIMER])
    w.writerow([])
    if hits:
        w.writerow(REPORT_HEADER)
        w.writerow(_report_row(a, hits[0]))
        w.writerow([])
    w.writerow([
        "Rank", "Gene", "BacMet ID", "Description", "Family", "Organism", "Identity %", "Align length",
        "Reference coverage %", "Query start", "Query end", "E-value", "Bitscore", "Passes gate",
    ])
    for h in hits:
        w.writerow([
            h.rank, h.gene_name, h.bacmet_id, h.description, h.family, h.organism, f"{h.identity:.2f}",
            h.align_length, _pct(h.ref_coverage), h.query_start, h.query_end,
            f"{h.evalue:.2e}", f"{h.bitscore:.1f}", "yes" if h.passes_gate else "no",
        ])
    return resp


def history(request):
    return render(request, "screening/history.html", {
        "batches": Batch.objects.all(),
        "analyses": Analysis.objects.filter(batch__isnull=True),
    })
