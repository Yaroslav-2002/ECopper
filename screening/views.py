import csv

from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render

from .forms import UploadForm
from .models import Analysis, Batch, Hit
from .services.blast import BlastError
from .services.scoring import explain, screen_many


def upload(request):
    form = UploadForm(request.POST, request.FILES) if request.method == "POST" else UploadForm()
    if form.is_bound and form.is_valid():
        data = form.cleaned_data["input"]
        records = data["records"]
        try:
            results = screen_many([r["sequence"] for r in records], data["seq_type"])
        except BlastError as exc:
            form.add_error(None, f"BLAST error: {exc}")
        else:
            with transaction.atomic():
                if not data["batch"]:
                    analysis = _save_analysis(None, data, records[0], *results[0])
                    return redirect("result", pk=analysis.pk)
                batch = Batch.objects.create(
                    filename=data["filename"],
                    seq_type=data["seq_type"],
                    sequence_count=len(records) + len(data["skipped"]),
                    screened_count=len(records),
                    skipped=data["skipped"],
                )
                # Only sequences with at least one copper hit are stored.
                for record, (best, verdict, hits) in zip(records, results):
                    if hits:
                        _save_analysis(batch, data, record, best, verdict, hits)
            return redirect("batch", pk=batch.pk)
    return render(request, "screening/upload.html", {"form": form})


def _save_analysis(batch, data, record, best, verdict, hits):
    analysis = Analysis.objects.create(
        batch=batch, filename=data["filename"], seq_type=data["seq_type"], query_id=record["query_id"], sequence=record["sequence"],
        length=len(record["sequence"]), best_score=best, verdict=verdict,
    )
    Hit.objects.bulk_create([Hit(analysis=analysis, **h) for h in hits])
    return analysis


def _batch_rows(batch):
    """Analyses of a batch with their top hit, best candidates first."""
    rows = [(a, a.hits.all()[0]) for a in batch.analyses.prefetch_related("hits")]
    rows.sort(key=lambda r: (-r[0].best_score, -r[1].identity, -r[1].bitscore))
    return rows


def batch_detail(request, pk):
    batch = get_object_or_404(Batch, pk=pk)
    rows = _batch_rows(batch)
    levels = [a.level for a, _ in rows]
    return render(request, "screening/batch.html", {
        "b": batch,
        "rows": rows,
        "counts": {lvl: levels.count(lvl) for lvl in ("strong", "possible", "weak")},
        "no_hit_count": batch.screened_count - len(rows),
    })


def batch_csv(request, pk):
    b = get_object_or_404(Batch, pk=pk)
    rows = _batch_rows(b)
    resp = HttpResponse(content_type="text/csv")
    resp["Content-Disposition"] = f'attachment; filename="curemediate_batch_{b.pk}.csv"'
    w = csv.writer(resp)
    nt = b.seq_type == "nucleotide"
    w.writerow(["File", b.filename])
    w.writerow(["Sequence type", b.get_seq_type_display()])
    w.writerow(["Date", b.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow(["Sequences in file", b.sequence_count])
    w.writerow(["Screened", b.screened_count])
    w.writerow(["Skipped", len(b.skipped)])
    w.writerow(["With Cu hits", len(rows)])
    w.writerow([])
    w.writerow([
        "Query ID", f"Length ({b.unit})", "Best score (/8)", "Verdict", "Top gene", "BacMet ID", "Organism",
        "Identity %", *(["Gene coverage %", "Query start", "Query end"] if nt else []),
        "E-value", "Bitscore", "C1", "C2", "C3", "C4",
    ])
    for a, h in rows:
        w.writerow([
            a.query_id, a.length, a.best_score, a.verdict, h.gene_name, h.bacmet_id, h.organism,
            f"{h.identity:.2f}", *(_nt_cols(h) if nt else []),
            f"{h.evalue:.2e}", f"{h.bitscore:.1f}", h.c1, h.c2, h.c3, h.c4,
        ])
    return resp


def _nt_cols(h):
    return ["" if h.coverage is None else f"{h.coverage:.1f}", h.query_start, h.query_end]


def result(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    hits = list(analysis.hits.all())
    top = hits[0] if hits else None
    return render(request, "screening/result.html", {
        "a": analysis,
        "hits": hits,
        "top": top,
        "breakdown": explain(top, analysis.seq_type) if top else [],
    })


def result_csv(request, pk):
    a = get_object_or_404(Analysis, pk=pk)
    resp = HttpResponse(content_type="text/csv")
    resp["Content-Disposition"] = f'attachment; filename="curemediate_analysis_{a.pk}.csv"'
    w = csv.writer(resp)
    nt = a.seq_type == "nucleotide"
    w.writerow(["Query ID", a.query_id])
    w.writerow(["File", a.filename])
    w.writerow(["Sequence type", a.get_seq_type_display()])
    w.writerow(["Date", a.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow([f"Length ({a.unit})", a.length])
    w.writerow(["Best score (/8)", a.best_score])
    w.writerow(["Verdict", a.verdict])
    w.writerow([])
    w.writerow([
        "Rank", "Gene", "BacMet ID", "Organism", "Identity %", "Align length",
        *(["Gene coverage %", "Query start", "Query end", "Description", "NCBI Gene"] if nt else []),
        "E-value", "Bitscore", "C1", "C2", "C3", "C4", "Total/8", "Europe PMC hits",
    ])
    for h in a.hits.all():
        w.writerow([
            h.rank, h.gene_name, h.bacmet_id, h.organism, f"{h.identity:.2f}", h.align_length,
            *([*_nt_cols(h), h.description, h.ncbi_gene_url] if nt else []),
            f"{h.evalue:.2e}", f"{h.bitscore:.1f}", h.c1, h.c2, h.c3, h.c4, h.total,
            "" if h.pmc_hits is None else h.pmc_hits,
        ])
    return resp


def history(request):
    return render(request, "screening/history.html", {
        "batches": Batch.objects.all(),
        "analyses": Analysis.objects.filter(batch__isnull=True),
    })
