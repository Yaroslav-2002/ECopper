import csv

from django.db import transaction
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render

from .forms import UploadForm
from .models import Analysis, Hit
from .services.blast import BlastError
from .services.scoring import explain, screen


def upload(request):
    form = UploadForm(request.POST, request.FILES) if request.method == "POST" else UploadForm()
    if form.is_bound and form.is_valid():
        data = form.cleaned_data["file"]
        try:
            best, verdict, hits = screen(data["sequence"])
        except BlastError as exc:
            form.add_error(None, f"BLAST error: {exc}")
        else:
            with transaction.atomic():
                analysis = Analysis.objects.create(
                    **data, length=len(data["sequence"]), best_score=best, verdict=verdict
                )
                Hit.objects.bulk_create([Hit(analysis=analysis, **h) for h in hits])
            return redirect("result", pk=analysis.pk)
    return render(request, "screening/upload.html", {"form": form})


def result(request, pk):
    analysis = get_object_or_404(Analysis, pk=pk)
    hits = list(analysis.hits.all())
    top = hits[0] if hits else None
    return render(request, "screening/result.html", {
        "a": analysis,
        "hits": hits,
        "top": top,
        "breakdown": explain(top) if top else [],
    })


def result_csv(request, pk):
    a = get_object_or_404(Analysis, pk=pk)
    resp = HttpResponse(content_type="text/csv")
    resp["Content-Disposition"] = f'attachment; filename="curemediate_analysis_{a.pk}.csv"'
    w = csv.writer(resp)
    w.writerow(["Query ID", a.query_id])
    w.writerow(["File", a.filename])
    w.writerow(["Date", a.created_at.strftime("%Y-%m-%d %H:%M")])
    w.writerow(["Length (aa)", a.length])
    w.writerow(["Best score (/8)", a.best_score])
    w.writerow(["Verdict", a.verdict])
    w.writerow([])
    w.writerow([
        "Rank", "Gene", "BacMet ID", "Organism", "Identity %", "Align length", "E-value",
        "Bitscore", "C1", "C2", "C3", "C4", "Total/8", "Europe PMC hits",
    ])
    for h in a.hits.all():
        w.writerow([
            h.rank, h.gene_name, h.bacmet_id, h.organism, f"{h.identity:.2f}", h.align_length,
            f"{h.evalue:.2e}", f"{h.bitscore:.1f}", h.c1, h.c2, h.c3, h.c4, h.total,
            "" if h.pmc_hits is None else h.pmc_hits,
        ])
    return resp


def history(request):
    return render(request, "screening/history.html", {"analyses": Analysis.objects.all()})
