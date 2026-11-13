from django.db import migrations, models


class Migration(migrations.Migration):
    """Ranking specification: per-sequence category/score/uncertainty replace the old per-hit C1-C4 scores."""

    dependencies = [
        ("screening", "0003_nucleotide"),
    ]

    operations = [
        migrations.RemoveField(model_name="analysis", name="best_score"),
        migrations.RemoveField(model_name="analysis", name="verdict"),
        migrations.AddField(model_name="analysis", name="category", field=models.CharField(blank=True, max_length=32)),
        migrations.AddField(model_name="analysis", name="score", field=models.PositiveSmallIntegerField(blank=True, null=True)),
        migrations.AddField(model_name="analysis", name="max_score", field=models.PositiveSmallIntegerField(blank=True, null=True)),
        migrations.AddField(model_name="analysis", name="score_pct", field=models.FloatField(blank=True, null=True)),
        migrations.AddField(model_name="analysis", name="uncertainty", field=models.CharField(blank=True, max_length=16)),
        migrations.AddField(model_name="analysis", name="flags", field=models.JSONField(blank=True, default=list)),
        migrations.AddField(model_name="analysis", name="breakdown", field=models.JSONField(blank=True, default=list)),
        migrations.AddField(model_name="analysis", name="conclusion", field=models.TextField(blank=True)),
        *[migrations.RemoveField(model_name="hit", name=f) for f in ("c1", "c2", "c3", "c4", "total", "pmc_hits")],
        migrations.RenameField(model_name="hit", old_name="coverage", new_name="ref_coverage"),
        migrations.AddField(model_name="hit", name="query_coverage", field=models.FloatField(blank=True, null=True)),
        migrations.AddField(model_name="hit", name="passes_gate", field=models.BooleanField(default=False)),
        migrations.AddField(model_name="hit", name="family", field=models.CharField(blank=True, max_length=32)),
    ]
