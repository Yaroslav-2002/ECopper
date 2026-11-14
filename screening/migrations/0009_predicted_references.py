from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [("screening", "0008_bacmet_exp_link")]

    operations = [
        migrations.AddField(model_name="bacmetgene", name="annotation_note", field=models.TextField(blank=True)),
        migrations.AddField(model_name="hit", name="dataset", field=models.CharField(blank=True, max_length=16)),
    ]
