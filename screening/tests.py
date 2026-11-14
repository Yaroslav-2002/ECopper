from pathlib import Path
import gzip
from tempfile import TemporaryDirectory
from unittest import mock, skipUnless

from django.conf import settings
from django.core.files.uploadedfile import SimpleUploadedFile
from django.test import TestCase
from django.urls import reverse

from .management.commands.build_cu_db import predicted_genes, read_annotation, read_mapping
from .models import Analysis, BacMetGene, Batch
from .views import _save_analysis
from .services import ranking, scoring
from .services.blast import db_exists
from .services import blast
from .services.scoring import NO_HIT, Ref, evaluate

EXP_URL = "http://bacmet.biomedicine.gu.se/general_search_get_tables.pl?keywords="
CUEO = Ref("cueO", "Escherichia coli", "Copper (Cu)", references="PMID: 11080161",
           description="Multicopper oxidase CueO", bacmet_url=EXP_URL + "BAC0103")
COPA = Ref("copA", "Escherichia coli", "Copper (Cu), Silver (Ag)",
           references="10639134; 10.1128/JB.182.5.1323-1327.2000", description="Copper-exporting P-type ATPase")


def hit(bacmet_id, identity=90.0, evalue=1e-50, bitscore=300.0, rcov=100):
    """A raw BLAST hit covering `rcov` residues of a 100-residue reference."""
    return {
        "bacmet_id": bacmet_id, "identity": identity, "align_length": 100, "evalue": evalue, "bitscore": bitscore,
        "query_start": 1, "query_end": 100, "query_length": 100,
        "subject_start": 1, "subject_end": rcov, "subject_length": 100,
    }


def score(**kw):
    args = {"dataset": "experimental", "compounds": "Copper (Cu)", "identity": 90, "coverage": 100,
            "evalue": 1e-50, "references": "12345"}
    return ranking.score(**{**args, **kw})


class RankingRuleTests(TestCase):
    def test_blast_bins(self):
        cases = {(100, 80): 5, (80, 100): 5, (79.9, 100): 4, (60, 100): 4, (59, 100): 3, (50, 80): 3,
                 (39.9, 100): 2, (20, 100): 2, (19.9, 100): 0, (100, 10): 0}
        for (identity, coverage), points in cases.items():
            with self.subTest(identity=identity, coverage=coverage):
                self.assertEqual(score(identity=identity, coverage=coverage)["score_blast"], points)

    def test_evalue_cap(self):
        r = score(evalue=1e-3)
        self.assertEqual(r["score_blast"], 1)
        self.assertIn("e-value > 1e-05 (BLAST capped)", r["notes"])
        self.assertEqual(score(evalue=1e-5)["score_blast"], 5)  # the cutoff itself is not capped
        self.assertEqual(score(evalue=1e-3, identity=10)["score_blast"], 0)  # a cap never raises points

    def test_bacmet_evidence(self):
        self.assertEqual(score(dataset="experimental")["score_bacmet"], 3)
        self.assertEqual(score(dataset="Predicted")["score_bacmet"], 1)
        r = score(dataset="")
        self.assertEqual(r["score_bacmet"], 0)
        self.assertIn("BacMet dataset missing", r["notes"])

    def test_other_metals(self):
        self.assertEqual(score(compounds="Copper (Cu), Zinc (Zn)")["score_metals"], 2)
        self.assertEqual(score(compounds="Copper (Cu)")["score_metals"], 1)
        # Acids and other non-metals are not metals.
        self.assertEqual(score(compounds="Copper (Cu), Hydrochloric acid (HCl) [class: Acid]")["score_metals"], 1)
        r = score(compounds=None)
        self.assertEqual(r["score_metals"], 0)
        self.assertIn("compounds missing", r["notes"])
        for compounds in ("Triclosan", "Hydrochloric acid (HCl) [class: Acid]", "Zinc (Zn)", "unknown"):
            with self.subTest(compounds=compounds):
                self.assertEqual(score(compounds=compounds)["score_metals"], 0)
        self.assertEqual(score(compounds="Cu, Zn")["score_metals"], 2)

    def test_metals_in(self):
        self.assertEqual(ranking.metals_in("Copper (Cu), Zinc (Zn), Sodium Deoxycholate (SDC) [class: Acid]"),
                         {"Cu", "Zn"})
        self.assertEqual(ranking.metals_in("copper; silver"), {"Cu", "Ag"})
        self.assertEqual(ranking.metals_in(""), set())

    def test_missing_blast_still_scored(self):
        r = score(identity=None)
        self.assertEqual((r["score_blast"], r["total_score"], r["tier"]), (0, 4, "Low"))
        self.assertIn("BLAST missing", r["notes"])

    def test_total_and_tiers(self):
        self.assertEqual((score()["total_score"], score()["tier"]), (9, "High"))
        r = score(identity=50, compounds="Copper (Cu), Zinc (Zn)")
        self.assertEqual((r["total_score"], r["tier"]), (8, "High"))
        r = score(identity=45)
        self.assertEqual((r["total_score"], r["tier"]), (7, "Medium"))
        r = score(identity=45, dataset="predicted")
        self.assertEqual((r["total_score"], r["tier"]), (5, "Medium"))
        self.assertEqual(score(identity=30, dataset="predicted")["tier"], "Low")
        self.assertEqual(ranking.MAX_SCORE, 10)

    def test_references(self):
        n, links = ranking.parse_references(
            "PMID: 10639134; 11283292 ; https://pubmed.ncbi.nlm.nih.gov/7783648/; doi:10.1016/j.femsre.2003.08.001.; "
            "CueO is a multi-copper oxidase; -")
        self.assertEqual(n, 5)
        self.assertEqual(links, [
            "https://pubmed.ncbi.nlm.nih.gov/10639134/", "https://pubmed.ncbi.nlm.nih.gov/11283292/",
            "https://pubmed.ncbi.nlm.nih.gov/7783648/", "https://doi.org/10.1016/j.femsre.2003.08.001",
        ])
        r = score(references="Some title only")
        self.assertEqual((r["n_references"], r["article_links"]), (1, []))
        self.assertIn("1 reference(s) without PMID/DOI", r["notes"])
        self.assertIn("no references", score(references="")["notes"])

    def test_references_not_scored(self):
        self.assertEqual(score(references="")["total_score"], score(references="1; 2; 3")["total_score"])

    def test_ranking_key(self):
        rows = [(9, 70, 1e-30), (10, 85, 1e-50), (9, 90, 1e-40), (9, 90, 1e-60), (None, None, None)]
        rows.sort(key=lambda r: ranking.ranking_key(*r))
        self.assertEqual(rows, [(10, 85, 1e-50), (9, 90, 1e-60), (9, 90, 1e-40), (9, 70, 1e-30), (None, None, None)])


class EvaluateTests(TestCase):
    def test_no_hits(self):
        r = evaluate("Q", [], {})
        self.assertEqual((r["tier"], r["total_score"], r["hits"]), (NO_HIT, 0, []))

    def test_cueo_example(self):
        r = evaluate("Query X", [hit("B1", identity=87)], {"B1": CUEO})
        top = r["hits"][0]
        self.assertEqual((top["score_blast"], top["score_bacmet"], top["score_metals"]), (5, 3, 1))
        self.assertEqual((r["total_score"], r["tier"]), (9, "High"))
        self.assertEqual(top["article_links"], ["https://pubmed.ncbi.nlm.nih.gov/11080161/"])
        self.assertEqual(top["notes"], "")
        self.assertEqual(r["conclusion"], "Query X: top candidate cueO (B1), 87% identity over 100% of the "
                                          "reference, scores 9/10 (High tier).")

    def test_similarity_uses_reference_coverage(self):
        top = evaluate("Q", [hit("B1", identity=90, rcov=50)], {"B1": CUEO})["hits"][0]
        self.assertEqual((top["ref_coverage"], top["similarity"], top["score_blast"]), (50, 45, 3))

    def test_hits_ranked_by_score_not_bitscore(self):
        # cueO has the higher bitscore, but copA scores more (other metal) at equal similarity.
        r = evaluate("Q", [hit("B1", bitscore=400), hit("B2", bitscore=300)], {"B1": CUEO, "B2": COPA})
        self.assertEqual([h["bacmet_id"] for h in r["hits"]], ["B2", "B1"])
        self.assertEqual([h["rank"] for h in r["hits"]], [1, 2])
        self.assertEqual(r["total_score"], 10)

    def test_unknown_reference_flagged(self):
        r = evaluate("Q", [hit("B9")], {})
        self.assertEqual((r["hits"], r["total_score"], r["tier"]), ([], 0, NO_HIT))

    def test_predicted_only_and_experimental_priority(self):
        pred = Ref("copA", compounds="Cu, Ag", dataset="predicted")
        r = evaluate("Q", [hit("PRE1", identity=100)], {"PRE1": pred})
        self.assertEqual(r["hits"][0]["score_bacmet"], 1)
        r = evaluate("Q", [hit("PRE1", identity=100), hit("B1", identity=30)], {"PRE1": pred, "B1": CUEO})
        self.assertEqual(r["hits"][0]["dataset"], "experimental")
        self.assertEqual(r["hits"][0]["score_bacmet"], 3)
        r = evaluate("Q", [hit("B1", evalue=0.1)], {"B1": CUEO})
        self.assertEqual((r["hits"], r["total_score"]), ([], 0))


class PredictedImportTests(TestCase):
    @mock.patch("screening.services.blast.db_exists", return_value=True)
    @mock.patch("screening.services.blast._run")
    def test_both_databases_searched_without_losing_predicted_identifiers(self, run, exists):
        run.side_effect = [
            "q0\tBAC0001\t90\t100\t1e-50\t300\t1\t100\t100\t1\t100\t100\tgene\n",
            "q0\tPRE123\t95\t100\t1e-60\t350\t1\t100\t100\t1\t100\t100\tBAC0001 family annotation\n",
        ]
        results = blast.run_blast_many(["MKVLAAGIVG" * 10])
        self.assertEqual([h["bacmet_id"] for h in results[0]], ["BAC0001", "PRE123"])
        self.assertEqual([c.args[0][c.args[0].index("-db") + 1] for c in run.call_args_list],
                         [settings.CU_DB_NAME, settings.CU_PRED_DB_NAME])

    def test_family_annotation_is_conservative_and_keeps_predicted_identity(self):
        rows = [
            {"bacmet_id": "BAC1", "gene_name": "copA", "compounds": "Copper (Cu), Silver (Ag)"},
            {"bacmet_id": "BAC2", "gene_name": "copA", "compounds": "Copper (Cu)"},
            {"bacmet_id": "BAC3", "gene_name": "cueR/ybbI", "compounds": "Copper (Cu), Zinc (Zn)"},
            {"bacmet_id": "BAC4", "gene_name": "other", "compounds": "Zinc (Zn)"},
        ]
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "pred.fasta.gz"
            with gzip.open(path, "wt") as f:
                for name in ("copA", "cueR/ybbI", "other", "unmapped"):
                    f.write(f">gi|1|ref|WP_{name}| protein [Predicted organism] [{name}]\nMKVLAAGIVG\n")
                f.write(">gi|2|sp|P12345|PROTEIN description [copA]\nMKVLAAGIVG\n")
            genes = list(predicted_genes(path, rows))
        self.assertEqual(len(genes), 3)
        self.assertEqual([g.compounds for g in genes], ["Cu", "Cu, Zn", "Cu"])
        self.assertEqual(genes[-1].organism, "")
        for gene in genes:
            self.assertEqual(gene.dataset, "predicted")
            self.assertTrue(gene.bacmet_id.startswith("PRE"))
            self.assertEqual(gene.references, "")
            self.assertEqual(gene.ncbi_gene_url, "")
            self.assertIn("inferred", gene.annotation_note)


def fake_screen(sequences, seq_type, query_ids):
    """Stand-in for BLAST: 'good' is a close copA homolog, anything else a 45% cueO homolog."""
    return [evaluate(q, [hit("B2", identity=95)], {"B2": COPA}) if q == "good"
            else evaluate(q, [hit("B1", identity=45)], {"B1": CUEO}) for q in query_ids]


class ViewTests(TestCase):
    FASTA = f">weak\n{'MKVLAAGIVG' * 10}\n>good\n{'MKVLAAGIVG' * 10}\n"

    @mock.patch("screening.views.screen_many")
    def test_multiple_sequences_rejected_before_screening(self, screen):
        for seq_type in ("protein", "nucleotide"):
            for source in ("sequence_text", "file"):
                with self.subTest(seq_type=seq_type, source=source):
                    value = (SimpleUploadedFile("input.fasta", self.FASTA.encode())
                             if source == "file" else self.FASTA)
                    resp = self.client.post(reverse("upload"), {"seq_type": seq_type, source: value})
                    self.assertContains(resp, "Only one sequence can be analyzed at a time.")
                    self.assertContains(resp, "Found 2 sequences.")
        screen.assert_not_called()
        self.assertFalse(Analysis.objects.exists())
        self.assertFalse(Batch.objects.exists())

    @mock.patch("screening.views.screen_many", side_effect=fake_screen)
    def test_single_sequence_upload(self, screen):
        for seq_type, sequence in (("protein", "MKVLAAGIVG" * 10), ("nucleotide", "ACGT" * 10)):
            for source in ("sequence_text", "file"):
                with self.subTest(seq_type=seq_type, source=source):
                    text = f">good\n{sequence}\n"
                    value = SimpleUploadedFile("input.fasta", text.encode()) if source == "file" else text
                    resp = self.client.post(reverse("upload"), {"seq_type": seq_type, source: value}, follow=True)
                    self.assertContains(resp, "10/10")
                    analysis = Analysis.objects.latest("pk")
                    self.assertEqual(analysis.sequence, sequence)
                    self.assertIsNone(analysis.batch_id)
                    screen.assert_called_with([sequence], seq_type, ["good"])

    @mock.patch("screening.views.screen_many")
    def test_predicted_report_labels_evidence_and_compound_source(self, screen):
        url = "https://www.ncbi.nlm.nih.gov/protein/WP_011730513.1"
        pred = Ref("copA", compounds="Cu", dataset="predicted",
                   bacmet_url=url,
                   annotation_note="Predicted family annotation inferred from EXP mapping: BAC0077")
        screen.return_value = [evaluate("good", [hit("PRE123")], {"PRE123": pred})]
        resp = self.client.post(reverse("upload"), {
            "seq_type": "protein", "sequence_text": ">good\n" + "MKVLAAGIVG" * 10,
        }, follow=True)
        self.assertContains(resp, "is a predicted BacMet homolog, not experimentally confirmed")
        self.assertContains(resp, "Predicted compound annotations are inferred")
        self.assertNotContains(resp, "its BacMet EXP entry PRE123")
        self.assertContains(resp, f'href="{url}"')
        self.assertContains(resp, "Source protein record (NCBI)")
        analysis = Analysis.objects.get()
        self.assertEqual(analysis.hits.get().score_bacmet, 1)
        self.assertEqual(analysis.hits.get().score_metals, 1)
        self.assertContains(self.client.get(reverse("result_csv", args=[analysis.pk])), "BacMet dataset")
        self.assertContains(self.client.get(reverse("result_csv", args=[analysis.pk])), url)

    def test_saved_predicted_hit_resolves_reference_accession(self):
        reference = BacMetGene.objects.create(bacmet_id="PREold", gene_name="copA", dataset="predicted",
                                              accession="WP_011730513.1", sequence="MK")
        a = Analysis.objects.create(filename="x", query_id="old", sequence="MK", length=2)
        h = a.hits.create(rank=1, bacmet_id=reference.bacmet_id, gene_name="copA", dataset="predicted",
                          identity=90, align_length=100, evalue=1e-50, bitscore=300)
        url = "https://www.ncbi.nlm.nih.gov/protein/WP_011730513.1"
        self.assertEqual(h.bacmet_link, url)
        self.assertEqual(scoring._references({reference.bacmet_id})[reference.bacmet_id].bacmet_url, url)
        self.assertContains(self.client.get(reverse("result", args=[a.pk])), f'href="{url}"')
        self.assertContains(self.client.get(reverse("result_csv", args=[a.pk])), url)

    @mock.patch("screening.views.screen_many", side_effect=fake_screen)
    def test_batch_report(self, _):
        # Existing batch reports remain readable, although new uploads require one sequence.
        batch = Batch.objects.create(filename="old.fasta", sequence_count=2, screened_count=2)
        data = {"filename": "old.fasta", "seq_type": "protein"}
        for query_id in ("weak", "good"):
            record = {"query_id": query_id, "sequence": "MKVLAAGIVG" * 10}
            _save_analysis(batch, data, record, fake_screen([record["sequence"]], "protein", [query_id])[0])
        resp = self.client.get(reverse("batch", args=[batch.pk]))
        self.assertEqual(resp.status_code, 200)
        body = resp.content.decode()
        self.assertIn("10/10", body)
        self.assertIn("Medium", body)
        self.assertIn("Confirm results experimentally", body)
        self.assertLess(body.index(">good<"), body.index(">weak<"))  # higher score first
        csv_text = self.client.get(reverse("batch_csv", args=[resp.context["b"].pk])).content.decode()
        self.assertIn("score_blast (/5),score_bacmet (/3),score_metals (/2),total_score (/10),Tier", csv_text)
        self.assertIn("https://pubmed.ncbi.nlm.nih.gov/10639134/ ; https://doi.org/10.1128/JB.182.5.1323-1327.2000",
                      csv_text)
        self.assertIn(f"{EXP_URL}BAC0103", csv_text)  # stored BacMet EXP entry link
        self.assertIn(f"{EXP_URL}copA", csv_text)  # no entry link: BacMet EXP search by gene name
        self.assertNotIn("connect_to_ALL", csv_text)  # never the predicted database
        for a in Analysis.objects.all():
            page = self.client.get(reverse("result", args=[a.pk]))
            self.assertContains(page, f"{a.tier} tier")
            self.assertContains(page, "Score breakdown")
            self.assertContains(page, "is experimentally confirmed")
            self.assertContains(self.client.get(reverse("result_csv", args=[a.pk])), "n_references")
        self.assertContains(self.client.get(reverse("history")), "Batch screens")

    def test_old_results_still_render(self):
        a = Analysis.objects.create(filename="x", query_id="old", sequence="MK", length=2)
        a.hits.create(rank=1, bacmet_id="B1", gene_name="cueO", identity=90, align_length=100, evalue=1e-50,
                      bitscore=300)
        self.assertContains(self.client.get(reverse("result", args=[a.pk])), "Old scoring scheme")
        self.assertContains(self.client.get(reverse("result_csv", args=[a.pk])), "Old scoring")
        self.assertContains(self.client.get(reverse("history")), "Old scoring")


EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@skipUnless(db_exists(), "BLAST protein database not built")
class BlastIntegrationTests(TestCase):
    @classmethod
    def setUpTestData(cls):
        # The test database starts empty: load the copper reference genes as build_cu_db does (without BLAST).
        rows = [r for r in read_mapping(settings.BACMET_MAPPING) if scoring.has_copper(r["compounds"])]
        annotation = read_annotation(settings.CU_ANNOTATION)
        BacMetGene.objects.bulk_create(
            BacMetGene(sequence="", **r, **annotation.get(r["bacmet_id"], {})) for r in rows)
        BacMetGene.objects.bulk_create(predicted_genes(settings.BACMET_PREDICTED_FASTA,
                                                     read_mapping(settings.BACMET_MAPPING)), batch_size=500)

    def read(self, name):
        return "".join(EXAMPLES.joinpath(name).read_text().splitlines()[1:])

    def test_cueo_example_scores_high(self):
        r = scoring.screen_many([self.read("cueO_ecoli.fasta")], "protein", ["cueO"])[0]
        best = r["hits"][0]
        self.assertGreaterEqual(best["similarity"], 80)
        self.assertEqual(best["score_blast"], 5)
        self.assertEqual(r["tier"], "High")
        self.assertEqual(best["dataset"], "experimental")
        self.assertTrue(any(h["dataset"] == "predicted" for h in r["hits"]))

    def test_dna_gene_found_with_blastx_on_both_strands(self):
        protein = scoring.screen_many([self.read("copA_ecoli.fasta")], "protein", ["copA protein"])[0]
        dna = self.read("copA_ecoli_gene.fna")
        reverse = dna[::-1].translate(str.maketrans("ACGT", "TGCA"))
        for label, seq in (("forward", dna), ("reverse complement", reverse)):
            with self.subTest(label):
                r = scoring.screen_many([seq], "nucleotide", [label])[0]
                best = r["hits"][0]
                self.assertEqual(best["bacmet_id"], protein["hits"][0]["bacmet_id"])
                self.assertGreaterEqual(best["identity"], 99)
                self.assertEqual((r["total_score"], r["tier"]), (protein["total_score"], protein["tier"]))
