from pathlib import Path
from unittest import mock, skipUnless

from django.test import TestCase
from django.urls import reverse

from .models import Analysis
from .services import scoring
from .services.blast import db_exists
from .services.scoring import (
    FLAG_INDIRECT, FLAG_SPECIFICITY, NO_EVIDENCE, NOT_TRANSFERABLE, Ref, classify, evaluate, metals, ranking_key,
)

MCO = Ref("cueO", "Escherichia coli", "Copper (Cu)", ("Multicopper oxidase CueO", "multicopper oxidase CueO"))
ATPASE = Ref("copA", "Escherichia coli", "Copper (Cu), Silver (Ag)", ("Copper-exporting P-type ATPase", ""))
REGULATOR = Ref("cueR", "Escherichia coli", "Copper (Cu)", ("HTH-type transcriptional regulator CueR", ""))
NO_DESC = Ref("copB", "Achromobacter sp.", "Copper (Cu)", ())


def hit(bacmet_id, identity=90.0, evalue=1e-50, bitscore=300.0, qlen=100, rcov=100):
    """A raw BLAST hit on a 100-residue query and a 100-residue reference."""
    return {
        "bacmet_id": bacmet_id, "identity": identity, "align_length": 100, "evalue": evalue, "bitscore": bitscore,
        "query_start": 1, "query_end": 100, "query_length": qlen,
        "subject_start": 1, "subject_end": rcov, "subject_length": 100,
    }


class GateTests(TestCase):
    def test_no_hits(self):
        self.assertEqual(evaluate("Q", [], {})["category"], NO_EVIDENCE)

    def test_each_gate_condition(self):
        refs = {"B1": MCO}
        for bad in ({"evalue": 1e-9}, {"bitscore": 49.9}, {"rcov": 79}):
            with self.subTest(bad):
                r = evaluate("Q", [hit("B1", **bad)], refs)
                self.assertEqual(r["category"], NO_EVIDENCE)
                self.assertIsNone(r["score"])

    def test_query_coverage_not_checked(self):
        # A long query (e.g. a contig) with the reference fully covered passes, however little of it aligns.
        r = evaluate("Q", [hit("B1", qlen=10000)], {"B1": MCO})
        self.assertTrue(r["hits"][0]["passes_gate"])
        self.assertNotIn("query_coverage", r["hits"][0])
        self.assertNotEqual(r["category"], NO_EVIDENCE)

    def test_gate_uses_best_bitscore_hit(self):
        # A weaker hit listed first must not decide the outcome.
        r = evaluate("Q", [hit("B2", bitscore=60, rcov=50), hit("B1", bitscore=300)], {"B1": MCO, "B2": MCO})
        self.assertEqual(r["hits"][0]["bacmet_id"], "B1")
        self.assertEqual(r["category"], "High potential")

    def test_low_identity_not_transferable(self):
        r = evaluate("Q", [hit("B1", identity=39.9)], {"B1": MCO})
        self.assertEqual(r["category"], NOT_TRANSFERABLE)
        self.assertIsNone(r["score"])


class ScoringTests(TestCase):
    def points(self, r):
        return {row["key"]: row["points"] for row in r["breakdown"]}

    def test_cueo_example(self):
        r = evaluate("Query X", [hit("B1", identity=87)], {"B1": MCO})
        self.assertEqual(self.points(r), {"A": 3, "B": 2, "C": 2, "D": None, "E": 1, "F": None})
        self.assertEqual((r["score"], r["max_score"], r["score_pct"]), (8, 8, 100))
        self.assertEqual((r["category"], r["uncertainty"], r["flags"]), ("High potential", "Low", []))
        self.assertEqual(r["conclusion"], (
            "Query X is 87% identical to Multicopper oxidase CueO: high Cu bioremediation potential, low "
            "uncertainty. Recommended lab check: growth at elevated Cu (MIC) and Cu removal measured by ICP-MS."))

    def test_regulator_is_indirect(self):
        r = evaluate("Q", [hit("B1", identity=85)], {"B1": REGULATOR})
        self.assertEqual(self.points(r)["C"], 0)
        self.assertEqual(r["flags"], [FLAG_INDIRECT])
        self.assertEqual(r["uncertainty"], "Moderate")  # exactly one flag

    def test_atpase_below_60_percent(self):
        r = evaluate("Q", [hit("B1", identity=55), hit("B2", identity=50, bitscore=250)], {"B1": ATPASE, "B2": ATPASE})
        p = self.points(r)
        self.assertEqual((p["A"], p["B"], p["C"], p["E"]), (2, 1, 1, 0))
        self.assertEqual(r["flags"], [FLAG_SPECIFICITY])
        self.assertEqual((r["score"], r["max_score"], r["category"], r["uncertainty"]), (4, 8, "Moderate", "High"))

    def test_mixed_families(self):
        r = evaluate("Q", [hit("B1", identity=82), hit("B2", bitscore=200)], {"B1": MCO, "B2": ATPASE})
        self.assertEqual(self.points(r)["E"], 0)
        self.assertEqual(r["uncertainty"], "High")

    def test_failed_hits_ignored_for_consistency(self):
        r = evaluate("Q", [hit("B1"), hit("B2", bitscore=40)], {"B1": MCO, "B2": ATPASE})
        self.assertEqual(self.points(r)["E"], 1)

    def test_missing_description_skips_c_and_e(self):
        r = evaluate("Q", [hit("B1", identity=45)], {"B1": NO_DESC})
        p = self.points(r)
        self.assertEqual((p["A"], p["B"], p["C"], p["E"]), (1, 2, None, None))
        self.assertEqual((r["score"], r["max_score"], r["category"]), (3, 5, "Moderate"))  # 60%
        self.assertEqual(r["uncertainty"], "High")  # identity < 50%

    def test_unknown_reference_skips_b(self):
        r = evaluate("Q", [hit("B9", identity=45)], {})
        self.assertEqual((r["score"], r["max_score"], r["category"]), (1, 3, "Low"))


class RuleTests(TestCase):
    def test_classify_by_description_not_gene_name(self):
        self.assertEqual(classify(("Multicopper oxidase", ""))[0], "multicopper oxidase")
        self.assertEqual(classify(("Copper-exporting P-type ATPase",))[0], "P-type ATPase")
        self.assertEqual(classify(("CopA",))[0], None)
        self.assertEqual(classify(("Copper resistance protein A",))[0], None)

    def test_classify_rules(self):
        cases = {
            "Copper-sensing transcriptional repressor CsoR": ("regulator/sensor", 0),
            "Sensor histidine kinase CusS": ("regulator/sensor", 0),
            "Copper chaperone CopZ": ("Cu chaperone", 2),
            "periplasmic Cu(I)/Cu(II)-binding protein CopK": ("Cu-binding protein", 2),
            "metallothionein": ("metallothionein", 2),
            "Cation efflux system protein CusA": ("RND/CBA efflux", 1),
            "Zinc transporter ZupT": ("transporter", 1),
            "copper homeostasis membrane protein CopD": ("homeostasis protein", 0),
            "Periplasmic or exported protein": (None, None),
        }
        for text, expected in cases.items():
            with self.subTest(text):
                self.assertEqual(classify((text,))[:2], expected)

    def test_fallback_to_second_description(self):
        self.assertEqual(classify(("CopB", "copper-binding protein CopB"))[:2], ("Cu-binding protein", 2))

    def test_specific_family_beats_generic_one(self):
        golt = ("Mercuric transport protein periplasmic component", "putative cation transport ATPase")
        self.assertEqual(classify(golt)[0], "P-type ATPase")
        self.assertEqual(classify(("Copper-exporting P-type ATPase", "P-type ATPase;soluble Cu(+) chaperone"))[0],
                         "P-type ATPase")

    def test_metals(self):
        self.assertEqual(metals("Copper (Cu), Hydrochloric acid (HCl) [class: Acid]"), {"Cu"})
        self.assertEqual(metals("Copper (Cu), Zinc (Zn), Sodium Deoxycholate (SDC) [class: Acid]"), {"Cu", "Zn"})
        self.assertEqual(metals(""), set())

    def test_ordering(self):
        rows = [
            ("Low", 40, "Low", 90), ("High potential", 80, "High", 70), ("High potential", 80, "Low", 60),
            ("High potential", 90, "High", 50), (NO_EVIDENCE, None, "", 99), (NOT_TRANSFERABLE, None, "", 30),
            ("High potential", 80, "Low", 95),
        ]
        rows.sort(key=lambda r: ranking_key(*r))
        self.assertEqual(rows, [
            ("High potential", 90, "High", 50), ("High potential", 80, "Low", 95), ("High potential", 80, "Low", 60),
            ("High potential", 80, "High", 70), ("Low", 40, "Low", 90), (NOT_TRANSFERABLE, None, "", 30),
            (NO_EVIDENCE, None, "", 99),
        ])


def fake_screen(sequences, seq_type, query_ids):
    """Stand-in for BLAST: 'good' is a close CueO homolog, anything else a 55% CopA homolog."""
    return [evaluate(q, [hit("B1", identity=87)], {"B1": MCO}) if q == "good"
            else evaluate(q, [hit("B2", identity=55)], {"B2": ATPASE}) for q in query_ids]


class ViewTests(TestCase):
    FASTA = f">weak\n{'MKVLAAGIVG' * 10}\n>good\n{'MKVLAAGIVG' * 10}\n"

    @mock.patch("screening.views.screen_many", side_effect=fake_screen)
    def test_batch_report(self, _):
        resp = self.client.post(reverse("upload"), {"seq_type": "protein", "sequence_text": self.FASTA}, follow=True)
        self.assertEqual(resp.status_code, 200)
        body = resp.content.decode()
        self.assertIn("High potential", body)
        self.assertIn("The best-hit method can miss novel genes", body)
        self.assertLess(body.index(">good<"), body.index(">weak<"))  # Step 5 ordering
        csv_text = self.client.get(reverse("batch_csv", args=[resp.context["b"].pk])).content.decode()
        self.assertIn("Identity %,Reference coverage %", csv_text)
        self.assertNotIn("Query coverage", csv_text)
        self.assertIn("A Identity (/3)", csv_text)
        self.assertIn("not in database", csv_text)
        for a in Analysis.objects.all():
            self.assertContains(self.client.get(reverse("result", args=[a.pk])), a.category)
            self.assertContains(self.client.get(reverse("result_csv", args=[a.pk])), "Uncertainty")


EXAMPLES = Path(__file__).resolve().parent.parent / "examples"


@skipUnless(db_exists("protein"), "BLAST protein database not built")
class BlastIntegrationTests(TestCase):
    def test_cueo_example_passes_gate(self):
        seq = "".join(EXAMPLES.joinpath("cueO_ecoli.fasta").read_text().splitlines()[1:])
        r = scoring.screen_many([seq], "protein", ["cueO"])[0]
        best = r["hits"][0]
        self.assertTrue(best["passes_gate"])
        self.assertGreaterEqual(best["ref_coverage"], 80)
        self.assertNotIn(r["category"], (NO_EVIDENCE, NOT_TRANSFERABLE))
