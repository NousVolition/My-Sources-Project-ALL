"""Mechanism, symmetry, censoring, randomization, and analysis validation."""
import csv
from itertools import product
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

from decentralized_coordination.simulation import (
    Parameters, crossed_mappings, exact_b, matched_schedule, network, seed_for, simulate, step_b)
from decentralized_coordination.analysis import influence_analysis, paired_contrasts, summarize, wilson
from decentralized_coordination.run_study import main, run_blocks
from decentralized_coordination.verify_results import verify

IDS = list(range(9))


class CoordinationTests(unittest.TestCase):
    def test_rule_switches_isolated_color_only(self):
        self.assertEqual(step_b([0, 0, 0, 0, 1, 0, 0, 0, 0]), (0,) * 9)
        self.assertEqual(step_b([0, 0, 0, 0, 1, 1, 1, 1, 1]), (0, 0, 0, 0, 1, 1, 1, 1, 1))

    def test_simultaneous_not_in_place(self):
        state = [0, 1, 0, 1, 0, 0, 0, 0, 0]
        self.assertEqual(step_b(state), (0, 0, 1, 0, 0, 0, 0, 0, 0))
        self.assertEqual(state, [0, 1, 0, 1, 0, 0, 0, 0, 0])

    def test_exhaustive_rotation_reflection_complement(self):
        for x in product((0, 1), repeat=9):
            y = step_b(x)
            self.assertEqual(step_b(x[1:] + x[:1]), y[1:] + y[:1])
            self.assertEqual(step_b(x[::-1]), y[::-1])
            self.assertEqual(step_b(tuple(1-v for v in x)), tuple(1-v for v in y))

    def test_exact_enumeration_independent_majority_oracle(self):
        counts = {}
        for x in product((0, 1), repeat=9):
            original = x
            t = 0
            while len(set(x)) > 1 and t < 120:
                x = tuple(int(x[i-1] + x[i] + x[(i+1) % 9] >= 2) for i in range(9))
                t += 10
            expected = t if len(set(x)) == 1 else None
            r = simulate("B", original, IDS, 3)
            self.assertEqual(r["consensus_time"], expected)
            counts[expected] = counts.get(expected, 0) + 1
        self.assertEqual(counts, {0: 2, 10: 60, 20: 54, 30: 18, 40: 18, None: 360})
        exact = exact_b()
        self.assertEqual(sum(r["period"] > 1 for r in exact), 0)
        self.assertEqual(max(r["transient_ticks"] for r in exact), 4)

    def test_adjacent_pairs_are_frozen(self):
        for state in product((0, 1), repeat=9):
            next_state = step_b(state)
            for i in range(9):
                if state[i] == state[(i+1) % 9]:
                    self.assertEqual(next_state[i], state[i])
                    self.assertEqual(next_state[(i+1) % 9], state[i])

    def test_initial_consensus_zero_actions(self):
        for c in "ABC":
            r = simulate(c, [1] * 9, IDS, 11)
            self.assertEqual((r["consensus_time"], r["switches"], r["refusals"]), (0, 0, 0))
            self.assertEqual(r["first_switch_ids"], [])
            self.assertIsNone(r["copy_hhi"])

    def test_hierarchy_follows_identity_in_every_seat(self):
        for mapping in crossed_mappings(33):
            r = simulate("A", [1, 0, 1, 0, 1, 0, 1, 0, 1], mapping, 22, leader_id=3)
            self.assertEqual(r["consensus_time"], 10)
            self.assertEqual(r["final"], [r["initial"][mapping.index(3)]] * 9)
            self.assertEqual(r["credits_by_id"][3], r["switches"])
            self.assertEqual(r["copy_hhi"], 1)

    def test_tied_first_switches_not_single_arbitrary_leader(self):
        r = simulate("A", [0, 1, 1, 0, 0, 0, 0, 0, 0], IDS, 1)
        self.assertEqual(r["first_switch_ids"], [1, 2])

    def test_failure_is_censored_not_event_at_120(self):
        r = simulate("B", [0, 0, 0, 1, 1, 1, 1, 1, 1], IDS, 1)
        self.assertFalse(r["consensus"])
        self.assertIsNone(r["consensus_time"])
        self.assertEqual(r["capped_time"], 120)

    def test_event_at_120_is_success(self):
        with patch("decentralized_coordination.simulation.random.Random") as rng:
            rng.return_value.random.side_effect = [0.9] * 11 + [0.1]
            r = simulate("A", [0, 1, 0, 0, 0, 0, 0, 0, 0], IDS, 1, Parameters(compliance=.5))
        self.assertTrue(r["consensus"])
        self.assertEqual(r["consensus_time"], 120)
        self.assertEqual(r["refusals"], 11)

    def test_all_refusal_no_initiative_prevents_switches(self):
        r = simulate("C", [0, 1, 0, 1, 0, 1, 0, 1, 0], IDS, 44,
                     Parameters(refusal=1, innovation=0, activation=1))
        self.assertEqual(r["switches"], 0)
        self.assertGreater(r["refusals"], 0)
        self.assertFalse(r["consensus"])
        self.assertEqual(r["refusals"], r["differing_proposals"])

    def test_initiatives_are_separate_from_refusals(self):
        r = simulate("C", [0, 1, 0, 1, 0, 1, 0, 1, 0], IDS, 44,
                     Parameters(refusal=1, innovation=1, update="sync10"))
        self.assertEqual(r["refusals"], 0)
        self.assertEqual(r["switches"], r["initiatives"])
        self.assertEqual(r["first_initiative_ids"], IDS)
        self.assertEqual(r["state_revisits"], 11)
        self.assertEqual(r["adjacent_tick_reversals"], 99)
        self.assertIsNone(r["copy_hhi"])

    def test_rng_reproduction_and_homogeneous_relabeling(self):
        initial = [0, 1, 1, 0, 1, 0, 0, 0, 1]
        for c in "ABC":
            r = simulate(c, initial, IDS, 77, keep_trace=True)
            self.assertEqual(r, simulate(c, initial, IDS, 77, keep_trace=True))
        r = simulate("C", initial, IDS, 77, keep_trace=True)
        alt = simulate("C", initial, IDS[::-1], 77, keep_trace=True)
        self.assertEqual(r["trajectory"], alt["trajectory"])
        self.assertEqual(r["credits_by_seat"], alt["credits_by_seat"])

    def test_no_hierarchy_memory_when_weight_is_one(self):
        initial = [0, 1, 0, 1, 0, 1, 0, 1, 0]
        a = simulate("C", initial, IDS, 2, former_leader=None, keep_trace=True)
        b = simulate("C", initial, IDS, 2, former_leader=0, keep_trace=True)
        self.assertEqual(a["trajectory"], b["trajectory"])

    def test_copy_credits_conserved_excluding_initiatives(self):
        for c in "ABC":
            r = simulate(c, [0, 1, 0, 1, 0, 1, 0, 1, 0], IDS, 14)
            self.assertAlmostEqual(sum(r["credits_by_seat"]), r["switches"] - r["initiatives"])

    def test_balanced_order_and_crossing(self):
        specs = matched_schedule(600, 1)
        orders = [r["order"] for r in specs]
        self.assertEqual(len(set(orders)), 6)
        self.assertEqual({orders.count(x) for x in set(orders)}, {100})
        mappings = crossed_mappings(13)
        for seat in range(9):
            self.assertEqual(sorted(m[seat] for m in mappings), IDS)

    def test_matched_blocks_and_pilot_size(self):
        _, rows, traces = run_blocks(4, 22, "pilot", trace=True)
        self.assertEqual(len(rows), 12)
        self.assertEqual(len(traces), 12)
        for block in range(4):
            group = [r for r in rows if r["block"] == block]
            self.assertEqual(len({str(r["initial"]) for r in group}), 1)
            self.assertEqual(len({str(r["seat_to_id"]) for r in group}), 1)

    def test_distinct_reproducible_streams(self):
        self.assertEqual(seed_for(1, "A", 3), seed_for(1, "A", 3))
        self.assertNotEqual(seed_for(1, "A", 3), seed_for(1, "B", 3))

    def test_invalid_inputs(self):
        for kwargs in ({"refusal": -1}, {"activation": 2}, {"update": "none"}, {"former_leader_weight": 0}):
            with self.assertRaises(ValueError):
                Parameters(**kwargs)
        for condition, initial, mapping in (("D", [0]*9, IDS), ("B", [0]*8, IDS), ("C", [0]*9, [0]*9)):
            with self.assertRaises(ValueError):
                simulate(condition, initial, mapping, 1)
        with self.assertRaises(ValueError):
            simulate("B", [0]*9, IDS, 1, topology="star")
        self.assertEqual(len(network("star")[0]), 8)

    def test_intervals_and_empty_success_summary(self):
        self.assertLess(wilson(0, 100)[1], .04)
        self.assertGreater(wilson(100, 100)[0], .96)
        r = simulate("B", [0, 0, 0, 1, 1, 1, 1, 1, 1], IDS, 1)
        s = summarize([r], "B")
        self.assertIsNone(s["median_seconds_successes_only"])
        self.assertEqual(s["mean_capped_seconds"], 120)

    def test_paired_zero_contrast(self):
        rows = [{"condition": c, "block": b, "consensus": b % 2, "capped_time": b, "switches": b}
                for b in range(10) for c in "ABC"]
        for row in paired_contrasts(rows, 1, draws=20):
            self.assertEqual(row["mean_difference"], 0)
            self.assertEqual(row["ci95"], [0, 0])

    def test_influence_recovers_known_identity_and_seat_effects(self):
        rows = [{"replicate": rep, "mapping": m, "seat": seat, "identity": identity,
                 "copied_switches": 2 * seat + identity}
                for rep in range(10) for m, mapping in enumerate(crossed_mappings(rep))
                for seat, identity in enumerate(mapping)]
        a = influence_analysis(rows)
        self.assertAlmostEqual(a["variance_fraction"]["position"], .8)
        self.assertAlmostEqual(a["variance_fraction"]["identity"], .2)
        self.assertAlmostEqual(a["heldout_mse_improvement"]["both"], 1)

    def test_tie_aware_persistence(self):
        rows = [{"replicate": rep, "mapping": m, "seat": seat, "identity": identity,
                 "copied_switches": 1}
                for rep in range(5) for m, mapping in enumerate(crossed_mappings(rep))
                for seat, identity in enumerate(mapping)]
        a = influence_analysis(rows)
        for v in a["top_copy_source_persistence"].values():
            self.assertAlmostEqual(v, 1 / 9)

    def test_cli_small_run_and_manifest(self):
        with tempfile.TemporaryDirectory() as d:
            out = Path(d) / "results"
            result = main(["--output", str(out), "--blocks", "6", "--sensitivity-blocks", "6",
                           "--influence-replicates", "5", "--no-plots"])
            self.assertEqual(result["B_exact"]["successes_120"], 152)
            manifest = json.loads((out / "manifest.json").read_text())
            self.assertEqual(manifest["total_simulated_rounds"], 12 + 18 + 126 + 180)
            with (out / "rounds.csv").open(newline="") as f:
                self.assertEqual(len(list(csv.DictReader(f))), 30)
            self.assertIn("rounds.csv", manifest["data_sha256"])
            self.assertNotIn(b"\r", (out / "rounds.csv").read_bytes())
            self.assertNotIn(b"\r", (out / "summary.json").read_bytes())
            self.assertEqual(verify(out, out), len(manifest["data_sha256"]))
            with (out / "rounds.csv").open("a") as f:
                f.write("corruption\n")
            with self.assertRaises(AssertionError):
                verify(out, out)


if __name__ == "__main__":
    unittest.main()
