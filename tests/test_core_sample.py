"""Targeted regression tests; also creates a small ROOT campaign for visual QA."""

import argparse
import json
import sys
import tempfile
import unittest
from pathlib import Path

import numpy as np
import uproot

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import core_sample as cs


def fixture(path, seed=12, run=666701):
    rng = np.random.default_rng(seed)
    positions, labels = [], []
    for x in range(3, 6):
        for y in range(3, 6):
            if (x, y) == (5, 5):  # Missing hole is never invented.
                continue
            n = 2400 if (x, y) == (4, 4) else 100
            center = np.array([(x - 4) * 2.54, (y - 4) * 1.524])
            sigma = np.array([.11 * 2.54, .09 * 1.524])
            p = rng.normal(size=(n, 2)) * sigma + center
            halo = rng.normal(size=(n // 5, 2)) * sigma * 3 + center
            positions.extend([p, halo])
            labels.extend([(x, y)] * (n + n // 5))
    p = np.concatenate(positions)
    labels = np.array(labels)
    n = len(p)
    data = {k: np.zeros(n, dtype=np.float64) for k in
            ("zfoil", "xfp", "xpfp", "yfp", "ypfp", "delta", "delta_low", "delta_high",
             "ytar", "xtar", "xptar", "yptar", "reactx", "reacty", "reactz", "xbpm_tar", "ybpm_tar")}
    data.update(run=np.full(n, run, np.int32), entry=np.arange(n, dtype=np.int64),
                foil=np.zeros(n, np.int32), ndel=np.full(n, 2, np.int32),
                xscol=labels[:, 0].astype(np.int32), yscol=labels[:, 1].astype(np.int32),
                xsieve=p[:, 0], ysieve=p[:, 1], delta=rng.uniform(-4.5, -.5, n),
                delta_low=np.full(n, -5.), delta_high=np.zeros(n))
    name = "rg01_test" if run == 666701 else "rg02_test"
    for suffix, axis, tree, column in (("a", "y", "TYCand", "xscol"), ("b", "x", "TXCand", "yscol")):
        directory = path / f"04{suffix}_candidate_trees_{axis}/root"
        directory.mkdir(parents=True, exist_ok=True)
        prefix = "Yscol" if axis == "y" else "Xscol"
        arrays = {k: v for k, v in data.items() if k != column}
        with uproot.recreate(directory / f"{prefix}Candidates_{name}.root") as root:
            root.mktree(tree, arrays)
    raw = path / f"source_{name}.root"
    mapping = {"gtr.y": "ytar", "gtr.x": "xtar", "react.x": "reactx", "react.y": "reacty",
               "react.z": "reactz", "gtr.dp": "delta", "gtr.ph": "yptar", "gtr.th": "xptar",
               "dc.y_fp": "yfp", "dc.yp_fp": "ypfp", "dc.x_fp": "xfp", "dc.xp_fp": "xpfp",
               "extcor.ysieve": "ysieve", "extcor.xsieve": "xsieve",
               "rb.raster.fr_xbpm_tar": "xbpm_tar", "rb.raster.fr_ybpm_tar": "ybpm_tar"}
    raw_data = {"H." + k: data[v] for k, v in mapping.items()}
    raw_data.update({"H.cer.npeSum": np.full(n, 10.), "H.cal.etottracknorm": np.full(n, 1.)})
    with uproot.recreate(raw) as root:
        root.mktree("T", raw_data)
    (path / "config").mkdir(exist_ok=True)
    table = path / "config/rungroups_test_inputs.tsv"
    cs.write_tsv(table, [dict(rungroup=name, optics_id=run, hms_angle_deg=12.490,
                             foils="0", nruns=1, runs="1544", rootfile=str(raw))])
    return data


class CoreTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = Path(self.tmp.name) / "HMS_test"
        self.data = fixture(self.path)
        self.cfg = cs.load_config(self.path)

    def tearDown(self):
        self.tmp.cleanup()

    def model(self, data=None):
        data = self.data if data is None else data
        holdout = cs.reserve(data, self.cfg, "rg01_test")
        spacing = cs.spacing_from_labels(data, np.flatnonzero(~holdout))
        model, rows = cs.select_slice(data, np.arange(len(holdout)), self.cfg, holdout, spacing)
        return holdout, model, rows

    def test_sparse_holes_and_missing_hole(self):
        holdout, model, rows = self.model()
        self.assertEqual(len(rows), 8)
        self.assertTrue(all(row["status"] == "accepted" for row in rows), rows)
        for hole, idx in cs.groups(self.data, ["xscol", "yscol"]).items():
            self.assertGreater(np.sum(model["quality"][idx] == 2), 10)
            self.assertGreater(np.sum((model["quality"][idx] == 2) & holdout[idx]), 0)

    def test_holdout_never_changes_density(self):
        holdout, first, _ = self.model()
        altered = {k: v.copy() for k, v in self.data.items()}
        altered["xsieve"][holdout] += 1000
        _, second, _ = self.model(altered)
        for key in ("density", "basin", "coremap", "xedges", "yedges"):
            np.testing.assert_array_equal(first[key], second[key])

    def test_joint_duplicate_and_large_entry(self):
        xfile = self.path / "04b_candidate_trees_x/root/XscolCandidates_rg01_test.root"
        yfile = self.path / "04a_candidate_trees_y/root/YscolCandidates_rg01_test.root"
        for path, tree in ((xfile, "TXCand"), (yfile, "TYCand")):
            with uproot.open(path) as root:
                data = root[tree].arrays(library="np")
            data["entry"] += 2 ** 33
            if tree == "TXCand":
                data = {k: np.append(v, v[:1]) for k, v in data.items()}
            with uproot.recreate(path) as root:
                root.mktree(tree, data)
        joint, excluded = cs.read_joint(xfile, yfile, 666701)
        self.assertEqual(excluded, [{"entry": 2 ** 33, "reason": "multiple_labels"}])
        self.assertEqual(int(joint["entry"][0]), 2 ** 33 + 1)

    def test_reordered_input_same_reserve(self):
        order = np.arange(len(self.data["entry"]))[::-1]
        reversed_data = {k: v[order] for k, v in self.data.items()}
        np.testing.assert_array_equal(cs.reserve(self.data, self.cfg, "rg01_test"),
                                      cs.reserve(reversed_data, self.cfg, "rg01_test")[order])

    def test_shared_budget_and_surplus(self):
        data = {}
        for name, n in (("busy", 2000), ("sparse", 200)):
            data[name] = dict(entry=np.arange(n), quality=np.full(n, 2), sample=np.zeros(n, np.int32),
                              zfoil=np.zeros(n), delta_low=np.full(n, -5.), delta_high=np.zeros(n),
                              xscol=np.full(n, 4), yscol=np.full(n, 4))
        cfg = {**self.cfg, "pools": {"busy": "shared", "sparse": "shared"}}
        report = cs.allocate(data, cfg)
        self.assertEqual(sum(row["fit"] for row in report), 400)
        self.assertEqual(sum(data["sparse"]["sample"] == 1), 160)
        self.assertEqual(sum(data["busy"]["sample"] == 1), 240)
        self.assertGreater(sum(data["sparse"]["sample"] == 3), 0)

    def test_campaign_limit_and_holdout_stays_fixed(self):
        protected, model, _ = self.model()
        data = {k: v.copy() for k, v in self.data.items()}
        data["quality"] = model["quality"]
        data["sample"] = np.where(protected, 2, 0).astype(np.int32)
        cfg = {**self.cfg, "fit_max": 80}
        cs.allocate({"rg01_test": data}, cfg)
        self.assertEqual(np.sum(data["sample"] == 1), 80)
        np.testing.assert_array_equal(data["sample"] == 2, protected)
        for idx in cs.groups(data, ["xscol", "yscol"]).values():
            self.assertGreater(np.sum(data["sample"][idx] == 1), 0)

    def test_low_statistics_not_forced_into_core(self):
        idx = np.flatnonzero((self.data["xscol"] == 5) & (self.data["yscol"] == 4))
        keep = np.ones(len(self.data["entry"]), bool)
        keep[idx[5:]] = False
        data = {k: v[keep] for k, v in self.data.items()}
        _, _, rows = self.model(data)
        row = next(r for r in rows if (r["xscol"], r["yscol"]) == (5, 4))
        self.assertEqual(row["status"], "low_stats")

    def test_mismatched_candidates_fail(self):
        xfile = self.path / "04b_candidate_trees_x/root/XscolCandidates_rg01_test.root"
        yfile = self.path / "04a_candidate_trees_y/root/YscolCandidates_rg01_test.root"
        with uproot.open(xfile) as root:
            data = root["TXCand"].arrays(library="np")
        data["xsieve"][0] += .1
        with uproot.recreate(xfile) as root:
            root.mktree("TXCand", data)
        with self.assertRaisesRegex(ValueError, "disagree in xsieve"):
            cs.read_joint(xfile, yfile, 666701)

    def test_campaign_outputs_repeat_and_partition(self):
        cs.run(self.path, "first")
        cs.run(self.path, "repeat")
        outputs = []
        for tag in ("first", "repeat"):
            with uproot.open(self.path / f"05c_core_sample/{tag}/root/CoreSample_rg01_test.root") as root:
                outputs.append(root["CoreSample"].arrays(library="np"))
        for key in outputs[0]:
            np.testing.assert_array_equal(outputs[0][key], outputs[1][key])
        a = outputs[0]
        self.assertTrue(np.all(a["core_keep"][np.isin(a["sample"], [1, 3])] == 1))
        self.assertEqual(len(np.unique(a["entry"])), len(a["entry"]))
        self.assertEqual(set(np.unique(a["sample"])), {0, 1, 2, 3})
        self.assertTrue((self.path / "05c_core_sample/first/plots/core_rg01_test_foil0_ndel2.png").exists())


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "--demo":
        parser = argparse.ArgumentParser()
        parser.add_argument("--demo", type=Path, required=True)
        args = parser.parse_args()
        fixture(args.demo.resolve())
        cs.run(args.demo.resolve(), "core")
    else:
        unittest.main()
