"""Viewer contract checks; these do not test the classification model."""
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np

from joint_surface.viewer import (counts, event_text, find_manifest, html_payload,
                                  load_campaign, sample_indices, trace_groups, write_html)


class ViewerTest(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.campaign = Path(self.temp.name)
        self.work = self.campaign / "00_joint_surface_work"
        self.work.mkdir()
        n = 120
        self.data = dict(entry=np.arange(n, dtype=np.int64) + 2**34,
                         xfp=np.linspace(-10, 10, n), xpfp=np.linspace(-.03, .02, n),
                         yfp=np.linspace(-8, 8, n), ypfp=np.linspace(-.02, .03, n),
                         delta=np.linspace(-10, 10, n, endpoint=False),
                         x_label=np.arange(n, dtype=np.int32) % 9,
                         y_label=np.arange(n, dtype=np.int32) % 7,
                         accepted_x=np.ones(n, bool), accepted_y=np.ones(n, bool),
                         reason_x=np.full(n, "accepted", dtype="U30"),
                         reason_y=np.full(n, "accepted", dtype="U30"),
                         proposed_x=np.arange(n) % 9, train_x=np.arange(n) % 2 == 0)
        self.data["accepted_x"][:4] = False
        self.data["reason_x"][:4] = "coverage"
        self.data["accepted_y"][5:7] = False
        self.data["reason_y"][5:7] = "occupancy"
        np.savez(self.work / "foil.npz", **self.data)
        self.manifest = {"schema_version": 1, "cases": [{"key": "rg01_foil0", "label": "foil 0 cm",
                         "events": "foil.npz", "predictions": {"main": "foil.npz"},
                         "delta_edges": [-10, -8, -5, 0, 5, 10]}]}
        (self.work / "viewer_manifest.json").write_text(json.dumps(self.manifest))

    def test_units_ids_counts_and_input_unchanged(self):
        source = (self.work / "foil.npz").read_bytes()
        case = load_campaign(self.campaign, 0)[0]
        self.assertEqual(find_manifest(self.campaign), self.work / "viewer_manifest.json")
        self.assertEqual(case.entry[0], 2**34)
        self.assertEqual(case.coordinates[0, 1], -30)
        self.assertEqual(case.coordinates[0, 3], -20)
        self.assertEqual(counts(case, "main"), dict(events=120, x=116, y=118, both=114))
        self.assertEqual(source, (self.work / "foil.npz").read_bytes())
        self.assertIn("rejected (coverage)", event_text(case, "main", 0))
        self.assertIn("train_x=True", event_text(case, "main", 0))

    def test_deterministic_bounded_sampling_includes_rare_rejections(self):
        case = load_campaign(self.campaign, 40)[0]
        np.testing.assert_array_equal(case.shown, sample_indices(case, 40))
        self.assertEqual(len(case.shown), 40)
        self.assertTrue(set(range(4)) <= set(case.shown))
        self.assertTrue({5, 6} <= set(case.shown))
        self.assertEqual(counts(case, "main", shown=True)["events"], 40)
        self.assertEqual(counts(case, "main")["events"], 120)

    def test_half_open_delta_and_color_dimension(self):
        case = load_campaign(self.campaign, 0)[0]
        case.coordinates[:3, 4] = [-10, -8, -7.999]
        self.assertEqual(counts(case, "main", 0)["events"], 10)
        groups = trace_groups(case, "main", "y", color_mode="x")
        self.assertTrue(all(name.startswith("X ") for name, *_ in groups))
        rejected = [ii for name, _, _, ii in groups if "rejected: coverage" in name][0]
        np.testing.assert_array_equal(rejected, np.arange(4))
        ii = [indices for *_, indices in trace_groups(case, "main", "x", event_filter="rejected")]
        self.assertEqual(sum(map(len, ii)), 6)

    def test_nonfinite_points_remain_in_counts_but_are_not_plotted(self):
        self.data["xfp"][0] = np.nan
        np.savez(self.work / "foil.npz", **self.data)
        case = load_campaign(self.campaign, 0)[0]
        self.assertNotIn(0, case.shown)
        self.assertEqual(counts(case, "main")["events"], 120)
        self.assertEqual(counts(case, "main", shown=True)["events"], 119)
        json.dumps(html_payload([case]), allow_nan=False)

    def test_inconsistent_prediction_order_rejected(self):
        self.data["entry"] = self.data["entry"][::-1]
        np.savez(self.work / "prediction.npz", **self.data)
        self.manifest["cases"][0]["predictions"]["main"] = "prediction.npz"
        (self.work / "viewer_manifest.json").write_text(json.dumps(self.manifest))
        with self.assertRaisesRegex(ValueError, "event order"):
            load_campaign(self.campaign)

    def test_offline_html_has_embedded_bundle_and_safe_json(self):
        case = load_campaign(self.campaign, 30)[0]
        case.label = "</script><img src=x>"
        local_js = self.campaign / "plotly.js"
        local_js.write_text("window.Plotly={version:'local-test'};")
        output = self.campaign / "view.html"
        write_html([case], output, plotly_js=local_js)
        html = output.read_text()
        self.assertIn("window.Plotly={version:'local-test'}", html)
        self.assertNotIn('<script src=', html)
        self.assertNotIn("__MEMBERSHIP_DATA__", html)
        self.assertIn("<\\/script><img src=x>", html)
        with self.assertRaisesRegex(ValueError, "Output exists"):
            write_html([case], output, plotly_js=local_js)

    def test_desktop_plot_can_render_without_browser_or_native_display(self):
        import matplotlib
        matplotlib.use("Agg", force=True)
        from matplotlib.figure import Figure
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        from joint_surface.desktop_viewer import draw_case, picked_index
        from mpl_toolkits.mplot3d import proj3d
        from types import SimpleNamespace
        case = load_campaign(self.campaign, 30)[0]
        fig = Figure(figsize=(12, 6))
        FigureCanvasAgg(fig)
        axes, artists, legends = draw_case(fig, case, "main")
        self.assertEqual(len(axes), 2)
        self.assertGreater(len(artists), 10)
        self.assertTrue(legends)
        fig.savefig(self.campaign / "native.png")
        self.assertGreater((self.campaign / "native.png").stat().st_size, 10000)
        # Event IDs remain correct after Matplotlib reorders points by depth.
        indices = case.shown
        expected = int(indices[10])
        p = case.coordinates[expected]
        x, y, _ = proj3d.proj_transform(p[0], p[1], p[4], axes[0].get_proj())
        pixel = axes[0].transData.transform((x, y))
        actual = picked_index(case, indices, axes[0], SimpleNamespace(x=pixel[0], y=pixel[1]))
        self.assertEqual(actual, expected)


if __name__ == "__main__":
    unittest.main()
