"""Compare one-pass samples with the pre-optimization ROOT reader."""
import contextlib
import hashlib
import io
from array import array
from pathlib import Path
import sys
import tempfile
import time

import ROOT

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import run_x_multifoils_fixedtheta as workflow
from spectrometer_config import from_campaign

quantile = workflow.quantile

# Frozen pre-optimization algorithm (e8a4410): independent read per foil/slice.
def legacy_split(
    rootfile: Path,
    ytar_cut_file: Path,
    foil: int,
    delta_min: float,
    delta_max: float,
    spec,
    verbose: bool = False,
) -> tuple[float | None, int]:
    started = time.monotonic()
    if verbose:
        print(f"SPLIT START foil={foil} delta=[{delta_min},{delta_max}) "
              f"opening {rootfile}", flush=True)
    input_file = ROOT.TFile.Open(str(rootfile), "READ")

    if not input_file or input_file.IsZombie():
        raise RuntimeError(f"cannot open {rootfile}")

    tree = input_file.Get("T") or input_file.Get("Tout")

    if not tree:
        raise RuntimeError(f"cannot find T or Tout in {rootfile}")

    cut_file = ROOT.TFile.Open(str(ytar_cut_file), "READ")

    if not cut_file or cut_file.IsZombie():
        raise RuntimeError(f"cannot open {ytar_cut_file}")

    ytar_cut = cut_file.Get(f"delta_vs_ytar_cut_foil{foil}")

    if not ytar_cut:
        raise RuntimeError(
            f"cannot find delta_vs_ytar_cut_foil{foil} "
            f"in {ytar_cut_file}"
        )

    values: list[float] = []

    total = tree.GetEntries()
    last_report = time.monotonic()
    if verbose:
        print(f"SPLIT READ total={total}", flush=True)
    for entry in range(total):
        if verbose and time.monotonic() - last_report >= 5.0:
            elapsed = time.monotonic() - started
            print(f"SPLIT PROGRESS foil={foil} delta=[{delta_min},{delta_max}) "
                  f"read={entry}/{total} selected={len(values)} "
                  f"elapsed={elapsed:.1f}s rate={entry / elapsed:.0f} events/s",
                  flush=True)
            last_report = time.monotonic()
        tree.GetEntry(entry)

        cer = float(getattr(tree, spec.cherenkov_branch))
        cal = float(getattr(tree, spec.branch("cal.etottracknorm")))
        delta = float(getattr(tree, spec.branch("gtr.dp")))
        ytar = float(getattr(tree, spec.branch("gtr.y")))
        xpfp = float(getattr(tree, spec.branch("dc.xp_fp")))

        if spec.name=="SHMS" and not spec.delta_min < delta < spec.delta_max:
            continue

        if cer <= 2.0 or cal <= 0.65:
            continue

        if not ytar_cut.IsInside(ytar, delta):
            continue

        if delta_min <= delta < delta_max:
            values.append(xpfp)

    if verbose:
        print(f"SPLIT DONE read={total}/{total} selected={len(values)} "
              f"elapsed={time.monotonic() - started:.1f}s", flush=True)
    input_file.Close()
    cut_file.Close()

    if len(values) < 1000:
        return None, len(values)

    q05 = quantile(values, 0.05)
    q95 = quantile(values, 0.95)

    if q05 is None or q95 is None:
        return None, len(values)

    return 0.5 * (q05 + q95), len(values)



def main():
    with tempfile.TemporaryDirectory(prefix="hallc-splits-") as folder:
        rootfile = Path(folder) / "events.root"
        cutfile = Path(folder) / "cuts.root"
        specs = [from_campaign("HMS_test"), from_campaign("SHMS_test")]
        names = [name for spec in specs for name in (
            spec.cherenkov_branch, spec.branch("cal.etottracknorm"),
            spec.branch("gtr.dp"), spec.branch("gtr.y"), spec.branch("dc.xp_fp"))]
        f = ROOT.TFile(str(rootfile), "RECREATE")
        t = ROOT.TTree("T", "synthetic events")
        data = {name: array("d", [0]) for name in names + ["unused.payload"]}
        for name, value in data.items():
            t.Branch(name, value, name + "/D")
        deltas = [-10, -8, -6, -4, -2, 0, 2, 4, 6, 8, 10, 16, 22]
        for i in range(26000):
            for spec in specs:
                data[spec.cherenkov_branch][0] = 2 if i % 17 == 0 else 7
                data[spec.branch("cal.etottracknorm")][0] = .65 if i % 19 == 0 else .9
                data[spec.branch("gtr.dp")][0] = deltas[i % len(deltas)]
                data[spec.branch("gtr.y")][0] = [-1, 0, 1][i % 3]
                data[spec.branch("dc.xp_fp")][0] = ((i * 37) % 1009) / 10000 - .05
            t.Fill()
        t.Write()
        f.Close()
        f = ROOT.TFile(str(cutfile), "RECREATE")
        # Overlapping foil cuts check independent membership; third is empty.
        for foil, (lo, hi) in enumerate([(-2, .5), (-.5, 2), (5, 6)]):
            cut = ROOT.TCutG(f"delta_vs_ytar_cut_foil{foil}", 5)
            for j, point in enumerate([(lo,-30),(hi,-30),(hi,30),(lo,30),(lo,-30)]):
                cut.SetPoint(j, *point)
            cut.Write()
        f.Close()
        checksums = [hashlib.sha256(p.read_bytes()).hexdigest() for p in (rootfile,cutfile)]
        edges = [-10, -6, -2, 2, 6, 10, 22]
        for spec in specs:
            expected = {
                (foil, ndel): legacy_split(rootfile, cutfile, foil, lo, hi, spec)
                for foil in range(3)
                for ndel, (lo, hi) in enumerate(zip(edges[:-1], edges[1:]))
            }
            output = io.StringIO()
            with contextlib.redirect_stdout(output):
                actual = workflow.compute_splits(rootfile, cutfile, 3, edges, spec, verbose=True)
            assert actual == expected, (spec.name, actual, expected)
            assert "passes=1" in output.getvalue()
            assert actual[(2, 0)] == (None, 0)
        assert checksums == [hashlib.sha256(p.read_bytes()).hexdigest() for p in (rootfile,cutfile)]
        print("PASS: HMS/SHMS counts and splits exactly match legacy; boundaries, PID thresholds, overlapping/empty foils; inputs unchanged")


if __name__ == "__main__":
    main()
