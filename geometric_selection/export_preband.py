#!/usr/bin/env python3
"""Export HMS/SHMS replay events after PID/foil selection, before X/Y bands.

Requires PyROOT and numpy. Inputs are opened READ; output must be a new directory.
Selections mirror make_yscol_candidate_tree.C, excluding its Y band selection.
"""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import uuid
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from spectrometer_config import from_campaign


BRANCHES = {
    "sumnpe": "H.cer.npeSum", "etracknorm": "H.cal.etottracknorm",
    "ytar": "H.gtr.y", "delta": "H.gtr.dp",
    "xfp": "H.dc.x_fp", "xpfp": "H.dc.xp_fp",
    "yfp": "H.dc.y_fp", "ypfp": "H.dc.yp_fp",
    "xsieve": "H.extcor.xsieve", "ysieve": "H.extcor.ysieve",
}


def metadata(path, run):
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        fields = [x.strip() for x in next(csv.reader([line]))]
        if not fields or fields[0] != str(run):
            continue
        if len(fields) < 7 or i + 2 >= len(lines):
            raise ValueError("Incomplete metadata for run %s" % run)
        nfoils, nedges = int(fields[3]), int(fields[5])
        foils = [float(x) for x in next(csv.reader([lines[i+1]])) if x.strip()]
        edges = [float(x) for x in next(csv.reader([lines[i+2]])) if x.strip()]
        if nfoils < 1 or len(foils) != nfoils or len(edges) != nedges or nedges < 2:
            raise ValueError("Metadata foil/delta counts disagree")
        if not all(math.isfinite(x) for x in foils + edges):
            raise ValueError("Nonfinite foil positions/delta boundaries")
        if any(a >= b for a, b in zip(edges, edges[1:])):
            raise ValueError("Delta boundaries must increase strictly")
        return dict(run=run, optics_id=fields[1], angle_deg=float(fields[2]),
                    sieve_flag=int(fields[4]), ymis=float(fields[6]),
                    foil_z_cm=foils, delta_edges=edges)
    raise ValueError("Run %s absent from %s" % (run, path))


def file_info(path, hash_content=True):
    path = path.resolve(strict=True)
    st = path.stat()
    info = dict(path=str(path), size_bytes=st.st_size, mtime_ns=st.st_mtime_ns)
    if hash_content:
        h = hashlib.sha256()
        with path.open("rb") as stream:
            for block in iter(lambda: stream.read(8*1024*1024), b""):
                h.update(block)
        info["sha256"] = h.hexdigest()
    else:
        info["sha256"] = None
        info["hash_note"] = "Not computed; enable --hash-replay for full-file SHA256"
    return info


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", type=Path, required=True, help="Uncut replay ROOT file")
    p.add_argument("--foil-cuts", type=Path, required=True, help="ROOT delta_vs_ytar_cut_foilN polygons")
    p.add_argument("--run-list", type=Path, required=True, help="DATfiles/list_of_optics_run.dat")
    p.add_argument("--metadata-run", type=int, required=True)
    p.add_argument("--campaign", required=True, help="Campaign path containing an HMS or HMS_... component")
    p.add_argument("--rungroup", help="Provenance tag; defaults to metadata optics ID")
    p.add_argument("--cal-cut", type=float, default=.65)
    p.add_argument("--cer-cut", type=float, default=6.)
    p.add_argument("--tree", choices=["auto", "T", "Tout"], default="auto",
                   help="auto prefers T, then Tout; original H.* branches are required")
    p.add_argument("--output", type=Path, required=True, help="Fresh directory; existing paths are refused")
    p.add_argument("--hash-replay", action="store_true", help="Hash complete replay file (may be slow)")
    return p.parse_args(argv)


def export(args):
    if os.path.lexists(args.output):
        raise FileExistsError("Output already exists; refusing to reuse: %s" % args.output)
    spec = from_campaign(args.campaign)
    branches = {key: spec.branch(value.split('.', 1)[1]) for key, value in BRANCHES.items()}
    branches['sumnpe'] = spec.cherenkov_branch
    if not all(math.isfinite(x) for x in (args.cal_cut, args.cer_cut)):
        raise ValueError('PID thresholds must be finite')
    info = metadata(args.run_list, args.metadata_run)
    if spec.name == 'SHMS' and info['sieve_flag'] != 1:
        raise ValueError('Only centered SHMS sieve is supported (SieveFlag=1)')
    import numpy as np
    import ROOT
    ROOT.gROOT.SetBatch(True)
    ROOT.DisableImplicitMT()  # deterministic source order; TCutG access remains serial
    replay = ROOT.TFile.Open(str(args.input), "READ")
    cuts_file = ROOT.TFile.Open(str(args.foil_cuts), "READ")
    try:
        if not replay or replay.IsZombie() or not cuts_file or cuts_file.IsZombie():
            raise OSError("Cannot open replay or foil-cut ROOT input")
        names = ["T", "Tout"] if args.tree == "auto" else [args.tree]
        tree = next((replay.Get(x) for x in names if replay.Get(x)), None)
        if not tree or not tree.InheritsFrom("TTree"):
            raise ValueError("Cannot find requested replay TTree")
        for branch in branches.values():
            leaf = tree.GetLeaf(branch)
            if not leaf or leaf.GetLenStatic() != 1 or leaf.GetLeafCount():
                raise ValueError("Missing or nonscalar replay branch: %s" % branch)
        cuts = []
        cut_records = []
        for foil in range(len(info["foil_z_cm"])):
            name = "delta_vs_ytar_cut_foil%d" % foil
            cut = cuts_file.Get(name)
            if not cut or not cut.InheritsFrom("TCutG") or cut.GetN() < 3:
                raise ValueError("Missing/invalid foil polygon: %s" % name)
            cuts.append(cut)
            cut_records.append(dict(name=name, points=[[float(cut.GetPointX(j)), float(cut.GetPointY(j))]
                                                       for j in range(cut.GetN())]))
        # C++ event loop avoids a slow Python/PyROOT call for each replay event.
        ns = "HmsPreband_" + uuid.uuid4().hex
        edge_text = ",".join(repr(x) for x in info["delta_edges"])
        ROOT.gInterpreter.Declare("""
          #include <TCutG.h>
          #include <vector>
          namespace %s {
            std::vector<TCutG*> cuts;
            const std::vector<double> edges = {%s};
            int foil(double y, double d) {
              int found = -1;
              for (unsigned i=0; i<cuts.size(); ++i)
                if (cuts[i]->IsInside(y,d)) found = i;
              return found;
            }
            int slice(double d) {
              for (unsigned i=0; i+1<edges.size(); ++i)
                if (d>=edges[i] && d<edges[i+1]) return i;
              return -1;
            }
          }
        """ % (ns, edge_text))
        for cut in cuts:
            getattr(ROOT, ns).cuts.push_back(cut)
        frame = ROOT.RDataFrame(tree)
        for label, branch in branches.items():
            frame = frame.Alias("pb_"+label, branch)
        baseline = frame.Filter('pb_etracknorm > %r && pb_sumnpe > %r && pb_delta > %r && pb_delta < %r' %
                                (args.cal_cut, args.cer_cut, spec.delta_min, spec.delta_max))
        selected = (baseline.Define("foil", ns+"::foil(pb_ytar,pb_delta)")
                    .Define("ndel", ns+"::slice(pb_delta)")
                    .Filter("foil >= 0 && ndel >= 0")
                    .Define("entry", "static_cast<Long64_t>(rdfentry_)"))
        count_baseline = baseline.Count()
        count_selected = selected.Count()
        args.output.mkdir(parents=True, exist_ok=False)
        keys = ["entry", "foil", "ndel"] + ["pb_"+x for x in branches]
        raw = selected.AsNumpy(keys)
        arrays = {x.removeprefix("pb_"): np.asarray(a).copy() for x, a in raw.items()}
        arrays["entry"] = arrays["entry"].astype(np.int64)
        arrays["foil"] = arrays["foil"].astype(np.int32)
        arrays["ndel"] = arrays["ndel"].astype(np.int32)
        edges = np.asarray(info["delta_edges"], dtype=np.float64)
        arrays["delta_low"] = edges[arrays["ndel"]]
        arrays["delta_high"] = edges[arrays["ndel"]+1]
        np.savez_compressed(args.output / "preband_events.npz", **arrays)
        provenance = dict(status="complete", arm=spec.name, tree=tree.GetName(),
                          campaign=args.campaign, rungroup=args.rungroup or info["optics_id"],
                          metadata=info, replay=file_info(args.input, args.hash_replay),
                          foil_cut_file=file_info(args.foil_cuts), run_list=file_info(args.run_list),
                          exporter=file_info(Path(__file__)), root_version=ROOT.gROOT.GetVersion(),
                          numpy_version=np.__version__, branch_mapping=branches,
                          selection=dict(pid='etracknorm > %g && sumnpe > %g' % (args.cal_cut, args.cer_cut),
                                         acceptance='%g < delta < %g' % (spec.delta_min, spec.delta_max), slices="low <= delta < high",
                                         foil="TCutG.IsInside(ytar,delta); last/highest matching foil wins",
                                         band_cuts="none", polygons=cut_records),
                          counts=dict(input=int(tree.GetEntries()), baseline=int(count_baseline.GetValue()),
                                      exported=int(count_selected.GetValue())),
                          output=file_info(args.output / "preband_events.npz"))
        with (args.output / "provenance.json").open("x") as stream:
            json.dump(provenance, stream, indent=2, allow_nan=False)
            stream.write("\n")
        print(json.dumps(provenance["counts"]))
        print(args.output / "preband_events.npz")
    finally:
        if cuts_file:
            cuts_file.Close()
        if replay:
            replay.Close()


def main():
    args = parse_args()
    try:
        export(args)
    except Exception as exc:
        print("ERROR: %s" % exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main())
