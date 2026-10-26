#!/usr/bin/env python3
"""Write conventional cut/diagnostic ROOT files using only PyROOT and NumPy.

Invoked as a file so a separate PyROOT interpreter need not import SciPy.
The numerical runner supplies a JSON payload and per-slice NPZ assignments.
All output files must be new. TCutG membership is checked before publication.
"""
import argparse
import json
import os
import re
from pathlib import Path
from array import array
import sys
import uuid
import numpy as np


def fresh(path):
    path = Path(path)
    if os.path.lexists(path):
        raise FileExistsError('Refusing to overwrite '+str(path))
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def write(payload):
    import ROOT
    ROOT.gROOT.SetBatch(True)
    axis, arm = payload['axis'], payload['spectrometer']
    prefix = 'H' if arm == 'HMS' else 'P'
    cut_path = fresh(payload['cut_output'])
    for row in payload['slices']:
        fresh(row['diagnostic_root'])
    namespace = 'GeometryCutQA_'+uuid.uuid4().hex
    ROOT.gInterpreter.Declare('''
      #include <TCutG.h>
      #include <vector>
      namespace %s {
        std::vector<TCutG*> cuts;
        std::vector<int> qa(const std::vector<double>& x, const std::vector<double>& y) {
          std::vector<int> out(2*x.size(),0);
          for (size_t k=0;k<x.size();++k) {
            out[k]=-1;
            for (size_t c=0;c<cuts.size();++c) if(cuts[c]->IsInside(x[k],y[k])) {
              out[k]=c; ++out[k+x.size()];
            }
          }
          return out;
        }
      }
    ''' % namespace)
    engine = getattr(ROOT, namespace)
    all_cuts, checks = [], []
    for row in payload['slices']:
        with np.load(row['assignments']) as source:
            data = {k: source[k].copy() for k in source.files}
        raw = np.c_[data[axis+'pfp'], data[axis+'fp']]
        model_label = data['model_label'].astype(int)
        engine.cuts.clear()
        slice_cuts, indices = [], []
        names = set()
        for component, polygon in enumerate(row['polygons']):
            # Y's provisional IDs identify disconnected components individually.
            # Physical yscol is supplied later by the unchanged sieve relabeler.
            index = int(polygon['index'])
            if axis == 'x':
                name = 'hXpFpXFp_cut_candidate_%02d_zone_global_nfoil_%d_ndel_%d' % (index, row['foil'], row['ndel'])
            else:
                name = 'hYpFpYFp_cut_yscol_%d_nfoil_%d_ndel_%d' % (component, row['foil'], row['ndel'])
            if name in names:
                raise ValueError('Duplicate X component ID; use distinct IDs for manual pieces: '+name)
            names.add(name)
            points = np.asarray(polygon['points'], float)
            if points.ndim != 2 or points.shape[1] != 2 or len(points) < 4 or not np.isfinite(points).all():
                raise ValueError('Invalid polygon '+name)
            cut = ROOT.TCutG(name, len(points))
            cut.SetVarX(prefix+'.dc.'+axis+'p_fp')
            cut.SetVarY(prefix+'.dc.'+axis+'_fp')
            cut.SetTitle('%s provisional band %d component %d' % (axis.upper(), index, component))
            for k, (x, y) in enumerate(points):
                cut.SetPoint(k, float(x), float(y))
            slice_cuts.append(cut)
            indices.append(index)
            engine.cuts.push_back(cut)
        values = np.asarray(engine.qa(ROOT.std.vector('double')(np.ascontiguousarray(raw[:, 0])),
                                      ROOT.std.vector('double')(np.ascontiguousarray(raw[:, 1]))), dtype=int)
        region, matches = values[:len(raw)], values[len(raw):]
        predicted = np.full(len(raw), -1, int)
        valid = region >= 0
        predicted[valid] = np.asarray(indices)[region[valid]]
        mismatch = float(np.mean(predicted != model_label))
        check = dict(foil=row['foil'], ndel=row['ndel'], events=len(raw),
                     cut_components=len(slice_cuts), multiple_matches=int((matches > 1).sum()),
                     unassigned=int((matches == 0).sum()), classifier_disagreements=int((predicted != model_label).sum()),
                     classifier_disagreement_fraction=mismatch,
                     provisional_component_map={c.GetName(): indices[i] for i, c in enumerate(slice_cuts)})
        checks.append(check)
        if check['multiple_matches'] or mismatch > payload.get('maximum_cut_disagreement', .002):
            raise ValueError('Cut approximation failed membership QA: '+json.dumps(check))
        all_cuts.extend(slice_cuts)
        _diagnostics(ROOT, row, data, raw, predicted, slice_cuts, check, arm)
    # Validation above cannot leave a partial canonical cut file behind.
    out = ROOT.TFile.Open(str(cut_path), 'NEW')
    if not out or out.IsZombie():
        raise OSError('Cannot exclusively create '+str(cut_path))
    for cut in all_cuts:
        out.cd()
        cut.Write()
    ROOT.TNamed('hallc_spectrometer', arm).Write()
    ROOT.TNamed('hallc_geometry_profile', arm+'_centered_v1').Write()
    ROOT.TNamed('geometric_selection_manifest', json.dumps(payload)).Write()
    ROOT.TNamed('geometric_cut_membership_checks', json.dumps(checks)).Write()
    out.Close()
    return checks


def _diagnostics(ROOT, row, data, raw, predicted, cuts, check, arm):
    out = ROOT.TFile.Open(str(row['diagnostic_root']), 'NEW')
    if not out or out.IsZombie():
        raise OSError('Cannot create diagnostics '+str(row['diagnostic_root']))
    keep = []
    def histogram(name, title, x, y=None):
        def limits(values):
            lo, hi = float(np.min(values)), float(np.max(values))
            pad = max(1e-8, (hi-lo)*.02)
            return lo-pad, hi+pad
        if y is None:
            hist = ROOT.TH1D(name, title, 240, *limits(x))
            for v in x: hist.Fill(float(v))
        else:
            hist = ROOT.TH2D(name, title, 180, *limits(x), 180, *limits(y))
            for a, b in zip(x, y): hist.Fill(float(a), float(b))
        hist.Write(); keep.append(hist)
        return hist
    h_orig = histogram('hOrig', 'Observed focal plane;fp slope;fp position', raw[:, 0], raw[:, 1])
    h_qp = histogram('hQP', 'Standardized orthogonal rotation;q;p', data['q'], data['p'])
    h_q = histogram('hQ', 'Observed q projection;q;events', data['q'])
    histogram('hP', 'Observed p projection;p;events', data['p'])
    if 'u' in data and np.isfinite(data['u']).any():
        histogram('hU', 'Observed unfolded X coordinate;u;events', data['u'][np.isfinite(data['u'])])
    for index in np.unique(predicted[predicted >= 0]):
        take = predicted == index
        histogram('hBandOrig_%d' % index, 'Selected component;fp slope;fp position', raw[take, 0], raw[take, 1])
        histogram('hBandQP_%d' % index, 'Selected component;q;p', data['q'][take], data['p'][take])
        histogram('hBandQ_%d' % index, 'Selected q;q;events', data['q'][take])
        if 'xsieve' in data and 'ysieve' in data:
            histogram('hBandSieve_%d' % index, 'Preliminary component at sieve;ysieve (cm);xsieve (cm)', data['ysieve'][take], data['xsieve'][take])
    tree = ROOT.TTree('GeometricAssignments', 'Provisional component IDs, before physical sieve relabeling')
    fields = {'entry': array('q', [0]), 'band': array('i', [0]), 'cut_band': array('i', [0])}
    tree.Branch('entry', fields['entry'], 'entry/L')
    tree.Branch('band', fields['band'], 'band/I')
    tree.Branch('cut_band', fields['cut_band'], 'cut_band/I')
    for i in range(len(raw)):
        fields['entry'][0] = int(data['entry'][i]); fields['band'][0] = int(data['model_label'][i]); fields['cut_band'][0] = int(predicted[i]); tree.Fill()
    tree.Write()
    canvas = ROOT.TCanvas('c_geometric_overview', 'Geometric selection overview', 1200, 900)
    canvas.Divide(2, 2); canvas.cd(1); h_orig.Draw('colz')
    for cut in cuts: cut.Draw('L same')
    canvas.cd(2); h_qp.Draw('colz'); canvas.cd(3); h_q.Draw()
    canvas.cd(4)
    note = ROOT.TPaveText(.05, .2, .95, .8, 'NDC')
    note.AddText('PROVISIONAL '+arm+' geometric selection')
    note.AddText('Physical sieve labels assigned in step 3')
    note.AddText('TCutG / classifier disagreements: %d' % check['classifier_disagreements'])
    note.Draw(); canvas.Write(); canvas.Close()
    ROOT.TNamed('hallc_spectrometer', arm).Write()
    ROOT.TNamed('geometric_membership_check', json.dumps(check)).Write()
    ROOT.TNamed('geometric_model', Path(row['model']).read_text()).Write()
    out.Close()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument('--payload', type=Path)
    mode.add_argument('--extract-cuts', type=Path, help='Export a reviewed automatic/edited ROOT slice to manual polygon JSON')
    parser.add_argument('--report', type=Path)
    parser.add_argument('--axis', choices=['x', 'y'])
    parser.add_argument('--foil', type=int)
    parser.add_argument('--ndel', type=int)
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    try:
        if args.extract_cuts:
            if args.axis is None or args.foil is None or args.ndel is None or args.output is None:
                raise ValueError('--extract-cuts requires --axis, --foil, --ndel and --output')
            import ROOT
            fresh(args.output)
            source = ROOT.TFile.Open(str(args.extract_cuts), 'READ')
            if not source or source.IsZombie(): raise OSError('Cannot open '+str(args.extract_cuts))
            pattern = (r'^hXpFpXFp_cut_candidate_(\d+)_zone_global_nfoil_%d_ndel_%d$' if args.axis == 'x'
                       else r'^hYpFpYFp_cut_yscol_(\d+)_nfoil_%d_ndel_%d$') % (args.foil, args.ndel)
            polygons, seen = [], set()
            for key in source.GetListOfKeys():
                match = re.fullmatch(pattern, key.GetName())
                if not match or key.GetName() in seen: continue
                seen.add(key.GetName()); cut = source.Get(key.GetName())
                if not cut.InheritsFrom('TCutG'): continue
                polygons.append(dict(index=int(match.group(1)), points=[[cut.GetPointX(i), cut.GetPointY(i)] for i in range(cut.GetN())]))
            source.Close()
            if not polygons: raise ValueError('No provisional cuts match the requested slice')
            with args.output.open('x') as stream: json.dump(polygons, stream, indent=2)
            print(args.output); return 0
        if args.report is None: raise ValueError('--payload requires --report')
        fresh(args.report)
        checks = write(json.loads(args.payload.read_text()))
        with args.report.open('x') as stream:
            json.dump(checks, stream, indent=2)
        return 0
    except Exception as exc:
        print('ERROR: '+str(exc), file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
