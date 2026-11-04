#!/usr/bin/env python3
"""Bundle the zero-constant YTAR candidate and matched sieve views for slide updates."""
import csv
import hashlib
import json
from pathlib import Path
import subprocess
import sys
import zipfile

PROJECT=Path(__file__).resolve().parent
CAMPAIGNS=[('HMS_5p878GeV','inner10',16.500),('HMS_6p117GeV','inner20',12.385),
           ('HMS_6p667GeV_joint_min5','inner20',12.490)]


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def table(path):
    with path.open(newline='') as stream:return list(csv.DictReader(stream,delimiter='\t'))


def matrix(path):
    result={}
    for line in path.read_text().splitlines():
        fields=line.split()
        if len(fields)==5 and len(fields[4])==5 and fields[4].isdigit():
            if fields[4] in result:raise ValueError(f'Duplicate term in {path}')
            result[fields[4]]=[float(v.replace('D','E').replace('d','e')) for v in fields[:4]]
    if not result:raise ValueError(f'No matrix rows: {path}')
    return result


def main():
    output=PROJECT/'collaboration_exports/nps_ytar_zero_slide_update_20261104.zip'
    if output.exists():raise FileExistsError(f'Export already exists: {output}')
    files={};items=[]
    def add(path):
        if not path.is_file():raise FileNotFoundError(path)
        files[str(path.relative_to(PROJECT))]=path
    for name,tag,angle in CAMPAIGNS:
        c=PROJECT/name
        refit=c/'07_diagnostics/ytar_zero_refit_20261103'
        candidate=refit/'ytar_zero.dat'
        new=matrix(candidate)
        if new.get('00000',[None,None])[1]!=0.:raise ValueError(f'{name}: YTAR constant is not exactly zero')
        add(candidate)
        for path in refit.rglob('*'):
            if path.is_file() and path.suffix.lower() in ('.tsv','.json','.md','.txt'):add(path)
        for report in ('terms.tsv','manifest.json','metrics.tsv'):
            if not any(path.name==report for path in refit.rglob(report)):
                raise FileNotFoundError(f'{refit}: missing {report}')
        original=c/'06e_beam_search'/f'beam20_{tag}'/'matrices/beam.dat'
        add(original)
        old=matrix(original)
        angular_same=all(all(new.get(term,[0.]*4)[j]==old.get(term,[0.]*4)[j] for j in (0,2)) for term in set(new)|set(old))
        view=c/'07_diagnostics/ytar_zero_sieve_review_20261104'
        status=json.loads((view/'manifest.json').read_text())
        if status.get('status')!='complete':raise ValueError(f'{view}: sieve generation has not completed')
        expected=dict(angle=angle,delta=[-10,10],cer=6.,cal=.65,foil_width=2.,matrix_offsets=[0.,0.,0.])
        def matches(saved):return saved.get('status')=='complete' and all(saved['settings'].get(k)==v for k,v in expected.items())
        if not matches(status):
            # The source command omitted PID/window flags; preserve its output and make the slide-matched view.
            view=c/'07_diagnostics/ytar_zero_sieve_slides_20261104'
            if not view.exists():
                print(f'{name}: generating a view with the original slide cuts (cer>6, cal>0.65, foil width 2 cm)',flush=True)
                subprocess.run([sys.executable,str(PROJECT/'sieve_afterburner.py'),str(c),view.name,
                    '--matrix',str(candidate),'--offsets','0','0','0','--angle',str(angle),
                    '--delta','-10','10','--cer','6','--cal','0.65','--foil-width','2'],cwd=PROJECT,check=True)
            status=json.loads((view/'manifest.json').read_text())
        if not matches(status):raise ValueError(f'{view}: settings differ from the original slides')
        if sha(view/'matrices/matrix.dat')!=sha(candidate):raise ValueError(f'{view}: matrix differs from the constrained refit')
        for relative in ('manifest.json','histograms.npz','tsv/counts.tsv','tsv/closure.tsv','matrices/matrix.dat'):add(view/relative)
        reference=c/'07_diagnostics'/f'sieve_2024flags_theta{angle:.3f}_all_delta'.replace('.','p')
        baseline=json.loads((reference/'manifest.json').read_text())
        if status['replays']!=baseline['replays']:raise ValueError(f'{name}: replay inputs differ from the saved 2024 view')
        for key in ('angle','delta','cer','cal','foil_width','xmis'):
            if status['settings'].get(key)!=baseline['settings'].get(key):raise ValueError(f'{name}: {key} differs from the saved 2024 view')
        counts=lambda path:{(r['rungroup'],r['zfoil']):int(r['n']) for r in table(path/'tsv/counts.tsv') if r['model']=='matrix'}
        if counts(view)!=counts(reference):raise ValueError(f'{name}: event counts differ from the saved 2024 view')
        add(reference/'manifest.json');add(reference/'tsv/counts.tsv')
        add(c/'config/oldfit.dat')
        items.append(dict(campaign=name,tag=tag,angle=angle,refit=str(refit.relative_to(PROJECT)),
                          sieve=str(view.relative_to(PROJECT)),angular_coefficients_unchanged=angular_same))
        print(f'{name}: zero YTAR constant and matched sieve events verified; angular coefficients unchanged: {angular_same}',flush=True)
    output.parent.mkdir(exist_ok=True)
    metadata=dict(schema='ytar_zero_slide_update_v1',campaigns=items,
                  files={name:dict(sha256=sha(path),bytes=path.stat().st_size) for name,path in files.items()})
    temporary=output.with_suffix('.partial')
    try:
        with zipfile.ZipFile(temporary,'w',compression=zipfile.ZIP_DEFLATED) as archive:
            archive.writestr('export_manifest.json',json.dumps(metadata,indent=2)+'\n')
            for name,path in sorted(files.items()):archive.write(path,name)
        temporary.rename(output)
    finally:temporary.unlink(missing_ok=True)
    print(f'Exported {output} ({output.stat().st_size/1024**2:.1f} MiB). No fits were run.')


if __name__=='__main__':
    try:main()
    except (OSError,ValueError,KeyError,subprocess.CalledProcessError) as exc:raise SystemExit(f'ERROR: {exc}')
