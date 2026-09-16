"""Manual campaign-size statistics benchmark; synthetic residuals, no ROOT input.

Run from the repository root: python3 tests/benchmark_comparison.py
Uses archived hole counts to reproduce event/group sizes. Plots are skipped;
tables and compressed synthetic residuals are written only to a temporary folder.
"""
import argparse
import json
import sys
import tempfile
from pathlib import Path
from time import perf_counter
from unittest.mock import patch

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from core_sample import read_tsv,PROJECT
from comparison_plots import report,MODELS


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('comparison',type=Path,nargs='?',
                        default=PROJECT/'HMS_6p667GeV/07_diagnostics/compare')
    args=parser.parse_args()
    rows=read_tsv(args.comparison/'tsv/residuals.tsv')
    holes=[r for r in rows if r['level']=='hole' and r['model']=='old' and r['target']=='xptar']
    columns={'pool':str,'zfoil':float,'ndel':int,'rungroup':str,'xscol':int,'yscol':int}
    data={k:np.repeat(np.array([cast(r[k]) for r in holes]),[int(r['n']) for r in holes])
          for k,cast in columns.items()}
    blocked=json.loads((args.comparison/'manifest.json').read_text())['counts'].get('blocked',0)
    for k in data:
        data[k]=np.concatenate((data[k],np.full(blocked,'blocked' if k=='pool' else data[k][0],dtype=data[k].dtype)))
    pool=data.pop('pool');data['ztarT']=data.pop('zfoil');n=len(pool)
    data['entry']=np.arange(n)
    rng=np.random.default_rng(667)
    # Shuffle so gathering resembles interleaved replay populations.
    order=rng.permutation(n);pool=pool[order];data={k:v[order] for k,v in data.items()}
    rr={m:rng.normal(size=(n,4)) for m in MODELS}
    known=rng.random(n)>.05;seen=rng.random(n)>.5
    mask=pool==pool[0]
    start=perf_counter()
    for _ in range(20):sum(mask)
    projected=(perf_counter()-start)/20*len(rows)
    print(f'{n:,} events; expected {len(rows):,} detail rows')
    print(f'Legacy repeated-mask counting alone: {projected:.1f}s projected from 20 calls')
    with tempfile.TemporaryDirectory() as tmp:
        out=Path(tmp);(out/'tsv').mkdir()
        start=perf_counter()
        with patch('comparison_plots.plots'):
            report(out,data,pool,rr,known,seen,{'method':'synthetic benchmark'})
        print(f'New report excluding plots: {perf_counter()-start:.1f}s measured')
        actual=read_tsv(out/'tsv/residuals.tsv')
        assert len(actual)==len(rows),(len(actual),len(rows))
        # Check every group, target, model, count and flag against the real layout.
        keys=('pool','level','zfoil','ndel','rungroup','xscol','yscol','model','target','n','qa_low')
        assert [tuple(r[k] for k in keys) for r in actual]==[tuple(r[k] for k in keys) for r in rows]
        print('All group identities, counts, flags and row order match the archived campaign.')


if __name__=='__main__':main()
