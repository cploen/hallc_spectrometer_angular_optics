#!/usr/bin/env python3
"""Export compact exact-holdout tracks for the isolated YTAR coupling calculation."""
import io
import json
from pathlib import Path
import zipfile
import numpy as np
import uproot
from beam_holdout import prepare
from core_sample import PROJECT, digest, read_tsv

CAMPAIGNS = [('HMS_5p878GeV','inner10'), ('HMS_6p117GeV','inner20'),
             ('HMS_6p667GeV_joint_min5','inner20')]


def main():
    output = PROJECT/'collaboration_exports/ytar_coupling_inputs_20261103.zip'
    if output.exists():
        raise FileExistsError(f'Refusing to overwrite {output}')
    output.parent.mkdir(exist_ok=True)
    temporary = output.with_suffix('.partial')
    manifest = dict(schema='ytar_coupling_inputs_v1', campaigns=[],
                    selection='Exact protected cores plus all supported shoulders; unchanged event IDs',
                    units=dict(xfp='cm',yfp='cm',xpfp='rad',ypfp='rad',rx='cm',ry='cm',rz='cm',angle='rad'))
    try:
        with zipfile.ZipFile(temporary, 'w', compression=zipfile.ZIP_STORED) as archive:
            for name,tag in CAMPAIGNS:
                campaign = PROJECT/name
                saved, beam, a, matrices, pool, counts, hashes = prepare(campaign,tag,f'beam20_{tag}')
                fields = ('entry','rungroup','angle','xfp','xpfp','yfp','ypfp','ztarT','ndel')
                packed = {key:a[key] for key in fields}
                for key in ('rx','ry','rz'):
                    packed[key] = np.empty(len(pool))
                source = campaign/'05c_core_sample'/tag
                table, = source.glob('rungroups_*_inputs.tsv')
                settings = {r['rungroup']:r for r in read_tsv(table)}
                for name_rg in np.unique(a['rungroup']):
                    chosen = np.flatnonzero(a['rungroup']==name_rg)
                    path = campaign/'06c_core_ntuple'/tag/'holdout/root'/f'Optics_{settings[name_rg]["optics_id"]}_-1_fit_tree_gmm.root'
                    with uproot.open(path) as root:
                        raw = root['TFit'].arrays(['entry','reactxcalc','reactycalc','ztar'],library='np')
                    order = np.argsort(raw['entry'])
                    if not np.array_equal(raw['entry'][order], a['entry'][chosen]):
                        raise ValueError(f'{name_rg}: event IDs disagree')
                    for destination,original in [('rx','reactxcalc'),('ry','reactycalc'),('rz','ztar')]:
                        packed[destination][chosen] = raw[original][order]
                keep = np.isin(pool,['core','noncore'])
                packed = {key:value[keep] for key,value in packed.items()}
                if not all(np.isfinite(value).all() for key,value in packed.items() if value.dtype.kind in 'fiu'):
                    raise ValueError(f'{name}: nonfinite input')
                buffer=io.BytesIO()
                np.savez_compressed(buffer,**packed)
                archive.writestr(name+'.npz',buffer.getvalue())
                for label,path in [('beam',saved/'matrices/beam.dat'),('old',saved/'seed.dat')]:
                    archive.write(path,f'{name}/{label}.dat')
                manifest['campaigns'].append(dict(campaign=name,tag=tag,n=int(keep.sum()),
                    beam_sha256=digest(saved/'matrices/beam.dat'),old_sha256=digest(saved/'seed.dat'),inputs=hashes))
                print(f'{name}: exported {keep.sum():,} identical protected core + shoulder tracks',flush=True)
            archive.writestr('manifest.json',json.dumps(manifest,indent=2)+'\n')
        temporary.rename(output)
    finally:
        temporary.unlink(missing_ok=True)
    print(f'Wrote {output} ({output.stat().st_size/1024**2:.1f} MiB). No fit or replay was run.')


if __name__=='__main__':
    main()
