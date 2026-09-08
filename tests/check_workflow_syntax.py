"""Load each ROOT macro in a fresh interpreter; check Python and shell syntax."""
import ast
from concurrent.futures import ThreadPoolExecutor
import json
from pathlib import Path
import subprocess

REPO=Path(__file__).resolve().parents[1]

def load(path):
    r=subprocess.run(['root','-l','-b','-q','-e',f'gSystem->Exit(gROOT->LoadMacro({json.dumps(str(path))})<0 ? 1 : 0);'],cwd=REPO,capture_output=True,text=True,timeout=90)
    assert r.returncode==0 and 'error:' not in r.stderr.lower(),str(path)+'\n'+r.stdout+r.stderr
    return path.relative_to(REPO)

if __name__=='__main__':
    macros=[*REPO.glob('*.C'),*(REPO/'diagnostics').rglob('*.C')]
    with ThreadPoolExecutor(max_workers=3) as pool:
        for path in pool.map(load,macros):print('PASS',path)
    for path in [*REPO.glob('*.py'),*(REPO/'diagnostics').rglob('*.py'),*(REPO/'tests').glob('*.py')]:ast.parse(path.read_text(),filename=str(path))
    for path in [*REPO.glob('*.sh'),*(REPO/'diagnostics').glob('*.sh')]:
        subprocess.run(['bash','-n',str(path)],check=True,cwd=REPO)
    print(f'PASS: {len(macros)} ROOT macro loads; Python and shell syntax')
