#!/usr/bin/env python3
"""Regenerate every Stage-B presentation figure, and optionally the videos.

    python make_all.py                 # all twelve figures
    python make_all.py --videos        # figures and both videos (~4 min)
    python make_all.py --only fig05    # one figure

Each script is run as a separate process so a failure in one does not take the
rest down, and so the matplotlib state of one never leaks into another.
"""
import argparse
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

FIGURES = [
    ('fig01_stageb_map.py', 'Stage B at a glance: blocks, gates, what B.3 inherits'),
    ('fig02_configuration.py', 'the box, the initial state, the deck'),
    ('fig03_eos_polar.py', 'polar EOS re-fit and the haline scaling (needs gsw)'),
    ('fig04_release_mechanism.py', 'hold-in-ice release, and the clip leak it exposed'),
    ('fig05_lifecycle.py', 'lock, release, settle, land — the motion gates'),
    ('fig06_conservation.py', 'budgets, ice rigidity, negative controls'),
    ('fig07_reproducibility.py', 'bit-identity ladder and the non-invasiveness contract'),
    ('fig08_conditioning.py', 'the escape timing does not converge, and why'),
    ('fig09_resolution.py', 'settling refinement and what is NOT claimed'),
    ('fig10_fields.py', 'mid-plane montage of temperature and meltwater'),
    ('fig11_mixing.py', 'profiles, plume, and one flag for B.3'),
    ('fig12_cost.py', 'what Stage B cost and what B.3 needs'),
]

VIDEOS = [
    ('video1_3d_lifecycle.py', ['--sub', '4', '--fps', '12', '--dpi', '130'],
     '3-D melt front + grain, with the mid-plane meltwater anomaly'),
    ('video2_slices.py', ['--sub', '4', '--fps', '12', '--dpi', '130'],
     'four fields on the plane through the grain'),
]


def run(script, extra=()):
    t0 = time.time()
    print(f'--- {script}', flush=True)
    r = subprocess.run([sys.executable, str(HERE / script), *extra],
                       cwd=HERE, capture_output=True, text=True)
    out = (r.stdout + r.stderr).strip()
    print(f'    {out.splitlines()[-1] if out else "(no output)"}'
          f'   [{time.time()-t0:.0f} s]', flush=True)
    if r.returncode:
        print(out, file=sys.stderr)
    return r.returncode


if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument('--videos', action='store_true', help='also render the movies')
    ap.add_argument('--only', help='substring of one script name')
    args = ap.parse_args()

    todo = [(s, ()) for s, _ in FIGURES]
    if args.videos:
        todo += [(s, a) for s, a, _ in VIDEOS]
    if args.only:
        todo = [(s, a) for s, a in todo if args.only in s]
        if not todo:
            sys.exit(f'nothing matches {args.only!r}')

    fails = [s for s, a in todo if run(s, a)]
    print()
    print(f'{len(todo) - len(fails)} of {len(todo)} succeeded')
    sys.exit(1 if fails else 0)
