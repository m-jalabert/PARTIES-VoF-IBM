#!/usr/bin/env python3
"""Build the checked-in solver in isolation; never submits a simulation."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

HERE = Path(__file__).resolve().parent
# Solver root = nearest ancestor holding src/ and the Makefile (robust to folder moves).
SOLVER = next(p for p in HERE.parents if (p/'src').is_dir() and (p/'Makefile').is_file())

p = argparse.ArgumentParser()
p.add_argument('--configuration', choices=['sediment', 'default', 'remap', 'smooth'], default='sediment')
p.add_argument('--destination', type=Path, required=True,
               help='New directory: existing directories are rejected')
p.add_argument('--jobs', type=int, default=4)
args = p.parse_args()
dest = args.destination.resolve()
dest.mkdir(parents=True, exist_ok=False)
shutil.copytree(SOLVER/'src', dest/'src', ignore=shutil.ignore_patterns('*.o', '*.d', '*.bak'))
for name in ['Makefile', 'make.def.generic']:
    shutil.copy2(SOLVER/name, dest/name)
if args.configuration in ['sediment', 'remap']:
    shutil.copy2(HERE/'Boundary.validation.h', dest/'src/Include/Boundary.h')
# 'smooth' is a DIAGNOSTIC build: the validated sediment configuration with the
# hard ice-penalization threshold disabled.  It exists to test whether that
# discontinuity is what destroys release-escape timestep convergence.  It is
# NOT a validated configuration and must not be used for a physical scenario.
if args.configuration == 'smooth':
    shutil.copy2(HERE/'Boundary.smooth_penalization.h', dest/'src/Include/Boundary.h')
if args.configuration == 'remap':
    source = dest/'src/Eulerian/VOF_DIFFUSE.c'
    text = source.read_text()
    start = text.index('void VOF_DIFFUSE_init(Cart3d_bag *db)')
    end = text.index('\nvoid VOF_DIFFUSE_set_boundary_values', start)
    closing = text.rfind('\n}', start, end)
    assert closing > start
    text = text[:closing] + '\n' + (HERE/'remap_fixture.inc').read_text() + text[closing:]
    source.write_text(text)
manifest = {'configuration': args.configuration, 'source_sha256': {}}
for file in sorted(dest.rglob('*')):
    if file.is_file():
        manifest['source_sha256'][str(file.relative_to(dest))] = hashlib.sha256(file.read_bytes()).hexdigest()
with (dest/'build.log').open('w') as log:
    subprocess.run(['make', '-C', str(dest), '-j'+str(args.jobs), 'parties'],
                   stdout=log, stderr=subprocess.STDOUT, check=True)
manifest['binary_sha256'] = hashlib.sha256((dest/'parties').read_bytes()).hexdigest()
(dest/'manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')
print(dest/'parties')
