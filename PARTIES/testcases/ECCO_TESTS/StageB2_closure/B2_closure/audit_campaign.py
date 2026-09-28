#!/usr/bin/env python3
"""Record Slurm accounting and immutable scratch artifacts for B2-named jobs."""
import csv
from datetime import datetime, timezone
import io
import json
from pathlib import Path
import subprocess

HERE = Path(__file__).resolve().parent
SCRATCH = Path('/anvil/scratch/x-mjalabert/ECCO_StageB')
raw = subprocess.check_output([
    'sacct', '-u', 'x-mjalabert', '-S', '2026-09-06', '-X', '--parsable2',
    '-o', 'JobID,JobName%40,State,ElapsedRaw,AllocCPUS,ExitCode,Submit,Start,End'
], text=True)
rows = [r for r in csv.DictReader(io.StringIO(raw), delimiter='|')
        if r['JobName'].startswith('B2')]
for row in rows:
    row['cpu_hours'] = int(row['AllocCPUS'])*int(row['ElapsedRaw'])/3600
    runs = list(SCRATCH.glob('B2closure_*_'+row['JobID']))
    row['scratch_runs'] = [str(p) for p in runs]
    row['binary_sha256'] = [(p/'binary.sha256').read_text().split()[0]
                            for p in runs if (p/'binary.sha256').exists()]
out = dict(recorded_utc=datetime.now(timezone.utc).isoformat(),
           actual_cpu_hours=sum(r['cpu_hours'] for r in rows), jobs=rows)
(HERE/'campaign_jobs.json').write_text(json.dumps(out, indent=2)+'\n')
print(json.dumps({k:v for k,v in out.items() if k != 'jobs'}, indent=2))
