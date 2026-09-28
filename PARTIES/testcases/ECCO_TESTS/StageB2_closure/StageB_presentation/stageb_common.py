#!/usr/bin/env python3
"""Shared paths, loaders and plot style for the Stage-B presentation figures.

Every number plotted by the scripts in this directory is read from a Stage-B run
directory on scratch or from a recorded audit JSON in ``ECCO_TESTS/StageB2_closure/B2_closure``.
Values that exist only in the prose of the roadmap / closure report are kept in
``ROADMAP_NUMBERS`` below with the section they come from, so a figure never
invents a number and the provenance of every annotated value is one grep away.

Field arrays in the PARTIES HDF5 output are indexed ``[z, y, x]`` and carry one
duplicated high boundary plane per direction; :func:`physical` strips it, which
is the same convention as ``B2_closure/analyze_fields.py``.
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np

# --------------------------------------------------------------------------
# paths
# --------------------------------------------------------------------------
HERE = Path(__file__).resolve().parent
FIGDIR = HERE / 'figures'
VIDDIR = HERE / 'videos'
DATADIR = HERE / 'data'
B2 = HERE.parent / 'B2_closure'
SCRATCH = Path('/anvil/scratch/x-mjalabert/ECCO_StageB')

for _d in (FIGDIR, VIDDIR, DATADIR):
    _d.mkdir(exist_ok=True)

#: Stage-B runs referenced by the figures.  Key -> scratch directory name.
RUNS = {
    # accepted B.2 coupled gate, final solver (release-ramp fix)
    'coupled16': 'B2closure_coupled16_rampfix_20492853',
    'coupled16_halfdt': 'B2closure_coupled16_rampfix_halfdt_20492854',
    'coupled16_quarterdt': 'B2closure_coupled16_rampfix_quarterdt_20492876',
    # same deck before the ramp fix (kept for the pre-release bit-identity check)
    'coupled16_old': 'B2closure_coupled16_20478206',
    'coupled16_old_halfdt': 'B2closure_coupled16_halfdt_20487690',
    'coupled16_old_quarterdt': 'B2closure_coupled16_quarterdt_20492786',
    # 32-rank twin of coupled16_old: roundoff-perturbation floor
    'coupled16_mpi32': 'B2closure_coupled16_mpi32_20492787',
    # hard ice-penalization threshold disabled (diagnostic only)
    'coupled16_smooth': 'B2closure_coupled16_smooth_20500720',
    'coupled16_smooth_halfdt': 'B2closure_coupled16_smooth_halfdt_20500721',
    # ice-free settling controls
    'settle16': 'B2closure_settle16_20460632',
    'settle24': 'B2closure_settle24_20477406',
    'settle32': 'B2closure_settle32_20477407',
    'settle24_dt10': 'B2closure_settle24_dt10_20500848',
    'settle24_dt05': 'B2closure_settle24_dt05_20500849',
    # coarse d/dx = 8 MPI pair under the ramp fix
    'v9_rampfix': 'B2closure_coupled_v9_rampfix_20492915',
    'v9_rampfix_mpi8': 'B2closure_coupled_v9_rampfix_mpi8_20492916',
    # negative control: warm water, phase change on, rock must not melt
    'hot_water': 'B2closure_hot_water_20477870',
    # closed-box budget seed
    'seed_v9': 'B2closure_seed_v9_20477869',
    # the first long 3-D shakeout (superseded physics, retained history)
    'release3d': 'StageB_release3d_20397146',
}


def run(key: str) -> Path:
    """Scratch directory for a named Stage-B run."""
    return SCRATCH / RUNS[key]


# --------------------------------------------------------------------------
# numbers quoted in the roadmap / closure report (prose-only, with provenance)
# --------------------------------------------------------------------------
ROADMAP_NUMBERS = {
    # B.1.1 EOS polar re-fit, roadmap  "B.1.1 Re-fit the EOS to polar conditions"
    'eos_yang': dict(Cb=0.011, T0=4.0, cS=-0.25, b0=0.77, rms=1.10053, maxerr=1.39148),
    'eos_polar': dict(Cb=0.006935, T0=3.6966, cS=-0.2151, b0=0.8116,
                      rms=0.00143, maxerr=0.00452),
    # B.1.1a  density-anomaly error vs TEOS-10 outside the fit box (kg/m^3),
    # each model referenced to S = 34 -- roadmap table "Error in the density
    # anomaly against TEOS-10"
    'dilution': dict(
        S=[31.5, 28, 24, 20, 10, 1],
        quadratic=[0.003, -0.000, -0.015, -0.040, -0.148, -0.267],
        linear=[0.006, 0.009, 0.008, 0.002, -0.033, -0.055]),
    # B.1.1a  rms/std of each model INSIDE the ambient box (per cent)
    'eos_inbox': {'linear  c0 + c1θ + c2s': 0.222,
                  'Roquet quadratic, vertex pinned': 0.088,
                  'full quadratic (ceiling)': 0.001},
    # B.1.3-fn  F_release calibration table
    'F_release_cal': {0.5: 0.565, 1.0: 0.709, 1.5: 0.867, 1.67: 0.900},
    'phi_liq_ceiling': 0.977,
    # B.1.4 / B.1.4-fix  meltwater-tracer budget error, 2-D release deck
    'tracer_leak': {2: -7.20e-2, 10: -6.75e-3, 15: 2.29e-2, 30: 9.30e-2},
    'tracer_noparticle': {2: -1.08e-7, 10: -4.77e-7, 15: -5.35e-7},
    # B.0 non-invasiveness, part 2
    'noninvasive': {2: 1.1e-12, 6: 7.0e-11, 10: 9.2e-9},
    'lineage_bar': {2: 1.7e-12, 6: 5.8e-5},
    'discriminating_size': 7.2e-3,
    # closure report, "The cause: the ice escape is an ill-conditioned residual"
    'escape': dict(
        buoyant_weight=0.785398,
        dt=[0.01, 0.005, 0.0025],
        F_IBM=[0.76887, 0.77784, 0.78457],
        net=[0.01653, 0.00756, 0.00083],
        tau_y330=[5.75, 6.63, 7.42],
        observed_order=0.15),
    # closure report, "The hard threshold costs a factor 8 in settling accuracy"
    'dt_convergence': {
        'ice-free, threshold never fires': 0.28,
        'coupled, threshold disabled': 0.41,
        'coupled, threshold active': 3.29},
    # closure report, B.2 acceptance table (job 20492853)
    'gates': [
        ('1. Locked grain stationary',
         'max |ΔX| = 0, max |U| = 0 for every locked step', 'exact', 'pass'),
        ('2. Impulse-free release',
         'peak |a| = 0.324 a_free  at t = 33.728', '< 2 a_free', 'pass'),
        ('3. Meltwater injected and conserved',
         'max |∫C_mw − ∫dF| / ∫dF = 7.31e-7', '< 1e-3', 'pass'),
        ('4. Settling and wall collision',
         'min U_y = −0.671, wall at t = 44.65, rest at y = 0.4997', 'reaches wall', 'pass'),
        ('5. Heat / salt / tracer budgets',
         'heat 1.01e-6, salt 4.45e-16, tracer 7.31e-7', 'machine / 1e-3', 'pass'),
        ('6. Ice rigidity (Darcy τ = 1e-3)',
         'worst bulk-ice / domain speed = 0.583 %', '< 1 %', 'pass'),
    ],
    # closure report, "MPI, restart, diagnostics and default-off regression"
    'reproducibility': [
        ('Stage-A no-impact gate', 0.0, 'λ identical to 16 digits'),
        ('Periodic remap, 8 vs 16 ranks', 0.0, 'all phase fields'),
        ('64 vs 32 ranks, full fall', 0.0, 'every sampled height'),
        ('Pre-release fields, ramp fix', 0.0, '17 snapshots, max diff 0'),
        ('Clean rebuild vs candidate', 0.0, '3 seed snapshots'),
        ('Restart vs continuous', 0.0, 'fluid, phase, scalars, particle'),
        ('Coupled MPI d/Δx = 8', 8.96e-6, 'max |Δy| / d before contact'),
        ('Profile reconstruction', 6.99e-15, 'scaled error, 32 snapshots'),
        ('Minimum-image $C_S$', 9.98e-13, 'vs an independent formula'),
    ],
    # roadmap, "RESOURCE STATUS AND FORWARD REQUIREMENT"
    'allocation': dict(limit=386853, used=309843, remaining=77010),
    'forward': [
        ('Finish Stage B only', 8000, 8000),
        ('+ 1 baseline production run', 17488, 36463),
        ('+ 3-run reduced sweep', 36463, 93389),
        ('+ 7-run one-at-a-time sweep', 74414, 207241),
        ('+ 7-run sweep and 1 hero run', 226216, 662648),
    ],
    'request': dict(minimum=110000, recommended=300000),
}

#: deck constants of the accepted B.2 coupled gate (read back from parties.inp)
COUPLED16 = dict(d=1.0, R=0.5, rho_s=2.5, g=1.0, y_ice0=3.6, y_grain0=4.1,
                 Lx=4.0, Ly=6.0, Lz=4.0, nx=64, ny=96, nz=64)

# --------------------------------------------------------------------------
# palette  (dataviz reference instance, light surface)
# --------------------------------------------------------------------------
SURFACE = '#fcfcfb'
PAGE = '#f9f9f7'
INK = '#0b0b0b'
INK2 = '#52514e'
MUTED = '#898781'
GRID = '#e1e0d9'
AXIS = '#c3c2b7'

SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
          '#e87ba4', '#008300', '#4a3aa7', '#e34948']
BLUE, ORANGE, AQUA, YELLOW, MAGENTA, GREEN, VIOLET, RED = SERIES

STATUS = dict(good='#0ca30c', warning='#fab219', serious='#ec835a', critical='#d03b3b')

#: single-hue sequential ramps (blue = default, orange = second context)
BLUE_RAMP = ['#cde2fb', '#b7d3f6', '#9ec5f4', '#86b6ef', '#6da7ec', '#5598e7',
             '#3987e5', '#2a78d6', '#256abf', '#1c5cab', '#184f95', '#104281', '#0d366b']
ORANGE_RAMP = ['#fde3d5', '#fbcdb5', '#f8b794', '#f5a074', '#f18a56', '#eb6834',
               '#d3562a', '#b64620', '#963718', '#752911', '#571d0b']


def _cmap(name, colors):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list(name, colors)


def cmaps():
    """Named sequential / diverging colormaps built from the palette ramps."""
    return dict(
        blue=_cmap('sb_blue', BLUE_RAMP),
        orange=_cmap('sb_orange', ORANGE_RAMP),
        diverging=_cmap('sb_div', ['#0d366b', '#2a78d6', '#9ec5f4',
                                   '#f0efec', '#f5a074', '#d03b3b', '#7a1d1d']),
    )


def use_style():
    """Apply the presentation style to matplotlib globally."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        'figure.facecolor': PAGE,
        'savefig.facecolor': PAGE,
        'axes.facecolor': SURFACE,
        'axes.edgecolor': AXIS,
        'axes.linewidth': 0.8,
        'axes.labelcolor': INK2,
        'axes.titlecolor': INK,
        'axes.titlesize': 11,
        'axes.titleweight': 'semibold',
        'axes.titlelocation': 'left',
        'axes.titlepad': 8,
        'axes.labelsize': 9.5,
        'axes.grid': True,
        'axes.axisbelow': True,
        'grid.color': GRID,
        'grid.linewidth': 0.7,
        'xtick.color': MUTED,
        'ytick.color': MUTED,
        'xtick.labelsize': 8.5,
        'ytick.labelsize': 8.5,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'text.color': INK,
        'font.family': 'sans-serif',
        'font.sans-serif': ['DejaVu Sans'],
        'font.size': 9.5,
        'legend.frameon': False,
        'legend.fontsize': 8.5,
        'legend.labelcolor': INK2,
        'lines.linewidth': 2.0,
        'lines.markersize': 4.5,
        'lines.solid_capstyle': 'round',
        'figure.dpi': 130,
        'savefig.dpi': 200,
        'savefig.bbox': 'tight',
        'savefig.pad_inches': 0.25,
        'mathtext.default': 'regular',
    })
    return plt


def use_ffmpeg():
    """Point matplotlib at a usable ffmpeg (Anvil has none on PATH by default).

    Order: $FFMPEG_BINARY, then the binary bundled with imageio-ffmpeg, then
    whatever is already on PATH (e.g. after `module load ffmpeg/4.2.2`).
    """
    import os
    import shutil
    import matplotlib
    cand = os.environ.get('FFMPEG_BINARY')
    if not cand:
        try:
            import imageio_ffmpeg
            cand = imageio_ffmpeg.get_ffmpeg_exe()
        except Exception:
            cand = shutil.which('ffmpeg')
    if cand:
        matplotlib.rcParams['animation.ffmpeg_path'] = cand
    from matplotlib.animation import FFMpegWriter
    if not FFMpegWriter.isAvailable():
        raise RuntimeError(
            'no usable ffmpeg: set $FFMPEG_BINARY, pip install imageio-ffmpeg, '
            'or module load ffmpeg/4.2.2')
    return cand


def finish(ax, spines=('top', 'right')):
    """Recede the chrome: drop the named spines, keep a hairline baseline."""
    for s in spines:
        ax.spines[s].set_visible(False)
    return ax


def save(fig, stem, caption=None):
    """Write <stem>.png and <stem>.pdf into figures/, print the paths.

    The caption is hard-wrapped rather than relying on ``wrap=True``: with
    ``savefig.bbox = 'tight'`` an unwrapped line silently widens the whole
    figure, which is how a 12.6-inch figure ends up 16 inches wide.
    """
    if caption:
        import textwrap
        width = int(fig.get_size_inches()[0] * 15)
        fig.text(0.004, 0.004, textwrap.fill(caption, width), color=MUTED,
                 fontsize=7.2, ha='left', va='bottom', linespacing=1.45)
    for ext in ('png', 'pdf'):
        out = FIGDIR / f'{stem}.{ext}'
        fig.savefig(out)
    print(f'wrote {FIGDIR/(stem + ".png")}  and .pdf')


# --------------------------------------------------------------------------
# loaders
# --------------------------------------------------------------------------
def deck(rundir) -> dict:
    """Parse a parties.inp into {key: [floats]} (same parser as analyze_profiles)."""
    out = {}
    for line in (Path(rundir) / 'parties.inp').read_text().splitlines():
        line = line.split('#')[0].split('//')[0].strip()
        if '=' not in line:
            continue
        key, val = line.split('=', 1)
        val = val.strip().strip('{}')
        try:
            out[key.strip()] = [float(x) for x in val.split(',')]
        except ValueError:
            pass
    return out


def mobile(rundir):
    """particle history: t, X(N,3), U(N,3) from mobile.dat."""
    a = np.loadtxt(Path(rundir) / 'mobile.dat', delimiter=',')
    return a[:, 0], a[:, 2:5], a[:, 5:8]


def release(rundir):
    """release.dat as a structured array (t, phi_liq, t_released, ramp, ...)."""
    return np.genfromtxt(Path(rundir) / 'release.dat', delimiter=',', names=True)


def release_time(rundir):
    """Recorded release instant, or None if the grain never released."""
    r = release(rundir)
    got = r[r['t_released'] >= 0]
    return float(got['t_released'][0]) if len(got) else None


def profiles(rundir):
    """ecco_profiles.csv as a structured array."""
    return np.genfromtxt(Path(rundir) / 'ecco_profiles.csv', delimiter=',', names=True)


def profile_blocks(rundir):
    """Split ecco_profiles.csv into (times, y, {field: array[time, y]})."""
    a = profiles(rundir)
    ny = len(np.unique(a['y']))
    times = np.unique(a['time'])
    y = a['y'][:ny]
    fields = {}
    for name in a.dtype.names:
        if name in ('time', 'step', 'y'):
            continue
        m = np.full((len(times), ny), np.nan)
        for i, t in enumerate(times):
            b = a[a['time'] == t][-ny:]
            if len(b) == ny:
                m[i] = b[name]
        fields[name] = m
    return times, y, fields


def snapshots(rundir):
    """Sorted list of (index, Data_*.h5 path) for a run."""
    d = Path(rundir)
    files = sorted(d.glob('Data_*.h5'),
                   key=lambda p: int(re.findall(r'_(\d+)\.h5$', p.name)[0]))
    return [(int(re.findall(r'_(\d+)\.h5$', p.name)[0]), p) for p in files]


def physical(a):
    """Strip the duplicated high boundary plane from a [z, y, x] field."""
    return a[:-1, :-1, :-1]


def load_fields(path, names=('VOF/C_L', 'VOF/C_S', 'Conc/0', 'Conc/1', 'Conc/2',
                             'u', 'v', 'w')):
    """Read named datasets plus the cell-centre grid and time from a Data_*.h5."""
    import h5py
    out = {}
    with h5py.File(path, 'r') as f:
        out['t'] = float(f['time'][0])
        for n in names:
            out[n] = physical(f[n][...])
        out['x'] = f['grid/xc'][...].ravel()[:-1]
        out['y'] = f['grid/yc'][...].ravel()[:-1]
        out['z'] = f['grid/zc'][...].ravel()[:-1]
    out['ice'] = np.clip(1.0 - out['VOF/C_L'] - out['VOF/C_S'], 0.0, 1.0)
    return out


def particle_state(path):
    """(X, U, R, t_released) from a Particle_*.h5."""
    import h5py
    with h5py.File(path, 'r') as p:
        X = p['mobile/X'][...].reshape(-1, 3)[0]
        U = p['mobile/U'][...].reshape(-1, 3)[0]
        R = float(p['mobile/R'][...].ravel()[0])
        tr = float(p['mobile/t_released'][...].ravel()[0])
    return X, U, R, tr


def audit(name):
    """Load a recorded JSON audit from ECCO_TESTS/StageB2_closure/B2_closure."""
    return json.loads((B2 / name).read_text())


def campaign_jobs():
    """Slurm accounting for the closure campaign (B2_closure/campaign_jobs.json)."""
    return audit('campaign_jobs.json')
