#!/usr/bin/env python3
"""Analyze a PARTIES melting Rayleigh-Benard run against Favier et al. (2019).

Benchmark: Favier, Purseed & Duchemin, JFM 858 (2019), case D:
Ra = 1e7, Pr = 1, thetaM = 0.05, lambda = 6, h0 = 0.05, horizontal periodic,
no-slip top/bottom. Their units: thermal-diffusion time on domain height.
PARTIES runs in free-fall units: t_diff = t_ff / Pe_T, hdot_diff = hdot_ff * Pe_T.

Quantitative targets (paper Fig. 7; eq. 48-51 carry the (1/2+St) storage factor):
  1. diffusive phase: h(t) ~ sqrt(h0^2 + 4*lam^2*t_diff), lam solving the
     one-phase Neumann relation lam*e^{lam^2}*erf(lam) = stefan*(1-thetaM)/sqrt(pi)
  2. onset of convection when Ra_e = Ra*(1-thetaM)*h^3 ~ 1708  => h_c ~ 0.0564
  3. post-onset: h ~ t with melting velocity (1/2+St)*hdot ~ gamma*Ra^(1/3)*(1-thetaM)^(4/3),
     gamma ~ 0.115 (beta = 1/3)  => hdot_diff ~ 23.14/(0.5+St) for Ra = 1e7.

Usage: analyze_melting_rb.py --run LABEL=DIR [--run ...] --output-dir OUT
"""
import argparse, csv, glob, re, math
from pathlib import Path
import h5py
import numpy as np

RA, PR, THETAM, H0 = 1.0e7, 1.0, 0.05, 0.05
PE = np.sqrt(RA / PR) * PR          # = Re*Pr = 3162.28
GAMMA, RAC = 0.115, 1707.76
HC = (RAC / (RA * (1 - THETAM))) ** (1.0 / 3.0)


def neumann_lambda(st_eff):
    """Solve lam*exp(lam^2)*erf(lam) = st_eff/sqrt(pi) by bisection."""
    lo, hi = 1e-8, 3.0
    f = lambda l: l * math.exp(l * l) * math.erf(l) - st_eff / math.sqrt(math.pi)
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if f(mid) > 0 else (mid, hi)
    return 0.5 * (lo + hi)


def snapshots(d):
    fs = glob.glob(str(Path(d) / "Data_*.h5"))
    return sorted(fs, key=lambda s: int(re.search(r"Data_(\d+)\.h5", s).group(1)))


def read_snap(fn, stefan):
    with h5py.File(fn, "r") as f:
        t = float(f["time"][0])
        ny = int(f["grid/NY"][0]) - 1
        nx = int(f["grid/NX"][0]) - 1
        F = np.array(f["VOF/C_L"][0])[:ny, :nx]
        th = np.array(f["Conc/0"][0])[:ny, :nx]
        u = np.array(f["u"][0]); v = np.array(f["v"][0])
        yc = np.array(f["grid/yc"][:ny])
    dy = yc[1] - yc[0]
    hbar = F.mean() * (ny * dy)                     # mean fluid height
    # bottom-wall heat flux, 2nd-order one-sided with Dirichlet theta=1 at y=0
    q = -(-8.0 * 1.0 + 9.0 * th[0] - th[1]).mean() / (3.0 * dy)
    nu = q * hbar / (1.0 - THETAM)
    return dict(t=t, hbar=hbar, qw=q, nu=nu,
                umax=float(max(np.abs(u).max(), np.abs(v).max())),
                Fmin=float(F.min()), Fmax=float(F.max()),
                finite=int(np.isfinite(F).all() and np.isfinite(th).all()))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", action="append", required=True)
    ap.add_argument("--stefan", type=float, default=1.0,
                    help="Favier Stefan number St = L/(cp dT) of the run(s)")
    ap.add_argument("--output-dir", required=True)
    a = ap.parse_args()
    out = Path(a.output_dir); out.mkdir(parents=True, exist_ok=True)

    for spec in a.run:
        label, d = spec.split("=", 1)
        rows = [read_snap(fn, a.stefan) for fn in snapshots(d)]
        t_ff = np.array([r["t"] for r in rows])
        h = np.array([r["hbar"] for r in rows])
        t_d = t_ff / PE
        # diffusive reference & onset (one-phase Neumann asymptote from h0)
        lam = neumann_lambda((1.0 / a.stefan) * (1 - THETAM))
        h_diff = np.sqrt(H0**2 + 4 * lam**2 * t_d)
        i_on = np.argmax(h > HC) if (h > HC).any() else None
        # post-onset linear fit (h in [0.15, 0.6])
        m = (h > 0.15) & (h < 0.6)
        hdot_d = np.nan
        if m.sum() > 3:
            hdot_d = np.polyfit(t_d[m], h[m], 1)[0]
        target = GAMMA * RA ** (1 / 3.0) * (1 - THETAM) ** (4 / 3.0) / (0.5 + a.stefan)
        with open(out / f"melting_rb_{label}_timeseries.csv", "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["t_ff", "t_diff", "hbar", "h_diffusive_ref", "qw", "nu", "umax"])
            for r, hd in zip(rows, h_diff):
                w.writerow([r["t"], r["t"] / PE, r["hbar"], hd, r["qw"], r["nu"], r["umax"]])
        print(f"[{label}] snapshots={len(rows)} t_ff={t_ff[-1]:.1f} h={h[-1]:.4f}")
        print(f"  diffusive-phase max |h/h_ref - 1| (h<h_c): "
              f"{np.nanmax(np.abs(h[h < HC] / h_diff[h < HC] - 1)) if (h < HC).any() else np.nan:.3e}")
        print(f"  onset height h_c target {HC:.4f}; first h>h_c at t_ff="
              f"{t_ff[i_on] if i_on is not None else np.nan}")
        print(f"  post-onset hdot (diff units) = {hdot_d:.2f}  "
              f"[Favier estimate gamma*Ra^1/3*(1-thetaM)^4/3/St = {target:.2f}]")
        print(f"  Fmin={min(r['Fmin'] for r in rows):.2e} "
              f"1-Fmax={1 - max(r['Fmax'] for r in rows):.2e} "
              f"finite={all(r['finite'] for r in rows)}")


if __name__ == "__main__":
    main()
