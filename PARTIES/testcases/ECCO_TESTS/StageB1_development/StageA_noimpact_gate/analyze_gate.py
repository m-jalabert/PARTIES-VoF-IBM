#!/usr/bin/env python3
"""Analyze the coarse Le=100 legacy-vs-Yang salt-transport gate."""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re

import h5py
import numpy as np


DATA_RE = re.compile(r"Data_(\d+)\.h5$")
PE_T = 100.0
PE_S = 10000.0
LE = PE_S / PE_T
STEFAN = 0.5
LIQUIDUS = 0.5
S_INF = 1.0
THETA_INF = 1.0
THETA_ICE = 0.0
X0 = 0.75
DELTA = 1.0e-6


def reference_state() -> tuple[float, float, float]:
    def values(lam: float) -> tuple[float, float, float]:
        lam_s = lam * math.sqrt(LE)
        denom = 1.0 + math.sqrt(math.pi) * lam_s * math.exp(lam_s**2) * (
            1.0 + math.erf(lam_s)
        )
        s_i = S_INF / denom
        theta_i = -LIQUIDUS * s_i
        left = lam * math.sqrt(math.pi) * math.exp(lam**2) / STEFAN
        right = (THETA_INF - theta_i) / (1.0 + math.erf(lam))
        right += (THETA_ICE - theta_i) / math.erfc(lam)
        return left - right, s_i, theta_i

    lo, hi = 1.0e-8, 1.0
    for _ in range(100):
        mid = 0.5 * (lo + hi)
        if values(lo)[0] * values(mid)[0] <= 0.0:
            hi = mid
        else:
            lo = mid
    lam = 0.5 * (lo + hi)
    _, s_i, theta_i = values(lam)
    return lam, s_i, theta_i


def numbered_files(run: Path) -> list[Path]:
    files: list[tuple[int, Path]] = []
    for path in run.glob("Data_*.h5"):
        match = DATA_RE.match(path.name)
        if match:
            files.append((int(match.group(1)), path))
    return [path for _, path in sorted(files)]


def physical(h5: h5py.File, name: str) -> np.ndarray:
    nx = len(h5["grid/xc"]) - 1
    ny = len(h5["grid/yc"]) - 1
    nz = len(h5["grid/zc"]) - 1
    return np.asarray(h5[name][:nz, :ny, :nx], dtype=np.float64)


def front_position(x: np.ndarray, f: np.ndarray) -> float:
    row = f.mean(axis=(0, 1))
    return float(np.interp(0.5, row[::-1], x[::-1]))


def analytic_salt(x: np.ndarray, time: float, lam: float, s_i: float) -> np.ndarray:
    eta_s = (x - X0) / (2.0 * math.sqrt(time / PE_S))
    lam_s = lam * math.sqrt(LE)
    return S_INF + (s_i - S_INF) * (1.0 + np.vectorize(math.erf)(eta_s)) / (
        1.0 + math.erf(lam_s)
    )


def analyze(run: Path, label: str) -> dict[str, float | str]:
    lam_ref, s_i_ref, theta_i_ref = reference_state()
    records: list[dict[str, float]] = []

    for path in numbered_files(run):
        with h5py.File(path, "r") as h5:
            nx = len(h5["grid/xc"]) - 1
            x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
            time = float(np.asarray(h5["time"]).flat[0])
            f = physical(h5, "VOF/C_L")
            theta = physical(h5, "Conc/0")
            salt = physical(h5, "Conc/1")
            u = physical(h5, "u")
            v = physical(h5, "v")
            w = physical(h5, "w")

        front = front_position(x, f)
        salt_x = salt.mean(axis=(0, 1))
        f_x = f.mean(axis=(0, 1))
        dx = float(np.mean(np.diff(x)))
        liquid = (f_x > 0.9) & (x < front - 4.0 * dx)
        if time > 0.0 and np.any(liquid):
            s_exact = analytic_salt(x, time, lam_ref, s_i_ref)
            profile_linf = float(np.max(np.abs(salt_x[liquid] - s_exact[liquid])))
        else:
            profile_linf = math.nan

        records.append(
            {
                "time": time,
                "front": front,
                "plain_salt": float(np.mean(salt)),
                "physical_salt": float(np.mean((f + DELTA) * salt)),
                "enthalpy": float(np.mean(theta + f / STEFAN)),
                "profile_linf": profile_linf,
                "max_speed": float(np.sqrt(u*u + v*v + w*w).max()),
                "finite": float(
                    all(np.isfinite(a).all() for a in (f, theta, salt, u, v, w))
                ),
                "salt_min": float(salt.min()),
                "salt_max": float(salt.max()),
            }
        )

    if len(records) < 2:
        raise RuntimeError(f"{run}: fewer than two snapshots")

    fit = [r for r in records if r["time"] >= 0.2]
    root_t = np.sqrt([r["time"] for r in fit])
    fronts = np.asarray([r["front"] for r in fit])
    slope, intercept = np.polyfit(root_t, fronts, 1)
    lam_fit = float(slope * math.sqrt(PE_T) / 2.0)
    pred = intercept + slope * root_t
    fit_rms = float(np.sqrt(np.mean((fronts - pred) ** 2)))

    first, last = records[0], records[-1]
    phys_scale = max(abs(first["physical_salt"]), 1.0e-30)
    plain_scale = max(abs(first["plain_salt"]), 1.0e-30)
    enth_scale = max(abs(first["enthalpy"]), 1.0e-30)
    profile_errors = [r["profile_linf"] for r in fit if math.isfinite(r["profile_linf"])]

    return {
        "branch": label,
        "outputs": float(len(records)),
        "final_time": last["time"],
        "finite": min(r["finite"] for r in records),
        "lambda_ref": lam_ref,
        "lambda_fit": lam_fit,
        "lambda_rel_error": abs(lam_fit / lam_ref - 1.0),
        "fit_rms": fit_rms,
        "s_interface_ref": s_i_ref,
        "theta_interface_ref": theta_i_ref,
        "profile_linf_final": last["profile_linf"],
        "profile_linf_max_fit": max(profile_errors) if profile_errors else math.nan,
        "physical_salt_rel_drift": abs(last["physical_salt"] - first["physical_salt"]) / phys_scale,
        "plain_salt_rel_drift": abs(last["plain_salt"] - first["plain_salt"]) / plain_scale,
        "enthalpy_rel_drift": abs(last["enthalpy"] - first["enthalpy"]) / enth_scale,
        "max_speed": max(r["max_speed"] for r in records),
        "salt_min": min(r["salt_min"] for r in records),
        "salt_max": max(r["salt_max"] for r in records),
    }


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--legacy", type=Path, required=True)
    parser.add_argument("--yang", type=Path, required=True)
    parser.add_argument("--output-prefix", default="gate")
    args = parser.parse_args()

    rows = [analyze(args.legacy.resolve(), "legacy"), analyze(args.yang.resolve(), "yang")]
    out_csv = Path(f"{args.output_prefix}_comparison.csv")
    with out_csv.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    legacy, yang = rows
    stable = bool(yang["finite"] == 1.0 and yang["max_speed"] < 1.0e-10)
    profile_improved = bool(yang["profile_linf_final"] < legacy["profile_linf_final"])
    front_improved = bool(yang["lambda_rel_error"] < legacy["lambda_rel_error"])
    physical_salt_ok = bool(yang["physical_salt_rel_drift"] < 1.0e-2)
    proceed_refined = stable and physical_salt_ok and (profile_improved or front_improved)

    lines = [
        "# Yang Le=100 coarse salt-transport gate",
        "",
        f"Reference: lambda={yang['lambda_ref']:.8g}, s_i={yang['s_interface_ref']:.8g}, "
        f"theta_i={yang['theta_interface_ref']:.8g}.",
        "",
        "| branch | lambda | rel. error | final salt-profile Linf | physical-salt drift | max speed |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['branch']} | {row['lambda_fit']:.8g} | {row['lambda_rel_error']:.3e} | "
            f"{row['profile_linf_final']:.3e} | {row['physical_salt_rel_drift']:.3e} | "
            f"{row['max_speed']:.3e} |"
        )
    lines.extend(
        [
            "",
            f"Stable/quiescent: **{stable}**.",
            f"Physical-salt drift below 1%: **{physical_salt_ok}**.",
            f"Liquid salinity profile improved: **{profile_improved}**.",
            f"Front constant improved: **{front_improved}**.",
            "",
            f"## Decision: {'RUN REFINED GATE' if proceed_refined else 'STOP — DO NOT RUN REFINED GATE'}",
            "",
            "This is a predeclared direction/stability gate. It is not a production acceptance test.",
        ]
    )
    Path(f"{args.output_prefix}_summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    # A scientific STOP decision is a successful gate outcome, not a batch
    # failure. Numerical or I/O failures still raise above and return nonzero.
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
