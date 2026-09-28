#!/usr/bin/env python3
"""Analyze the Yang production run and restart sensitivity branches.

Examples
--------
Analyze the completed production run::

    python3 analyze_yang_production.py \
      --run production=/bigscratch/mjalabert314/Yang_production \
      --output-dir /tmp/yang_analysis

Compare branches resumed from the t=20 output::

    python3 analyze_yang_production.py \
      --run tau_1e-4=/bigscratch/mjalabert314/Yang_penalty_branches_t20/tau_1e-4 \
      --run tau_2p5e-4=/bigscratch/mjalabert314/Yang_penalty_branches_t20/tau_2p5e-4 \
      --run tau_5e-4=/bigscratch/mjalabert314/Yang_penalty_branches_t20/tau_5e-4 \
      --output-dir /bigscratch/mjalabert314/Yang_penalty_branches_t20/analysis \
      --baseline-file /bigscratch/mjalabert314/Yang_production/Data_0.h5

PARTIES writes one duplicated high-side storage plane per coordinate.  All
integrals and extrema in this script use only ``[:nz, :ny, :nx]``, where the
physical dimensions are obtained from ``grid/{zc,yc,xc}[:-1]``.  This avoids
the volume-normalization error caused by including the duplicated planes.

For restart branches that do not contain Data_0.h5, the default normalization
is the exact configured initial ice volume, 0.1.  Pass ``--baseline-file`` to
derive it directly from the original Data_0.h5 instead.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
from typing import Iterable

import h5py
import numpy as np


DEFAULT_REFERENCE = Path(__file__).with_name("yang_fig4_sm5_reference.csv")
DATA_RE = re.compile(r"Data_(\d+)\.h5$")
FLOAT_RE = r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[Ee][+-]?\d+)?"
TIME_FILTER_TOL = 1.0e-4


def parse_run(value: str) -> tuple[str, Path]:
    """Parse LABEL=PATH supplied to --run."""
    if "=" not in value:
        raise argparse.ArgumentTypeError("--run must have the form LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("--run needs a non-empty label and path")
    return label, Path(raw_path).expanduser().resolve()


def parse_compare_times(value: str) -> list[float]:
    try:
        times = [float(item) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("--compare-times must be comma-separated numbers") from exc
    if not times:
        raise argparse.ArgumentTypeError("--compare-times cannot be empty")
    return times


def data_files(directory: Path) -> list[Path]:
    """Return numbered snapshots, ordered by their output index."""
    found = []
    for path in directory.glob("Data_*.h5"):
        match = DATA_RE.search(path.name)
        if match:
            found.append((int(match.group(1)), path))
    return [path for _, path in sorted(found)]


def snapshot_time(path: Path) -> float:
    """Read only the scalar time, avoiding full-field I/O during filtering."""
    with h5py.File(path, "r") as h5:
        return float(np.asarray(h5["time"]).flat[0])


def physical_shape(h5: h5py.File) -> tuple[int, int, int]:
    """Return the (nz, ny, nx) physical cell counts."""
    nx = len(h5["grid/xc"]) - 1
    ny = len(h5["grid/yc"]) - 1
    nz = len(h5["grid/zc"]) - 1
    if min(nx, ny, nz) < 1:
        raise ValueError(f"invalid physical grid shape {(nz, ny, nx)}")
    return nz, ny, nx


def physical(dataset: h5py.Dataset, shape: tuple[int, int, int]) -> np.ndarray:
    """Load only physical cells, excluding all duplicated high-side planes."""
    nz, ny, nx = shape
    expected_storage = (nz + 1, ny + 1, nx + 1)
    if dataset.shape != expected_storage:
        raise ValueError(
            f"{dataset.name} shape {dataset.shape} is not the expected "
            f"PARTIES storage shape {expected_storage}"
        )
    return np.asarray(dataset[:nz, :ny, :nx], dtype=np.float64)


def front_positions(x: np.ndarray, phase_yx: np.ndarray) -> np.ndarray:
    """Interpolate the x position of the F=0.5 front for every y row."""
    front = np.full(phase_yx.shape[0], np.nan)
    for j, row in enumerate(phase_yx):
        # The Yang slab is liquid on the low-x side, so reversing gives the
        # increasing 0 -> 1 phase coordinate expected by np.interp.
        reverse = row[::-1]
        if np.nanmin(reverse) <= 0.5 <= np.nanmax(reverse):
            front[j] = np.interp(0.5, reverse, x[::-1])
    return front


def read_snapshot(path: Path, stefan: float) -> dict[str, float]:
    with h5py.File(path, "r") as h5:
        shape = physical_shape(h5)
        nz, ny, nx = shape
        x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
        y = np.asarray(h5["grid/yc"][:ny], dtype=np.float64)
        phase = physical(h5["VOF/C_L"], shape)
        theta = physical(h5["Conc/0"], shape)
        salt = physical(h5["Conc/1"], shape)
        u = physical(h5["u"], shape)
        v = physical(h5["v"], shape)
        w = physical(h5["w"], shape)
        time = float(np.asarray(h5["time"]).flat[0])

    finite = all(np.isfinite(field).all() for field in (phase, theta, salt, u, v, w))
    speed2 = u*u + v*v + w*w
    speed = np.sqrt(speed2)
    solid = phase < 0.1  # Yang's direct-forcing threshold is solid fraction > 0.9.
    deep_ice = phase < 1.0e-3
    phase_yx = phase.mean(axis=0)
    front = front_positions(x, phase_yx)
    valid_front = front[np.isfinite(front)]

    if valid_front.size:
        front_mean = float(valid_front.mean())
        front_std = float(valid_front.std())
        front_min = float(valid_front.min())
        front_max = float(valid_front.max())
        # Average 2% of the wall-normal rows at each horizontal boundary.  The
        # sign is useful for the fresh-lid signature: negative means more ice
        # remains at the top (a smaller front x position).
        edge = max(1, int(round(0.02*len(front))))
        front_top_minus_bottom = float(np.nanmean(front[-edge:]) - np.nanmean(front[:edge]))
    else:
        front_mean = front_std = front_min = front_max = math.nan
        front_top_minus_bottom = math.nan

    ice_volume = float(np.mean(1.0 - phase))
    enthalpy = float(np.mean(theta + phase/stefan))
    return {
        "output_index": float(DATA_RE.search(path.name).group(1)),
        "time": time,
        "ice_volume": ice_volume,
        "liquid_fraction_mean": float(phase.mean()),
        "salt_mean": float(salt.mean()),
        "enthalpy_mean": enthalpy,
        "F_min": float(phase.min()),
        "F_max": float(phase.max()),
        "theta_min": float(theta.min()),
        "theta_max": float(theta.max()),
        "salt_min": float(salt.min()),
        "salt_max": float(salt.max()),
        "max_speed": float(speed.max()),
        "solid_speed_rms": float(np.sqrt(speed2[solid].mean())) if np.any(solid) else math.nan,
        "solid_speed_max": float(speed[solid].max()) if np.any(solid) else math.nan,
        "deep_ice_salt_max_abs": float(np.abs(salt[deep_ice]).max()) if np.any(deep_ice) else math.nan,
        "front_mean": front_mean,
        "front_std": front_std,
        "front_min": front_min,
        "front_max": front_max,
        "front_peak_to_peak": front_max - front_min,
        "front_top_minus_bottom": front_top_minus_bottom,
        "finite": float(finite),
        "nx": float(nx),
        "ny": float(ny),
        "nz": float(nz),
        "y_extent": float(y[-1] - y[0]) if len(y) > 1 else 0.0,
    }


def scan_log(path: Path) -> dict[str, float]:
    """Stream run.log so the multi-gigabyte production log is never loaded."""
    result = {
        "log_present": float(path.exists()),
        "simulation_complete": 0.0,
        "log_failure_token": 0.0,
        "max_accepted_projection_divergence": math.nan,
        "max_ch_iterations": 0.0,
        "max_velocity_iterations": 0.0,
    }
    if not path.exists():
        return result

    divergence_re = re.compile(rf"Maximum Velocity Divergence:\s*({FLOAT_RE})")
    ch_re = re.compile(r"biharmonic converged to .*? after (\d+) iterations", re.I)
    velocity_re = re.compile(r"velocity converged to .*? after (\d+) iterations", re.I)
    failure_re = re.compile(r"\b(?:nan|inf)\b|did not converge|Error called by|MPI_ABORT", re.I)
    max_accepted = math.nan
    with path.open(errors="replace") as stream:
        for line in stream:
            if "Simulation complete" in line:
                result["simulation_complete"] = 1.0
            if failure_re.search(line):
                result["log_failure_token"] = 1.0
            match = divergence_re.search(line)
            if match:
                value = float(match.group(1))
                if value <= 1.0e-6:
                    max_accepted = value if math.isnan(max_accepted) else max(max_accepted, value)
            match = ch_re.search(line)
            if match:
                result["max_ch_iterations"] = max(
                    result["max_ch_iterations"], float(match.group(1))
                )
            match = velocity_re.search(line)
            if match:
                result["max_velocity_iterations"] = max(
                    result["max_velocity_iterations"], float(match.group(1))
                )
    result["max_accepted_projection_divergence"] = max_accepted
    return result


def scan_logs(directory: Path) -> dict[str, float]:
    """Aggregate the production log and branch logs named by SLURM job ID."""
    paths = []
    conventional = directory / "run.log"
    if conventional.exists():
        paths.append(conventional)
    paths.extend(sorted(directory.glob("run_*.log")))
    if not paths:
        return scan_log(conventional)

    summaries = [scan_log(path) for path in paths]
    divergences = [
        summary["max_accepted_projection_divergence"]
        for summary in summaries
        if math.isfinite(summary["max_accepted_projection_divergence"])
    ]
    return {
        "log_present": 1.0,
        "simulation_complete": max(summary["simulation_complete"] for summary in summaries),
        "log_failure_token": max(summary["log_failure_token"] for summary in summaries),
        "max_accepted_projection_divergence": max(divergences) if divergences else math.nan,
        "max_ch_iterations": max(summary["max_ch_iterations"] for summary in summaries),
        "max_velocity_iterations": max(
            summary["max_velocity_iterations"] for summary in summaries
        ),
    }


def load_reference(path: Path) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    times, volumes, uncertainties = [], [], []
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            times.append(float(row["time_t_ff"]))
            volumes.append(float(row["ice_volume_ratio"]))
            uncertainties.append(float(row["digitization_uncertainty"]))
    order = np.argsort(times)
    return (
        np.asarray(times)[order],
        np.asarray(volumes)[order],
        np.asarray(uncertainties)[order],
    )


def interpolate_bounded(x: float, xp: np.ndarray, fp: np.ndarray) -> float:
    if x < xp[0] or x > xp[-1]:
        return math.nan
    return float(np.interp(x, xp, fp))


def half_volume_time(times: np.ndarray, ratios: np.ndarray) -> float:
    """Linearly interpolate the first downward crossing of V/V0=0.5."""
    for i in range(1, len(times)):
        if ratios[i - 1] > 0.5 >= ratios[i]:
            return float(np.interp(0.5, ratios[i - 1:i + 1][::-1], times[i - 1:i + 1][::-1]))
        if ratios[i] == 0.5:
            return float(times[i])
    return math.nan


def write_rows(path: Path, fieldnames: Iterable[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", action="append", type=parse_run, required=True, metavar="LABEL=PATH",
        help="run directory; repeat to compare production/restart branches",
    )
    parser.add_argument("--output-dir", type=Path, default=Path.cwd())
    parser.add_argument("--prefix", default="yang_production")
    parser.add_argument("--reference-csv", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument(
        "--baseline-file", type=Path,
        help="optional original Data_0.h5 used to measure V0 instead of nominal 0.1",
    )
    parser.add_argument(
        "--initial-ice-volume", type=float, default=0.1,
        help="V0 used when --baseline-file is absent (default: exact configured 0.1)",
    )
    parser.add_argument("--stefan", type=float, default=0.25)
    parser.add_argument(
        "--compare-times", type=parse_compare_times,
        default=parse_compare_times("40,60,80,100,118,120,140,160,180,198"),
        help="comma-separated times for interpolated run/reference comparisons",
    )
    parser.add_argument("--start-time", type=float, default=-math.inf)
    parser.add_argument("--end-time", type=float, default=math.inf)
    parser.add_argument(
        "--skip-log", action="store_true",
        help="do not scan run.log (useful for a quick HDF-only analysis)",
    )
    args = parser.parse_args()

    if args.stefan <= 0.0 or args.initial_ice_volume <= 0.0:
        parser.error("--stefan and --initial-ice-volume must be positive")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    initial_ice_volume = args.initial_ice_volume
    normalization_source = "configured_initial_ice_volume"
    if args.baseline_file:
        baseline = read_snapshot(args.baseline_file.resolve(), args.stefan)
        initial_ice_volume = baseline["ice_volume"]
        normalization_source = str(args.baseline_file.resolve())

    ref_t, ref_v, ref_u = load_reference(args.reference_csv.resolve())
    timeseries_rows: list[dict] = []
    comparison_rows: list[dict] = []
    summary_rows: list[dict] = []

    for label, directory in args.run:
        files = data_files(directory)
        if not files:
            raise FileNotFoundError(f"no Data_*.h5 files in {directory}")

        snapshots = []
        for path in files:
            time = snapshot_time(path)
            # Restarted output times retain a few microseconds of accumulated
            # timestep offset (for example, t=80 is stored as 80.0000064).
            # Treat requested diagnostic endpoints as physical times rather
            # than excluding those snapshots through exact float comparison.
            if (
                args.start_time - TIME_FILTER_TOL
                <= time
                <= args.end_time + TIME_FILTER_TOL
            ):
                snapshots.append(read_snapshot(path, args.stefan))
        if not snapshots:
            raise ValueError(f"no snapshots in requested time range for {directory}")
        snapshots.sort(key=lambda row: row["time"])

        salt0 = snapshots[0]["salt_mean"]
        enthalpy0 = snapshots[0]["enthalpy_mean"]
        for snap in snapshots:
            row = {"run": label, "directory": str(directory), **snap}
            row["ice_volume_ratio"] = snap["ice_volume"] / initial_ice_volume
            row["salt_relative_drift_from_first_analyzed"] = (
                (snap["salt_mean"] - salt0) / max(abs(salt0), 1.0e-30)
            )
            row["enthalpy_relative_drift_from_first_analyzed"] = (
                (snap["enthalpy_mean"] - enthalpy0) / max(abs(enthalpy0), 1.0e-30)
            )
            timeseries_rows.append(row)

        sim_t = np.asarray([snap["time"] for snap in snapshots])
        sim_v = np.asarray([snap["ice_volume"] / initial_ice_volume for snap in snapshots])
        for compare_time in args.compare_times:
            sim_value = interpolate_bounded(compare_time, sim_t, sim_v)
            reference_value = interpolate_bounded(compare_time, ref_t, ref_v)
            uncertainty = interpolate_bounded(compare_time, ref_t, ref_u)
            if math.isfinite(sim_value) and math.isfinite(reference_value):
                comparison_rows.append({
                    "run": label,
                    "time_t_ff": compare_time,
                    "simulated_ice_volume_ratio": sim_value,
                    "reference_ice_volume_ratio": reference_value,
                    "sim_minus_reference": sim_value - reference_value,
                    "relative_error": (sim_value - reference_value) / reference_value,
                    "reference_digitization_uncertainty": uncertainty,
                })

        log = scan_logs(directory) if not args.skip_log else {}
        metrics = {
            "normalization_initial_ice_volume": initial_ice_volume,
            "time_first": float(sim_t[0]),
            "time_final": float(sim_t[-1]),
            "outputs_analyzed": float(len(snapshots)),
            "half_volume_time": half_volume_time(sim_t, sim_v),
            "final_ice_volume_ratio": float(sim_v[-1]),
            "max_abs_salt_relative_drift_from_first_analyzed": max(
                abs(row["salt_relative_drift_from_first_analyzed"])
                for row in timeseries_rows if row["run"] == label
            ),
            "max_abs_enthalpy_relative_drift_from_first_analyzed": max(
                abs(row["enthalpy_relative_drift_from_first_analyzed"])
                for row in timeseries_rows if row["run"] == label
            ),
            "F_min_all_outputs": min(snap["F_min"] for snap in snapshots),
            "F_max_all_outputs": max(snap["F_max"] for snap in snapshots),
            "max_speed_all_outputs": max(snap["max_speed"] for snap in snapshots),
            "max_solid_speed_all_outputs": max(snap["solid_speed_max"] for snap in snapshots),
            **log,
        }
        for compare_time in (40.0, 60.0, 80.0):
            value = interpolate_bounded(compare_time, sim_t, sim_v)
            if math.isfinite(value):
                metrics[f"ice_volume_ratio_t{int(compare_time)}"] = value
        v40 = interpolate_bounded(40.0, sim_t, sim_v)
        v80 = interpolate_bounded(80.0, sim_t, sim_v)
        if math.isfinite(v40) and math.isfinite(v80):
            metrics["ice_volume_slope_t40_t80"] = (v80 - v40) / 40.0

        for metric, value in metrics.items():
            summary_rows.append({"run": label, "metric": metric, "value": value})

    timeseries_fields = [
        "run", "directory", "output_index", "time", "ice_volume",
        "ice_volume_ratio", "liquid_fraction_mean", "salt_mean",
        "salt_relative_drift_from_first_analyzed", "enthalpy_mean",
        "enthalpy_relative_drift_from_first_analyzed", "F_min", "F_max",
        "theta_min", "theta_max", "salt_min", "salt_max", "max_speed",
        "solid_speed_rms", "solid_speed_max", "deep_ice_salt_max_abs",
        "front_mean", "front_std", "front_min", "front_max",
        "front_peak_to_peak", "front_top_minus_bottom", "finite", "nx", "ny",
        "nz", "y_extent",
    ]
    comparison_fields = [
        "run", "time_t_ff", "simulated_ice_volume_ratio",
        "reference_ice_volume_ratio", "sim_minus_reference", "relative_error",
        "reference_digitization_uncertainty",
    ]
    write_rows(args.output_dir / f"{args.prefix}_timeseries.csv", timeseries_fields, timeseries_rows)
    write_rows(args.output_dir / f"{args.prefix}_reference_comparison.csv", comparison_fields, comparison_rows)
    write_rows(args.output_dir / f"{args.prefix}_summary.csv", ("run", "metric", "value"), summary_rows)

    print(f"Analyzed {len(args.run)} run(s); normalization V0={initial_ice_volume:.12g}")
    print(f"Normalization source: {normalization_source}")
    print(f"Reference: {args.reference_csv.resolve()}")
    print(f"Outputs: {args.output_dir.resolve() / (args.prefix + '_*.csv')}")


if __name__ == "__main__":
    main()
