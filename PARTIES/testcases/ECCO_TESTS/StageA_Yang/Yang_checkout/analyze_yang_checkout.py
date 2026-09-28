#!/usr/bin/env python3
"""Post-process the reduced integrated Yang checkout matrix.

Expected runs below --data-root:
  run_nz4_1rank, run_nz4_4rank, run_nz8_4rank, run_static_control

Writes summary/time-series/interface CSV files and three PNG figures.
"""

from pathlib import Path
import argparse
import csv
import glob
import re

import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

COUPLED_RUNS = ("run_nz4_1rank", "run_nz4_4rank", "run_nz8_4rank")
CONTROL_RUN = "run_static_control"
RUNS = COUPLED_RUNS + (CONTROL_RUN,)
STEFAN_BY_RUN = {run: 0.25 for run in COUPLED_RUNS} | {CONTROL_RUN: 0.0}


def files(directory):
    found = glob.glob(str(directory / "Data_*.h5"))
    return sorted(found, key=lambda s: int(re.search(r"Data_(\d+)", s).group(1)))


def physical(dataset):
    """Remove the high-side storage/ghost plane in each direction."""
    return np.asarray(dataset[:-1, :-1, :-1], dtype=float)


def interface_positions(x, Fmean):
    positions = []
    for row in Fmean:
        positions.append(np.interp(0.5, row[::-1], x[::-1]))
    return np.asarray(positions)


def read_snapshot(filename, stefan):
    with h5py.File(filename, "r") as h5:
        x = np.asarray(h5["grid/xc"][:-1])
        y = np.asarray(h5["grid/yc"][:-1])
        F = physical(h5["VOF/C_L"])
        theta = physical(h5["Conc/0"])
        salt = physical(h5["Conc/1"])
        u = physical(h5["u"])
        v = physical(h5["v"])
        w = physical(h5["w"])
        time = float(h5["time"][0])
    speed = np.sqrt(u*u + v*v + w*w)
    Fxy = F.mean(axis=0)
    front = interface_positions(x, Fxy)
    deep = F < 1.0e-3
    liquid = F > 0.5
    return {
        "time": time, "x": x, "y": y, "F": F, "theta": theta,
        "salt": salt, "u": u, "v": v, "w": w, "speed": speed,
        "liquid_volume": float(F.mean()),
        "salt_total": float(salt.mean()),
        "theta_total": float(theta.mean()),
        # With phase change disabled (stefan=0), theta itself is the conserved
        # thermal quantity used by the stationary control.
        "enthalpy": float((theta + F/stefan).mean()) if stefan else float(theta.mean()),
        "F_min": float(F.min()), "F_max": float(F.max()),
        "theta_min": float(theta.min()), "theta_max": float(theta.max()),
        "salt_min": float(salt.min()), "salt_max": float(salt.max()),
        "max_speed": float(speed.max()),
        "max_speed_liquid": float(speed[liquid].max()) if np.any(liquid) else 0.0,
        "max_speed_deep_ice": float(speed[deep].max()) if np.any(deep) else 0.0,
        "max_salt_deep_ice": float(np.abs(salt[deep]).max()) if np.any(deep) else 0.0,
        "front_mean": float(front.mean()), "front_std": float(front.std()),
        "front_min": float(front.min()), "front_max": float(front.max()),
        "span_F": float(np.max(np.abs(F-F.mean(axis=0, keepdims=True)))),
        "span_theta": float(np.max(np.abs(theta-theta.mean(axis=0, keepdims=True)))),
        "span_salt": float(np.max(np.abs(salt-salt.mean(axis=0, keepdims=True)))),
        "span_speed": float(np.max(np.abs(speed-speed.mean(axis=0, keepdims=True)))),
        "finite": bool(all(np.isfinite(a).all() for a in (F, theta, salt, u, v, w))),
    }


def log_metrics(directory):
    path = directory / "run.log"
    text = path.read_text(errors="replace") if path.exists() else ""
    divergence = [float(v) for v in re.findall(r"Maximum Velocity Divergence:\s*([0-9.eE+-]+)", text)]
    # The log contains every pressure-correction outer iteration.  Values at or
    # below the solver's 1e-6 stopping criterion are the accepted projection
    # residuals; larger values are intermediate iterations, not failed steps.
    accepted_divergence = [v for v in divergence if v <= 1.0e-6]
    velocity_iterations = [int(v) for v in re.findall(r"velocity converged to .*? after (\d+) iterations", text)]
    ch_iterations = [int(v) for v in re.findall(r"biharmonic converged to .*? after (\d+) iterations", text)]
    bad = bool(re.search(r"\b(?:nan|inf)\b|did not converge|Error called by|MPI_ABORT", text, re.I))
    return {
        "completed": "Simulation complete" in text,
        "log_failure_token": bad,
        "max_projection_divergence": max(divergence, default=np.nan),
        "max_accepted_divergence": max(accepted_divergence, default=np.nan),
        "max_velocity_iterations": max(velocity_iterations, default=0),
        "max_ch_iterations": max(ch_iterations, default=0),
    }


def max_field_difference(a, b):
    return max(float(np.max(np.abs(a[key]-b[key])))
               for key in ("F", "theta", "salt", "u", "v", "w"))


def max_mean_difference(a, b):
    return max(float(np.max(np.abs(a[key].mean(axis=0)-b[key].mean(axis=0))))
               for key in ("F", "theta", "salt", "u", "v", "w"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path,
                        default=Path("/bigscratch/mjalabert314/Yang_checkout"))
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)

    snapshots, logs = {}, {}
    timeseries, interfaces, summary = [], [], []
    for run in RUNS:
        directory = args.data_root/run
        fs = files(directory)
        if not fs:
            raise FileNotFoundError(f"No Data_*.h5 in {directory}")
        snaps = [read_snapshot(f, STEFAN_BY_RUN[run]) for f in fs]
        snapshots[run] = snaps
        logs[run] = log_metrics(directory)
        initial = snaps[0]
        for snap in snaps:
            salt_drift = ((snap["salt_total"]-initial["salt_total"])/
                          max(abs(initial["salt_total"]), 1e-30))
            enthalpy_drift = ((snap["enthalpy"]-initial["enthalpy"])/
                              max(abs(initial["enthalpy"]), 1e-30))
            timeseries.append((run, snap["time"], snap["liquid_volume"],
                               snap["liquid_volume"]-initial["liquid_volume"],
                               snap["salt_total"], salt_drift, snap["enthalpy"],
                               enthalpy_drift, snap["max_speed"],
                               snap["max_speed_deep_ice"], snap["max_salt_deep_ice"]))
            Fxy = snap["F"].mean(axis=0)
            front = interface_positions(snap["x"], Fxy)
            for j, pos in enumerate(front):
                interfaces.append((run, snap["time"], snap["y"][j], pos))

        final = snaps[-1]
        initial_deep = initial["F"] < 1e-3
        deep_salt_change = float(np.max(np.abs(final["salt"][initial_deep]-
                                               initial["salt"][initial_deep])))
        max_F_change = float(np.max(np.abs(final["F"]-initial["F"])))
        vals = {
            "time_final": final["time"], "outputs": len(snaps),
            "completed": int(logs[run]["completed"]),
            "log_failure_token": int(logs[run]["log_failure_token"]),
            "all_fields_finite": int(final["finite"]),
            "F_min": final["F_min"], "F_max": final["F_max"],
            "theta_min": final["theta_min"], "theta_max": final["theta_max"],
            "salt_min": final["salt_min"], "salt_max": final["salt_max"],
            "liquid_volume_change": final["liquid_volume"]-initial["liquid_volume"],
            "salt_relative_drift": abs((final["salt_total"]-initial["salt_total"])/
                                       max(abs(initial["salt_total"]), 1e-30)),
            "enthalpy_relative_drift": abs((final["enthalpy"]-initial["enthalpy"])/
                                           max(abs(initial["enthalpy"]), 1e-30)),
            "deep_ice_salt_change": deep_salt_change,
            "max_F_change": max_F_change,
            "max_speed": final["max_speed"],
            "max_speed_liquid": final["max_speed_liquid"],
            "max_speed_deep_ice": final["max_speed_deep_ice"],
            "front_mean": final["front_mean"], "front_std": final["front_std"],
            "front_displacement": final["front_mean"]-initial["front_mean"],
            "span_F": final["span_F"], "span_theta": final["span_theta"],
            "span_salt": final["span_salt"], "span_speed": final["span_speed"],
            **logs[run],
        }
        for metric, value in vals.items():
            summary.append((run, metric, float(value)))

    one = snapshots["run_nz4_1rank"][-1]
    four = snapshots["run_nz4_4rank"][-1]
    eight = snapshots["run_nz8_4rank"][-1]
    mpi_difference = max_field_difference(one, four)
    nz_difference = max_mean_difference(four, eight)
    summary.extend([
        ("cross_run", "nz4_1_vs_4_worst_difference", mpi_difference),
        ("cross_run", "nz4_vs_nz8_mean_worst_difference", nz_difference),
    ])

    summary_map = {(run, metric): value for run, metric, value in summary}
    run_passes = []
    for run in RUNS:
        common_pass = (
            summary_map[(run, "completed")] == 1 and
            summary_map[(run, "log_failure_token")] == 0 and
            summary_map[(run, "all_fields_finite")] == 1 and
            summary_map[(run, "F_min")] >= -1e-10 and
            summary_map[(run, "F_max")] <= 1.0+1e-10 and
            summary_map[(run, "max_accepted_divergence")] <= 1e-6)
        if run == CONTROL_RUN:
            passed = (
                common_pass and
                summary_map[(run, "salt_relative_drift")] < 1e-12 and
                summary_map[(run, "enthalpy_relative_drift")] < 1e-12 and
                # The analytic tanh initialization undergoes a small CH
                # equilibration.  Require negligible global-volume change and
                # less than 0.01 local F change / 0.013 coarse grid cells of
                # front motion, rather than an inappropriate roundoff test.
                abs(summary_map[(run, "liquid_volume_change")]) < 1e-7 and
                summary_map[(run, "max_F_change")] < 1e-2 and
                abs(summary_map[(run, "front_displacement")]) < 1e-4 and
                summary_map[(run, "max_speed")] < 1e-12)
        else:
            passed = (
                common_pass and
                summary_map[(run, "salt_relative_drift")] < 1e-8 and
                summary_map[(run, "enthalpy_relative_drift")] < 1e-6 and
                summary_map[(run, "deep_ice_salt_change")] < 1e-6 and
                summary_map[(run, "liquid_volume_change")] > 0.0)
        run_passes.append(passed)
        summary.append((run, "run_pass", float(passed)))
    overall = all(run_passes) and mpi_difference < 1e-9 and nz_difference < 1e-7
    summary.extend([
        ("cross_run", "mpi_pass", float(mpi_difference < 1e-9)),
        ("cross_run", "nz_consistency_pass", float(nz_difference < 1e-7)),
        ("overall", "checkout_pass", float(overall)),
    ])

    with open(args.output_dir/"yang_checkout_summary.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(("run", "metric", "value")); w.writerows(summary)
    with open(args.output_dir/"yang_checkout_timeseries.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(("run", "time", "liquid_volume", "liquid_volume_change",
                    "salt_total", "salt_relative_drift", "enthalpy",
                    "enthalpy_relative_drift", "max_speed", "max_speed_deep_ice",
                    "max_salt_deep_ice"))
        w.writerows(timeseries)
    with open(args.output_dir/"yang_checkout_interface.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(("run", "time", "y", "x_F_0p5")); w.writerows(interfaces)

    main_run = snapshots["run_nz4_4rank"][-1]
    extent = (main_run["x"][0], main_run["x"][-1], main_run["y"][0], main_run["y"][-1])
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.2))
    for ax, key, title in zip(axes, ("F", "theta", "salt"),
                              ("liquid fraction F", "temperature theta", "salinity s")):
        field = main_run[key].mean(axis=0)
        im = ax.imshow(field, origin="lower", extent=extent, aspect="equal")
        ax.contour(main_run["x"], main_run["y"], main_run["F"].mean(axis=0),
                   levels=[.5], colors="w", linewidths=1)
        ax.set(title=title, xlabel="x", ylabel="y")
        fig.colorbar(im, ax=ax, shrink=.82)
    fig.suptitle("Reduced integrated Yang checkout, final state")
    fig.tight_layout(); fig.savefig(args.output_dir/"yang_checkout_fields.png", dpi=160); plt.close(fig)

    fig, axes = plt.subplots(2, 2, figsize=(11, 8))
    for run in RUNS:
        rows = [r for r in timeseries if r[0] == run]
        t = [r[1] for r in rows]
        axes[0, 0].plot(t, [r[3] for r in rows], "o-", label=run)
        axes[0, 1].semilogy(t, np.maximum(np.abs([r[5] for r in rows]), 1e-18), "o-", label=run)
        axes[1, 0].semilogy(t, np.maximum(np.abs([r[7] for r in rows]), 1e-18), "o-", label=run)
        axes[1, 1].semilogy(t, np.maximum([r[8] for r in rows], 1e-18), "o-", label=run)
    axes[0, 0].set(title="Liquid-volume increase", ylabel="Delta <F>")
    axes[0, 1].set(title="Salt conservation", ylabel="relative drift")
    axes[1, 0].set(title="Thermal/enthalpy conservation", ylabel="relative drift")
    axes[1, 1].set(title="Flow growth", ylabel="max speed")
    for ax in axes.flat:
        ax.set_xlabel("time"); ax.grid(alpha=.3); ax.legend(fontsize=8)
    fig.tight_layout(); fig.savefig(args.output_dir/"yang_checkout_timeseries.png", dpi=160); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7, 4.5))
    labels = ["1 vs 4 ranks\nNZ=4", "NZ=4 vs 8\nz-mean"]
    values = [mpi_difference, nz_difference]
    ax.bar(labels, values)
    ax.set_yscale("log"); ax.axhline(1e-9, color="C1", ls="--", label="MPI tolerance")
    ax.axhline(1e-7, color="C3", ls=":", label="NZ tolerance")
    ax.set(ylabel="worst field difference", title="Checkout consistency checks")
    ax.grid(axis="y", alpha=.3); ax.legend(); fig.tight_layout()
    fig.savefig(args.output_dir/"yang_checkout_consistency.png", dpi=160); plt.close(fig)

    print("Yang checkout analysis written to", args.output_dir)
    print("overall pass:", overall)
    print("MPI difference:", mpi_difference)
    print("NZ-mean difference:", nz_difference)


if __name__ == "__main__":
    main()
