#!/usr/bin/env python3
"""Matched-ice-fraction transport audit for the Yang benchmark.

The expensive Yang simulations need not be rerun for this audit.  The script
scans existing ``Data_*.h5`` snapshots, selects the nearest states to requested
ice-volume ratios, and reconstructs the instantaneous Stefan source used by
PARTIES:

    m = St/(Pe_T*eps) * (theta - theta_L) * F*(1-F)/(sqrt(2)*Cn)
    theta_L = T_melt - liquidus_slope*s.

It also measures liquid-side thermal and salinity recovery lengths, the
phase-weighted salt diffusive flux at the F=0.5 surface, near-interface flow
and buoyancy, and vertically binned melt delivery.  All HDF5 calculations
exclude PARTIES' duplicated high-side storage planes.

Example (fresh and salty archives visible on the same filesystem)::

    python3 analyze_yang_matched_transport.py \
      --run sm0=/anvil/scratch/x-mjalabert/Yang_Sm0_07182026 \
      --run sm5=/bigscratch/mjalabert314/Yang_production \
      --output-dir matched_transport_audit

If only the four matched salty snapshots were transferred, use the completed
history to select them and compute the local finite-difference rate::

    python3 analyze_yang_matched_transport.py \
      --run sm5=/anvil/scratch/x-mjalabert/Yang_production \
      --history-csv sm5=/anvil/scratch/x-mjalabert/Yang_production/final_timeseries.csv \
      --baseline-matched-csv matched_transport_audit_sm0/yang_matched_transport_matched.csv \
      --output-dir matched_transport_audit_sm5

The scan reads only ``time`` and ``VOF/C_L`` from every output.  Full scalar
and velocity fields are read only for the four selected matched states.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path
import re
from typing import Iterable, Sequence

import h5py
import numpy as np


DATA_RE = re.compile(r"Data_(\d+)\.h5$")


def parse_run(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("--run must have the form LABEL=PATH")
    label, raw_path = value.split("=", 1)
    if not label or not raw_path:
        raise argparse.ArgumentTypeError("--run needs a non-empty label and path")
    return label, Path(raw_path).expanduser().resolve()


def parse_float_list(value: str) -> list[float]:
    try:
        result = [float(item) for item in value.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("expected comma-separated numbers") from exc
    if not result:
        raise argparse.ArgumentTypeError("list cannot be empty")
    return result


def data_files(directory: Path) -> list[Path]:
    found: list[tuple[int, Path]] = []
    for path in directory.glob("Data_*.h5"):
        match = DATA_RE.search(path.name)
        if match:
            found.append((int(match.group(1)), path))
    return [path for _, path in sorted(found)]


def parse_input(path: Path) -> dict[str, float | list[float]]:
    """Parse the scalar/list values needed from a PARTIES input file."""
    values: dict[str, float | list[float]] = {}
    for raw_line in path.read_text().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if not line or line.startswith("[") or "=" not in line:
            continue
        key, raw_value = (item.strip() for item in line.split("=", 1))
        try:
            if raw_value.startswith("{") and raw_value.endswith("}"):
                values[key] = [
                    float(item.strip())
                    for item in raw_value[1:-1].split(",")
                    if item.strip()
                ]
            else:
                values[key] = float(raw_value)
        except ValueError:
            continue
    return values


def scalar_parameter(
    values: dict[str, float | list[float]], key: str, default: float | None = None
) -> float:
    value = values.get(key, default)
    if value is None or isinstance(value, list):
        raise ValueError(f"missing scalar parameter {key}")
    return float(value)


def list_parameter(
    values: dict[str, float | list[float]], key: str, minimum: int
) -> list[float]:
    value = values.get(key)
    if not isinstance(value, list) or len(value) < minimum:
        raise ValueError(f"parameter {key} must contain at least {minimum} values")
    return [float(item) for item in value]


def physical_shape(h5: h5py.File) -> tuple[int, int, int]:
    nx = len(h5["grid/xc"]) - 1
    ny = len(h5["grid/yc"]) - 1
    nz = len(h5["grid/zc"]) - 1
    if min(nx, ny, nz) < 1:
        raise ValueError(f"invalid physical grid shape {(nz, ny, nx)}")
    return nz, ny, nx


def scalar_plane(h5: h5py.File, name: str, ny: int, nx: int) -> np.ndarray:
    dataset = h5[name]
    if dataset.shape[0] < 1 or dataset.shape[1] < ny or dataset.shape[2] < nx:
        raise ValueError(f"unexpected shape {dataset.shape} for {name}")
    return np.asarray(dataset[0, :ny, :nx], dtype=np.float64)


def centered_velocity_planes(
    h5: h5py.File, ny: int, nx: int
) -> tuple[np.ndarray, np.ndarray]:
    """Interpolate face-centered u and v to scalar cell centers."""
    u_face = np.asarray(h5["u"][0, :ny, : nx + 1], dtype=np.float64)
    v_face = np.asarray(h5["v"][0, : ny + 1, :nx], dtype=np.float64)
    return 0.5 * (u_face[:, :-1] + u_face[:, 1:]), 0.5 * (
        v_face[:-1, :] + v_face[1:, :]
    )


def scan_snapshot(path: Path) -> dict[str, float | str]:
    with h5py.File(path, "r") as h5:
        _, ny, nx = physical_shape(h5)
        phase = scalar_plane(h5, "VOF/C_L", ny, nx)
        time = float(np.asarray(h5["time"]).flat[0])
    match = DATA_RE.search(path.name)
    assert match is not None
    return {
        "output_index": float(match.group(1)),
        "path": str(path),
        "time": time,
        "ice_volume": float(np.mean(1.0 - phase)),
    }


def scan_history(
    path: Path, directory: Path, initial_ice_volume: float
) -> list[dict[str, float | str]]:
    """Use a saved volume history when only selected HDF snapshots are local."""
    rows = read_csv(path)
    if not rows:
        raise ValueError(f"empty history CSV {path}")
    time_key = "t" if "t" in rows[0] else "time"
    ratio_key = "VV0" if "VV0" in rows[0] else "ice_volume_ratio"
    if time_key not in rows[0] or ratio_key not in rows[0]:
        raise ValueError(
            f"{path} needs t/VV0 or time/ice_volume_ratio columns"
        )
    result: list[dict[str, float | str]] = []
    for sequence_index, row in enumerate(rows):
        output_index = int(float(row.get("output_index", sequence_index)))
        ratio = float(row[ratio_key])
        result.append(
            {
                "output_index": float(output_index),
                "path": str(directory / f"Data_{output_index}.h5"),
                "time": float(row[time_key]),
                "ice_volume": ratio * initial_ice_volume,
            }
        )
    return result


def inferred_initial_salt_mean(parameters: dict[str, float | list[float]]) -> float:
    """Infer the Yang slab's domain-mean initial salt if Data_0 is not local."""
    xmin = scalar_parameter(parameters, "xmin", 0.0)
    xmax = scalar_parameter(parameters, "xmax", 1.0)
    slab = scalar_parameter(parameters, "vof_slab_x0")
    salt_top = scalar_parameter(parameters, "cbd2", 0.0)
    salt_bottom = scalar_parameter(parameters, "cbd5", 0.0)
    liquid_width_fraction = min(max((slab - xmin) / (xmax - xmin), 0.0), 1.0)
    return liquid_width_fraction * 0.5 * (salt_top + salt_bottom)


def front_locations(x: np.ndarray, phase: np.ndarray) -> np.ndarray:
    """Find the east, liquid-to-ice F=0.5 crossing in each horizontal row."""
    ny = phase.shape[0]
    front = np.full(ny, np.nan)
    for j in range(ny):
        row = phase[j]
        crossings = np.flatnonzero((row[:-1] >= 0.5) & (row[1:] < 0.5))
        if crossings.size == 0:
            continue
        i = int(crossings[-1])
        denominator = row[i + 1] - row[i]
        fraction = 0.5 if denominator == 0.0 else (0.5 - row[i]) / denominator
        front[j] = x[i] + fraction * (x[i + 1] - x[i])
    return front


def sample_rows_at_x(
    field: np.ndarray, x: np.ndarray, positions: np.ndarray
) -> np.ndarray:
    """Linearly sample one x position per y row; invalid positions remain NaN."""
    ny, nx = field.shape
    result = np.full(ny, np.nan)
    dx = float(x[1] - x[0])
    coordinate = (positions - x[0]) / dx
    finite_position = np.isfinite(positions)
    lower = np.zeros(ny, dtype=np.int64)
    lower[finite_position] = np.floor(coordinate[finite_position]).astype(np.int64)
    valid = finite_position & (lower >= 0) & (lower < nx - 1)
    rows = np.flatnonzero(valid)
    indices = lower[valid]
    fraction = coordinate[valid] - indices
    result[valid] = (
        (1.0 - fraction) * field[rows, indices]
        + fraction * field[rows, indices + 1]
    )
    return result


def first_recovery_distance(
    samples: np.ndarray,
    interface_value: np.ndarray,
    bulk_value: np.ndarray,
    distances_cells: np.ndarray,
    contrast_threshold: float,
) -> np.ndarray:
    """First distance where 90% of the interface-to-bulk contrast is recovered."""
    ny = samples.shape[0]
    result = np.full(ny, np.nan)
    contrast = bulk_value - interface_value
    usable = np.isfinite(contrast) & (np.abs(contrast) > contrast_threshold)
    normalized = np.full_like(samples, np.nan)
    np.divide(
        samples - interface_value[:, None],
        contrast[:, None],
        out=normalized,
        where=usable[:, None],
    )
    extended_distance = np.concatenate(([0.0], distances_cells))
    for j in np.flatnonzero(usable):
        profile = np.concatenate(([0.0], normalized[j]))
        crossing = np.flatnonzero(profile >= 0.9)
        if crossing.size == 0:
            continue
        k = int(crossing[0])
        if k == 0:
            result[j] = 0.0
            continue
        a, b = profile[k - 1], profile[k]
        if not (math.isfinite(a) and math.isfinite(b)) or b == a:
            continue
        result[j] = extended_distance[k - 1] + (0.9 - a) * (
            extended_distance[k] - extended_distance[k - 1]
        ) / (b - a)
    return result


def weighted_mean(values: np.ndarray, weights: np.ndarray) -> float:
    valid = np.isfinite(values) & np.isfinite(weights) & (weights > 0.0)
    denominator = float(np.sum(weights[valid]))
    if denominator == 0.0:
        return math.nan
    return float(np.sum(values[valid] * weights[valid]) / denominator)


def finite_percentile(values: np.ndarray, percentile: float) -> float:
    finite = values[np.isfinite(values)]
    return float(np.percentile(finite, percentile)) if finite.size else math.nan


def summarize_distribution(prefix: str, values: np.ndarray) -> dict[str, float]:
    return {
        f"{prefix}_p10": finite_percentile(values, 10.0),
        f"{prefix}_median": finite_percentile(values, 50.0),
        f"{prefix}_p90": finite_percentile(values, 90.0),
    }


def selected_snapshot(
    path: Path,
    parameters: dict[str, float | list[float]],
    target_ratio: float,
    actual_ratio: float,
    finite_difference_rate: float,
    initial_salt_mean: float,
    recovery_cells: np.ndarray,
    vertical_bins: int,
) -> tuple[dict[str, float | str], list[dict[str, float]]]:
    stefan = scalar_parameter(parameters, "stefan")
    cn = scalar_parameter(parameters, "Cn")
    band_eps = scalar_parameter(parameters, "melt_band_eps")
    t_melt = scalar_parameter(parameters, "T_melt", 0.0)
    liquidus_slope = scalar_parameter(parameters, "liquidus_slope", 0.0)
    pe = list_parameter(parameters, "Pe", 2)
    ice_diffusivity_ratio = list_parameter(parameters, "kappa_ice_ratio", 2)
    eos_q = scalar_parameter(parameters, "eos_q", 2.0)
    eos_beta_t = scalar_parameter(parameters, "eos_betaT", 1.0)
    eos_beta_s = scalar_parameter(parameters, "eos_betaS", 0.0)
    eos_tmd0 = scalar_parameter(parameters, "eos_Tmd0", 0.0)
    eos_tmd_slope = scalar_parameter(parameters, "eos_Tmd_slope", 0.0)

    with h5py.File(path, "r") as h5:
        nz, ny, nx = physical_shape(h5)
        if nz != 1:
            raise ValueError(f"the matched transport audit currently requires nz=1, got {nz}")
        x = np.asarray(h5["grid/xc"][:nx], dtype=np.float64)
        y = np.asarray(h5["grid/yc"][:ny], dtype=np.float64)
        phase = scalar_plane(h5, "VOF/C_L", ny, nx)
        theta = scalar_plane(h5, "Conc/0", ny, nx)
        salt = scalar_plane(h5, "Conc/1", ny, nx)
        u, v = centered_velocity_planes(h5, ny, nx)
        time = float(np.asarray(h5["time"]).flat[0])

    dx = float(np.median(np.diff(x)))
    dy = float(np.median(np.diff(y)))
    front = front_locations(x, phase)
    valid_front = np.isfinite(front)
    front_slope = np.full(ny, np.nan)
    if np.count_nonzero(valid_front) >= 3:
        valid_indices = np.flatnonzero(valid_front)
        front_slope[valid_front] = np.gradient(front[valid_front], y[valid_front])
        # Do not bridge large gaps of fully melted/frozen rows.
        gaps = np.flatnonzero(np.diff(valid_indices) > 1)
        for gap in gaps:
            front_slope[valid_indices[gap]] = math.nan
            front_slope[valid_indices[gap + 1]] = math.nan

    theta_liquidus = t_melt - liquidus_slope * salt
    band_weight = phase * (1.0 - phase) / (math.sqrt(2.0) * cn)
    superheat = theta - theta_liquidus
    melt_source = stefan / (pe[0] * band_eps) * superheat * band_weight
    buoyancy = -eos_beta_t * np.abs(
        theta - eos_tmd0 - eos_tmd_slope * salt
    ) ** eos_q + eos_beta_s * salt

    row_weight = np.mean(band_weight, axis=1)
    row_melt = np.mean(melt_source, axis=1)
    row_abs_v = np.divide(
        np.sum(np.abs(v) * band_weight, axis=1),
        np.sum(band_weight, axis=1),
        out=np.full(ny, np.nan),
        where=np.sum(band_weight, axis=1) > 0.0,
    )
    row_buoyancy = np.divide(
        np.sum(buoyancy * band_weight, axis=1),
        np.sum(band_weight, axis=1),
        out=np.full(ny, np.nan),
        where=np.sum(band_weight, axis=1) > 0.0,
    )

    sample_theta = np.column_stack(
        [sample_rows_at_x(theta, x, front - cells * dx) for cells in recovery_cells]
    )
    sample_salt = np.column_stack(
        [sample_rows_at_x(salt, x, front - cells * dx) for cells in recovery_cells]
    )
    interface_theta = sample_rows_at_x(theta, x, front)
    interface_salt = sample_rows_at_x(salt, x, front)
    interface_liquidus = t_melt - liquidus_slope * interface_salt
    bulk_theta = sample_theta[:, -1]
    bulk_salt = sample_salt[:, -1]

    delta_t90_cells = first_recovery_distance(
        sample_theta,
        interface_liquidus,
        bulk_theta,
        recovery_cells,
        contrast_threshold=1.0e-4,
    )
    delta_s90_cells = first_recovery_distance(
        sample_salt,
        interface_salt,
        bulk_salt,
        recovery_cells,
        contrast_threshold=1.0e-3,
    )

    flux_cells = min(4.0, float(recovery_cells[-1]))
    flux_index = int(np.argmin(np.abs(recovery_cells - flux_cells)))
    flux_cells = float(recovery_cells[flux_index])
    theta_at_flux = sample_theta[:, flux_index]
    salt_at_flux = sample_salt[:, flux_index]
    heat_flux_proxy = (theta_at_flux - interface_liquidus) / (
        flux_cells * dx * pe[0]
    )
    # At F=0.5, kappa(F)/kappa_liq = kr + (1-kr)*F.  Positive is toward ice.
    salt_interface_diffusivity = ice_diffusivity_ratio[1] + 0.5 * (
        1.0 - ice_diffusivity_ratio[1]
    )
    salt_flux_proxy = salt_interface_diffusivity * (salt_at_flux - interface_salt) / (
        flux_cells * dx * pe[1]
    )

    valid_rows = valid_front & (row_weight > 0.0)
    positive_row_melt = row_melt[valid_rows & (row_melt > 0.0)]
    median_positive_melt = (
        float(np.median(positive_row_melt)) if positive_row_melt.size else math.nan
    )
    starved = valid_rows & (
        row_melt < 0.25 * median_positive_melt
        if math.isfinite(median_positive_melt)
        else False
    )
    refreezing = valid_rows & (row_melt < 0.0)

    if np.count_nonzero(valid_front) >= 2:
        fy = front[valid_front]
        yy = y[valid_front]
        integrand = np.sqrt(1.0 + np.gradient(fy, yy) ** 2)
        front_length = float(
            np.sum(0.5 * (integrand[:-1] + integrand[1:]) * np.diff(yy))
        )
    else:
        front_length = math.nan

    salt_denominator = float(np.sum(np.abs(salt)))
    salt_solid_leakage_index = (
        float(np.sum((1.0 - phase) * np.abs(salt)) / salt_denominator)
        if salt_denominator > 0.0
        else 0.0
    )
    deep_ice = phase < 0.01
    salt_mean = float(np.mean(salt))

    match = DATA_RE.search(path.name)
    assert match is not None
    metrics: dict[str, float | str] = {
        "target_ice_volume_ratio": target_ratio,
        "actual_ice_volume_ratio": actual_ratio,
        "absolute_match_error": abs(actual_ratio - target_ratio),
        "output_index": float(match.group(1)),
        "time": time,
        "finite_difference_melt_rate": finite_difference_rate,
        "reconstructed_melt_source_mean": float(np.mean(melt_source)),
        "melt_rate_reconstruction_ratio": float(np.mean(melt_source))
        / finite_difference_rate
        if finite_difference_rate != 0.0
        else math.nan,
        "positive_melt_source_mean": float(np.mean(np.maximum(melt_source, 0.0))),
        "negative_melt_source_mean": float(np.mean(np.minimum(melt_source, 0.0))),
        "interface_weighted_superheat": weighted_mean(superheat, band_weight),
        "interface_condition_residual_rms": math.sqrt(
            weighted_mean(superheat * superheat, band_weight)
        ),
        "interface_weighted_theta": weighted_mean(theta, band_weight),
        "interface_weighted_salt": weighted_mean(salt, band_weight),
        "interface_weighted_buoyancy": weighted_mean(buoyancy, band_weight),
        "interface_weighted_vertical_speed_abs": weighted_mean(np.abs(v), band_weight),
        "interface_weighted_vertical_speed_rms": math.sqrt(
            weighted_mean(v * v, band_weight)
        ),
        "interface_weighted_horizontal_speed_abs": weighted_mean(np.abs(u), band_weight),
        "liquid_side_heat_flux_proxy_mean": float(np.nanmean(heat_flux_proxy)),
        "liquid_side_salt_flux_proxy_mean": float(np.nanmean(salt_flux_proxy)),
        "front_valid_row_fraction": float(np.mean(valid_front)),
        "front_mean": float(np.nanmean(front)),
        "front_std": float(np.nanstd(front)),
        "front_length": front_length,
        "starved_interface_row_fraction": float(np.count_nonzero(starved))
        / max(np.count_nonzero(valid_rows), 1),
        "refreezing_interface_row_fraction": float(np.count_nonzero(refreezing))
        / max(np.count_nonzero(valid_rows), 1),
        "row_melt_coefficient_of_variation": float(np.std(row_melt[valid_rows]))
        / abs(float(np.mean(row_melt[valid_rows])))
        if np.any(valid_rows) and float(np.mean(row_melt[valid_rows])) != 0.0
        else math.nan,
        "salt_mean": salt_mean,
        "salt_relative_drift_from_initial": (salt_mean - initial_salt_mean)
        / abs(initial_salt_mean)
        if initial_salt_mean != 0.0
        else 0.0,
        "salt_solid_leakage_index": salt_solid_leakage_index,
        "deep_ice_salt_mean_abs": float(np.mean(np.abs(salt[deep_ice])))
        if np.any(deep_ice)
        else math.nan,
        "enthalpy_mean": float(np.mean(theta + phase / stefan)),
        "dx": dx,
        "dy": dy,
        "Pe_T": pe[0],
        "Pe_S": pe[1],
        "Cn_over_dx": cn / dx,
        "recovery_reference_cells": float(recovery_cells[-1]),
        **summarize_distribution("thermal_delta90_cells", delta_t90_cells),
        **summarize_distribution("salinity_delta90_cells", delta_s90_cells),
        **summarize_distribution("row_melt", row_melt[valid_rows]),
        **summarize_distribution("liquid_side_heat_flux_proxy", heat_flux_proxy),
        **summarize_distribution("liquid_side_salt_flux_proxy", salt_flux_proxy),
    }
    valid_front_count = max(np.count_nonzero(valid_front), 1)
    metrics["thermal_delta90_profile_fraction"] = float(
        np.count_nonzero(np.isfinite(delta_t90_cells) & valid_front)
    ) / valid_front_count
    metrics["salinity_delta90_profile_fraction"] = float(
        np.count_nonzero(np.isfinite(delta_s90_cells) & valid_front)
    ) / valid_front_count
    resolved_salt = np.isfinite(delta_s90_cells)
    metrics["salinity_delta90_resolved_fraction_ge_4dx"] = (
        float(np.mean(delta_s90_cells[resolved_salt] >= 4.0))
        if np.any(resolved_salt)
        else math.nan
    )

    vertical_rows: list[dict[str, float]] = []
    edges = np.linspace(float(y[0] - 0.5 * dy), float(y[-1] + 0.5 * dy), vertical_bins + 1)
    for bin_index in range(vertical_bins):
        in_bin = (
            (y >= edges[bin_index])
            & (y < edges[bin_index + 1] if bin_index < vertical_bins - 1 else y <= edges[bin_index + 1])
            & valid_rows
        )
        vertical_rows.append(
            {
                "target_ice_volume_ratio": target_ratio,
                "actual_ice_volume_ratio": actual_ratio,
                "time": time,
                "vertical_bin": float(bin_index),
                "y_low": edges[bin_index],
                "y_high": edges[bin_index + 1],
                "active_rows": float(np.count_nonzero(in_bin)),
                "front_mean": float(np.nanmean(front[in_bin])) if np.any(in_bin) else math.nan,
                "row_melt_mean": float(np.nanmean(row_melt[in_bin])) if np.any(in_bin) else math.nan,
                "row_melt_std": float(np.nanstd(row_melt[in_bin])) if np.any(in_bin) else math.nan,
                "interface_vertical_speed_abs_mean": float(np.nanmean(row_abs_v[in_bin]))
                if np.any(in_bin)
                else math.nan,
                "interface_buoyancy_mean": float(np.nanmean(row_buoyancy[in_bin]))
                if np.any(in_bin)
                else math.nan,
                "thermal_delta90_cells_median": finite_percentile(delta_t90_cells[in_bin], 50.0),
                "salinity_delta90_cells_median": finite_percentile(delta_s90_cells[in_bin], 50.0),
            }
        )
    return metrics, vertical_rows


def write_rows(path: Path, fieldnames: Sequence[str], rows: Iterable[dict]) -> None:
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def paired_comparison_rows(
    baseline_rows: Sequence[dict], comparison_rows: Sequence[dict]
) -> list[dict[str, float | str]]:
    """Build a long-form metric comparison at each requested volume target."""
    baseline_by_target = {
        float(row["target_ice_volume_ratio"]): row for row in baseline_rows
    }
    excluded = {
        "run",
        "directory",
        "input_file",
        "target_ice_volume_ratio",
        "output_index",
    }
    result: list[dict[str, float | str]] = []
    for comparison in comparison_rows:
        target = float(comparison["target_ice_volume_ratio"])
        if target not in baseline_by_target:
            continue
        baseline = baseline_by_target[target]
        metrics = sorted((set(baseline) & set(comparison)) - excluded)
        for metric in metrics:
            try:
                baseline_value = float(baseline[metric])
                comparison_value = float(comparison[metric])
            except (TypeError, ValueError):
                continue
            result.append(
                {
                    "target_ice_volume_ratio": target,
                    "baseline_run": baseline.get("run", "baseline"),
                    "comparison_run": comparison.get("run", "comparison"),
                    "metric": metric,
                    "baseline_value": baseline_value,
                    "comparison_value": comparison_value,
                    "comparison_minus_baseline": comparison_value - baseline_value,
                    "comparison_over_baseline": comparison_value / baseline_value
                    if baseline_value != 0.0
                    else math.nan,
                }
            )
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--run", action="append", type=parse_run, required=True, metavar="LABEL=PATH"
    )
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--prefix", default="yang_matched_transport")
    parser.add_argument(
        "--targets", type=parse_float_list, default=parse_float_list("0.8,0.7,0.6,0.5")
    )
    parser.add_argument(
        "--recovery-cells",
        type=parse_float_list,
        default=parse_float_list("2,4,8,16,32,64,128"),
        help="liquid-side sampling distances in grid cells",
    )
    parser.add_argument("--vertical-bins", type=int, default=12)
    parser.add_argument("--input-name", default="parties.inp")
    parser.add_argument(
        "--baseline-matched-csv",
        type=Path,
        help="optional prior matched CSV; writes a long-form paired comparison",
    )
    parser.add_argument(
        "--history-csv",
        action="append",
        type=parse_run,
        default=[],
        metavar="LABEL=PATH",
        help="saved volume history for a run when only selected Data HDFs are local",
    )
    parser.add_argument(
        "--initial-ice-volume",
        type=float,
        default=0.1,
        help="V0 for --history-csv data (default: exact configured Yang value 0.1)",
    )
    args = parser.parse_args()

    if args.vertical_bins < 1:
        parser.error("--vertical-bins must be positive")
    if args.initial_ice_volume <= 0.0:
        parser.error("--initial-ice-volume must be positive")
    recovery_cells = np.asarray(sorted(set(args.recovery_cells)), dtype=np.float64)
    if recovery_cells[0] <= 0.0:
        parser.error("--recovery-cells values must be positive")
    if any(not 0.0 < target <= 1.0 for target in args.targets):
        parser.error("--targets values must lie in (0,1]")
    args.output_dir.mkdir(parents=True, exist_ok=True)

    scan_rows: list[dict] = []
    matched_rows: list[dict] = []
    vertical_rows: list[dict] = []
    history_by_label = dict(args.history_csv)

    for label, directory in args.run:
        files = data_files(directory)
        if not files and label not in history_by_label:
            raise FileNotFoundError(f"no Data_*.h5 files in {directory}")
        input_path = directory / args.input_name
        if not input_path.exists():
            raise FileNotFoundError(f"missing input file {input_path}")
        parameters = parse_input(input_path)

        if label in history_by_label:
            run_scan = scan_history(
                history_by_label[label], directory, args.initial_ice_volume
            )
            initial_ice_volume = args.initial_ice_volume
        else:
            run_scan = [scan_snapshot(path) for path in files]
            initial_ice_volume = float(run_scan[0]["ice_volume"])
        run_scan.sort(key=lambda row: float(row["time"]))
        for row in run_scan:
            row["run"] = label
            row["directory"] = str(directory)
            row["ice_volume_ratio"] = float(row["ice_volume"]) / initial_ice_volume
            scan_rows.append(row)

        initial_path = directory / "Data_0.h5"
        if initial_path.exists():
            with h5py.File(initial_path, "r") as h5:
                _, ny, nx = physical_shape(h5)
                initial_salt_mean = float(
                    np.mean(scalar_plane(h5, "Conc/1", ny, nx))
                )
        else:
            initial_salt_mean = inferred_initial_salt_mean(parameters)

        times = np.asarray([float(row["time"]) for row in run_scan])
        volumes = np.asarray([float(row["ice_volume"]) for row in run_scan])
        ratios = volumes / initial_ice_volume
        selections: list[tuple[float, int]] = []
        for target in args.targets:
            selected = int(np.argmin(np.abs(ratios - target)))
            selections.append((target, selected))
        missing_paths = [
            Path(str(run_scan[selected]["path"]))
            for _, selected in selections
            if not Path(str(run_scan[selected]["path"])).exists()
        ]
        if missing_paths:
            names = ", ".join(path.name for path in missing_paths)
            raise FileNotFoundError(
                f"missing matched snapshots for {label}: {names}"
            )

        for target, selected in selections:
            lo = max(0, selected - 1)
            hi = min(len(run_scan) - 1, selected + 1)
            if lo == hi:
                finite_difference_rate = math.nan
            else:
                finite_difference_rate = -(volumes[hi] - volumes[lo]) / (
                    times[hi] - times[lo]
                )
            metrics, bins = selected_snapshot(
                Path(str(run_scan[selected]["path"])),
                parameters,
                target,
                float(ratios[selected]),
                finite_difference_rate,
                initial_salt_mean,
                recovery_cells,
                args.vertical_bins,
            )
            matched_rows.append(
                {
                    "run": label,
                    "directory": str(directory),
                    "input_file": str(input_path),
                    **metrics,
                }
            )
            vertical_rows.extend({"run": label, **row} for row in bins)
            print(
                f"{label}: target V/V0={target:.3f}, selected "
                f"{ratios[selected]:.6f} at t={times[selected]:.6f} "
                f"({Path(str(run_scan[selected]['path'])).name})"
            )

    scan_fields = [
        "run", "directory", "output_index", "path", "time", "ice_volume", "ice_volume_ratio"
    ]
    matched_fields = list(matched_rows[0])
    vertical_fields = list(vertical_rows[0])
    write_rows(args.output_dir / f"{args.prefix}_scan.csv", scan_fields, scan_rows)
    write_rows(args.output_dir / f"{args.prefix}_matched.csv", matched_fields, matched_rows)
    write_rows(
        args.output_dir / f"{args.prefix}_vertical_bins.csv", vertical_fields, vertical_rows
    )
    if args.baseline_matched_csv:
        baseline_rows = read_csv(args.baseline_matched_csv.resolve())
        paired_rows = paired_comparison_rows(baseline_rows, matched_rows)
        if not paired_rows:
            raise ValueError(
                "baseline/current matched CSVs have no common ice-volume targets"
            )
        write_rows(
            args.output_dir / f"{args.prefix}_paired.csv",
            list(paired_rows[0]),
            paired_rows,
        )
    print(f"Wrote {args.output_dir.resolve() / (args.prefix + '_*.csv')}")


if __name__ == "__main__":
    main()
