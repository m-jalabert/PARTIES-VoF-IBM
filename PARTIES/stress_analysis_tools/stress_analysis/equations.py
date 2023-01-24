"""Stress analysis equations."""
from datetime import date
from typing import Dict, Union

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid

from stress_analysis import data_ingestion as di


def stress_x(sim: di.Simulation, idx: int) -> pd.DataFrame:
    """Horizontally averaged x-direction stress balance.

    Dimensionless formulation of equation 7.12 in E. Biegert's 2018 PhD
    disertation (p. 110). Local ibm output is needed. All quantities are
    horizontally averaged and reported with respect to `yv` coordinates,
    defined at the cell edges in y.

    Args:
        sim: A parties simulation.
        idx: Output file index. e.g. `Data_{idx}.h5`.

    Returns:
        The horizontally averaged terms of equation 7.12 in E. Biegert's
        disertation. The profiles are reported with respect to the staggered y
        coordinates defined at the vertical cell edges.
    """
    stress_params = {
        "Re": sim.Re,
        "x_1": add_ghost_coords(sim.fluid.yc, sim.h),
        "x_2": sim.fluid.xc,
        "u_1": interpolate_xu_to_xc(sim.fluid.u(idx)),
        "u_2": interpolate_yv_to_yc(sim.fluid.v(idx)),
        "ax_1": 1,
        "ax_2": 2,
        "ind": sim.fluid.vfv(idx),
    }
    stress_params["u_1"] = add_ghost_nodes(
        stress_params["u_1"], top=sim.ubulk, bot=-sim.ubulk
    )
    stress_params["u_2"] = add_ghost_nodes(stress_params["u_2"])

    # Clculate stresses
    stress = {"idx": idx, "yv": sim.fluid.yv}
    stress["visc_f"], stress["visc_p"] = visc_stress(stress_params)
    stress["adv_f"], stress["adv_p"] = adv_stress(stress_params)
    stress_params["u_1"] = stress_params["u_2"] = None

    stress_params["f_ibm"] = interpolate_xu_to_xc(sim.fluid.fx_ibm(idx))
    stress["ibm"] = ibm_stress(stress_params)
    stress_params["f_ibm"] = None

    stress["fluid"] = stress["visc_f"] + stress["adv_f"]
    stress["particle"] = stress["visc_p"] + stress["adv_p"] + stress["ibm"]
    stress["total"] = stress["fluid"] + stress["particle"]

    stress["external"] = stress["visc_f"][-1] + stress["visc_p"][-1]
    stress["external"] += sim.dp_dx * sim.Re * (sim.Ly - stress["yv"])

    return pd.DataFrame(stress)


def stress_y(sim: di.Simulation, idx: int) -> pd.DataFrame:
    """Horizontally averaged stress balance in the y-direction.

    Dimensionless formulation of equation 7.14 in E. Biegert's 2018 PhD
    disertation (p. 110). Local ibm output is needed. All quantities are
    horizontally averaged and reported with respect to `yv` coordinates, defined
    at the cell edges in y.

    Args:
        sim: A parties simulation.
        idx: Output file index. e.g. `Data_{idx}.h5`.

    Returns:
        The horizontally averaged terms of equation 7.14 in E. Biegert's
        disertation. The profiles are reported with respect to the staggered y
        coordinates defined at the vertical cell edges.
    """
    stress_params = {
        "Re": sim.Re,
        "x_1": add_ghost_coords(sim.fluid.yc, sim.h),
        "u_1": interpolate_yv_to_yc(sim.fluid.v(idx)),
        "ax_1": 1,
        "ax_2": 1,
        "ind": sim.fluid.vfv(idx),
    }
    stress_params["x_2"] = stress_params["x_1"]
    stress_params["u_1"] = add_ghost_nodes(stress_params["u_1"])
    stress_params["u_2"] = stress_params["u_1"]

    # Calculate stresses
    stress = {"idx": idx, "yv": sim.fluid.yv}
    stress["visc_f"], stress["visc_p"] = visc_stress(stress_params)
    stress["adv_f"], stress["adv_p"] = adv_stress(stress_params)
    stress_params["u_1"] = stress_params["u_2"] = None

    stress_params["f_ibm"] = interpolate_yv_to_yc(sim.fluid.fy_ibm(idx))
    stress["ibm"] = ibm_stress(stress_params)
    stress_params["f_ibm"] = None

    stress_params["p"] = add_set_bounds(sim.fluid.p(idx), top=1, bot=1)
    stress["press_f"], stress["press_p"] = pressure_stress(stress_params)
    stress_params["p"] = None

    top_wall_visc = stress["visc_f"][-1] + stress["visc_p"][-1]
    top_wall_press = stress["press_f"][-1] + stress["press_p"][-1]
    stress["external"] = np.full(sim.Ny + 1, (top_wall_visc + top_wall_press))

    stress["fluid"] = stress["visc_f"] + stress["adv_f"] + stress["press_f"]
    stress["particle"] = (
        stress["visc_p"] + stress["adv_p"] + stress["ibm"] + stress["press_p"]
    )
    stress["total"] = stress["fluid"] + stress["particle"]
    return pd.DataFrame(stress)


def particle_pressure(sim: di.Simulation, idx: int) -> pd.DataFrame:
    """Compute particle pressure on `walls` made of fixed partilcles.

    Args:
        sim: A parties simulation.
        idx: Output file index. e.g. `Particle_{idx}.h5`.
    Returns:
        The particle pressure reported on the top and bottom fixed particle
        `walls`. The average of the two is also given.
    """
    fy_coll = sim.p_fixed.F_coll(idx)[:, 1]
    y_part = sim.p_fixed.X(idx)[:, 1]
    area = sim.Lx * sim.Lz
    half_height = sim.Ly / 2
    part_press = {
        "idx": idx,
        "top": fy_coll[y_part > half_height].sum() / area,
        "bot": fy_coll[y_part < half_height].sum() / area,
    }
    part_press["avg"] = (part_press["top"] + part_press["bot"]) / 2
    return pd.DataFrame(part_press, index=[0])


# Stress term calculations
def visc_stress(
    stress_params: Dict[str, Union[np.ndarray, int]]
) -> np.ndarray:
    """Compute horizontally avaraged viscous stresses."""
    du1dx1 = np.gradient(
        stress_params["u_1"],
        stress_params["x_1"],
        axis=stress_params["ax_1"],
        edge_order=2,
    )
    du2dx2 = np.gradient(
        stress_params["u_2"],
        stress_params["x_2"],
        axis=stress_params["ax_2"],
        edge_order=2,
    )
    visc = du1dx1 + du2dx2
    visc = interpolate_yv_to_yc(visc)
    visc_f = plane_average(visc, stress_params["ind"])
    visc_p = plane_average(visc, 1 - stress_params["ind"])
    return np.array([visc_f, visc_p])


def adv_stress(stress_params: Dict[str, Union[np.ndarray, int]]) -> np.ndarray:
    """Compute horizontally avaraged advective stresses."""
    adv = stress_params["u_1"] * stress_params["u_2"]
    adv = interpolate_yv_to_yc(adv)
    adv_f = plane_average(adv, stress_params["ind"])
    adv_p = plane_average(adv, 1 - stress_params["ind"])
    return -stress_params["Re"] * np.array([adv_f, adv_p])


def ibm_stress(stress_params: Dict[str, Union[np.ndarray, int]]) -> np.ndarray:
    """Compute horizontally averaged stress due to local ibm forcing."""
    ibm = stress_params["f_ibm"].mean(axis=(0, 2))
    ibm = np.append(ibm, 0)
    ibm = np.insert(ibm, 0, 0)
    ibm = integrate_backwards(ibm, stress_params["x_1"])
    ibm = (ibm[1:] + ibm[:-1]) / 2
    return -stress_params["Re"] * ibm


def pressure_stress(
    stress_params: Dict[str, Union[np.ndarray, int]]
) -> np.ndarray:
    """Compute horizontally averaged stresses due to pressure."""
    press = interpolate_yv_to_yc(stress_params["p"])
    press_f = plane_average(press, stress_params["ind"])
    press_p = plane_average(press, 1 - stress_params["ind"])
    return -stress_params["Re"] * np.array([press_f, press_p])


def plane_average(quantity: np.ndarray, indicator: np.ndarray) -> np.ndarray:
    """Conditionally plane averge a quantity using an indicator function."""
    average_quantity = quantity * (1 - indicator)
    average_quantity = average_quantity.mean(axis=(0, 2))
    return average_quantity


def integrate_backwards(quantity: np.ndarray, coord: np.ndarray) -> np.ndarray:
    """Cumulatively integrate quantity from last to first value."""
    return np.flip(cumulative_trapezoid(np.flip(quantity), coord, initial=0))


def interpolate_xu_to_xc(field: np.ndarray) -> np.ndarray:
    """Interpolate between stagered `u` to cell centered coordinates."""
    return (field[:, :, 1:] + field[:, :, :-1]) / 2


def interpolate_yv_to_yc(field: np.ndarray) -> np.ndarray:
    """Interpolate betwen stagered `v` to cell centered coordinates."""
    return (field[:, 1:, :] + field[:, :-1, :]) / 2


def add_ghost_coords(coords: np.ndarray, h: float) -> np.ndarray:
    """Add ghost node coordinates to cordinate array."""
    coords = np.append(coords, coords[-1] + h)
    coords = np.insert(coords, 0, coords[0] - h)
    return coords


def add_ghost_nodes(
    field: np.ndarray, top: float = 0.0, bot: float = 0.0
) -> np.ndarray:
    """Add ghost gridpoints to satisfy boundary conditions at domain walls."""
    field = np.append(field, -field[:, -1:, :] + 2 * top, axis=1)
    field = np.insert(field, [0], -field[:, :1, :] + 2 * bot, axis=1)
    return field


def add_set_bounds(
    field: np.ndarray, top: float = 0, bot: float = 0
) -> np.ndarray:
    """Add specified quantities to the top and bottom of the domain in y."""
    field = np.append(field, top * field[:, -1:, :], axis=1)
    field = np.insert(field, [0], bot * field[:, 0:1, :], axis=1)
    return field


def save_stress(profiles: pd.DataFrame, component: str) -> None:
    """Save stress profiles to file.

    Args:
        profiles: Profiles for each term in the stress balance.
        component: The direction of the stress balance.
    """
    filename = f"stress_{component}_{profiles['idx'][0]}.csv"
    desc_short = f"Profiles of {component}-stress components"
    desc_long = (
        f"Horizontally averaged profiles of the {component}-stress "
        "component terms in equation 7.12/7.14 in E. Biegert's 2018 "
        "PhD Disertation. Reported quantities are normalized as in "
        "the simulation output."
    )
    save_to_file(profiles, filename, desc_short, desc_long=desc_long)


def save_to_file(
    dataframe: pd.DataFrame,
    filename: str,
    desc_short: str,
    desc_long: str = "",
) -> None:
    """Save DataFrame to file with formatted descriptions."""
    today = date.today()
    desc_long = split_input(desc_long, 78)
    with open(filename, "w", encoding="utf8") as file:
        file.write(f"% Filename: {filename}\n%\n")
        file.write(f"% Description: {desc_short}\n%\n")
        for des in desc_long:
            file.write(f"% {des}\n")
        file.write(f"%\n% Last modified: {today.strftime(r'%d/%m/%Y')}\n%\n")
        file.write("%" + "-" * 79 + "\n")
    cols = dataframe.keys()[dataframe.keys() != "idx"]
    dataframe.to_csv(filename, mode="a", columns=cols)


def split_input(user_string: str, line_length: int) -> int:
    """Format string to lines respecting the specified line length."""
    output = []
    words = user_string.split(" ")
    total_length = 0
    while total_length < len(user_string) and len(words) > 0:
        line = []
        next_word = words[0]
        line_len = len(next_word) + 1
        while (line_len < line_length) and len(words) > 0:
            words.pop(0)
            line.append(next_word)
            if len(words) > 0:
                next_word = words[0]
                line_len += len(next_word) + 1
        line = " ".join(line)
        output.append(line)
        total_length += len(line)

    return output
