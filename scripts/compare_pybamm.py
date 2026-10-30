#!/usr/bin/env python
"""Compare a dfn-sim discharge curve with PyBaMM's DFN model.

Both solvers use the parameters of configs/demo.yaml, the same OCP fits and
the same isothermal, constant-property assumptions. Remaining differences are
the concentrated-solution term (PyBaMM keeps kappa_D, dfn-sim sets it to
zero), the discretisation and the time tolerances.

Example:
    .venv/Scripts/python.exe scripts/compare_pybamm.py results/cc60.csv
"""

import argparse
import csv
import sys


def build_parameter_values(area, current_density):
    import pybamm

    faraday = pybamm.constants.F

    def neg_ocp(sto):
        sto = pybamm.maximum(sto, 1e-6)
        return 0.05 + 0.8 * pybamm.exp(-10.0 * sto) + 0.1 * (1.0 - sto) ** 2

    def pos_ocp(sto):
        return 3.9 + 0.2 * (1.0 - sto)

    def exchange_current(ce, cs, cmax, temperature):
        ce = pybamm.maximum(ce, 1e-3)
        cs = pybamm.maximum(pybamm.minimum(cs, cmax - 1e-3), 1e-3)
        return (
            faraday
            * 1.0e-6
            * pybamm.sqrt(ce)
            * pybamm.sqrt(cs)
            * pybamm.sqrt(cmax - cs)
        )

    values = pybamm.ParameterValues("Chen2020")
    values.update(
        {
            "Electrode height [m]": 0.1,
            "Electrode width [m]": area / 0.1,
            "Nominal cell capacity [A.h]": 0.595,
            "Current function [A]": current_density * area,
            "Ambient temperature [K]": 298.15,
            "Initial temperature [K]": 298.15,
            "Lower voltage cut-off [V]": 3.5,
            "Upper voltage cut-off [V]": 4.2,
            "Negative electrode thickness [m]": 85e-6,
            "Positive electrode thickness [m]": 75e-6,
            "Separator thickness [m]": 12e-6,
            "Negative particle radius [m]": 5e-6,
            "Positive particle radius [m]": 5e-6,
            "Negative electrode porosity": 0.25,
            "Positive electrode porosity": 0.33,
            "Separator porosity": 0.47,
            "Negative electrode Bruggeman coefficient (electrode)": 1.5,
            "Positive electrode Bruggeman coefficient (electrode)": 1.5,
            "Negative electrode conductivity [S.m-1]": 100.0,
            "Positive electrode conductivity [S.m-1]": 10.0,
            "Maximum concentration in negative electrode [mol.m-3]": 30000.0,
            "Maximum concentration in positive electrode [mol.m-3]": 60000.0,
            "Initial concentration in negative electrode [mol.m-3]": 24000.0,
            "Initial concentration in positive electrode [mol.m-3]": 18000.0,
            "Initial concentration in electrolyte [mol.m-3]": 1000.0,
            "Negative particle diffusivity [m2.s-1]": 3e-14,
            "Positive particle diffusivity [m2.s-1]": 4e-15,
            "Electrolyte diffusivity [m2.s-1]": 5e-10,
            "Electrolyte conductivity [S.m-1]": 1.0,
            "Cation transference number": 0.26,
            "Negative electrode OCP [V]": neg_ocp,
            "Positive electrode OCP [V]": pos_ocp,
            "Negative electrode exchange-current density [A.m-2]": exchange_current,
            "Positive electrode exchange-current density [A.m-2]": exchange_current,
        }
    )
    return values


def solve_reference(area, current_density, t_end, points):
    import numpy as np
    import pybamm

    model = pybamm.lithium_ion.DFN()
    simulation = pybamm.Simulation(
        model, parameter_values=build_parameter_values(area, current_density)
    )
    solution = simulation.solve(t_eval=np.linspace(0.0, t_end, points))
    return (
        solution["Time [s]"].entries,
        solution["Voltage [V]"].entries,
    )


def load_curve(path):
    times, voltages, currents = [], [], []
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            times.append(float(row["t_s"]))
            voltages.append(float(row["voltage_v"]))
            currents.append(float(row["current_a_m2"]))
    return times, voltages, currents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv_path", help="CSV produced by dfn-sim --out")
    parser.add_argument("--area", type=float, default=0.0165)
    parser.add_argument(
        "--current-density",
        type=float,
        default=None,
        help="A/m^2; default: read from the CSV",
    )
    parser.add_argument("--points", type=int, default=801)
    parser.add_argument("--save", help="save a comparison plot to this file")
    args = parser.parse_args()

    try:
        import numpy as np
    except ImportError:
        print("error: numpy is required (pip install numpy)", file=sys.stderr)
        return 1

    ours_t, ours_v, ours_i = load_curve(args.csv_path)
    current_density = (
        args.current_density
        if args.current_density is not None
        else abs(ours_i[0])
    )
    ref_t, ref_v = solve_reference(
        args.area, current_density, ours_t[-1] + 1.0, args.points
    )

    limit = min(ours_t[-1], ref_t[-1])
    mask = np.asarray(ours_t) <= limit
    t = np.asarray(ours_t)[mask]
    v = np.asarray(ours_v)[mask]
    v_ref = np.interp(t, ref_t, ref_v)
    error = v - v_ref

    rmse = float(np.sqrt(np.mean(error**2)))
    worst = float(np.max(np.abs(error)))
    worst_t = float(t[int(np.argmax(np.abs(error)))])

    print(f"dfn-sim samples       : {len(ours_t)}")
    print(f"dfn-sim t_end         : {ours_t[-1]:.1f} s")
    print(f"PyBaMM  t_end         : {ref_t[-1]:.1f} s")
    print(f"compared window       : 0 .. {limit:.1f} s ({len(t)} points)")
    print(f"V(0)  dfn-sim/PyBaMM  : {v[0]:.6f} / {v_ref[0]:.6f} V")
    print(f"RMSE                  : {1000.0 * rmse:.3f} mV")
    print(f"max |dV|              : {1000.0 * worst:.3f} mV at t = {worst_t:.0f} s")

    if args.save:
        import matplotlib.pyplot as plt

        fig, (ax_v, ax_e) = plt.subplots(
            2, 1, figsize=(7, 5), sharex=True, gridspec_kw={"height_ratios": [3, 1]}
        )
        ax_v.plot(ours_t, ours_v, label="dfn-sim")
        ax_v.plot(ref_t, ref_v, "--", label="PyBaMM DFN")
        ax_v.set_ylabel("voltage [V]")
        ax_v.grid(alpha=0.3)
        ax_v.legend()
        ax_e.plot(t, 1000.0 * error, color="tab:red")
        ax_e.set_ylabel("dfn-sim - PyBaMM [mV]")
        ax_e.set_xlabel("time [s]")
        ax_e.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(args.save, dpi=150)
        print(f"plot written to {args.save}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
