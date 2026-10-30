#!/usr/bin/env python
"""Normalized sensitivity of the 1C discharge to model parameters.

For each parameter p the script runs dfn-sim at p*(1-delta) and p*(1+delta),
and reports

    S = (result(+delta) - result(-delta)) / (2 * delta * result(0))

for the discharged capacity and the cutoff time. A value of +1 means the
result changes by the same relative amount as the parameter; 0 means the
result is insensitive to it.

Example:
    .venv/Scripts/python.exe scripts/sensitivity.py --binary build-klu/Release/dfn-sim
"""

import argparse
import csv
import os
import re
import subprocess
import sys
import tempfile


TARGETS = [
    ("negative", "d_s"),
    ("positive", "d_s"),
    ("negative", "k"),
    ("positive", "k"),
    ("electrolyte", "d_e"),
    ("electrolyte", "t_plus"),
    ("negative", "thickness"),
    ("positive", "thickness"),
    ("negative", "porosity"),
]


def parse_parameters(path):
    sections = {}
    current = None
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            stripped = line.split("#", 1)[0].rstrip()
            if not stripped.strip():
                continue
            if ":" not in stripped:
                continue
            key, value = stripped.split(":", 1)
            key = key.strip()
            value = value.strip()
            if not value:
                current = key
                sections[current] = {}
            elif current is not None:
                sections[current][key] = float(value)
    return sections


def write_parameters(path, sections):
    lines = []
    for section, values in sections.items():
        lines.append(f"{section}:")
        for key, value in values.items():
            lines.append(f"  {key}: {value:.12g}")
        lines.append("")
    with open(path, "w", encoding="utf-8") as handle:
        handle.write("\n".join(lines))


def run_case(binary, config_path, workdir, label, crate=1.0):
    csv_path = os.path.join(workdir, f"{label}.csv")
    completed = subprocess.run(
        [
            binary,
            "--config",
            config_path,
            "--mode",
            "cc",
            "--crate",
            f"{crate:.6g}",
            "--dt",
            "600",
            "--out",
            csv_path,
        ],
        capture_output=True,
        text=True,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"{label}: {completed.stderr.strip()}")
    summary = completed.stderr
    capacity = re.search(r"cap=([0-9.eE+-]+) Ah/m\^2", summary)
    t_end = re.search(r"t_end=([0-9.eE+-]+) s", summary)
    if capacity is None or t_end is None:
        raise RuntimeError(f"{label}: cannot parse summary '{summary.strip()}'")

    voltage = None
    with open(csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            if abs(float(row["t_s"]) - 600.0) < 1e-6:
                voltage = float(row["voltage_v"])
                break
    if voltage is None:
        raise RuntimeError(f"{label}: no sample at t = 600 s")

    return float(capacity.group(1)), float(t_end.group(1)), voltage


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="build-klu/Release/dfn-sim")
    parser.add_argument("--config", default="configs/demo.yaml")
    parser.add_argument("--delta", type=float, default=0.10)
    parser.add_argument("--crate", type=float, default=1.0)
    args = parser.parse_args()

    base = parse_parameters(args.config)
    with tempfile.TemporaryDirectory() as workdir:
        base_path = os.path.join(workdir, "base.yaml")
        write_parameters(base_path, base)
        q0, t0, v0 = run_case(args.binary, base_path, workdir, "base", args.crate)
        print(
            f"reference: cap = {q0:.4f} Ah/m^2, t_cut = {t0:.1f} s, "
            f"V(600 s) = {v0:.5f} V, delta = {args.delta:.0%}"
        )
        print()
        print("parameter\tS(cap)\tS(t_cut)\tS(V600)")

        results = []
        for section, key in TARGETS:
            if key not in base.get(section, {}):
                print(f"skip {section}.{key} (not in the config)", file=sys.stderr)
                continue
            original = base[section][key]
            sensitivities = []
            for sign in (-1.0, 1.0):
                perturbed = {s: dict(v) for s, v in base.items()}
                perturbed[section][key] = original * (1.0 + sign * args.delta)
                path = os.path.join(workdir, f"{section}_{key}_{sign:+.0f}.yaml")
                write_parameters(path, perturbed)
                label = f"{section}_{key}_{sign:+.0f}"
                sensitivities.append(
                    run_case(args.binary, path, workdir, label, args.crate)
                )
            (q_minus, t_minus, v_minus), (q_plus, t_plus, v_plus) = sensitivities
            s_cap = (q_plus - q_minus) / (2.0 * args.delta * q0)
            s_t = (t_plus - t_minus) / (2.0 * args.delta * t0)
            s_v = (v_plus - v_minus) / (2.0 * args.delta * v0)
            results.append((section, key, s_cap, s_t, s_v))

        results.sort(key=lambda item: -abs(item[2]))
        for section, key, s_cap, s_t, s_v in results:
            print(f"{section}.{key}\t{s_cap:+.3f}\t{s_t:+.3f}\t{s_v:+.4f}")


if __name__ == "__main__":
    raise SystemExit(main())
