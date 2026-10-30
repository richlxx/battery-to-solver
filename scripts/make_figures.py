#!/usr/bin/env python
"""Regenerate the figures under book/assets from the built executables.

Run from the repository root, after a KLU build:

    .venv/Scripts/python.exe scripts/make_figures.py
"""

import argparse
import csv
import pathlib
import re
import subprocess
import sys
import tempfile

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import (
    Circle,
    FancyArrowPatch,
    FancyBboxPatch,
    Patch,
    Rectangle,
)

ROOT = pathlib.Path(__file__).resolve().parents[1]
ASSETS = ROOT / "book" / "assets"
PREVIEW = False

plt.rcParams.update(
    {
        "font.family": "Microsoft YaHei",
        "axes.unicode_minus": False,
        "svg.fonttype": "path",
        "figure.dpi": 110,
        "savefig.bbox": "tight",
        "axes.edgecolor": "#5b616b",
        "axes.labelcolor": "#30343b",
        "text.color": "#30343b",
        "xtick.color": "#5b616b",
        "ytick.color": "#5b616b",
    }
)

BLUE = "#3b6ea5"
ORANGE = "#d9822b"
GREEN = "#4f8a5b"
GRAY = "#8a8f98"
INK = "#30343b"
FILL = "#eef2f7"
PALE_ORANGE = "#f8e7d2"
PALE_BLUE = "#dfe9f5"
PALE_GRAY = "#e3e5e8"


def exe(name):
    path = ROOT / "build-klu" / "Release" / (name + ".exe")
    if not path.exists():
        path = ROOT / "build-klu" / "Release" / name
    if not path.exists():
        raise SystemExit(f"missing executable {path}; build with DFN_USE_KLU=ON")
    return str(path)


def cuda_exe(name):
    path = ROOT / "build-cuda" / "Release" / (name + ".exe")
    if not path.exists():
        path = ROOT / "build-cuda" / "Release" / name
    if not path.exists():
        raise SystemExit(f"missing executable {path}; build with DFN_USE_CUDA=ON")
    return str(path)


def run(command):
    completed = subprocess.run(
        command, capture_output=True, text=True, cwd=ROOT, encoding="utf-8"
    )
    if completed.returncode != 0:
        raise SystemExit(
            f"command failed: {' '.join(command)}\n{completed.stderr.strip()}"
        )
    return completed.stdout, completed.stderr


def save(fig, name):
    ASSETS.mkdir(parents=True, exist_ok=True)
    path = ASSETS / name
    fig.savefig(path)
    if PREVIEW:
        fig.savefig(path.with_suffix(".png"), dpi=110)
    plt.close(fig)
    print(f"wrote {path.relative_to(ROOT)}")


def read_csv(text):
    return list(csv.DictReader(text.splitlines()))


def box(ax, x, y, w, h, text, fc=FILL, fs=10.5):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.2,
        edgecolor="#5b616b",
        facecolor=fc,
    )
    ax.add_patch(patch)
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs)


def arrow(ax, start, end, color=INK, style="-|>", ls="-"):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle=style,
            mutation_scale=13,
            linewidth=1.2,
            color=color,
            linestyle=ls,
            shrinkA=0,
            shrinkB=0,
        )
    )


def ocv_figure():
    stdout, _ = run([exe("ocv_curve"), "--n", "101"])
    rows = read_csv(stdout)
    soc = np.array([float(r["soc"]) for r in rows])
    voltage = np.array([float(r["ocv_v"]) for r in rows])

    fig, ax = plt.subplots(figsize=(6.4, 3.8))
    ax.plot(soc, voltage, color=BLUE, lw=2)
    ax.set_xlabel("SOC（荷电状态）")
    ax.set_ylabel("开路电压 OCV [V]")
    ax.grid(alpha=0.3)
    ax.annotate(
        "锂越少，负极电位越高",
        xy=(0.08, voltage[8]),
        xytext=(0.18, 3.72),
        arrowprops={"arrowstyle": "->", "color": GRAY},
        fontsize=9,
        color=INK,
    )
    save(fig, "fig00-ocv.svg")


def tangent_figure():
    x = np.linspace(0.6, 2.3, 300)
    f = x * x
    x0, x1 = 1.0, 1.8
    secant = 1.0 + (x1 + x0) * (x - x0)
    tangent = 1.0 + 2.0 * (x - x0)

    fig, ax = plt.subplots(figsize=(6.0, 4.0))
    ax.plot(x, f, color=BLUE, lw=2, label="$f(x)=x^2$")
    ax.plot(x, secant, color=GRAY, ls="--", lw=1.5, label="割线（差商）")
    ax.plot(x, tangent, color=ORANGE, lw=1.8, label="切线（导数）")
    ax.plot([x0], [x0 * x0], "o", color=INK)
    ax.plot([x1], [x1 * x1], "o", color=INK)
    ax.annotate("$x_0$", (x0, 0), xytext=(x0, -0.35), ha="center", fontsize=10)
    ax.annotate("$x_0+h$", (x1, 0), xytext=(x1, -0.35), ha="center", fontsize=10)
    ax.annotate(
        "斜率 $f'(x_0)=2x_0$",
        xy=(1.85, 2.7),
        xytext=(1.52, 4.0),
        arrowprops={"arrowstyle": "->", "color": ORANGE},
        color=ORANGE,
        fontsize=10,
    )
    ax.set_xlim(0.6, 2.3)
    ax.set_ylim(-0.4, 5.2)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$f(x)$")
    ax.grid(alpha=0.3)
    ax.legend(loc="upper left", frameon=False)
    save(fig, "fig01-tangent.svg")


def error_figure():
    stdout, _ = run([exe("diff_error")])
    rows = read_csv(stdout)
    h = np.array([float(r["h"]) for r in rows])
    forward = np.array([float(r["forward_err"]) for r in rows])
    central = np.array([float(r["central_err"]) for r in rows])

    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.loglog(h, forward, "o-", ms=3, color=BLUE, label="前向差分（一阶）")
    ax.loglog(h, central, "s-", ms=3, color=ORANGE, label="中心差分（二阶）")
    ax.loglog(h, 5e-2 * h, color=BLUE, ls=":", lw=1)
    ax.loglog(h, 5e-2 * h**2, color=ORANGE, ls=":", lw=1)
    ax.set_xlabel("步长 $h$")
    ax.set_ylabel("绝对误差")
    ax.set_xlim(1e-1, 1e-15)
    ax.set_ylim(1e-17, 1e-1)
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False)
    save(fig, "fig02-error.svg")


def dfn_figure():
    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(9.4, 3.9), gridspec_kw={"width_ratios": [1.35, 1]}
    )

    layers = [
        (0.0, 0.45, PALE_GRAY, "负极\n集流体"),
        (0.45, 3.05, PALE_BLUE, "负极"),
        (3.05, 4.25, "#f6f2df", "隔膜"),
        (4.25, 6.85, PALE_ORANGE, "正极"),
        (6.85, 7.30, PALE_GRAY, "正极\n集流体"),
    ]
    for left, right, color, label in layers:
        ax.add_patch(
            Rectangle((left, 0.0), right - left, 2.0, facecolor=color,
                      edgecolor="#8a8f98", lw=0.8)
        )
        ax.text(
            (left + right) / 2,
            1.0 if label != "隔膜" else 0.35,
            label,
            ha="center",
            va="center",
            fontsize=9.5,
        )

    rng = np.random.default_rng(3)
    for left, right in [(0.55, 3.0), (4.3, 6.8)]:
        for _ in range(26):
            cx = rng.uniform(left + 0.12, right - 0.12)
            cy = rng.uniform(0.2, 1.8)
            ax.add_patch(Circle((cx, cy), 0.09, facecolor="white",
                                edgecolor=GRAY, lw=0.7, alpha=0.95))

    for xpos in (1.4, 2.2, 3.6, 5.0, 5.9):
        arrow(ax, (xpos, 2.2), (xpos + 0.55, 2.2), color=ORANGE)
    ax.text(3.9, 2.55, "放电：Li$^+$ 从负极迁移到正极", ha="center",
            fontsize=9.5, color=ORANGE)

    ax.text(3.65, -0.3, "$x$（穿过电极、隔膜）", ha="center", fontsize=9.5)
    ax.text(0.05, 2.12, "$\\phi_s=0$", fontsize=9)
    ax.text(7.05, 2.12, "$\\phi_s=V$", fontsize=9)
    ax.set_xlim(-0.15, 7.45)
    ax.set_ylim(-0.5, 2.7)
    ax.set_axis_off()

    ax2.add_patch(Circle((0.5, 0.5), 0.42, facecolor=PALE_BLUE,
                         edgecolor=BLUE, lw=1.4))
    for radius in (0.14, 0.28):
        ax2.add_patch(Circle((0.5, 0.5), radius, facecolor="none",
                             edgecolor=BLUE, lw=0.8, ls=":"))
    ax2.add_patch(Circle((0.5, 0.5), 0.02, facecolor=BLUE))
    arrow(ax2, (0.5, 0.5), (0.92, 0.5), color=INK)
    ax2.text(0.71, 0.57, "$r$", fontsize=11)
    ax2.text(1.02, 0.5, "$R$", fontsize=10, va="center")
    arrow(ax2, (1.03, 0.82), (1.03, 1.02), color=ORANGE)
    ax2.text(1.08, 0.95, "$j$（反应通量）", fontsize=10, color=ORANGE)
    ax2.text(0.5, 0.0, "$-D_s\\,\\partial c_s/\\partial r = j$",
             ha="center", fontsize=10)
    ax2.text(0.5, -0.12, "颗粒半径约 5 μm", ha="center", fontsize=9,
             color=GRAY)
    ax2.set_xlim(0.0, 1.75)
    ax2.set_ylim(-0.28, 1.15)
    ax2.set_aspect("equal")
    ax2.set_axis_off()

    save(fig, "fig03-dfn.svg")


def grid_figure():
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.4))

    for i in range(5):
        ax.add_patch(Rectangle((i, 0.0), 1.0, 1.0, facecolor=PALE_BLUE,
                               edgecolor=BLUE, lw=1.1))
        ax.plot([i + 0.5], [0.5], "o", color=ORANGE, ms=6)
        ax.text(i + 0.5, -0.22, f"$i={i}$", ha="center", fontsize=9)
    for i in range(6):
        ax.plot([i, i], [-0.05, 1.05], color=GRAY, lw=0.8)
    ax.annotate("", xy=(3.0, 1.25), xytext=(2.0, 1.25),
                arrowprops={"arrowstyle": "<->", "color": GRAY})
    ax.text(2.5, 1.32, "控制体 $V_i$", ha="center", fontsize=9.5)
    ax.annotate("", xy=(2.0, -0.52), xytext=(3.0, -0.52),
                arrowprops={"arrowstyle": "<->", "color": ORANGE})
    ax.text(2.5, -0.75, "界面通量 $A_i F_i$", ha="center", fontsize=9.5,
            color=ORANGE)
    ax.set_xlim(-0.2, 5.2)
    ax.set_ylim(-0.9, 1.45)
    ax.set_axis_off()
    ax.set_title("一维网格：格子中心与界面")

    for radius, color in [(1.0, PALE_GRAY), (0.72, PALE_BLUE),
                          (0.45, "#cfe0f2")]:
        ax2.add_patch(Circle((0, 0), radius, facecolor=color,
                             edgecolor=BLUE, lw=1.0))
    ax2.plot([0], [0], "o", color=ORANGE, ms=5)
    arrow(ax2, (0, 0), (1.0, 0), color=INK)
    ax2.text(0.5, 0.08, "$r$", fontsize=10)
    ax2.text(1.05, 0.0, "$R$", fontsize=10, va="center")
    ax2.text(0.0, -1.22, "球壳：$A_i=4\\pi r_i^2$，$V_i=\\frac{4\\pi}{3}(r_{i+1}^3-r_i^3)$",
             ha="center", fontsize=9)
    ax2.set_xlim(-1.2, 1.35)
    ax2.set_ylim(-1.35, 1.2)
    ax2.set_aspect("equal")
    ax2.set_axis_off()
    ax2.set_title("球坐标：洋葱模型")

    save(fig, "fig04-grid.svg")


def sparsity_figure():
    stdout, _ = run([exe("sparsity_map")])
    blocks = []
    rows = []
    in_rows = False
    n = 0
    nnz = 0
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("n,"):
            n = int(line.split(",")[1])
        elif line.startswith("nnz,"):
            nnz = int(line.split(",")[1])
        elif line.startswith("block,"):
            _, name, lo, hi = line.split(",")
            blocks.append((name, int(lo), int(hi)))
        elif line == "row,col":
            in_rows = True
        elif in_rows and line:
            row, col = line.split(",")
            rows.append((int(col), int(row)))

    cols = np.array([c for c, _ in rows])
    rws = np.array([r for _, r in rows])
    labels = {
        "cs_neg": "负极固相",
        "cs_pos": "正极固相",
        "ce": "电解液 $c_e$",
        "phis_neg": "负极 $\\phi_s$",
        "phis_pos": "正极 $\\phi_s$",
        "phie": "电解液 $\\phi_e$",
        "cs_surf": "表面浓度",
    }

    fig, ax = plt.subplots(figsize=(6.4, 5.8))
    for name, lo, hi in blocks:
        ax.axvline(lo, color="#c8ccd2", lw=0.6, zorder=0)
        ax.axhline(lo, color="#c8ccd2", lw=0.6, zorder=0)

    ax.scatter(cols, rws, s=2.2, color=BLUE, marker="s", alpha=0.85)

    legend = "\n".join(
        f"{lo}–{hi - 1}  {labels[name]}" for name, lo, hi in blocks
    )
    ax.text(
        0.03,
        0.03,
        legend,
        transform=ax.transAxes,
        fontsize=7.8,
        va="bottom",
        linespacing=1.45,
        bbox={"boxstyle": "round,pad=0.35", "facecolor": "white",
              "edgecolor": "#c8ccd2", "alpha": 0.92},
    )

    ax.set_xlim(-2, n + 2)
    ax.set_ylim(n + 2, -2)
    ax.set_aspect("equal")
    ax.set_xlabel("列（未知量）")
    ax.set_ylabel("行（方程）")
    ax.set_title(f"DFN Jacobian 结构：$n={n}$，非零元 {nnz}（{nnz / n:.1f}$n$）",
                 fontsize=10.5)
    ax.tick_params(labelsize=8)
    save(fig, "fig04-sparsity.svg")


def newton_figure():
    f = lambda x: x * x - 2.0
    df = lambda x: 2.0 * x
    xs = [2.5, 1.65, 1.424]
    x = np.linspace(1.05, 2.6, 300)

    fig, ax = plt.subplots(figsize=(6.2, 4.0))
    ax.plot(x, f(x), color=BLUE, lw=2, label="$f(x)=x^2-2$")
    ax.axhline(0.0, color=GRAY, lw=0.8)
    for k, xk in enumerate(xs):
        tangent = f(xk) + df(xk) * (x - xk)
        ax.plot(x, tangent, color=ORANGE, lw=1.0, alpha=0.65)
        ax.plot([xk], [0.0], "o", color=ORANGE, ms=5)
        ax.plot([xk], [f(xk)], "o", color=INK, ms=4)
        ax.vlines(xk, 0.0, f(xk), color=GRAY, ls=":", lw=0.8)
        ax.text(xk, -0.42, f"$x_{k}$", ha="center", fontsize=9.5)
    ax.plot([np.sqrt(2)], [0], "*", color=GREEN, ms=12)
    ax.annotate("根 $\\sqrt{2}$", xy=(np.sqrt(2), 0.0), xytext=(1.48, 0.72),
                color=GREEN, fontsize=9.5,
                arrowprops=dict(arrowstyle="-", color=GREEN, lw=0.7))
    ax.set_xlim(1.05, 2.6)
    ax.set_ylim(-1.2, 6.2)
    ax.set_xlabel("$x$")
    ax.set_ylabel("$f(x)$")
    ax.grid(alpha=0.3)
    ax.legend(frameon=False, loc="upper left")
    save(fig, "fig06-newton.svg")


def ida_figure():
    fig, ax = plt.subplots(figsize=(9.0, 4.3))
    box(ax, 0.0, 1.45, 2.5, 1.5,
        "DaeModel\n\nresidual(t, y, y')\njacobian(t, y, y', cj)\nvariable_types()",
        fc=PALE_BLUE, fs=9.5)
    box(ax, 3.25, 1.45, 2.5, 1.5,
        "IDA\n\nBDF 步长与阶数\n局部误差控制\n事件检测", fc=FILL, fs=10)
    box(ax, 6.5, 1.65, 2.2, 1.1,
        "KLU\n\n符号分析（一次）\n数值分解（每步）", fc=PALE_ORANGE, fs=9.5)
    box(ax, 3.25, 0.15, 2.5, 0.85,
        "IDACalcIC\n初值一致性", fc="#e7f0e8", fs=9.5)
    box(ax, 3.25, 3.25, 2.5, 0.75,
        "事件函数 g(t, y, y')", fc="#f6f2df", fs=9.5)

    arrow(ax, (2.5, 2.2), (3.25, 2.2))
    ax.text(2.88, 2.3, "$F$, $J$", ha="center", fontsize=9.5)
    arrow(ax, (5.75, 2.2), (6.5, 2.2))
    ax.text(6.12, 2.3, "$J$, $c_j$", ha="center", fontsize=9.5)
    arrow(ax, (6.5, 1.85), (5.75, 1.85), color=GRAY, ls="--")
    ax.text(6.12, 1.7, "$\\Delta y$", ha="center", fontsize=9, color=GRAY)
    arrow(ax, (5.75, 1.85), (5.75, 1.0), color=GRAY, ls="--")
    arrow(ax, (4.5, 1.0), (4.5, 1.45), color=GRAY, ls="--")
    arrow(ax, (4.5, 3.25), (4.5, 2.95), color=GRAY, ls="--")
    ax.text(4.62, 3.05, "到达截止电压", fontsize=9, color=GRAY)

    arrow(ax, (3.25, 1.6), (2.5, 1.6), color=GREEN)
    ax.text(2.87, 1.66, "下一时间步", ha="center", fontsize=9, color=GREEN)

    ax.set_xlim(-0.2, 8.9)
    ax.set_ylim(0.0, 4.2)
    ax.set_axis_off()
    save(fig, "fig08-ida.svg")


def protocol_figure():
    t = np.linspace(0.0, 1.0, 400)
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 5.0))

    panels = [
        ("恒流 CC", np.ones_like(t), 3.95 - 0.45 * t**1.6),
        ("恒压 CV", np.exp(-4.5 * t), 3.50 * np.ones_like(t)),
        (
            "恒流恒压 CC-CV",
            np.where(t < 0.55, 1.0, np.exp(-9.0 * (t - 0.55))),
            np.where(t < 0.55, 3.95 - 0.45 * (t / 0.55) ** 1.5, 3.5),
        ),
        (
            "脉冲",
            ((t % 0.2) < 0.1).astype(float),
            3.95 - 0.22 * t - 0.06 * ((t % 0.2) < 0.1),
        ),
    ]

    for k, (ax, (title, current, voltage)) in enumerate(zip(axes.flat, panels)):
        ax2 = ax.twinx()
        line_i, = ax.plot(t, current, color=BLUE, ls="--", lw=1.4)
        line_v, = ax2.plot(t, voltage, color=ORANGE, lw=1.8)
        if k == 0:
            ax.legend([line_i, line_v], ["电流（左轴）", "电压（右轴）"],
                      fontsize=8, frameon=False, loc="lower left")
        ax.set_title(title, fontsize=10.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.05, 1.15)
        ax2.set_ylim(3.4, 4.05)
        ax.set_yticks([0, 0.5, 1.0])
        ax2.set_yticks([3.5, 3.7, 3.9])
        ax.tick_params(labelsize=8, colors=BLUE)
        ax2.tick_params(labelsize=8, colors=ORANGE)
        ax.set_xlabel("时间（归一化）", fontsize=9)
        ax.set_ylabel("电流", fontsize=9, color=BLUE)
        ax2.set_ylabel("电压 [V]", fontsize=9, color=ORANGE)
        ax.grid(alpha=0.25)

    fig.tight_layout()
    save(fig, "fig09-protocol.svg")


def cccv_figure(workdir):
    path = pathlib.Path(workdir) / "cccv.csv"
    run(
        [
            exe("dfn-sim"),
            "--mode",
            "cccv",
            "--crate",
            "1",
            "--dt",
            "30",
            "--hold",
            "900",
            "--out",
            str(path),
        ]
    )
    with open(path, newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    t = np.array([float(r["t_s"]) for r in rows])
    v = np.array([float(r["voltage_v"]) for r in rows])
    i = np.array([float(r["current_a_m2"]) for r in rows])
    switch = next(float(r["t_s"]) for r in rows if r["phase"] == "cv")

    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    ax2 = ax.twinx()
    ax.plot(t, v, color=BLUE, lw=1.8, label="端电压")
    ax2.plot(t, i, color=ORANGE, lw=1.4, label="电流密度")
    ax.axvline(switch, color=GRAY, ls="--", lw=1.0)
    ax.text(switch - 80, 3.90, "切换到恒压", ha="right", fontsize=9, color=GRAY)
    ax.set_xlabel("时间 [s]")
    ax.set_ylabel("电压 [V]", color=BLUE)
    ax2.set_ylabel("电流密度 [A/m²]", color=ORANGE)
    ax.tick_params(axis="y", colors=BLUE)
    ax2.tick_params(axis="y", colors=ORANGE)
    ax.set_ylim(3.45, 4.02)
    ax2.set_ylim(0, 40)
    ax.grid(alpha=0.25)
    save(fig, "fig09-cccv.svg")


def compare_figure(workdir):
    sys.path.insert(0, str(ROOT / "scripts"))
    import compare_pybamm

    path = pathlib.Path(workdir) / "cc.csv"
    run(
        [
            exe("dfn-sim"),
            "--mode",
            "cc",
            "--crate",
            "1",
            "--dt",
            "60",
            "--out",
            str(path),
        ]
    )
    t_ours, v_ours, current = compare_pybamm.load_curve(str(path))
    current_density = abs(current[0])
    t_ref, v_ref = compare_pybamm.solve_reference(
        0.0165, current_density, t_ours[-1] + 1.0, 801
    )

    t_ours = np.array(t_ours)
    v_ours = np.array(v_ours)
    limit = min(t_ours[-1], t_ref[-1])
    mask = t_ours <= limit
    error = (v_ours[mask] - np.interp(t_ours[mask], t_ref, v_ref)) * 1000.0
    rmse = float(np.sqrt(np.mean(error**2)))

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7.2, 5.2), sharex=True,
        gridspec_kw={"height_ratios": [3, 1.2]},
    )
    ax.plot(t_ours, v_ours, color=BLUE, lw=1.8, label="dfn-sim")
    ax.plot(t_ref, v_ref, color=ORANGE, ls="--", lw=1.6, label="PyBaMM DFN")
    ax.set_ylabel("电压 [V]")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, loc="lower left")
    ax2.plot(t_ours[mask], error, color=GREEN, lw=1.2)
    ax2.axhline(0.0, color=GRAY, lw=0.8)
    ax2.set_ylabel("差值 [mV]")
    ax2.set_xlabel("时间 [s]")
    ax2.grid(alpha=0.25)
    ax2.text(0.02, 0.85, f"RMSE = {rmse:.2f} mV", transform=ax2.transAxes,
             fontsize=9.5, color=GREEN)
    fig.tight_layout()
    save(fig, "fig09-compare.svg")


def perf_figure(workdir):
    specs = [
        (4, 2, 4, 4, 4),
        (8, 4, 8, 6, 6),
        (16, 8, 16, 8, 8),
        (32, 16, 32, 12, 12),
        (64, 32, 64, 16, 16),
    ]
    ns, walls, steps, nnz = [], [], [], []
    for spec in specs:
        path = pathlib.Path(workdir) / f"perf_{spec[0]}.csv"
        _, stderr = run(
            [
                exe("dfn-sim"),
                "--mode",
                "cc",
                "--crate",
                "1",
                "--dt",
                "600",
                "--nx",
                str(spec[0]),
                str(spec[1]),
                str(spec[2]),
                "--nr",
                str(spec[3]),
                str(spec[4]),
                "--out",
                str(path),
            ]
        )
        ns.append(int(re.search(r" n=(\d+)", stderr).group(1)))
        nnz.append(int(re.search(r" nnz=(\d+)", stderr).group(1)))
        walls.append(float(re.search(r"wall=([0-9.]+) ms", stderr).group(1)))
        steps.append(int(re.search(r"steps=(\d+)", stderr).group(1)))

    ns = np.array(ns)
    walls = np.array(walls)
    slope = np.polyfit(np.log(ns), np.log(walls), 1)[0]

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))
    ax.loglog(ns, walls, "o-", color=BLUE, lw=1.6, ms=5)
    ax.loglog(ns, walls[0] * (ns / ns[0]) ** slope, color=GRAY, ls=":",
              label=f"拟合斜率 {slope:.2f}")
    ax.set_xlabel("未知量数 $n$")
    ax.set_ylabel("墙钟时间 [ms]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("1C 放电耗时", fontsize=10.5)

    ax2.plot(ns, steps, "s-", color=ORANGE, lw=1.5, ms=5)
    ax2.set_xscale("log")
    ax2.set_ylim(0, max(steps) * 1.25)
    ax2.set_xlabel("未知量数 $n$")
    ax2.set_ylabel("IDA 步数")
    ax2.grid(alpha=0.25)
    ax2.set_title("步数基本不变（容差驱动）", fontsize=10.5)

    fig.tight_layout()
    save(fig, "fig09-perf.svg")


def spm_figure(workdir):
    stdout, _ = run([exe("spm_compare")])
    currents, rmse, dfn_ms, spm_ms = [], [], [], []
    for line in stdout.splitlines():
        parts = line.split("\t")
        if len(parts) == 10 and parts[0].replace(".", "", 1).isdigit():
            currents.append(float(parts[0]))
            rmse.append(float(parts[4]))
            dfn_ms.append(float(parts[6]))
            spm_ms.append(float(parts[7]))

    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(9.0, 3.8))
    ax.plot(currents, rmse, "o-", color=BLUE, lw=1.6, ms=5)
    ax.set_xlabel("电流密度 [A/m²]")
    ax.set_ylabel("SPM − DFN RMSE [mV]")
    ax.set_ylim(0, 34)
    ax.set_title("电压误差随倍率增长", fontsize=10.5)
    ax.grid(alpha=0.25)

    fine_path = pathlib.Path(workdir) / "fine.csv"
    _, stderr = run(
        [
            exe("dfn-sim"),
            "--mode",
            "cc",
            "--crate",
            "1",
            "--dt",
            "600",
            "--nx",
            "32",
            "16",
            "32",
            "--nr",
            "12",
            "12",
            "--out",
            str(fine_path),
        ]
    )
    fine_ms = float(re.search(r"wall=([0-9.]+) ms", stderr).group(1))

    labels = ["DFN\n8/4/8", "SPM", "DFN\n32/16/32"]
    values = [dfn_ms[3], spm_ms[3], fine_ms]
    colors = [BLUE, ORANGE, BLUE]
    ax2.bar(labels, values, color=colors, alpha=0.85, width=0.55)
    ax2.set_yscale("log")
    ax2.set_ylabel("1C 放电耗时 [ms]")
    ax2.set_title("同样的精度，不同的规模", fontsize=10.5)
    ax2.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    save(fig, "fig10-spm.svg")


def fill_figure():
    n = 5
    before = np.zeros((n, n), dtype=bool)
    before[0, :] = True
    before[:, 0] = True
    np.fill_diagonal(before, True)
    after = before.copy()
    after[1:, 1:] = True
    fill = after & ~before

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 3.5))
    titles = ["消元前：箭头结构", "消去节点 0 后：右下角被填满"]
    for ax, pattern, title, show_fill in zip(
        axes, [before, after], titles, [False, True]
    ):
        for i in range(n):
            for j in range(n):
                if show_fill and fill[i, j]:
                    color = PALE_ORANGE
                elif pattern[i, j]:
                    color = PALE_BLUE
                else:
                    color = "white"
                ax.add_patch(
                    Rectangle((j - 0.5, i - 0.5), 1, 1, facecolor=color,
                              edgecolor="#c9ced6", lw=0.8)
                )
        ax.set_xlim(-0.5, n - 0.5)
        ax.set_ylim(n - 0.5, -0.5)
        ax.set_xticks(range(n))
        ax.set_yticks(range(n))
        ax.set_xlabel("列")
        ax.set_ylabel("行")
        ax.set_title(title, fontsize=10)
        ax.tick_params(length=0)
        for spine in ax.spines.values():
            spine.set_visible(False)

    axes[1].legend(
        handles=[
            Patch(facecolor=PALE_BLUE, edgecolor="#c9ced6", label="原有非零元"),
            Patch(facecolor=PALE_ORANGE, edgecolor="#c9ced6", label="消元产生的填充"),
        ],
        loc="upper center",
        bbox_to_anchor=(0.5, -0.22),
        ncol=2,
        frameon=False,
        fontsize=9,
    )
    fig.tight_layout()
    save(fig, "fig05-fill.svg")


def stability_figure():
    x = np.linspace(-8.0, 0.9, 500)
    explicit = 1.0 + x
    implicit = 1.0 / (1.0 - x)

    fig, ax = plt.subplots(figsize=(6.4, 3.9))
    ax.axhspan(1.0, 8.0, color=PALE_ORANGE, alpha=0.55, lw=0)
    ax.axhspan(-8.0, -1.0, color=PALE_ORANGE, alpha=0.55, lw=0)
    ax.axhline(1.0, color=GRAY, ls=":", lw=1.0)
    ax.axhline(-1.0, color=GRAY, ls=":", lw=1.0)
    ax.axvline(-2.0, color=GRAY, ls="--", lw=1.0)
    ax.plot(x, explicit, color=BLUE, lw=2, label="显式欧拉 $g=1+h\\lambda$")
    ax.plot(x, implicit, color=ORANGE, lw=2, label="隐式欧拉 $g=1/(1-h\\lambda)$")
    ax.set_xlim(-8.0, 0.9)
    ax.set_ylim(-8.0, 7.5)
    ax.set_xlabel("$h\\lambda$")
    ax.set_ylabel("放大因子 $g$")
    ax.grid(alpha=0.25)
    ax.annotate(
        "显式：$h\\lambda<-2$ 时 $|g|>1$",
        xy=(-4.0, -3.0),
        xytext=(-7.6, 2.4),
        color=BLUE,
        fontsize=9,
        arrowprops={"arrowstyle": "->", "color": BLUE, "lw": 0.8},
    )
    ax.annotate(
        "$h\\lambda=-2$",
        xy=(-2.0, -1.0),
        xytext=(-2.35, -4.6),
        color=GRAY,
        fontsize=9,
    )
    ax.legend(loc="upper left", frameon=False, fontsize=9)
    fig.tight_layout()
    save(fig, "fig07-stability.svg")


def fit_figures(workdir):
    wd = pathlib.Path(workdir)
    kappas = {}
    for p in ("a", "b", "c"):
        data = wd / f"fit_{p}.csv"
        sensitivity = wd / f"sens_{p}.csv"
        run([exe("param_fit"), "--protocol", p, "--generate", "--out",
             str(data), "--noise", "0.001", "--seed", "7"])
        stdout, _ = run([exe("param_fit"), "--protocol", p, "--sensitivity",
                         "--data", str(data), "--sensitivity-csv",
                         str(sensitivity)])
        kappas[p] = float(re.search(r"kappa ([0-9.eE+-]+)", stdout).group(1))

    curve_path = wd / "fit_c_curve.csv"
    run([exe("param_fit"), "--protocol", "c", "--fit", "--data",
         str(wd / "fit_c.csv"), "--starts", "8", "--curve", str(curve_path)])

    rows = read_csv((wd / "fit_c.csv").read_text())
    t = np.array([float(r["t_s"]) for r in rows]) / 3600.0
    current = np.array([float(r["current_a_m2"]) for r in rows])
    voltage = np.array([float(r["voltage_v"]) for r in rows])

    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7.6, 4.2), sharex=True,
        gridspec_kw={"height_ratios": [1.0, 1.4]}
    )
    ax.step(t, current, where="post", color=BLUE, lw=1.6)
    ax.axhline(0.0, color=GRAY, lw=0.8, ls=":")
    ax.set_ylabel("电流密度 [A/m²]", color=BLUE)
    ax.tick_params(axis="y", colors=BLUE)
    ax.grid(alpha=0.25)
    ax.set_title("协议 C：静置 / 0.5C 放电 / 静置 / 0.5C 充电 / 1C / 2C",
                 fontsize=10)
    ax2.plot(t, voltage, color=ORANGE, lw=1.0)
    ax2.set_xlabel("时间 [h]")
    ax2.set_ylabel("端电压 [V]", color=ORANGE)
    ax2.tick_params(axis="y", colors=ORANGE)
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "fig11-protocol.svg")

    fig, (ax, ax2) = plt.subplots(
        1, 2, figsize=(9.0, 3.7), gridspec_kw={"width_ratios": [1.7, 1]}
    )
    sens = read_csv((wd / "sens_c.csv").read_text())
    ts = np.array([float(r["t_s"]) for r in sens]) / 3600.0
    ax.plot(ts, 1000.0 * np.array([float(r["s_k"]) for r in sens]),
            color=BLUE, lw=1.3, label="$\\partial V/\\partial\\ln k$")
    ax.plot(ts, 1000.0 * np.array([float(r["s_dneg"]) for r in sens]),
            color=ORANGE, lw=1.3, label="$\\partial V/\\partial\\ln D_{s,-}$")
    ax.plot(ts, 1000.0 * np.array([float(r["s_dpos"]) for r in sens]),
            color=GREEN, lw=1.3, label="$\\partial V/\\partial\\ln D_{s,+}$")
    ax.set_xlabel("时间 [h]")
    ax.set_ylabel("敏感度 [mV / 单位对数]")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=8.5)
    ax.set_title("协议 C 的敏感度", fontsize=10)

    labels = ["A\n单次 0.5C 放电", "B\n+ 静置与充电", "C\n+ 1C/2C"]
    values = [kappas[p] for p in "abc"]
    bars = ax2.bar(labels, values, color=[GRAY, ORANGE, BLUE], alpha=0.85,
                   width=0.55)
    ax2.set_yscale("log")
    ax2.set_ylabel("$\\kappa(J^\\top W J)$")
    ax2.set_title("条件数：实验设计决定可辨识性", fontsize=10)
    ax2.grid(axis="y", alpha=0.25)
    for bar, value in zip(bars, values):
        ax2.text(bar.get_x() + bar.get_width() / 2, value * 1.15,
                 f"{value:.0f}", ha="center", fontsize=9)
    fig.tight_layout()
    save(fig, "fig11-identifiability.svg")

    fit_rows = read_csv(curve_path.read_text())
    ft = np.array([float(r["t_s"]) for r in fit_rows]) / 3600.0
    data_v = np.array([float(r["data_v"]) for r in fit_rows])
    fit_v = np.array([float(r["fit_v"]) for r in fit_rows])
    residual = 1000.0 * (fit_v - data_v)
    rmse = np.sqrt(np.mean(residual ** 2))
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(7.6, 4.4), sharex=True,
        gridspec_kw={"height_ratios": [2.4, 1]}
    )
    ax.plot(ft, data_v, ".", ms=2.5, color=GRAY, alpha=0.6, label="合成数据")
    ax.plot(ft, fit_v, color=BLUE, lw=1.4, label="拟合")
    ax.set_ylabel("端电压 [V]")
    ax.grid(alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title(f"协议 C 的拟合：RMSE = {rmse:.2f} mV", fontsize=10)
    ax2.plot(ft, residual, color=ORANGE, lw=0.8)
    ax2.axhline(0.0, color=GRAY, lw=0.8, ls=":")
    ax2.set_xlabel("时间 [h]")
    ax2.set_ylabel("拟合 − 数据 [mV]")
    ax2.grid(alpha=0.25)
    fig.tight_layout()
    save(fig, "fig11-fit.svg")

    benches = []
    for n, threads in ((1000, 1), (4000, 16)):
        stdout, _ = run([cuda_exe("spm_batch"), "--cpu", "--samples", str(n),
                         "--threads", str(threads)])
        per = float(re.search(r"per_sample_ms=([0-9.eE+-]+)", stdout).group(1))
        benches.append((n, threads, 1000.0 * per))
    gpu = []
    for n in (128, 512, 2048, 8192, 65536):
        stdout, _ = run([cuda_exe("spm_batch"), "--samples", str(n),
                         "--repeat", "2"])
        per = float(re.search(r"per_sample_us=([0-9.eE+-]+)", stdout).group(1))
        gpu.append((n, per))

    fig, ax = plt.subplots(figsize=(6.8, 4.0))
    ax.loglog([b[0] for b in benches[:1]], [b[2] for b in benches[:1]], "o",
              color=GRAY, ms=8, label="CPU 串行")
    ax.loglog([b[0] for b in benches[1:]], [b[2] for b in benches[1:]], "s",
              color=ORANGE, ms=8, label="CPU 16 线程")
    ax.axhline(benches[0][2], color=GRAY, ls="--", lw=0.8, alpha=0.7)
    ax.axhline(benches[1][2], color=ORANGE, ls="--", lw=0.8, alpha=0.7)
    ax.loglog([g[0] for g in gpu], [g[1] for g in gpu], "o-", color=BLUE,
              lw=1.6, ms=6, label="GPU（RTX 3060 Ti）")
    ax.text(1.4e2, benches[0][2] * 0.82,
            f"串行 {benches[0][2] / 1000.0:.2f} ms/样本", color=GRAY,
            fontsize=8.5)
    ax.text(1.4e2, benches[1][2] * 1.18,
            f"16 线程 {benches[1][2] / 1000.0:.2f} ms/样本", color=ORANGE,
            fontsize=8.5)
    ax.set_xlabel("参数样本数 $N$")
    ax.set_ylabel("每个样本耗时 [µs]")
    ax.grid(True, which="both", alpha=0.25)
    ax.legend(frameon=False, fontsize=9)
    ax.set_title("GPU 只在批量足够大时占优", fontsize=10)
    fig.tight_layout()
    save(fig, "fig11-speedup.svg")
    print("fit kappa:", kappas, "cpu:", benches, "gpu:", gpu)


def main():
    global PREVIEW
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", help="figure names to regenerate")
    parser.add_argument("--preview", action="store_true",
                        help="also write PNG previews")
    args = parser.parse_args()
    PREVIEW = args.preview

    figures = {
        "ocv": ocv_figure,
        "tangent": tangent_figure,
        "error": error_figure,
        "dfn": dfn_figure,
        "grid": grid_figure,
        "sparsity": sparsity_figure,
        "fill": fill_figure,
        "newton": newton_figure,
        "ida": ida_figure,
        "stability": stability_figure,
        "protocol": protocol_figure,
    }

    with tempfile.TemporaryDirectory() as workdir:
        for name, function in figures.items():
            if args.only and name not in args.only:
                continue
            function()
        if not args.only or "cccv" in args.only:
            cccv_figure(workdir)
        if not args.only or "compare" in args.only:
            compare_figure(workdir)
        if not args.only or "perf" in args.only:
            perf_figure(workdir)
        if not args.only or "spm" in args.only:
            spm_figure(workdir)
        if not args.only or "fit11" in args.only:
            fit_figures(workdir)


if __name__ == "__main__":
    main()
