import argparse
import csv

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot forward/central difference errors.")
    parser.add_argument("csv_path", help="CSV produced by diff_error")
    parser.add_argument("--save", help="save to a PNG instead of showing")
    args = parser.parse_args()

    h, forward, central = [], [], []
    with open(args.csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            h.append(float(row["h"]))
            forward.append(float(row["forward_err"]))
            central.append(float(row["central_err"]))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.loglog(h, forward, "o-", label="forward (order 1)")
    ax.loglog(h, central, "s-", label="central (order 2)")
    ax.invert_xaxis()
    ax.set_xlabel("step size h")
    ax.set_ylabel("absolute error")
    ax.grid(True, which="both", alpha=0.3)
    ax.legend()
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
