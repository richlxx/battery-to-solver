import argparse
import csv

import matplotlib.pyplot as plt


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot the OCV curve lab output.")
    parser.add_argument("csv_path", help="CSV produced by ocv_curve")
    parser.add_argument("--save", help="save to a PNG instead of showing")
    args = parser.parse_args()

    soc, ocv = [], []
    with open(args.csv_path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            soc.append(float(row["soc"]))
            ocv.append(float(row["ocv_v"]))

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(soc, ocv, lw=2)
    ax.set_xlabel("SOC")
    ax.set_ylabel("OCV [V]")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()

    if args.save:
        fig.savefig(args.save, dpi=150)
    else:
        plt.show()


if __name__ == "__main__":
    main()
