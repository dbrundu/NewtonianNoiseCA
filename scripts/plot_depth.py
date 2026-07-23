#!/usr/bin/env python3
"""
plot_depth.py -- plot sqrt(S_h) vs detector depth r0 from a tidy scan CSV.

Expects a scan over r0 at a single frequency (see run_scan.py --scan r0=...).
Overlays the model-independent large-depth asymptote sqrt(S_h) ~ 1/r0.

Example:
  ./plot_depth.py results/pbl_vs_r0.csv --out results/pbl_vs_r0.png \
      --title "PBL: depth scaling at f=2 Hz" --etd data/ET-D.txt --etd-f 2
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nnlib
import asymptotics as asy


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--errbars", action="store_true")
    ap.add_argument("--etd", help="ET-D file (f  sqrt_Sh) to read a horizontal line from")
    ap.add_argument("--etd-f", type=float, help="frequency [Hz] at which to read ET-D")
    args = ap.parse_args()

    scan_key, data = nnlib.read_csv(args.csv)
    if scan_key != "r0":
        print("warning: expected an r0 scan, got scan key '{}'".format(scan_key))

    pts = sorted(data, key=lambda d: float(d["scan_val"]))
    r0 = [float(d["scan_val"]) for d in pts]
    sh = [d["sqrt_Sh"] for d in pts]

    fig, ax = plt.subplots(figsize=(7, 5))
    if args.errbars:
        err = [d["sqrt_Sh"] * d["rel_err"] for d in pts]
        ax.errorbar(r0, sh, yerr=err, marker="o", capsize=3, label="numerical")
    else:
        ax.plot(r0, sh, marker="o", label="numerical")

    # 1/r0 asymptote anchored to the largest-depth point
    ax.plot(r0, asy.inv_r0(r0, (r0[-1], sh[-1])), "k--", lw=1,
            label=r"$\sqrt{S_h}\propto 1/r_0$")

    if args.etd and args.etd_f:
        ef, es = [], []
        for line in open(args.etd):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c = line.split(); ef.append(float(c[0])); es.append(float(c[1]))
        # nearest ET-D value to the requested frequency
        j = min(range(len(ef)), key=lambda k: abs(ef[k] - args.etd_f))
        ax.axhline(es[j], color="grey", lw=2, label="ET-D @ {} Hz".format(args.etd_f))

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"detector depth  $r_0$  [m]")
    ax.set_ylabel(r"$\sqrt{S_h}$  [$1/\sqrt{\mathrm{Hz}}$]")
    ax.set_title(args.title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
