#!/usr/bin/env python3
"""
plot_spectra.py -- plot sqrt(S_h) vs frequency from a tidy scan CSV.

One curve per scanned value (e.g. per h_model, or per wind speed). Optional
asymptotic cross-checks and an ET-D sensitivity overlay.

Example:
  ./plot_spectra.py results/hi_band_vs_f.csv --out results/hi_band_vs_f.png \
      --title "HI turbulence: correlation-model band" \
      --tail-slope --frozen U=10,r0=5 --etd data/ET-D.txt

--etd FILE expects two whitespace columns: f[Hz]  sqrt_Sh[1/sqrtHz].
"""

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import nnlib
import asymptotics as asy


def group_by_scan(data):
    groups = {}
    for d in data:
        groups.setdefault(d["scan_val"], []).append(d)
    for v in groups.values():
        v.sort(key=lambda d: d["f"])
    return groups


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv")
    ap.add_argument("--out", required=True)
    ap.add_argument("--title", default="")
    ap.add_argument("--errbars", action="store_true", help="draw VEGAS error bars")
    ap.add_argument("--tail-slope", action="store_true",
                    help="overlay the f^-6 eddy-decay tail reference")
    ap.add_argument("--frozen", help="overlay frozen wind envelope, arg 'U=..,r0=..'")
    ap.add_argument("--etd", help="ET-D sensitivity file (f  sqrt_Sh)")
    args = ap.parse_args()

    scan_key, data = nnlib.read_csv(args.csv)
    groups = group_by_scan(data)

    fig, ax = plt.subplots(figsize=(7, 5))
    for i, (val, pts) in enumerate(sorted(groups.items(), key=str)):
        f  = [p["f"] for p in pts]
        sh = [p["sqrt_Sh"] for p in pts]
        label = "{}={}".format(scan_key, val)
        if args.errbars:
            err = [p["sqrt_Sh"] * p["rel_err"] for p in pts]
            ax.errorbar(f, sh, yerr=err, marker="o", ms=3, capsize=2, label=label)
        else:
            ax.plot(f, sh, marker="o", ms=3, label=label)

    # asymptotic cross-checks anchored to the first curve
    any_pts = sorted(groups.items(), key=str)[0][1]
    fa = any_pts[0]["f"]; ya = any_pts[0]["sqrt_Sh"]
    fs = [p["f"] for p in any_pts]
    if args.tail_slope:
        fb = any_pts[-1]["f"]; yb = any_pts[-1]["sqrt_Sh"]
        ax.plot(fs, asy.powerlaw(fs, asy.TAIL_SLOPE, (fb, yb)),
                "k--", lw=1, label=r"$f^{-6}$ tail")
    if args.frozen:
        kv = dict(tok.split("=") for tok in args.frozen.split(","))
        U = float(kv["U"]); r0 = float(kv["r0"])
        ax.plot(fs, asy.frozen_envelope(fs, U, r0, (fa, ya)),
                "k:", lw=1, label=r"frozen $\propto f^{-49/12}e^{-\omega r_0/U}$")
    if args.etd:
        ef, es = [], []
        for line in open(args.etd):
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            c = line.split()
            ef.append(float(c[0])); es.append(float(c[1]))
        ax.plot(ef, es, color="grey", lw=2, label="ET-D")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("frequency  $f$  [Hz]")
    ax.set_ylabel(r"$\sqrt{S_h}$  [$1/\sqrt{\mathrm{Hz}}$]")
    ax.set_title(args.title)
    ax.grid(True, which="both", alpha=0.3)
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(args.out, dpi=130)
    print("wrote", args.out)


if __name__ == "__main__":
    main()
