#!/usr/bin/env python3
"""
tune_vegas.py -- VEGAS convergence study for a single frequency point.

Runs the same (model, parameters, frequency) at a growing number of VEGAS calls
and records the estimate and its reported relative error. Use it to pick the
Ncalls needed for smooth, non-noisy scan curves before launching a full scan.

Example (PBL, gaussian h, f=2 Hz, r0=5 m):
  ./tune_vegas.py --model anisotropic --f 2 --set r0=5 u_star=1.5 h_model=gaussian \
      --ncalls 2e5,5e5,1e6,2e6,4e6,8e6 --niter 20 \
      --out results/tune_pbl_f2.csv --plot results/tune_pbl_f2.png

VEGAS is deterministic here, so the run-to-run scatter is zero; the reported
``rel_err`` (VEGAS' own error estimate) is the convergence indicator, and the
drift of the estimate as Ncalls grows shows the residual bias.
"""

import argparse
import nnlib


def parse_sets(items):
    out = {}
    for it in items or []:
        k, _, v = it.partition("=")
        try:
            out[k] = float(v)
        except ValueError:
            out[k] = v
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--model", required=True, choices=["isotropic", "anisotropic"])
    ap.add_argument("--f", type=float, required=True, help="frequency [Hz]")
    ap.add_argument("--set", nargs="*", help="fixed overrides KEY=VAL ...")
    ap.add_argument("--ncalls", required=True, help="comma list, e.g. 2e5,1e6,4e6")
    ap.add_argument("--niter", type=int, default=20)
    ap.add_argument("--backend", default="tbb")
    ap.add_argument("--out", help="output CSV path")
    ap.add_argument("--plot", help="output PNG path (needs matplotlib)")
    args = ap.parse_args()

    base = dict(nnlib.DEFAULTS[args.model])
    base.update(parse_sets(args.set))
    nnlib.freq_sweep(base, args.f, args.f, 1)          # single-frequency run
    base["Niterations"] = args.niter

    ncalls_list = [int(float(x)) for x in args.ncalls.split(",")]

    print("model={} f={} Hz  fixed={}".format(args.model, args.f,
          {k: base[k] for k in base if k not in nnlib.DEFAULTS[args.model]
           or base[k] != nnlib.DEFAULTS[args.model][k]}))
    print("{:>12} {:>16} {:>12}".format("Ncalls", "sqrt_Sh", "rel_err"))
    results = []
    for nc in ncalls_list:
        base["Ncalls"] = nc
        (_, sh, err), = nnlib.run(args.model, base, args.backend)
        results.append((nc, sh, err))
        print("{:>12d} {:>16.6e} {:>12.3g}".format(nc, sh, err))

    if args.out:
        import csv, os
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        with open(args.out, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["Ncalls", "sqrt_Sh", "rel_err"])
            w.writerows(results)
        print("wrote", args.out)

    if args.plot:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        nc = [r[0] for r in results]
        sh = [r[1] for r in results]
        err = [r[1] * r[2] for r in results]      # absolute error on sqrt_Sh
        fig, (a1, a2) = plt.subplots(2, 1, figsize=(6, 6), sharex=True)
        a1.errorbar(nc, sh, yerr=err, marker="o", capsize=3)
        a1.set_xscale("log"); a1.set_yscale("log")
        a1.set_ylabel(r"$\sqrt{S_h}$  [$1/\sqrt{\mathrm{Hz}}$]")
        a1.set_title("VEGAS convergence, {}  f={} Hz".format(args.model, args.f))
        a1.grid(True, which="both", alpha=0.3)
        a2.plot(nc, [r[2] for r in results], marker="s", color="C3")
        a2.set_xscale("log"); a2.set_yscale("log")
        a2.axhline(0.01, ls="--", color="grey", label="1% target")
        a2.set_xlabel("Ncalls"); a2.set_ylabel("rel_err"); a2.legend()
        a2.grid(True, which="both", alpha=0.3)
        fig.tight_layout()
        fig.savefig(args.plot, dpi=130)
        print("wrote", args.plot)


if __name__ == "__main__":
    main()
