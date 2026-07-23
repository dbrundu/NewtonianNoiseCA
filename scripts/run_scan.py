#!/usr/bin/env python3
"""
run_scan.py -- run an integrator over a scan of one parameter and save a CSV.

Examples
--------
Correlation-model band (HI), sqrt(S_h) vs frequency, log-spaced:
  ./run_scan.py --model isotropic \
      --scan h_model=gaussian,sweeping,elliptic,exponential \
      --set r0=5 u=10 beta_sw=0.3 \
      --freq 0.5,20,40,log --ncalls 4000000 --niter 20 \
      --out results/hi_band_vs_f.csv

Wind-speed family (PBL), sqrt(S_h) vs frequency:
  ./run_scan.py --model anisotropic --scan u_star=0.5,1.0,1.5 \
      --set r0=5 --freq 1,10,30,log --ncalls 8000000 --niter 20 \
      --out results/pbl_vs_f_varU.csv

Depth scan at fixed frequency f=2 Hz (single-frequency runs, Npoints=1):
  ./run_scan.py --model anisotropic --scan r0=1,2,5,10,20,50,100,200 \
      --set u_star=1.5 --freq 2,2,1 --ncalls 8000000 --niter 20 \
      --out results/pbl_vs_r0.csv

The output CSV is tidy: columns  <scan_key>, f_Hz, sqrt_Sh, rel_err.
"""

import argparse
import nnlib


def parse_scan(spec):
    key, _, vals = spec.partition("=")
    values = []
    for tok in vals.split(","):
        tok = tok.strip()
        try:
            values.append(float(tok))
        except ValueError:
            values.append(tok)          # e.g. h_model strings
    return key, values


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
    ap.add_argument("--scan", required=True, help="KEY=v1,v2,...  (one parameter)")
    ap.add_argument("--set", nargs="*", help="fixed overrides KEY=VAL ...")
    ap.add_argument("--freq", help="FMIN,FMAX,N[,log] in Hz")
    ap.add_argument("--ncalls", type=int, help="VEGAS calls per iteration")
    ap.add_argument("--niter", type=int, help="VEGAS iterations")
    ap.add_argument("--backend", default="tbb", choices=["tbb", "cpp", "omp", "cuda"])
    ap.add_argument("--out", required=True, help="output CSV path")
    args = ap.parse_args()

    scan_key, scan_values = parse_scan(args.scan)
    overrides = parse_sets(args.set)

    if args.freq:
        toks = args.freq.split(",")
        fmin, fmax, npts = float(toks[0]), float(toks[1]), int(toks[2])
        log = len(toks) > 3 and toks[3].lower().startswith("log")
        nnlib.freq_sweep(overrides, fmin, fmax, npts, log)
    if args.ncalls:
        overrides["Ncalls"] = args.ncalls
    if args.niter:
        overrides["Niterations"] = args.niter

    print("model={} scan {}={} fixed={}".format(
        args.model, scan_key, scan_values, overrides))
    rows = nnlib.scan(args.model, scan_key, scan_values, overrides, args.backend)
    nnlib.write_csv(args.out, rows, scan_key)
    print("wrote {} ({} rows)".format(args.out, len(rows)))


if __name__ == "__main__":
    main()
