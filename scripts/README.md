# scripts/ — parameter scans, plotting, and VEGAS tuning

Helpers to run the integrators over parameter scans, tune the Monte-Carlo
statistics, and plot the results with asymptotic cross-checks.

## Requirements

- The integrators built in `../build/` (e.g. `make -C ../build Integrator_isotropic_tbb Integrator_anisotropic_tbb`).
- Python 3. The **runner/tuner** (`nnlib.py`, `run_scan.py`, `tune_vegas.py`) are
  pure standard library. The **plotters** need `matplotlib` and `numpy`
  (`pip install matplotlib numpy`).

## Files

| file | what it does |
|---|---|
| `nnlib.py` | library: writes libconfig files from dicts, runs a binary, parses `f  sqrt_Sh  rel_err` output. Holds the default parameter sets. |
| `run_scan.py` | run a scan over one parameter → tidy CSV (`<key>, f_Hz, sqrt_Sh, rel_err`). |
| `tune_vegas.py` | convergence study at one frequency: sweep `Ncalls`, report/plot the estimate and its `rel_err`. |
| `asymptotics.py` | 2022-paper asymptotic shapes for overlays (`f^-6` tail, frozen wind envelope, `1/r0` depth law). |
| `plot_spectra.py` | `sqrt(S_h)` vs frequency, one curve per scan value, optional asymptotes + ET-D. |
| `plot_depth.py` | `sqrt(S_h)` vs depth `r0`, with the `1/r0` overlay. |

The integrators now print three columns to stdout — `f[Hz]  sqrt_Sh  rel_err`
(the last is VEGAS' own relative-error estimate) — and accept a `--config` path,
so the scripts drive them from anywhere.

## Typical workflow

1. **Tune** the statistics until `rel_err` is acceptably small (say < 1%) at the
   hardest (highest-frequency) point you care about:

   ```bash
   ./tune_vegas.py --model anisotropic --f 10 --set r0=5 u_star=1.5 h_model=gaussian \
       --ncalls 5e5,1e6,2e6,4e6,8e6,1.6e7 --niter 20 \
       --out results/tune_f10.csv --plot results/tune_f10.png
   ```

2. **Scan** with the chosen `Ncalls`. Correlation-model band vs frequency:

   ```bash
   ./run_scan.py --model anisotropic \
       --scan h_model=analytic,gaussian,sweeping,exponential \
       --set r0=5 u_star=1.5 --freq 1,10,40,log --ncalls 8e6 --niter 20 \
       --out results/pbl_band_vs_f.csv
   ```

   Depth scaling at a fixed frequency (single-frequency runs, `Npoints=1`):

   ```bash
   ./run_scan.py --model anisotropic --scan r0=1,2,5,10,20,50,100,200 \
       --set u_star=1.5 h_model=gaussian --freq 2,2,1 --ncalls 8e6 --niter 20 \
       --out results/pbl_vs_r0.csv
   ```

3. **Plot**:

   ```bash
   ./plot_spectra.py results/pbl_band_vs_f.csv --out results/pbl_band_vs_f.png \
       --title "PBL: correlation-model band" --tail-slope --errbars
   ./plot_depth.py results/pbl_vs_r0.csv --out results/pbl_vs_r0.png \
       --title "PBL: depth scaling at f=2 Hz"
   ```

## Notes on statistics

- VEGAS here is deterministic (fixed seed), so repeated runs are identical; the
  reported `rel_err` is the convergence indicator, and the drift of the estimate
  as `Ncalls` grows shows the residual bias. Tune against both.
- The **top-hat** model in the 4D PBL integrator is discontinuous and
  Monte-Carlo–limited at high frequency; use `h_model=analytic` (fast closed
  form) for the top-hat, and the smooth models (`gaussian`, `sweeping`,
  `elliptic`, `exponential`) for the band.
- ET-D overlays require a local sensitivity file (`--etd`), two columns
  `f[Hz]  sqrt_Sh`; not shipped here.
