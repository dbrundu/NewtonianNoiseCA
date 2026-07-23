#!/usr/bin/env bash
###############################################################################
# NewtonianNoiseCA -- one-shot: build everything, run the parameter scans, and
# produce the plots. Clone the repo and just run:
#
#     ./run.sh
#
# It (1) clones the header-only HYDRA library if needed and builds the multicore
# (TBB) integrators, (2) runs three scans -- the HI and PBL correlation-model
# bands vs frequency and the PBL depth scaling at f = 2 Hz -- and (3) plots them.
# Each VEGAS point is TBB-parallel across all cores, so on a 128-core node one
# point uses all 128 cores; the scans run serially by design.
#
# Everything is overridable via environment variables (defaults in brackets):
#
#   NCALLS   [20000000]  VEGAS calls per iteration (raise for smoother curves)
#   NITER    [20]        VEGAS iterations
#   FMIN,FMAX,NFREQ [1,10,40]  frequency sweep (Hz), log-spaced
#   R0       [5]         detector depth for the band scans (m)
#   USTAR    [1.5]       friction velocity for the PBL model (m/s)
#   U        [10]        wind speed for the HI model (m/s)
#   R0LIST   [1,2,5,10,20,50,100,200]  depths for the depth scan (m)
#   HIMODELS  [gaussian,sweeping,exponential]
#   PBLMODELS [analytic,gaussian,sweeping,exponential]
#   OUT      [scripts/results]   output directory
#   HYDRA_INCLUDE_DIR  path to an existing Hydra checkout (else auto-cloned)
#   BUILD_JOBS [4]      parallel compile jobs
#
# On an HPC cluster, load your toolchain first (adapt to your module system):
#   module load gcc/13 cmake tbb libconfig python
###############################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$HERE"

PY="${PYTHON:-python3}"
NCALLS="${NCALLS:-20000000}"
NITER="${NITER:-20}"
FMIN="${FMIN:-1}"; FMAX="${FMAX:-10}"; NFREQ="${NFREQ:-40}"
R0="${R0:-5}"; USTAR="${USTAR:-1.5}"; U="${U:-10}"
R0LIST="${R0LIST:-1,2,5,10,20,50,100,200}"
HIMODELS="${HIMODELS:-gaussian,sweeping,exponential}"
PBLMODELS="${PBLMODELS:-analytic,gaussian,sweeping,exponential}"
BACKEND="${BACKEND:-tbb}"
OUT="${OUT:-$HERE/scripts/results}"

echo "=============================================================="
echo " NewtonianNoiseCA one-shot run"
echo "   Ncalls=$NCALLS Niter=$NITER  f=[$FMIN,$FMAX]x$NFREQ (log)"
echo "   r0=$R0 u_star=$USTAR U=$U  backend=$BACKEND"
echo "   out=$OUT"
echo "=============================================================="

# --- 1. build --------------------------------------------------------------
bash "$HERE/scripts/build.sh"

# --- 2. scans --------------------------------------------------------------
mkdir -p "$OUT"

echo "[run] HI correlation-model band vs frequency ..."
$PY "$HERE/scripts/run_scan.py" --model isotropic --scan "h_model=$HIMODELS" \
    --set "r0=$R0" "u=$U" --freq "$FMIN,$FMAX,$NFREQ,log" \
    --ncalls "$NCALLS" --niter "$NITER" --backend "$BACKEND" \
    --out "$OUT/hi_band_vs_f.csv"

echo "[run] PBL correlation-model band vs frequency ..."
$PY "$HERE/scripts/run_scan.py" --model anisotropic --scan "h_model=$PBLMODELS" \
    --set "r0=$R0" "u_star=$USTAR" --freq "$FMIN,$FMAX,$NFREQ,log" \
    --ncalls "$NCALLS" --niter "$NITER" --backend "$BACKEND" \
    --out "$OUT/pbl_band_vs_f.csv"

echo "[run] PBL depth scaling at f=2 Hz ..."
$PY "$HERE/scripts/run_scan.py" --model anisotropic --scan "r0=$R0LIST" \
    --set "u_star=$USTAR" "h_model=gaussian" --freq "2,2,1" \
    --ncalls "$NCALLS" --niter "$NITER" --backend "$BACKEND" \
    --out "$OUT/pbl_vs_r0.csv"

# --- 3. plots (skipped gracefully if matplotlib is absent) -----------------
if $PY -c "import matplotlib" >/dev/null 2>&1; then
    echo "[run] plotting ..."
    $PY "$HERE/scripts/plot_spectra.py" "$OUT/hi_band_vs_f.csv" \
        --out "$OUT/hi_band_vs_f.png" --title "HI: correlation-model band (r0=$R0 m)" --tail-slope
    $PY "$HERE/scripts/plot_spectra.py" "$OUT/pbl_band_vs_f.csv" \
        --out "$OUT/pbl_band_vs_f.png" --title "PBL: correlation-model band (r0=$R0 m)" --tail-slope
    $PY "$HERE/scripts/plot_depth.py" "$OUT/pbl_vs_r0.csv" \
        --out "$OUT/pbl_vs_r0.png" --title "PBL: depth scaling at f=2 Hz"
else
    echo "[run] matplotlib not found -> skipping plots (CSVs are in $OUT)"
fi

echo "=============================================================="
echo " done. results in: $OUT"
ls -1 "$OUT" | sed 's/^/   /'
echo "=============================================================="
