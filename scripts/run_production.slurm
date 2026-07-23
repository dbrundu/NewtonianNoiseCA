#!/usr/bin/env bash
###############################################################################
# SLURM template for a production scan on one multicore node.
#
# Each VEGAS point is TBB-parallel across all allocated cores, so we request a
# whole node and run the scan serially (one point at a time saturates the cores).
# Submit with:  sbatch scripts/run_production.slurm
#
# ADAPT the SBATCH lines, the module loads, and the scan parameters below.
###############################################################################
#SBATCH --job-name=nn-scan
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128        # 128-core node -> TBB uses all 128
#SBATCH --time=24:00:00
#SBATCH --output=nn-scan-%j.log
set -euo pipefail

# --- cluster toolchain (adapt to your module system) ------------------------
# module load gcc/13 cmake tbb libconfig python matplotlib
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-128}"   # only used by the _omp variant

REPO="${SLURM_SUBMIT_DIR:-$(cd "$(dirname "$0")/.." && pwd)}"
cd "$REPO"

# TBB follows the job's CPU affinity mask (SLURM cpuset). If your cluster does
# not bind cores, uncomment to pin explicitly:
#   RUN="taskset -c 0-$((${SLURM_CPUS_PER_TASK:-128}-1))"
RUN=""

# --- build once (safe to re-run) --------------------------------------------
bash scripts/build.sh

# --- optional: right-size the statistics at the hardest point first ---------
# python3 scripts/tune_vegas.py --model anisotropic --f 10 --set r0=5 u_star=1.5 \
#     h_model=gaussian --ncalls 2e6,8e6,3.2e7,1.28e8 --niter 20 \
#     --out scripts/results/tune_f10.csv --plot scripts/results/tune_f10.png

# --- production scan: PBL correlation-model band vs frequency ----------------
$RUN python3 scripts/run_scan.py --model anisotropic \
    --scan h_model=analytic,gaussian,sweeping,exponential \
    --set r0=5 u_star=1.5 --freq 1,10,40,log \
    --ncalls 40000000 --niter 20 \
    --out scripts/results/pbl_band_vs_f.csv

# --- production scan: depth scaling at f=2 Hz -------------------------------
$RUN python3 scripts/run_scan.py --model anisotropic \
    --scan r0=1,2,5,10,20,50,100,200 \
    --set u_star=1.5 h_model=gaussian --freq 2,2,1 \
    --ncalls 40000000 --niter 20 \
    --out scripts/results/pbl_vs_r0.csv

# --- plots (matplotlib Agg backend works headless) --------------------------
python3 scripts/plot_spectra.py scripts/results/pbl_band_vs_f.csv \
    --out scripts/results/pbl_band_vs_f.png \
    --title "PBL: correlation-model band, r0=5 m" --tail-slope
python3 scripts/plot_depth.py scripts/results/pbl_vs_r0.csv \
    --out scripts/results/pbl_vs_r0.png --title "PBL: depth scaling at f=2 Hz"

echo "done."
