#!/usr/bin/env bash
###############################################################################
# One-shot build of the multicore (TBB) integrators.
#
# Clones the header-only HYDRA library if HYDRA_INCLUDE_DIR is not already set,
# then configures and builds the _tbb targets. Safe to re-run.
#
#   HYDRA_INCLUDE_DIR   path to a Hydra checkout (optional; auto-cloned if unset)
#   BUILD_JOBS          parallel compile jobs (default 4; Hydra TUs are RAM-heavy)
#
# On an HPC cluster load your toolchain first, e.g.:
#   module load gcc/13 cmake tbb libconfig python
# (needs: C++20 GCC>=10, CMake>=3.16, TBB, libconfig++, TCLAP headers)
###############################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"   # repo root
cd "$HERE"

# --- HYDRA headers ----------------------------------------------------------
if [ -z "${HYDRA_INCLUDE_DIR:-}" ] || [ ! -f "${HYDRA_INCLUDE_DIR}/hydra/Hydra.h" ]; then
    HYDRA_INCLUDE_DIR="$HERE/.deps/Hydra"
    if [ ! -f "$HYDRA_INCLUDE_DIR/hydra/Hydra.h" ]; then
        echo "[build] cloning HYDRA into .deps/Hydra ..."
        mkdir -p "$HERE/.deps"
        git clone --depth 1 https://github.com/MultithreadCorner/Hydra.git "$HYDRA_INCLUDE_DIR"
    fi
fi
export HYDRA_INCLUDE_DIR
echo "[build] HYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR"

# --- configure + build the multicore targets --------------------------------
BUILD_DIR="${NN_BUILD_DIR:-$HERE/build}"
cmake -S "$HERE" -B "$BUILD_DIR" -DHYDRA_INCLUDE_DIR="$HYDRA_INCLUDE_DIR"
cmake --build "$BUILD_DIR" -j "${BUILD_JOBS:-4}" \
    --target Integrator_isotropic_tbb Integrator_anisotropic_tbb

echo "[build] done -> $BUILD_DIR/Integrator_{isotropic,anisotropic}_tbb"
