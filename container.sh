#!/usr/bin/env bash
###############################################################################
# Build-and-run NewtonianNoiseCA inside a container -- one command, no sudo.
#
#     ./container.sh                 # build image (once) + run ./run.sh
#     REBUILD=1 ./container.sh       # force-rebuild the image
#     NCALLS=1e8 ./container.sh      # forward run.sh knobs (see run.sh header)
#     ./container.sh bash            # drop into a shell in the container
#
# Uses rootless podman if available (no sudo, no daemon), else docker. The repo
# is bind-mounted at /work and built inside the container, so results land back
# in ./scripts/results on the host, owned by you.
###############################################################################
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IMAGE="${IMAGE:-newtoniannoise}"

if command -v podman >/dev/null 2>&1; then
    ENGINE=podman
    OWN=(--userns=keep-id)                     # map container root -> host user
elif command -v docker >/dev/null 2>&1; then
    ENGINE=docker
    OWN=(--user "$(id -u):$(id -g)")           # write host files as you, not root
else
    echo "error: neither podman nor docker found. Ask an admin, or build natively (./run.sh)." >&2
    exit 1
fi
echo "[container] engine: $ENGINE   image: $IMAGE"

if [ "${REBUILD:-0}" = "1" ] || ! "$ENGINE" image inspect "$IMAGE" >/dev/null 2>&1; then
    echo "[container] building image (first run downloads the base image)..."
    "$ENGINE" build -t "$IMAGE" "$HERE"
fi

# forward any run.sh knobs that are set in the environment
ENVS=()
for v in NCALLS NITER FMIN FMAX NFREQ R0 USTAR U R0LIST HIMODELS PBLMODELS BACKEND BUILD_JOBS; do
    [ -n "${!v:-}" ] && ENVS+=(-e "$v=${!v}")
done

# default: run the pipeline via bash (robust to the mount's exec bit); allow an
# override, e.g. `./container.sh bash` for an interactive shell
CMD=("$@")
[ ${#CMD[@]} -eq 0 ] && CMD=(bash run.sh)

# label=disable: let the container read the bind mount on SELinux hosts (Fedora)
exec "$ENGINE" run --rm -i --security-opt label=disable "${OWN[@]}" "${ENVS[@]}" \
    -v "$HERE":/work -w /work \
    "$IMAGE" "${CMD[@]}"
