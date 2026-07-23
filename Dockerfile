# Toolchain image for NewtonianNoiseCA: C++20 compiler, TBB, libconfig, TCLAP,
# Python+matplotlib, and the header-only HYDRA library. The source is NOT copied
# in -- it is bind-mounted at run time (see container.sh) and built by run.sh,
# so `-march=native` matches whatever host actually runs the container and code
# stays editable on the host.
#
# Works with rootless podman (no sudo) or docker. Build/run via ./container.sh.
FROM ubuntu:24.04

ENV DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt-get install -y --no-install-recommends \
        g++ cmake make git ca-certificates pkg-config \
        libtbb-dev libconfig++-dev libtclap-dev \
        python3 python3-matplotlib python3-numpy \
    && rm -rf /var/lib/apt/lists/*

# header-only HYDRA, baked into the image (rarely changes)
RUN git clone --depth 1 https://github.com/MultithreadCorner/Hydra.git /opt/Hydra
ENV HYDRA_INCLUDE_DIR=/opt/Hydra

# keep the container build separate from any native build/ on the host, and give
# HOME/matplotlib a writable location when running as an arbitrary --user
ENV NN_BUILD_DIR=/work/build-container \
    HOME=/tmp \
    MPLCONFIGDIR=/tmp/mpl

WORKDIR /work
CMD ["bash", "run.sh"]
