# Atmospheric Newtonian Noise Calculator

[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A C++20 tool to calculate the numerical integrals for the Newtonian noise induced by
temperature fluctuations in the atmosphere. This application is built on top of
[HYDRA](https://github.com/MultithreadCorner/Hydra), which lets the same code run on
CPU (single-thread, TBB, OpenMP) and NVidia GPU (CUDA) backends.

## Reference

This repository was used to carry out the demanding numerical computations presented in:

> D. Brundu, M. Cadoni, M. Oi, P. Olla, and A. P. Sanna,
> *Atmospheric Newtonian noise modeling for third-generation gravitational wave detectors*,
> Phys. Rev. D **106**, 064040 (2022).
> [doi:10.1103/PhysRevD.106.064040](https://doi.org/10.1103/PhysRevD.106.064040)

## Dependencies

The tool depends on:

- [HYDRA >= v.4.0](https://github.com/MultithreadCorner/Hydra) (header only)
- [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/)
- [TCLAP >= v1.2.1](http://tclap.sourceforge.net/) (header only)
- [GCC >= v.10](https://gcc.gnu.org/) (or another C++20-capable compiler)

For the best performance at least the TBB or OMP backend is recommended (TBB and OpenMP
are auto-detected at configure time). Optionally,
[CUDA >= 10.0](https://developer.nvidia.com/cuda-toolkit) is needed to target NVidia GPUs.

## Installation, Build and Run

After installing all the dependencies (except HYDRA, which is header only), check out the
[HYDRA v.4](https://github.com/MultithreadCorner/Hydra) library alongside this application:

```bash
mkdir <MyDevDirectory>
cd <MyDevDirectory>
git clone https://github.com/MultithreadCorner/Hydra.git Hydra
git clone https://github.com/dbrundu/NewtonianNoiseCA.git NewtonianNoiseCA
```

Then set up the proper environment variables (point `HYDRA_INCLUDE_DIR` at the directory
that contains the `hydra/` headers; set `CC`/`CXX` only if you need a compiler other than
the system default):

```bash
export CC=/usr/bin/gcc
export CXX=/usr/bin/g++
export HYDRA_INCLUDE_DIR=<path-to-hydra>
```

From the `NewtonianNoiseCA` folder, create a `build` directory and run `cmake`:

```bash
cd NewtonianNoiseCA
mkdir build
cd build
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../
```

To run the application, set up the desired parameters in the configuration file inside the
`etc/` folder (e.g. `configuration_isotropic.cfg`), then build and run one of the targets.
Each target is generated once per available backend, with the suffix `_cpp`, `_tbb`, `_omp`
or `_cuda`:

```bash
make Integrator_isotropic_tbb
./Integrator_isotropic_tbb
```

By default each executable reads `../etc/configuration_<model>.cfg` (i.e. run it from
inside `build/`). A different configuration file can be passed explicitly:

```bash
./Integrator_isotropic_tbb --config /path/to/my_config.cfg
```

The computed noise spectrum is printed to stdout, one value per frequency point.

## License

This project is licensed under the GNU General Public License v3.0 — see the [LICENSE](LICENSE) file.
