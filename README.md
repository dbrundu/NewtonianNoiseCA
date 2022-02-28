## Newtonian Noise 
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

A C++14 tool to calculate numerical integrals for newtonian noise induced by temperature fluctuations in atmosphere. This application is based on [HYDRA v.3](https://github.com/MultithreadCorner/Hydra).

## Dependencies
The tool depends on [HYDRA >= v.3.2.1](https://github.com/MultithreadCorner/Hydra), [ROOT >= v.6.14](https://github.com/root-project/root), [libconfig >= v1.5](https://hyperrealm.github.io/libconfig/), [TCLAP >= v1.2.1](http://tclap.sourceforge.net/), [GSL >= 2.7](https://www.gnu.org/software/gsl/), [Eigen3](https://gitlab.com/libeigen/eigen) and [FFTW >= 3.3](https://www.fftw.org/). For the best performances at least TBB or OMP backends are needed. Optionally  [CUDA >= 10.0](https://developer.nvidia.com/cuda-toolkit) is needed for nVidia GPUs. [GCC >= v.8](https://gcc.gnu.org/) is needed. 


## Installation, Build and Run
After the installation of all the dependencies (except HYDRA, which is header only) checkout the [HYDRA v.3](https://github.com/MultithreadCorner/Hydra) library and also this application:
```bash
mkdir <MyDevDirectory>
cd <MyDevDirectory>
git clone https://github.com/MultithreadCorner/Hydra.git Hydra
git clone https://github.com/dbrundu/NewtonianNoiseCA.git NewtonianNoiseCA
```

Then you have to setup the proper enveironment variables:
```bash
export CC=/usr/bin/gcc-8
export CXX=/usr/bin/g++-8
export ROOTSYS=<path-to-root-build>
export HYDRA_INCLUDE_DIR=<path-to-hydra>
...
```

Starting from the NewtonianNoiseCA folder, please create a `build` directory for convenience and run the cmake command:
```bash
cd NewtonianNoiseCA
mkdir build
cd build
cmake -DHYDRA_INCLUDE_DIR=$HYDRA_INCLUDE_DIR ../
```

At this point you can run the examples: please read at the examples section for a description on how to setup and run them. If you want to run the actual application please setup the proper configuration using the file `configuration.cgf` inside the `etc/` folder, then build and run the application as:
```bash
make make Integrator_isotropic_tbb
./Integrator_isotropic_tbb 
```


