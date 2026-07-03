/*
 *  Created on: 25/02/2022
 *      Author: Davide Brundu
 */

#ifndef INTEGRATOR_SRC_INL_
#define INTEGRATOR_SRC_INL_


// STD
#include <assert.h>
#include <time.h>
#include <cmath>
#include <math.h>
#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <chrono>
#include <limits>
#include <algorithm>

// OTHERs
#include <tclap/CmdLine.h>
#include <libconfig.h++>
#include <sys/stat.h>

// HYDRA
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/FunctorArithmetic.h>
#include <hydra/Plain.h>
#include <hydra/VegasState.h>
#include <hydra/Vegas.h>
#include <hydra/Lambda.h>
#include <hydra/Sobol.h>
#include <hydra/host/System.h>
#include <hydra/device/System.h>


// THIS LIB
#include <core/Utils.h>
#include <core/Types.h>
#include <core/functions/SgAnisotropic.h>
#include <core/functions/SgIsotropic.h>


declarg(Phi_t,  LDouble)
declarg(P_t,   LDouble)
declarg(K_t,   LDouble)

using namespace hydra::arguments;
namespace libconf = libconfig;


int main(int argc, char** argv)
{

    std::cout.precision(12);

    // command line: path to the configuration file (defaults to the etc/ copy)
    std::string cfg_path = "../etc/configuration_isotropic.cfg";
    try {
        TCLAP::CmdLine cmd("NewtonianNoiseCA isotropic integrator", ' ', "1.0");
        TCLAP::ValueArg<std::string> cfgArg("c", "config",
            "Path to the libconfig configuration file", false, cfg_path, "path", cmd);
        cmd.parse(argc, argv);
        cfg_path = cfgArg.getValue();
    }
    catch (TCLAP::ArgException& e) {
        std::cerr << "Command-line error: " << e.error() << " for arg " << e.argId() << std::endl;
        return 1;
    }

    // get configuration
    libconf::Config Cfg;
    try {
        Cfg.readFile(cfg_path.c_str());
    }
    catch (const libconf::FileIOException&) {
        std::cerr << "Could not read configuration file: " << cfg_path << std::endl;
        return 1;
    }
    catch (const libconf::ParseException& e) {
        std::cerr << "Parse error in " << e.getFile() << ":" << e.getLine()
                  << " - " << e.getError() << std::endl;
        return 1;
    }
    const libconf::Setting& cfg_root  = Cfg.getRoot();


    // constants
    constexpr size_t Ndim   = 3;
    constexpr double G      = 6.67e-11;
    constexpr double rho0   = 1.225;
    constexpr double L      = 1E4;
    constexpr double T0     = 300.;

    int Ncalls       = (int) cfg_root["Ncalls"];
    int Niterations  = (int) cfg_root["Niterations"];
    double MaxError  = (double) cfg_root["MaxError"];


    // physics parameters
    auto par = SgIsotropicParams();
    par.omega        = 2*M_PI;
    par.alpha        = (double) cfg_root["alpha"];
    par.psi          = (double) cfg_root["psi"];
    par.r0           = (double) cfg_root["r0"];
    par.u            = (double) cfg_root["u"];

    double min_phi    = (double) cfg_root["min_phi"];
    double max_phi    = (double) cfg_root["max_phi"];

    double min_p     = (double) cfg_root["min_p"];
    double max_p     = (double) cfg_root["max_p"];

    double min_k     = (double) cfg_root["min_k"];
    double max_k     = (double) cfg_root["max_k"];

    double omega_min = (double) cfg_root["omega_min"];
    double omega_max = (double) cfg_root["omega_max"];
    int    Npoints   = (int)    cfg_root["Npoints"];


    //integration region limits
    double  min[Ndim]   = { min_phi , min_p , min_k };
    double  max[Ndim]   = { max_phi,  max_p,  max_k };


    // Vegas integrator state
    auto integrator = hydra::VegasState<Ndim,  hydra::device::sys_t>(min,max);
    integrator.SetAlpha(1.5);
    integrator.SetIterations( Niterations );
    integrator.SetUseRelativeError(1);
    integrator.SetMaxError( MaxError );
    integrator.SetCalls( Ncalls );
    integrator.SetNDimensions(Ndim);
    integrator.SetMode(-1);
    integrator.SetTrainingCalls( 50000/10 );
    integrator.SetTrainingIterations( Niterations );


    for(int i=0; i<Npoints; ++i)
    {
        double omega_step = (omega_max-omega_min)/(Npoints-1.0);
        double omega      = omega_min + i*omega_step;

        par.omega = omega;

        double scale = sqrt(0.2) * G * rho0;
        scale /= (omega*omega* L * T0);
        scale *= scale;


        // Integrand
        auto integrand = SgIsotropic<Phi_t, P_t, K_t>(par);


        // Integral
        auto Vegas_d = hydra::Vegas< Ndim,  hydra::device::sys_t >(integrator);

        auto result  = Vegas_d.Integrate(integrand);

        std::cout << sqrt(scale*result.first) << std::endl;
    }

    return 0;


}


#endif
