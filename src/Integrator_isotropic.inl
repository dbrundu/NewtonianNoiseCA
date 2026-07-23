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
#include <core/functions/SgAnisotropic.h>
#include <core/functions/SgIsotropic.h>


declarg(Phi_t,  double)
declarg(P_t,   double)
declarg(K_t,   double)

using namespace hydra::arguments;
namespace libconf = libconfig;


// Selectable space-time correlation model (the function h in Eq. (int)).
enum HModel { H_GAUSSIAN, H_TOPHAT, H_SWEEPING, H_ELLIPTIC, H_EXPONENTIAL };


// Integrate one frequency point for a given correlation policy Corr. Returns the
// VEGAS estimate and its absolute error. Kept as a template so the pluggable h
// is a compile-time plug behind a run-time switch.
template<class Corr>
std::pair<double,double> integrate_isotropic(SgIsotropicParams const& par, Corr const& corr,
                           hydra::VegasState<3, hydra::device::sys_t>& state)
{
    auto integrand = SgIsotropic<Corr, Phi_t, P_t, K_t>(par, corr);
    auto vegas     = hydra::Vegas<3, hydra::device::sys_t>(state);
    auto r         = vegas.Integrate(integrand);
    return { r.first, r.second };
}


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


    // space-time correlation model (WP-A1) and virtual-temperature/humidity
    // enhancement c_Tv^2/c_T^2 (WP-A2). All optional -> old configs still work.
    std::string h_model = "gaussian";
    double sigma_U      = 1.5;   // sweeping model:  sigma_U ~ u_*
    double beta_sw      = 0.0;   // elliptic model:  V_c / U_c  (0 = Gaussian)
    double cTv2_ratio   = 1.0;   // humidity:  c_Tv^2 / c_T^2   (1 = dry air)
    bool   omega_log    = false; // log-spaced frequency sweep (default: linear)
    cfg_root.lookupValue("h_model",    h_model);
    cfg_root.lookupValue("sigma_U",    sigma_U);
    cfg_root.lookupValue("beta_sw",    beta_sw);
    cfg_root.lookupValue("cTv2_ratio", cTv2_ratio);
    cfg_root.lookupValue("omega_log",  omega_log);

    HModel hmod;
    if      (h_model == "gaussian")    hmod = H_GAUSSIAN;
    else if (h_model == "tophat")      hmod = H_TOPHAT;
    else if (h_model == "sweeping")    hmod = H_SWEEPING;
    else if (h_model == "elliptic")    hmod = H_ELLIPTIC;
    else if (h_model == "exponential") hmod = H_EXPONENTIAL;
    else {
        std::cerr << "Unknown h_model '" << h_model
                  << "' (expected gaussian|tophat|sweeping|elliptic|exponential)" << std::endl;
        return 1;
    }
    std::cerr << "# h_model=" << h_model << "  cTv2_ratio=" << cTv2_ratio << std::endl;


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


    std::cout << "# f[Hz]\tsqrt_Sh[1/sqrtHz]\trel_err\n";
    for(int i=0; i<Npoints; ++i)
    {
        double omega;
        if      (Npoints <= 1) omega = omega_min;
        else if (omega_log)    omega = omega_min * pow(omega_max/omega_min, i/(Npoints-1.0));
        else                   omega = omega_min + i*(omega_max-omega_min)/(Npoints-1.0);

        par.omega = omega;

        double scale = sqrt(0.2) * G * rho0;
        scale /= (omega*omega* L * T0);
        scale *= scale;
        scale *= cTv2_ratio;                 // virtual-temperature (humidity) enhancement


        // Integral with the selected space-time correlation model
        std::pair<double,double> res(0.0, 0.0);
        switch(hmod) {
            case H_GAUSSIAN:    res = integrate_isotropic(par, CorrGaussian{},       integrator); break;
            case H_TOPHAT:      res = integrate_isotropic(par, CorrTophat{},         integrator); break;
            case H_SWEEPING:    res = integrate_isotropic(par, CorrSweeping{sigma_U}, integrator); break;
            case H_ELLIPTIC:    res = integrate_isotropic(par, CorrElliptic{beta_sw}, integrator); break;
            case H_EXPONENTIAL: res = integrate_isotropic(par, CorrExponential{},     integrator); break;
        }

        double sqrtSh = sqrt(scale*res.first);
        double relerr = (res.first > 0.0) ? 0.5*res.second/res.first : 0.0;  // sqrt halves rel. error
        std::cout << omega/(2*M_PI) << '\t' << sqrtSh << '\t' << relerr << std::endl;
    }

    return 0;


}


#endif
