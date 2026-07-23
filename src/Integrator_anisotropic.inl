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
#include <hydra/host/System.h>
#include <hydra/device/System.h>


// THIS LIB
#include <core/Utils.h>
#include <core/functions/SgAnisotropic.h>
#include <core/functions/SgIsotropic.h>


declarg(X3_t,  double)
declarg(Z_t,   double)
declarg(K_t,   double)
declarg(Phi_t, double)

using namespace hydra::arguments;
namespace libconf = libconfig;


// Selectable space-time correlation model h in Eq. (S_g PBL). "analytic" uses
// the fast 3D closed-form top-hat azimuthal integral (SgAnisotropic); the other
// members keep the azimuthal integral explicit and evaluate a pluggable h in 4D
// (SgAnisotropicH). Running "tophat" (4D numerical) against "analytic" (3D
// closed form) is the correctness cross-check.
enum HModel { H_ANALYTIC, H_GAUSSIAN, H_TOPHAT, H_SWEEPING, H_ELLIPTIC, H_EXPONENTIAL };


// Integrate one frequency point with the 4D general-h PBL integrand. Returns the
// VEGAS estimate and its absolute error.
template<class Corr>
std::pair<double,double> integrate_aniso(SgAnisotropicParams const& par, Corr const& corr,
                       hydra::VegasState<4, hydra::device::sys_t>& state)
{
    auto integrand = SgAnisotropicH<Corr, X3_t, Z_t, K_t, Phi_t>(par, corr);
    auto vegas     = hydra::Vegas<4, hydra::device::sys_t>(state);
    auto r         = vegas.Integrate(integrand);
    return { r.first, r.second };
}


int main(int argc, char** argv)
{

    // command line: path to the configuration file (defaults to the etc/ copy)
    std::string cfg_path = "../etc/configuration_anisotropic.cfg";
    try {
        TCLAP::CmdLine cmd("NewtonianNoiseCA anisotropic integrator", ' ', "1.0");
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
    constexpr double G      = 6.67e-11;
    constexpr double rho0   = 1.225;
    constexpr double L      = 1E4;
    constexpr double T0     = 300.;

    int Ncalls       = (int) cfg_root["Ncalls"];
    int Niterations  = (int) cfg_root["Niterations"];
    double MaxError  = (double) cfg_root["MaxError"];


    // physics parameters
    auto par = SgAnisotropicParams();
    par.omega        = 2*M_PI;
    par.delta        = (double) cfg_root["delta"];
    par.T_star       = (double) cfg_root["T_star"];
    par.psi          = (double) cfg_root["psi"];
    par.r0           = (double) cfg_root["r0"];
    par.u_star       = (double) cfg_root["u_star"];

    double min_x3    = (double) cfg_root["min_x3"];
    double max_x3    = (double) cfg_root["max_x3"];

    double min_z     = (double) cfg_root["min_z"];
    double max_z     = (double) cfg_root["max_z"];

    double min_k     = (double) cfg_root["min_k"];
    double max_k     = (double) cfg_root["max_k"];

    double omega_min = (double) cfg_root["omega_min"];
    double omega_max = (double) cfg_root["omega_max"];
    int    Npoints   = (int)    cfg_root["Npoints"];


    // space-time correlation model (WP-A1) + humidity enhancement (WP-A2).
    // Optional -> old configs default to the analytic top-hat of the 2022 paper.
    std::string h_model = "analytic";
    double sigma_U      = 1.5;   // sweeping model:  sigma_U ~ u_*
    double beta_sw      = 0.0;   // elliptic model:  V_c / U_c
    double cTv2_ratio   = 1.0;   // humidity:  c_Tv^2 / c_T^2
    double min_phi      = 0.0;
    double max_phi      = 2.0*M_PI;
    bool   omega_log    = false; // log-spaced frequency sweep (default: linear)
    cfg_root.lookupValue("h_model",    h_model);
    cfg_root.lookupValue("sigma_U",    sigma_U);
    cfg_root.lookupValue("beta_sw",    beta_sw);
    cfg_root.lookupValue("cTv2_ratio", cTv2_ratio);
    cfg_root.lookupValue("min_phi",    min_phi);
    cfg_root.lookupValue("max_phi",    max_phi);
    cfg_root.lookupValue("omega_log",  omega_log);

    HModel hmod;
    if      (h_model == "analytic")    hmod = H_ANALYTIC;
    else if (h_model == "gaussian")    hmod = H_GAUSSIAN;
    else if (h_model == "tophat")      hmod = H_TOPHAT;
    else if (h_model == "sweeping")    hmod = H_SWEEPING;
    else if (h_model == "elliptic")    hmod = H_ELLIPTIC;
    else if (h_model == "exponential") hmod = H_EXPONENTIAL;
    else {
        std::cerr << "Unknown h_model '" << h_model
                  << "' (expected analytic|gaussian|tophat|sweeping|elliptic|exponential)" << std::endl;
        return 1;
    }
    std::cerr << "# h_model=" << h_model << "  cTv2_ratio=" << cTv2_ratio << std::endl;


    // Vegas integrator states: 3D for the analytic top-hat path, 4D (adds the
    // azimuthal angle phi) for the pluggable-h path.
    double  min3[3]   = { par.delta , min_z , min_k };
    double  max3[3]   = { max_x3,     max_z,  max_k  };
    double  min4[4]   = { par.delta , min_z , min_k , min_phi };
    double  max4[4]   = { max_x3,     max_z,  max_k , max_phi };

    auto integ3 = hydra::VegasState<3, hydra::device::sys_t>(min3, max3);
    auto integ4 = hydra::VegasState<4, hydra::device::sys_t>(min4, max4);

    auto configure = [&](auto& st, int nd){
        st.SetAlpha(1.5);
        st.SetIterations( Niterations );
        st.SetUseRelativeError(1);
        st.SetMaxError( MaxError );
        st.SetCalls( Ncalls );
        st.SetNDimensions(nd);
        st.SetMode(-1);
        st.SetTrainingCalls( 50000/10 );
        st.SetTrainingIterations( Niterations );
    };
    configure(integ3, 3);
    configure(integ4, 4);


    std::cout << "# f[Hz]\tsqrt_Sh[1/sqrtHz]\trel_err\n";
    for(int i=0; i<Npoints; ++i){

        double omega;
        if      (Npoints <= 1) omega = omega_min;
        else if (omega_log)    omega = omega_min * pow(omega_max/omega_min, i/(Npoints-1.0));
        else                   omega = omega_min + i*(omega_max-omega_min)/(Npoints-1.0);

        par.omega = omega;

        double scale  = G * rho0 / (omega * omega * L * T0);
        scale *= scale;
        scale *= cTv2_ratio;                 // virtual-temperature (humidity) enhancement

        std::pair<double,double> res(0.0, 0.0);
        switch(hmod) {
            case H_ANALYTIC: {
                auto integrand = SgAnisotropic<X3_t, Z_t, K_t>(par);
                auto vegas     = hydra::Vegas<3, hydra::device::sys_t>(integ3);
                auto rr        = vegas.Integrate(integrand);
                res            = { rr.first, rr.second };
                break;
            }
            case H_GAUSSIAN:    res = integrate_aniso(par, CorrGaussian{},       integ4); break;
            case H_TOPHAT:      res = integrate_aniso(par, CorrTophat{},         integ4); break;
            case H_SWEEPING:    res = integrate_aniso(par, CorrSweeping{sigma_U}, integ4); break;
            case H_ELLIPTIC:    res = integrate_aniso(par, CorrElliptic{beta_sw}, integ4); break;
            case H_EXPONENTIAL: res = integrate_aniso(par, CorrExponential{},     integ4); break;
        }

        double sqrtSh = sqrt(scale*res.first);
        double relerr = (res.first > 0.0) ? 0.5*res.second/res.first : 0.0;  // sqrt halves rel. error
        std::cout << omega/(2*M_PI) << '\t' << sqrtSh << '\t' << relerr << std::endl;

    }

    return 0;


}


#endif
