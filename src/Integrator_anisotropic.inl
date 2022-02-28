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
#include <hydra/Lambda.h>
#include <hydra/host/System.h>
#include <hydra/device/System.h>


// THIS LIB
#include <core/Utils.h>
#include <core/Types.h>
#include <core/functions/SgAnisotropic.h>
#include <core/functions/SgIsotropic.h>




declarg(X3_t,  LDouble)
declarg(Z_t,   LDouble)
declarg(K_t,   LDouble)

using namespace hydra::arguments;
namespace libconf = libconfig;





int main(int argv, char** argc)
{

    // get configuration
    libconf::Config Cfg;
    Cfg.readFile("../etc/configuration_anisotropic.cfg");
    const libconf::Setting& cfg_root  = Cfg.getRoot();
    
    
    // constants 
    constexpr size_t Ndim   = 3;
    constexpr double G      = 6.67e-11;
    constexpr double rho0   = 1.225;
    constexpr double L      = 1E4;
    constexpr double T0     = 300.;
    int Ncalls = (int) cfg_root["Ncalls"];
    
    
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
    
    
    //integration region limits
    double  min[Ndim]   = { par.delta , min_z , min_k };
    double  max[Ndim]   = { max_x3,     max_z,  max_k  };
    
    

    for(int i=0; i<Npoints; ++i){
    
        double omega = omega_min + i * (omega_max-omega_min)/(Npoints-1.0);
        par.omega = omega;
    
        double scale  = G * rho0 / (omega * omega * L * T0); 
        scale *= scale;

        // Integrand
        auto integrand = SgAnisotropic<X3_t, Z_t, K_t>(par);
        

        // Integral
        
        auto integrator = hydra::Plain< Ndim,  hydra::device::sys_t >(min, max, Ncalls);

        auto result     = integrator.Integrate(integrand);

        std::cout << "Result: " << sqrt(scale*result.first) << std::endl;
        
    }

    return 0;


}











#endif 
