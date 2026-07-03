/*
 *  Created on: 25/02/2022
 *      Author: Davide Brundu
 */


#ifndef SGANISOTROPIC_H_
#define SGANISOTROPIC_H_


#include <cmath>
#include <tuple>
#include <limits>
#include <stdexcept>
#include <assert.h>
#include <utility>
#include <ratio>


#include <hydra/detail/Config.h>
#include <hydra/detail/BackendPolicy.h>
#include <hydra/Types.h>
#include <hydra/Function.h>
#include <hydra/Pdf.h>
#include <hydra/Integrator.h>
#include <hydra/Tuple.h>
#include <hydra/functions/Utils.h>
#include <hydra/functions/Math.h>


struct SgAnisotropicParams
{

    double omega;
    double delta;
    double T_star;
    double psi;
    double r0;
    double u_star;

};


/*
 *
 *  @class Sg anisotropic
 *
 */
template< typename ArgTypeX3, typename ArgTypeZ, typename ArgTypeK, typename Signature = double(ArgTypeX3, ArgTypeZ, ArgTypeK)>
class SgAnisotropic: public hydra::BaseFunctor< SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK>, Signature, 6>
{

    using ThisBaseFunctor = hydra::BaseFunctor< SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK>, Signature, 6 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;


public:


    constexpr static int IDX  = 3;

    SgAnisotropic() : ThisBaseFunctor({(double)2*M_PI,
                                    (double)1.0,
                                    (double)1.0,
                                    (double)0.0,
                                    (double)5.0,
                                    (double)1.0})
    { }


    SgAnisotropic( libconfig::Setting const& cfg):
    ThisBaseFunctor({(double) cfg["omega"],
                     (double) cfg["delta"],
                     (double) cfg["T_star"],
                     (double) cfg["psi"],
                     (double) cfg["r0"],
                     (double) cfg["u_star"]})
    {}


    SgAnisotropic( double omega, double delta, double T_star, double psi, double r0, double u_star) :
    ThisBaseFunctor( omega, delta,  T_star,  psi, r0, u_star)
    {}


    SgAnisotropic( SgAnisotropicParams params):
    ThisBaseFunctor({params.omega,
                    params.delta,
                    params.T_star,
                    params.psi,
                    params.r0,
                    params.u_star})
    {}


    __hydra_dual__
    SgAnisotropic(SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK> const& other):
        ThisBaseFunctor(other)
    {}


    __hydra_dual__
    SgAnisotropic& operator=( SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }


    __hydra_dual__
    void SetParameters(SgAnisotropicParams const& params){

        _par[0] = params.omega;
        _par[1] = params.delta;
        _par[2] = params.T_star;
        _par[3] = params.psi;
        _par[4] = params.r0;
        _par[5] = params.u_star;

    }


    __hydra_dual__ inline
    double Evaluate( ArgTypeX3 x3, ArgTypeZ z, ArgTypeK k )  const
    {

        double omega  = _par[0];
        double delta  = _par[1];
        double T_star = _par[2];
        double psi    = _par[3];
        double r0     = _par[4];
        double u_star = _par[5];

        //double c = 200. ; // up to 200
        //double cutoff = sqrt( std::max((double)0.0 , (double) 1. - ( c/(2.*k*r0 )) ) );
        //if(z<cutoff) return 0.0;

        double p           = sqrt(1.0 - z*z);

        //double L_x3        = std::max(x3, delta);

        double U_x3        = u_star/0.4 * (log( M_E * x3 / delta ));

        // (k*x3)^(2/3) and (k*x3)^(13/3) via a single cube root instead of two std::pow:
        //   cbrt(k*x3) = (k*x3)^(1/3)  ->  ^(2/3) = kxc*kxc,  ^(13/3) = (k*x3)^4 * kxc
        double kx          = k * x3;
        double kxc         = cbrt(kx);

        double tau_k       = x3 / (u_star * std::max( (double)1.0 , kxc*kxc ) );

        double alpha_den   = p * tau_k * k * U_x3;
        double alpha_plus  = (omega * tau_k + 1.0) / alpha_den;
        double alpha_minus = (omega * tau_k - 1.0) / alpha_den;

        double phi_plus    = psi + acos( sgn(alpha_minus) * std::min( (double)1.0, std::abs(alpha_minus) ) ) ;
        double phi_minus   = psi + acos( sgn(alpha_plus)  * std::min( (double)1.0, std::abs(alpha_plus) ) ) ;

        double A           = phi_plus + phi_minus + 0.5 * ( sin(2.*phi_plus) - sin(2.*phi_minus) );

        double B           = 1.0 + exp(-2. * k * p * x3) - 2*exp(-1 * k * p * x3) * cos(k * z * x3);

        double r           = x3*x3*x3 * A * B * exp(-2. * k * p * r0) / ( std::max((double)1.0 , kx*kx*kx*kx*kxc ) );

        return T_star*T_star*r / (u_star*M_PI);

    }


};


#endif
