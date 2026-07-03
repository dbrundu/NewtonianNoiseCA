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

    LDouble omega;
    LDouble delta;
    LDouble T_star;
    LDouble psi;
    LDouble r0;
    LDouble u_star;

};


/*
 *
 *  @class Sg anisotropic
 *
 */
template< typename ArgTypeX3, typename ArgTypeZ, typename ArgTypeK, typename Signature = LDouble(ArgTypeX3, ArgTypeZ, ArgTypeK)>
class SgAnisotropic: public hydra::BaseFunctor< SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK>, Signature, 6>
{

    using ThisBaseFunctor = hydra::BaseFunctor< SgAnisotropic<ArgTypeX3, ArgTypeZ, ArgTypeK>, Signature, 6 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;


public:


    constexpr static int IDX  = 3;

    SgAnisotropic() : ThisBaseFunctor({(LDouble)2*M_PI,
                                    (LDouble)1.0,
                                    (LDouble)1.0,
                                    (LDouble)0.0,
                                    (LDouble)5.0,
                                    (LDouble)1.0})
    { }


    SgAnisotropic( libconfig::Setting const& cfg):
    ThisBaseFunctor({(LDouble) cfg["omega"],
                     (LDouble) cfg["delta"],
                     (LDouble) cfg["T_star"],
                     (LDouble) cfg["psi"],
                     (LDouble) cfg["r0"],
                     (LDouble) cfg["u_star"]})
    {}


    SgAnisotropic( LDouble omega, LDouble delta, LDouble T_star, LDouble psi, LDouble r0, LDouble u_star) :
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
    LDouble Evaluate( ArgTypeX3 x3, ArgTypeZ z, ArgTypeK k )  const
    {

        LDouble omega  = _par[0];
        LDouble delta  = _par[1];
        LDouble T_star = _par[2];
        LDouble psi    = _par[3];
        LDouble r0     = _par[4];
        LDouble u_star = _par[5];

        //LDouble c = 200. ; // up to 200
        //LDouble cutoff = sqrt( std::max((LDouble)0.0 , (LDouble) 1. - ( c/(2.*k*r0 )) ) );
        //if(z<cutoff) return 0.0;

        LDouble p           = sqrt(1.0 - z*z);

        //LDouble L_x3        = std::max(x3, delta);

        LDouble U_x3        = u_star/0.4 * (log( M_E * x3 / delta ));

        // (k*x3)^(2/3) and (k*x3)^(13/3) via a single cube root instead of two std::pow:
        //   cbrt(k*x3) = (k*x3)^(1/3)  ->  ^(2/3) = kxc*kxc,  ^(13/3) = (k*x3)^4 * kxc
        LDouble kx          = k * x3;
        LDouble kxc         = cbrt(kx);

        LDouble tau_k       = x3 / (u_star * std::max( (LDouble)1.0 , kxc*kxc ) );

        LDouble alpha_den   = p * tau_k * k * U_x3;
        LDouble alpha_plus  = (omega * tau_k + 1.0) / alpha_den;
        LDouble alpha_minus = (omega * tau_k - 1.0) / alpha_den;

        LDouble phi_plus    = psi + acos( sgn(alpha_minus) * std::min( (LDouble)1.0, std::abs(alpha_minus) ) ) ;
        LDouble phi_minus   = psi + acos( sgn(alpha_plus)  * std::min( (LDouble)1.0, std::abs(alpha_plus) ) ) ;

        LDouble A           = phi_plus + phi_minus + 0.5 * ( sin(2.*phi_plus) - sin(2.*phi_minus) );

        LDouble B           = 1.0 + exp(-2. * k * p * x3) - 2*exp(-1 * k * p * x3) * cos(k * z * x3);

        LDouble r           = x3*x3*x3 * A * B * exp(-2. * k * p * r0) / ( std::max((LDouble)1.0 , kx*kx*kx*kx*kxc ) );

        return T_star*T_star*r / (u_star*M_PI);

    }


};


#endif
