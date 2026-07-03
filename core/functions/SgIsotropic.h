/*
 *  Created on: 25/02/2022
 *      Author: Davide Brundu
 */


#ifndef SGISOTROPIC_H_
#define SGISOTROPIC_H_


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


struct SgIsotropicParams
{

    double omega;
    double alpha;
    double psi;
    double r0;
    double u;

};


/*
 *
 *  @class Sg isotropic
 *
 */
template< typename ArgTypePhi, typename ArgTypeP, typename ArgTypeK, typename Signature = double(ArgTypePhi, ArgTypeP, ArgTypeK)>
class SgIsotropic: public hydra::BaseFunctor< SgIsotropic<ArgTypePhi, ArgTypeP, ArgTypeK>, Signature, 5>
{

    using ThisBaseFunctor = hydra::BaseFunctor< SgIsotropic<ArgTypePhi, ArgTypeP, ArgTypeK>, Signature, 5 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;


public:


    constexpr static int IDX  = 2;

    SgIsotropic() : ThisBaseFunctor({(double)2*M_PI,
                                    (double)1.0,
                                    (double)1.0,
                                    (double)0.0,
                                    (double)5.0,
                                    (double)10.0})
    { }


    SgIsotropic( libconfig::Setting const& cfg):
    ThisBaseFunctor({cfg["omega"],
                     cfg["alpha"],
                     cfg["psi"],
                     cfg["r0"],
                     cfg["u"]})
    {}


    SgIsotropic( SgIsotropicParams const& params):
    ThisBaseFunctor({(double)params.omega,
                    (double)params.alpha,
                    (double)params.psi,
                    (double)params.r0,
                    (double)params.u})
    {}


    __hydra_dual__
    SgIsotropic(SgIsotropic<ArgTypePhi, ArgTypeP, ArgTypeK> const& other):
        ThisBaseFunctor(other)
    {}


    __hydra_dual__
    SgIsotropic& operator=( SgIsotropic<ArgTypePhi, ArgTypeP, ArgTypeK> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        return *this;
    }


    __hydra_dual__
    void SetParameters(SgIsotropicParams const& params){

        _par[0] = params.omega;
        _par[1] = params.alpha;
        _par[2] = params.psi;
        _par[3] = params.r0;
        _par[4] = params.u;

    }


    __hydra_dual__ inline
    double Evaluate( ArgTypePhi phi, ArgTypeP p, ArgTypeK k )  const
    {

        double omega  = _par[0];
        double alpha  = _par[1];
        double psi    = _par[2];
        double r0     = _par[3];
        double u      = _par[4];

        double a = omega - k * u * p * cos(phi-psi);

        // k^(-13/3) and k^(-4/3) via a single cube root instead of two std::pow:
        //   cbrt(k) = k^(1/3)  ->  k^(13/3) = k^4 * cbrt(k),  k^(4/3) = k * cbrt(k)
        double kc = cbrt(k);

        double r = p / sqrt(1.0 - p*p);

        r *= cos(phi)*cos(phi);
        r /= k*k*k*k * kc;
        r *= exp(-2.0 * k * p * r0);
        r *= exp( -alpha*alpha * a * a / (k * kc) );

        return r;


    }


};


#endif
