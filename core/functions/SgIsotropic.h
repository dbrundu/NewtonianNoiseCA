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

#include <core/functions/CorrelationModels.h>


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
 *  Integrand of Eq. (int) for homogeneous-isotropic turbulence, integrated over
 *  (phi, p = k_perp/k, k). The space-time temperature correlation h is supplied
 *  as a policy (see CorrelationModels.h); the default reproduces the Gaussian
 *  model of the 2022 paper.
 *
 */
template< typename Corr, typename ArgTypePhi, typename ArgTypeP, typename ArgTypeK,
          typename Signature = double(ArgTypePhi, ArgTypeP, ArgTypeK)>
class SgIsotropic: public hydra::BaseFunctor< SgIsotropic<Corr, ArgTypePhi, ArgTypeP, ArgTypeK>, Signature, 5>
{

    using ThisBaseFunctor = hydra::BaseFunctor< SgIsotropic<Corr, ArgTypePhi, ArgTypeP, ArgTypeK>, Signature, 5 >;
    using ThisBaseFunctor::_par;
    using Param = hydra::Parameter;


public:


    constexpr static int IDX  = 2;

    SgIsotropic() : ThisBaseFunctor({(double)2*M_PI,
                                    (double)1.0,
                                    (double)1.0,
                                    (double)0.0,
                                    (double)5.0}),
                    _corr()
    { }


    SgIsotropic( SgIsotropicParams const& params, Corr const& corr = Corr()):
    ThisBaseFunctor({(double)params.omega,
                    (double)params.alpha,
                    (double)params.psi,
                    (double)params.r0,
                    (double)params.u}),
    _corr(corr)
    {}


    __hydra_dual__
    SgIsotropic(SgIsotropic<Corr, ArgTypePhi, ArgTypeP, ArgTypeK> const& other):
        ThisBaseFunctor(other),
        _corr(other._corr)
    {}


    __hydra_dual__
    SgIsotropic& operator=( SgIsotropic<Corr, ArgTypePhi, ArgTypeP, ArgTypeK> const& other){
        if(this == &other) return *this;
        ThisBaseFunctor::operator=(other);
        _corr = other._corr;
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

        // k^(-13/3) and k^(-4/3) via a single cube root instead of two std::pow:
        //   cbrt(k) = k^(1/3)  ->  k^(13/3) = k^4 * cbrt(k),  k^(2/3) = cbrt(k)^2
        double kc    = cbrt(k);
        double tau_k = alpha / (kc*kc);              // eddy turnover time E^{-1/3} k^{-2/3}
        double kdotU = k * u * p * cos(phi-psi);     // k.U = k_perp U cos(phi-psi)

        double r = p / sqrt(1.0 - p*p);

        r *= cos(phi)*cos(phi);
        r /= k*k*k*k * kc;
        r *= exp(-2.0 * k * p * r0);
        r *= _corr(omega, kdotU, tau_k, k, u);       // pluggable space-time correlation h

        return r;


    }


private:

    Corr _corr;


};


#endif
