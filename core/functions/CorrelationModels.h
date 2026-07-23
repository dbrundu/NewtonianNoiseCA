/*
 *  Created on: 2026
 *      Author: Davide Brundu
 *
 *  Space-time temperature-correlation models: the function h in Eq. (int) and
 *  Eq. (S_g PBL) of Brundu et al., Phys. Rev. D 106, 064040 (2022). Each model
 *  is a small policy struct with a __hydra_dual__ operator() returning the
 *  temporal weight h[...] from the local turbulence quantities:
 *
 *      omega  : angular frequency
 *      kdotU  : k.U = k_perp U cos(phi-psi)   (mean-wind Doppler shift)
 *      tau_k  : eddy turnover time at wavenumber k
 *      k      : wavenumber magnitude
 *      U      : local mean-wind speed magnitude (constant U in the HI model,
 *               height-dependent U_{x3} in the PBL model)
 *
 *  The policy is a template parameter of the SgIsotropic / SgAnisotropicH
 *  integrands, so the choice of correlation model is a compile-time plug and a
 *  run-time config switch. This realises work package WP-A1 of the follow-up
 *  program: turning the "weak dependence on the form of h" remark of the 2022
 *  paper into a quantified family.
 *
 *  The frozen limit h -> delta(omega - k.U) is a singular case handled
 *  analytically (hypergeometric result, App. B of the 2022 paper) and is not
 *  part of this pointwise-evaluable family.
 */

#ifndef CORRELATIONMODELS_H_
#define CORRELATIONMODELS_H_

#include <cmath>
#include <hydra/detail/Config.h>


// Gaussian, finite correlation time -- Eq. (hgauss):
//     h = exp[ -(tau_k (omega - k.U))^2 ].
// The model used in the 2022 HI computation; the default.
struct CorrGaussian
{
    __hydra_dual__ inline
    double operator()(double omega, double kdotU, double tau_k, double /*k*/, double /*U*/) const
    {
        double z = tau_k * (omega - kdotU);
        return exp(-z*z);
    }
};


// Top-hat: h = 1 for |tau_k (omega - k.U)| < 1, else 0.
// The model whose analytic azimuthal integral is used by the fast SgAnisotropic
// path; provided here so the numerical PBL integrand can reproduce it as a
// regression check, and so the HI integrand can select it too.
struct CorrTophat
{
    __hydra_dual__ inline
    double operator()(double omega, double kdotU, double tau_k, double /*k*/, double /*U*/) const
    {
        double z = tau_k * (omega - kdotU);
        return (fabs(z) < 1.0) ? 1.0 : 0.0;
    }
};


// Random-sweeping (Kraichnan): the Eulerian decorrelation time at scale k^{-1}
// is the sweeping time 1/(k sigma_U) rather than the eddy-turnover time tau_k.
//     h = exp[ -((omega - k.U) / (k sigma_U))^2 ].
// sigma_U ~ u_* (surface layer). One-line change of the Gaussian width.
struct CorrSweeping
{
    double sigma_U;

    __hydra_dual__ inline
    double operator()(double omega, double kdotU, double /*tau_k*/, double k, double /*U*/) const
    {
        double z = (omega - kdotU) / (k * sigma_U);
        return exp(-z*z);
    }
};


// Elliptic (He-Zhang): a convection velocity U_c and a sweeping velocity
// V_c = beta_sw U_c combine into a single elliptic argument
//     h = exp[ -( (tau_k (omega - k.U_c))^2 + (k V_c tau_k)^2 ) ].
// U_c is the local mean wind (kdotU already carries it; U supplies its
// magnitude for the sweeping term). beta_sw = 0 recovers the Gaussian.
struct CorrElliptic
{
    double beta_sw;   // V_c / U_c

    __hydra_dual__ inline
    double operator()(double omega, double kdotU, double tau_k, double k, double U) const
    {
        double zc = tau_k * (omega - kdotU);
        double zs = k * beta_sw * U * tau_k;
        return exp(-(zc*zc + zs*zs));
    }
};


// Exponential / Lorentzian-spectrum stress test: h = exp(-|tau_k (omega - k.U)|).
// The associated spectrum decays only as omega^{-2}, which drives the small-k
// (IR) divergence discussed in the 2022 paper -- kept to map the boundary of
// validity, not as a physical model.
struct CorrExponential
{
    __hydra_dual__ inline
    double operator()(double omega, double kdotU, double tau_k, double /*k*/, double /*U*/) const
    {
        double z = tau_k * (omega - kdotU);
        return exp(-fabs(z));
    }
};


#endif
