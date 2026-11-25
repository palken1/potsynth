/*
 * potsynth_vsh.h
 *
 * Vector spherical harmonic front-end: vector fields that need not be
 * gradients of a potential. A model holds nexp fields, each
 *
 * V(r,theta,phi) = \sum_{j=0}^{nrad-1} f_j(r) \sum_{l=lmin}^{lmax} \sum_{m=-min(l,mmax)}^{min(l,mmax)}
 *                    [ a^j_lm Y_lm rhat + b^j_lm Psi_lm + c^j_lm Phi_lm ]
 *
 * with the real harmonics Y_lm = P_l^m(cos theta) cos(m phi) for m >= 0 and
 * P_l^|m|(cos theta) sin(|m| phi) for m < 0, and the vector harmonics of
 * Barrera et al. (1985),
 *
 * Psi_lm = dY_lm/dtheta thetahat + (1/sin theta) dY_lm/dphi phihat
 * Phi_lm = rhat x Psi_lm = -(1/sin theta) dY_lm/dphi thetahat + dY_lm/dtheta phihat
 *
 * The radial functions f_j are those of the caller (potsynth_radial in
 * potsynth.h): the coefficients a^j, b^j, c^j may be values at radial nodes
 * with any linear interpolation rule, spline coefficients, or anything else
 * that is linear in them. Without radial functions the fields are functions
 * on the sphere.
 *
 * The model is an ordinary workspace: it is evaluated with potsynth_eval() or
 * potsynth_eval_rows(), and all fields share one grid of the fast method, so
 * that one window sum per point gives every field.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#ifndef INCLUDED_potsynth_vsh_h
#define INCLUDED_potsynth_vsh_h

#include "potsynth.h"

#ifdef __cplusplus
extern "C" {
#endif

/* quantities computed (bit flags); the outputs of each field are ordered
 * V_r, V_theta, V_phi, divergence, curl */
#define POTSYNTH_VSH_VR              0x1 /* radial component V_r */
#define POTSYNTH_VSH_VH              0x2 /* horizontal components (V_theta, V_phi) */
#define POTSYNTH_VSH_V               0x3 /* field (V_r, V_theta, V_phi) */
#define POTSYNTH_VSH_DIV             0x4 /* surface divergence of the tangential part on the unit
                                            sphere, -\sum f_j l(l+1) b^j_lm Y_lm (divide by r) */
#define POTSYNTH_VSH_CURL            0x8 /* radial component of the surface curl of the tangential part
                                            on the unit sphere, -\sum f_j l(l+1) c^j_lm Y_lm (divide by r) */

/* normalization of the associated Legendre functions of the coefficients */
#define POTSYNTH_VSH_NORM_FULL       0   /* 4-pi fully normalized */
#define POTSYNTH_VSH_NORM_SCHMIDT    1   /* Schmidt semi-normalized */

/*
 * Prototypes
 */

POTSYNTH_API size_t potsynth_vsh_nout(const int quantities);
POTSYNTH_API potsynth_workspace * potsynth_vsh_alloc(const size_t lmin, const size_t lmax, const size_t mmax,
                                                     const size_t nexp,
                                                     const double * a, const double * b, const double * c,
                                                     const int norm, const int quantities,
                                                     const double rmin, const double rmax,
                                                     const potsynth_radial * radial,
                                                     const potsynth_parameters * params);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDED_potsynth_vsh_h */
