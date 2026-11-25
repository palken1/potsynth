/*
 * vsh.c
 *
 * Vector spherical harmonic front-end. Each of the nexp fields is
 *
 * V = \sum_j f_j(r) \sum_lm [ a^j_lm Y_lm rhat + b^j_lm Psi_lm + c^j_lm Phi_lm ]
 *
 * (potsynth_vsh.h). The fields are the outputs of a model with several
 * coefficient sets and the radial functions f_j of the caller
 * (potsynth_alloc_sets()):
 *
 * V_r     = [a] P
 * V_theta = [b] dP/dtheta - [c] (P/sin theta) d/dphi
 * V_phi   = [b] (P/sin theta) d/dphi + [c] dP/dtheta
 * div     = -l(l+1) [b] P
 * curl    = -l(l+1) [c] P
 *
 * where [x] T D stands for \sum_lm T[Pbar_l^m] D[x_lm cos(m phi) + x_l,-m sin(m phi)],
 * with the coefficient set of each radial function. The components of each
 * field for the accuracy are V_r, the horizontal vector (V_theta, V_phi), the
 * divergence and the curl.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "potsynth_internal.h"
#include "potsynth_vsh.h"

/*
potsynth_vsh_nout()
  Number of outputs per point of each field for a set of quantities
*/

size_t
potsynth_vsh_nout(const int quantities)
{
  if (quantities & ~(POTSYNTH_VSH_V | POTSYNTH_VSH_DIV | POTSYNTH_VSH_CURL))
    return 0; /* invalid flags, as potsynth_mag_nout() and potsynth_grav_nout() */

  return ((quantities & POTSYNTH_VSH_VR) ? 1 : 0) + ((quantities & POTSYNTH_VSH_VH) ? 2 : 0) +
         ((quantities & POTSYNTH_VSH_DIV) ? 1 : 0) + ((quantities & POTSYNTH_VSH_CURL) ? 1 : 0);
}

/*
potsynth_vsh_alloc()
  Allocate a model of nexp vector fields

Inputs: lmin       - minimum SH degree (>= 0)
        lmax       - maximum SH degree
        mmax       - maximum SH order (<= lmax)
        nexp       - number of fields
        a          - coefficients of Y rhat: field e and radial function j at
                     a + (e*nrad + j)*nlm, indexed by potsynth_idx(l,m,lmin,mmax),
                     with nlm = potsynth_nlm(lmin,lmax,mmax) and nrad the number of
                     radial functions (1 without); or NULL for zero
        b          - coefficients of Psi, as a; or NULL for zero
        c          - coefficients of Phi, as a; or NULL for zero
        norm       - normalization of the coefficients (POTSYNTH_VSH_NORM_xxx)
        rmin       - minimum radius of evaluation points (fast method)
        rmax       - maximum radius of evaluation points (fast method)
        radial     - radial functions f_j, or NULL: fields on the sphere (rmin and rmax
                     are then not used, and the radii of the points may be NULL)
        quantities - quantities computed (POTSYNTH_VSH_xxx)
        params     - parameters, or NULL for defaults

Return: pointer to workspace, or NULL on error; evaluate it with potsynth_eval()
or potsynth_eval_rows(): quantity q of field e is output e*nout + q, with
nout = potsynth_vsh_nout(quantities)

Notes:
1) The coefficients are copied once, into the model. An output of an array
that is NULL is zero but is still stored in the grid of the fast method: request
only the quantities that are needed (e.g. POTSYNTH_VSH_VH for fields without a
radial part)

2) The accuracy eps is relative to the largest magnitude of each component of
each field: of V_r, of the horizontal vector (V_theta, V_phi), of the divergence
and of the curl
*/

potsynth_workspace *
potsynth_vsh_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const size_t nexp,
                   const double * a, const double * b, const double * c,
                   const int norm, const double rmin, const double rmax,
                   const potsynth_radial * radial, const int quantities,
                   const potsynth_parameters * params)
{
  const size_t nq = potsynth_vsh_nout(quantities);
  const size_t nrad = (radial != NULL) ? radial->n : 1;
  const potsynth_radial sphere = { 1, NULL, NULL };
  const double * src[3] = { a, b, c };
  size_t off[3], nper = 0;
  potsynth_workspace * w;
  potsynth_output * out;
  potsynth_term * term;
  const double ** sets;
  double * sc, * rl;
  size_t nlm, nset, e, s, j, l, o, t;

  if (lmin > lmax || mmax > lmax)
    {
      POTSYNTH_ERROR_NULL ("degrees and orders must satisfy lmin <= lmax, mmax <= lmax", POTSYNTH_EINVAL);
    }
  else if (nexp == 0 || nrad == 0)
    {
      POTSYNTH_ERROR_NULL ("number of fields and of radial functions must be positive", POTSYNTH_EINVAL);
    }
  else if (a == NULL && b == NULL && c == NULL)
    {
      POTSYNTH_ERROR_NULL ("coefficients must be given", POTSYNTH_EINVAL);
    }
  else if (nq == 0 || (quantities & ~(POTSYNTH_VSH_V | POTSYNTH_VSH_DIV | POTSYNTH_VSH_CURL)))
    {
      POTSYNTH_ERROR_NULL ("invalid quantities", POTSYNTH_EINVAL);
    }
  else if (norm != POTSYNTH_VSH_NORM_FULL && norm != POTSYNTH_VSH_NORM_SCHMIDT)
    {
      POTSYNTH_ERROR_NULL ("invalid normalization", POTSYNTH_EINVAL);
    }

  /* coefficient sets of a field: those of the arrays given, each with nrad sets
     (one for each radial function) next to each other */
  for (s = 0; s < 3; ++s)
    {
      off[s] = nper;
      if (src[s] != NULL)
        ++nper;
    }

  nlm = potsynth_nlm(lmin, lmax, mmax);
  nset = nexp * nper * nrad;

  out = malloc(nexp * nq * sizeof(potsynth_output));
  term = malloc(nexp * 7 * sizeof(potsynth_term));
  sets = malloc(nset * sizeof(double *));
  sc = malloc((lmax + 1) * sizeof(double));
  rl = malloc((lmax + 1) * sizeof(double));
  if (out == NULL || term == NULL || sets == NULL || sc == NULL || rl == NULL)
    {
      free(out);
      free(term);
      free(sets);
      free(sc);
      free(rl);
      POTSYNTH_ERROR_NULL ("unable to allocate vsh model", POTSYNTH_ENOMEM);
    }

  /* the normalization is a factor per degree of the terms (Schmidt
     P_l^m = Pbar_l^m / sqrt(2l+1)), so that the coefficients of the caller are
     used as they are; rl: the same with the factor -l(l+1) of divergence and curl */
  for (l = 0; l <= lmax; ++l)
    {
      sc[l] = (norm == POTSYNTH_VSH_NORM_SCHMIDT) ? 1.0 / sqrt(2.0 * (double) l + 1.0) : 1.0;
      rl[l] = -(double) l * (double) (l + 1) * sc[l];
    }

  /* set (e, s, j): field e, array s, radial function j */
  for (e = 0; e < nexp; ++e)
    {
      for (s = 0; s < 3; ++s)
        {
          if (src[s] == NULL)
            continue;

          for (j = 0; j < nrad; ++j)
            sets[(e * nper + off[s]) * nrad + j] = src[s] + (e * nrad + j) * nlm;
        }
    }

  /* outputs of each field, with its components for the accuracy; a term of an
     array that is not given has factor 0 on a set of the field */
  for (e = 0, o = 0, t = 0; e < nexp; ++e)
    {
      const size_t base = e * nper;
      const size_t sa = (base + ((a != NULL) ? off[0] : 0)) * nrad;
      const size_t sb = (base + ((b != NULL) ? off[1] : 0)) * nrad;
      const size_t sc_ = (base + ((c != NULL) ? off[2] : 0)) * nrad;
      const double fa = (a != NULL) ? 1.0 : 0.0;
      const double fb = (b != NULL) ? 1.0 : 0.0;
      const double fc = (c != NULL) ? 1.0 : 0.0;
      const int comp = (int) (4 * e) + 1;

      if (quantities & POTSYNTH_VSH_VR)
        {
          term[t].set = sa; term[t].fac = fa; term[t].rfac = sc; term[t].theta = POTSYNTH_THETA_P; term[t].dphi = 0;
          out[o].rpow = 0; out[o].nterm = 1; out[o].term = &term[t]; out[o].component = comp; out[o].weight = 1.0;
          ++o; ++t;
        }

      if (quantities & POTSYNTH_VSH_VH)
        {
          potsynth_term * tt = &term[t], * tp = &term[t + 2];

          tt[0].set = sb; tt[0].fac = fb; tt[0].rfac = sc; tt[0].theta = POTSYNTH_THETA_DP; tt[0].dphi = 0;
          tt[1].set = sc_; tt[1].fac = -fc; tt[1].rfac = sc; tt[1].theta = POTSYNTH_THETA_PSIN; tt[1].dphi = 1;

          tp[0].set = sb; tp[0].fac = fb; tp[0].rfac = sc; tp[0].theta = POTSYNTH_THETA_PSIN; tp[0].dphi = 1;
          tp[1].set = sc_; tp[1].fac = fc; tp[1].rfac = sc; tp[1].theta = POTSYNTH_THETA_DP; tp[1].dphi = 0;

          out[o].rpow = 0; out[o].nterm = 2; out[o].term = tt; out[o].component = comp + 1; out[o].weight = 1.0; ++o;
          out[o].rpow = 0; out[o].nterm = 2; out[o].term = tp; out[o].component = comp + 1; out[o].weight = 1.0; ++o;
          t += 4;
        }

      if (quantities & POTSYNTH_VSH_DIV)
        {
          term[t].set = sb; term[t].fac = fb; term[t].rfac = rl; term[t].theta = POTSYNTH_THETA_P; term[t].dphi = 0;
          out[o].rpow = 0; out[o].nterm = 1; out[o].term = &term[t]; out[o].component = comp + 2; out[o].weight = 1.0;
          ++o; ++t;
        }

      if (quantities & POTSYNTH_VSH_CURL)
        {
          term[t].set = sc_; term[t].fac = fc; term[t].rfac = rl; term[t].theta = POTSYNTH_THETA_P; term[t].dphi = 0;
          out[o].rpow = 0; out[o].nterm = 1; out[o].term = &term[t]; out[o].component = comp + 3; out[o].weight = 1.0;
          ++o; ++t;
        }
    }

  /* fields on the sphere: one radial function without radial dependence, on a unit shell */
  if (radial != NULL)
    w = potsynth_alloc_setp(lmin, lmax, mmax, nset, sets, 1.0, rmin, rmax, radial, nexp * nq, out, params);
  else
    w = potsynth_alloc_setp(lmin, lmax, mmax, nset, sets, 1.0, 1.0, 1.0, &sphere, nexp * nq, out, params);

  free(out);
  free(term);
  free(sets);
  free(sc);
  free(rl);

  return w; /* an error is reported by potsynth_alloc_setp() */
}
