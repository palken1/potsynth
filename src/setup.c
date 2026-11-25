/*
 * setup.c
 *
 * Model set-up: parameters, coefficients, method selection
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "potsynth_internal.h"

static int auto_method(const size_t lmax, const size_t nfunc, const potsynth_parameters * params);

/*
potsynth_default_parameters()
  Return default parameters
*/

potsynth_parameters
potsynth_default_parameters(void)
{
  potsynth_parameters params;

  params.eps = 1.0e-10;
  params.method = POTSYNTH_METHOD_AUTO;
  params.sigma = 0.0;
  params.w = 0;
  params.rank_tol = 0.0;
  params.nthreads = 0;
  params.npts_hint = 0;
  params.theta_min = 0.0;
  params.theta_max = 0.0;

  return params;
}

/*
potsynth_nlm()
  Number of SH coefficients for degrees lmin..lmax and orders up to mmax,
i.e. the length of a coefficient vector indexed by potsynth_idx()

Inputs: lmin - minimum SH degree
        lmax - maximum SH degree
        mmax - maximum SH order
*/

size_t
potsynth_nlm(const size_t lmin, const size_t lmax, const size_t mmax)
{
  size_t l, sum = 0;

  for (l = lmin; l <= lmax; ++l)
    {
      const size_t M = (l < mmax) ? l : mmax;
      sum += 2 * M + 1;
    }

  return sum;
}

/*
potsynth_alloc()
  Allocate a model workspace with user-defined functionals

Inputs: lmin   - minimum SH degree (>= 0)
        lmax   - maximum SH degree (<= POTSYNTH_LMAX)
        mmax   - maximum SH order (<= lmax)
        coef   - 4-pi fully normalized coefficients, indexed by potsynth_idx(l,m,lmin,mmax),
                 length potsynth_nlm(lmin,lmax,mmax); copied
        a      - reference radius (same unit as r)
        rmin   - minimum radius of evaluation points (fast method; ignored by direct method)
        rmax   - maximum radius of evaluation points (fast method; ignored by direct method)
        nfunc  - number of functionals
        func   - functional definitions, length nfunc; rfac arrays are copied
        params - parameters, or NULL for defaults

Return: pointer to workspace, or NULL on error

Notes:
1) Choose [rmin,rmax] as tight as the data allow: a narrower shell needs fewer
radial functions, so evaluation is faster
*/

potsynth_workspace *
potsynth_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef,
               const double a, const double rmin, const double rmax,
               const size_t nfunc, const potsynth_functional * func, const potsynth_parameters * params)
{
  potsynth_workspace * w;
  potsynth_output * out;
  potsynth_term * term;
  size_t c;

  if (coef == NULL || func == NULL || nfunc == 0)
    {
      POTSYNTH_ERROR_NULL ("coef and functionals must be given", POTSYNTH_EINVAL);
    }

  out = malloc(nfunc * sizeof(potsynth_output));
  term = malloc(nfunc * sizeof(potsynth_term));
  if (out == NULL || term == NULL)
    {
      free(out);
      free(term);
      POTSYNTH_ERROR_NULL ("unable to allocate functional terms", POTSYNTH_ENOMEM);
    }

  /* each functional is one term of coefficient set 0 */
  for (c = 0; c < nfunc; ++c)
    {
      term[c].set = 0;
      term[c].fac = 1.0;
      term[c].rfac = func[c].rfac;
      term[c].theta = func[c].theta;
      term[c].dphi = func[c].dphi;

      out[c].rpow = func[c].rpow;
      out[c].nterm = 1;
      out[c].term = &term[c];
      out[c].component = 0;
      out[c].weight = 0.0;
    }

  w = potsynth_alloc_sets(lmin, lmax, mmax, 1, coef, a, rmin, rmax, NULL, nfunc, out, params);

  free(out);
  free(term);

  return w;
}

/*
potsynth_alloc_sets()
  Allocate a model workspace with several coefficient sets: each output is a
sum of terms, and a term applies one theta and phi operator to one set. All
outputs share one grid (fast method)

Inputs: lmin   - minimum SH degree (>= 0)
        lmax   - maximum SH degree (<= POTSYNTH_LMAX)
        mmax   - maximum SH order (<= lmax)
        nset   - number of coefficient sets
        coef   - 4-pi fully normalized coefficients of the sets, set-major: set s at
                 coef + s * potsynth_nlm(lmin,lmax,mmax), indexed by
                 potsynth_idx(l,m,lmin,mmax); copied
        a      - reference radius (same unit as r)
        rmin   - minimum radius of evaluation points (fast method; ignored by direct method)
        rmax   - maximum radius of evaluation points (fast method; ignored by direct method)
        radial - radial functions of the caller in place of the powers of a/r, or NULL:
                 each output is then sum_j f_j(r) x (its terms with the sets set + j), so
                 that the set of a term is the first of radial->n consecutive sets; a
                 and rpow are not used, and [rmin,rmax] must lie in the range of the
                 functions. radial->f = NULL with n = 1: no radial dependence
        nout   - number of outputs
        out    - output definitions, length nout; terms and rfac arrays are copied
        params - parameters, or NULL for defaults

Return: pointer to workspace, or NULL on error

Notes:
1) The terms of an output must form one field component: either all with the
same operators, or each one of (POTSYNTH_THETA_DP, dphi = 0) and
(POTSYNTH_THETA_PSIN, dphi = 1), the two operators of a horizontal vector, or,
for an output with a given component and weight, each one of
(POTSYNTH_THETA_S2A, dphi = 0) and (POTSYNTH_THETA_S2B, dphi = 1), the two
operators of a horizontal tensor

2) potsynth_alloc() is the case of one set and one term with factor 1 per output
*/

potsynth_workspace *
potsynth_alloc_sets(const size_t lmin, const size_t lmax, const size_t mmax,
                    const size_t nset, const double * coef,
                    const double a, const double rmin, const double rmax,
                    const potsynth_radial * radial,
                    const size_t nout, const potsynth_output * out, const potsynth_parameters * params)
{
  if (coef == NULL || nset == 0 || lmin > lmax || mmax > lmax)
    {
      POTSYNTH_ERROR_NULL ("coefficient sets must be given, with lmin <= lmax and mmax <= lmax", POTSYNTH_EINVAL);
    }
  else
    {
      const size_t nlm = potsynth_nlm(lmin, lmax, mmax);
      const double ** sets = malloc(nset * sizeof(double *));
      potsynth_workspace * w;
      size_t s;

      if (sets == NULL)
        {
          POTSYNTH_ERROR_NULL ("unable to allocate coefficient sets", POTSYNTH_ENOMEM);
        }

      for (s = 0; s < nset; ++s)
        sets[s] = coef + s * nlm;

      w = potsynth_alloc_setp(lmin, lmax, mmax, nset, sets, a, rmin, rmax, radial, nout, out, params);
      free(sets);

      return w;
    }
}

/*
potsynth_alloc_setp()
  potsynth_alloc_sets() with the coefficient sets given by pointers, so that a
front-end can use the arrays of its caller without copying them into one block

Inputs: sets - coefficients of each set, sets[s] indexed by
               potsynth_idx(l,m,lmin,mmax); length nset; copied
        (the other arguments as for potsynth_alloc_sets())

Return: pointer to workspace, or NULL on error
*/

potsynth_workspace *
potsynth_alloc_setp(const size_t lmin, const size_t lmax, const size_t mmax,
                    const size_t nset, const double * const * sets,
                    const double a, const double rmin, const double rmax,
                    const potsynth_radial * radial,
                    const size_t nout, const potsynth_output * out, const potsynth_parameters * params)
{
  const potsynth_parameters p = (params != NULL) ? *params : potsynth_default_parameters();
  const size_t L1 = lmax + 1;
  const size_t nlm = L1 * (lmax + 2) / 2;
  const size_t nfunc = nout;
  const size_t nrad = (radial != NULL) ? radial->n : 0;
  potsynth_workspace * w;
  size_t npart = 0;
  size_t c, l, m, s, t, k;
  int status;

  if (lmin > lmax)
    {
      POTSYNTH_ERROR_NULL ("lmin must be <= lmax", POTSYNTH_EINVAL);
    }
  else if (mmax > lmax)
    {
      POTSYNTH_ERROR_NULL ("mmax must be <= lmax", POTSYNTH_EINVAL);
    }
  else if (lmax > POTSYNTH_LMAX)
    {
      POTSYNTH_ERROR_NULL ("lmax exceeds POTSYNTH_LMAX", POTSYNTH_EUNSUP);
    }
  else if (sets == NULL || out == NULL || nout == 0 || nset == 0)
    {
      POTSYNTH_ERROR_NULL ("coefficient sets and outputs must be given", POTSYNTH_EINVAL);
    }
  else if (radial != NULL && (radial->n == 0 || (radial->f == NULL && radial->n != 1)))
    {
      POTSYNTH_ERROR_NULL ("radial functions must be given (or one function without radial dependence)",
                           POTSYNTH_EINVAL);
    }
  else if (!(a > 0.0) || !isfinite(a))
    {
      POTSYNTH_ERROR_NULL ("reference radius must be positive", POTSYNTH_EINVAL);
    }
  else if (!(p.eps > 0.0) || p.sigma < 0.0 || (p.sigma > 0.0 && p.sigma < 1.0) || p.w < 0 ||
           p.rank_tol < 0.0 || p.nthreads < 0)
    {
      /* sigma below 1 would put the spectrum beyond the grid rows */
      POTSYNTH_ERROR_NULL ("invalid parameters", POTSYNTH_EINVAL);
    }
  else if (p.method != POTSYNTH_METHOD_AUTO && p.method != POTSYNTH_METHOD_FAST &&
           p.method != POTSYNTH_METHOD_DIRECT)
    {
      POTSYNTH_ERROR_NULL ("invalid method", POTSYNTH_EINVAL);
    }
  else if (!(p.theta_min >= 0.0) || !(p.theta_max >= 0.0) || p.theta_max > M_PI ||
           (p.theta_max > 0.0 && p.theta_min > p.theta_max) || (p.theta_max == 0.0 && p.theta_min > 0.0))
    {
      POTSYNTH_ERROR_NULL ("colatitude range must satisfy 0 <= theta_min <= theta_max <= pi", POTSYNTH_EINVAL);
    }

  for (c = 0; c < nfunc; ++c)
    {
      if (out[c].nterm == 0 || out[c].term == NULL)
        {
          POTSYNTH_ERROR_NULL ("each output needs at least one term", POTSYNTH_EINVAL);
        }
      else if (radial == NULL && (long) lmin + out[c].rpow < 0)
        {
          POTSYNTH_ERROR_NULL ("lmin + rpow must be >= 0", POTSYNTH_EINVAL);
        }
      else if (out[c].component < 0 || !(out[c].weight >= 0.0) || !isfinite(out[c].weight))
        {
          POTSYNTH_ERROR_NULL ("component and weight of an output must be >= 0", POTSYNTH_EINVAL);
        }

      for (t = 0; t < out[c].nterm; ++t)
        {
          const potsynth_term * tm = &(out[c].term[t]);

          if (tm->theta < POTSYNTH_THETA_P || tm->theta > POTSYNTH_THETA_S2B)
            {
              POTSYNTH_ERROR_NULL ("invalid theta operator", POTSYNTH_EINVAL);
            }
          else if (tm->dphi != 0 && tm->dphi != 1)
            {
              POTSYNTH_ERROR_NULL ("dphi must be 0 or 1", POTSYNTH_EINVAL);
            }
          else if (tm->set >= nset || tm->set + ((nrad > 0) ? nrad : 1) > nset)
            {
              POTSYNTH_ERROR_NULL ("coefficient set of a term out of range", POTSYNTH_EINVAL);
            }
          else if (!isfinite(tm->fac))
            {
              POTSYNTH_ERROR_NULL ("factor of a term must be finite", POTSYNTH_EINVAL);
            }
          else if (t > 0 &&
                   !(tm->theta == out[c].term[0].theta && tm->dphi == out[c].term[0].dphi) &&
                   !(potsynth_op_family(tm->theta, tm->dphi) == 1 &&
                     potsynth_op_family(out[c].term[0].theta, out[c].term[0].dphi) == 1) &&
                   !(out[c].component > 0 && out[c].weight > 0.0 &&
                     potsynth_op_family(tm->theta, tm->dphi) == 2 &&
                     potsynth_op_family(out[c].term[0].theta, out[c].term[0].dphi) == 2))
            {
              POTSYNTH_ERROR_NULL ("terms of an output must form one field component", POTSYNTH_EINVAL);
            }
        }

      npart += out[c].nterm;
    }

  w = calloc(1, sizeof(potsynth_workspace));
  if (w == NULL)
    {
      POTSYNTH_ERROR_NULL ("unable to allocate workspace", POTSYNTH_ENOMEM);
    }

  w->lmin = lmin;
  w->lmax = lmax;
  w->mmax = mmax;
  w->a = a;
  w->rmin = rmin;
  w->rmax = rmax;
  w->nfunc = nfunc;
  w->nout = nfunc;
  w->params = p;
  w->check_err = -1.0;
  w->nthreads = potsynth_max_threads(p.nthreads);
  w->th_lo = p.theta_min;
  w->th_hi = (p.theta_max > 0.0) ? p.theta_max : M_PI;

  if (p.method == POTSYNTH_METHOD_AUTO)
    {
      /* the fast method needs the radial range of the points */
      if (rmin > 0.0 && rmin <= rmax && isfinite(rmax))
        w->method = auto_method(lmax, nfunc, &p);
      else
        w->method = POTSYNTH_METHOD_DIRECT;
    }
  else
    w->method = p.method;

  if (w->method == POTSYNTH_METHOD_FAST && !(rmin > 0.0 && rmin <= rmax && isfinite(rmax)))
    {
      potsynth_free(w);
      POTSYNTH_ERROR_NULL ("fast method needs 0 < rmin <= rmax", POTSYNTH_EINVAL);
    }

  w->nrad = nrad;
  w->rfn = (radial != NULL) ? radial->f : NULL;
  w->rfn_params = (radial != NULL) ? radial->params : NULL;
  w->nset = nset;
  w->nalf = nlm;
  w->npart = npart;
  w->rpow = malloc(nfunc * sizeof(int));
  w->ccomp = malloc(nfunc * sizeof(int));
  w->cwt = malloc(nfunc * sizeof(double));
  w->pbase = malloc((nfunc + 1) * sizeof(size_t));
  w->pset = malloc(npart * sizeof(size_t));
  w->ptheta = malloc(npart * sizeof(int));
  w->pdphi = malloc(npart * sizeof(int));
  w->prfac = potsynth_aligned_calloc(npart * L1 * sizeof(double));
  w->C = potsynth_aligned_calloc(nset * nlm * sizeof(double));
  w->S = potsynth_aligned_calloc(nset * nlm * sizeof(double));
  if (!w->rpow || !w->ccomp || !w->cwt || !w->pbase || !w->pset || !w->ptheta || !w->pdphi || !w->prfac || !w->C || !w->S)
    {
      potsynth_free(w);
      POTSYNTH_ERROR_NULL ("unable to allocate model arrays", POTSYNTH_ENOMEM);
    }

  for (c = 0, k = 0; c < nfunc; ++c)
    {
      w->rpow[c] = (nrad > 0) ? 0 : out[c].rpow;
      w->ccomp[c] = out[c].component;
      w->cwt[c] = out[c].weight;
      w->pbase[c] = k;

      for (t = 0; t < out[c].nterm; ++t, ++k)
        {
          const potsynth_term * tm = &(out[c].term[t]);

          w->pset[k] = tm->set;
          w->ptheta[k] = tm->theta;
          w->pdphi[k] = tm->dphi;

          /* the constant factor is applied only when it is not 1, so that the per-degree
             factors of models of potsynth_alloc() are those of version 1bb53c5 bit for bit
             (the radial weights and hence the fast results may differ in the last bits) */
          for (l = 0; l <= lmax; ++l)
            {
              const double rf = (tm->rfac != NULL) ? tm->rfac[l] : 1.0;
              w->prfac[k * L1 + l] = (tm->fac == 1.0) ? rf : tm->fac * rf;
            }
        }
    }
  w->pbase[nfunc] = k;

  /* field components; the outputs of a component given by the caller must have
     the same radial power (one radial basis bounds the error of the component) */
  {
    size_t * grp = malloc(nfunc * sizeof(size_t));
    size_t c2;
    int bad = 0;

    if (grp == NULL)
      {
        potsynth_free(w);
        POTSYNTH_ERROR_NULL ("unable to allocate model arrays", POTSYNTH_ENOMEM);
      }

    w->ncomp = potsynth_components(w, grp);

    for (c = 0; c < nfunc && !bad; ++c)
      {
        for (c2 = 0; c2 < c && !bad; ++c2)
          bad = (grp[c2] == grp[c] && w->rpow[c2] != w->rpow[c]);
      }

    free(grp);

    if (bad)
      {
        potsynth_free(w);
        POTSYNTH_ERROR_NULL ("outputs of a component must have the same radial power", POTSYNTH_EINVAL);
      }
  }

  /* public layout (degree-major, potsynth_idx) -> internal M-major */
  {
    for (s = 0; s < nset; ++s)
      {
        const double * cs = sets[s];
        double * C = w->C + s * nlm;
        double * S = w->S + s * nlm;

        for (l = lmin; l <= lmax; ++l)
          {
            const size_t M = (l < mmax) ? l : mmax;

            for (m = 0; m <= M; ++m)
              {
                const size_t idx = potsynth_alf_idx(l, m, lmax);
                const long ll = (long) l, lm = (long) m, llmin = (long) lmin, lmmax = (long) mmax;

                C[idx] = cs[potsynth_idx(ll, lm, llmin, lmmax)];
                if (m > 0)
                  S[idx] = cs[potsynth_idx(ll, -lm, llmin, lmmax)];
              }
          }
      }
  }

  w->alf_workspace_p = potsynth_alf_alloc(lmax);
  if (w->alf_workspace_p == NULL)
    {
      potsynth_free(w);
      return NULL; /* reported by potsynth_alf_alloc() */
    }

  if (w->method == POTSYNTH_METHOD_FAST)
    {
      status = potsynth_fast_init(w);

      /* check the chosen window against the direct sum at points where its error is
       * likely to be largest (potsynth_fast_check()); widen it while the error
       * exceeds eps (automatic width only) */
      if (status == POTSYNTH_SUCCESS && w->params.w <= 0)
        {
          int iter;
          for (iter = 0; ; ++iter)
            {
              const double ratio = potsynth_fast_check(w);
              if (ratio < 0.0)
                {
                  status = (int) -ratio;
                  potsynth_error("the window check after set-up failed (out of memory, or the radial "
                                 "functions of the caller failed)", __FILE__, __LINE__, status);
                  break;
                }
              w->check_err = ratio;
              if (ratio <= w->params.eps || w->w >= POTSYNTH_WMAX || iter == 4)
                break;
              potsynth_aligned_free(w->grid);
              w->grid = NULL;
              w->params.w = w->w + 1;
              status = potsynth_fast_init(w);
              w->params.w = 0;
              if (status)
                break;
            }
        }
    }
  else
    status = potsynth_direct_init(w);

  if (status)
    {
      potsynth_free(w);
      return NULL; /* reported above */
    }

  /* the fast method evaluates from the grid: its copy of the coefficients, which
     can be large for a model with many sets, is not needed any more */
  if (w->method == POTSYNTH_METHOD_FAST)
    {
      potsynth_aligned_free(w->C);
      potsynth_aligned_free(w->S);
      w->C = NULL;
      w->S = NULL;
    }

  return w;
}

/*
potsynth_term_factors()
  Factors per degree of one term of a part for order m, used by both methods to
build their weight rows: f[l - m] = beta[l] rfac[l] x (factor of the operator
term), for l from lstart to lmax, and 0 below lstart

Inputs: w    - workspace
        p    - part
        sub  - term of the theta operator of the part
        m    - order
        beta - factors per degree beta[l] (radial expansion, fast method), or NULL
        f    - (output) factors, f[l - m] for l = m..lmax

Return: lstart, the first degree with a nonzero factor: max(lmin, m), or lmax + 1
for P/sin and S2B at m = 0, where they are not defined
*/

size_t
potsynth_term_factors(const potsynth_workspace * w, const size_t p, const size_t sub,
                      const size_t m, const double * beta, double * f)
{
  const size_t L = w->lmax;
  const int op = w->ptheta[p];
  const double * rfac = w->prfac + p * (L + 1);
  size_t lstart = (w->lmin > m) ? w->lmin : m;
  size_t l;

  if ((op == POTSYNTH_THETA_PSIN || op == POTSYNTH_THETA_S2B) && m == 0)
    lstart = L + 1;

  for (l = m; l <= L; ++l)
    {
      if (l < lstart)
        f[l - m] = 0.0;
      else if (beta != NULL)
        f[l - m] = beta[l] * rfac[l] * potsynth_op_factor(op, sub, l, m);
      else
        f[l - m] = rfac[l] * potsynth_op_factor(op, sub, l, m);
    }

  return lstart;
}

/*
potsynth_radial_eval()
  Radial functions of the caller at a radius (models with w->nrad > 0)

Inputs: w     - workspace
        r     - radius
        first - (output) first function that is nonzero at r
        n     - (output) number of functions that are nonzero at r
        f     - (output) their values, f[0..n-1]; length nrad

Return: success, or POTSYNTH_EDOM if r is outside of the range of the functions
*/

int
potsynth_radial_eval(const potsynth_workspace * w, const double r, size_t * first, size_t * n, double * f)
{
  *first = 0;
  *n = 1;

  if (w->rfn == NULL)
    {
      f[0] = 1.0; /* no radial dependence */
      return POTSYNTH_SUCCESS;
    }

  if (w->rfn(r, first, n, f, NULL, w->rfn_params) != POTSYNTH_SUCCESS || *n == 0 || *first + *n > w->nrad)
    return POTSYNTH_EDOM;

  return POTSYNTH_SUCCESS;
}

/*
potsynth_free()
  Free model workspace
*/

void
potsynth_free(potsynth_workspace * w)
{
  if (w == NULL)
    return; /* as free(NULL) */

  if (w->rpow)
    free(w->rpow);

  if (w->ccomp)
    free(w->ccomp);

  if (w->cwt)
    free(w->cwt);

  if (w->pbase)
    free(w->pbase);

  if (w->pset)
    free(w->pset);

  if (w->ptheta)
    free(w->ptheta);

  if (w->pdphi)
    free(w->pdphi);

  if (w->prfac)
    potsynth_aligned_free(w->prfac);

  if (w->C)
    potsynth_aligned_free(w->C);

  if (w->S)
    potsynth_aligned_free(w->S);

  if (w->alf_workspace_p)
    potsynth_alf_free(w->alf_workspace_p);

  if (w->grid)
    potsynth_aligned_free(w->grid);

  if (w->radial_workspace_p)
    potsynth_radial_free(w->radial_workspace_p);

  potsynth_direct_free(w);

  if (w->omap)
    free(w->omap);

  free(w);
}

/*
potsynth_get_info()
  Return parameters chosen during set-up

Inputs: info - (output) parameters
        w    - workspace

Notes:
1) The fast method fields are 0 for the direct method
*/

int
potsynth_get_info(potsynth_info * info, const potsynth_workspace * w)
{
  memset(info, 0, sizeof(potsynth_info));

  info->lmin = w->lmin;
  info->lmax = w->lmax;
  info->mmax = w->mmax;
  info->nfunc = w->nfunc;
  info->nout = w->nout;
  info->method = w->method;

  if (w->method == POTSYNTH_METHOD_FAST)
    {
      info->w = w->w;
      info->sigma = w->sigma;
      info->rho = w->radial_workspace_p->rho;
      info->ncheb = w->radial_workspace_p->ncheb;
      info->nf1 = w->nf1;
      info->nf2 = w->nf2;
      info->nrows = w->nrows;
      info->grid_bytes = w->nrows * w->ncols * w->nchp * sizeof(double);
      info->err_radial = w->radial_workspace_p->err;
      info->err_check = w->check_err;
    }

  info->ncomp = w->ncomp;

  return POTSYNTH_SUCCESS;
}

/*
auto_method()
  Choose the evaluation method for POTSYNTH_METHOD_AUTO

Inputs: lmax   - maximum SH degree
        nfunc  - number of functionals
        params - parameters (eps, npts_hint)

Return: POTSYNTH_METHOD_FAST or POTSYNTH_METHOD_DIRECT

Notes:
1) Without a hint (npts_hint = 0), the fast method is used for
lmax >= POTSYNTH_AUTO_LMIN

2) With a hint, the method with the smaller estimated total time for
npts_hint points is chosen:

fast:   t_setup + n t_fast,   direct: n t_direct

with the empirical cost models (seconds; measured on 16 threads, consistent
with the measured break-even points to within about a factor 2, and near the
break-even both methods cost about the same):

t_setup  = 6.3e-12 nch L^3 + 7.8e-8 nch L^2
t_direct = 3.5e-11 L^2 (1 + (nfunc - 1)/2)
t_fast   = 1.5e-11 w^2 nch + 5e-8

where nch = nfunc * rho is the number of grid functionals; the radial rank is
estimated from eps as rho ~ 2 + 0.7 log10(1/eps)
*/

static int
auto_method(const size_t lmax, const size_t nfunc, const potsynth_parameters * params)
{
  if (params->npts_hint == 0)
    {
      return (lmax >= POTSYNTH_AUTO_LMIN) ? POTSYNTH_METHOD_FAST : POTSYNTH_METHOD_DIRECT;
    }
  else
    {
      const double L = (double) lmax;
      const double n = (double) params->npts_hint;
      const double sigma = (params->sigma > 0.0) ? params->sigma : 2.0;
      const double ww = (params->w > 0) ? (double) params->w :
                        (double) potsynth_window_choose(params->eps, sigma);
      const double rho = 2.0 + 0.7 * fmax(0.0, -log10(params->eps));
      const double nch = (double) nfunc * rho;
      const double t_setup = 6.3e-12 * nch * L * L * L + 7.8e-8 * nch * L * L;
      const double t_direct = 3.5e-11 * L * L * (1.0 + 0.5 * ((double) nfunc - 1.0));
      const double t_fast = 1.5e-11 * ww * ww * nch + 5.0e-8;

      return (t_setup + n * t_fast < n * t_direct) ? POTSYNTH_METHOD_FAST : POTSYNTH_METHOD_DIRECT;
    }
}

/*
potsynth_set_outputs()
  Define the outputs of a model as linear combinations of its functionals; used
by the front-ends (e.g. 6 tensor components from 5 functionals)

Inputs: w    - workspace
        nout - number of outputs per point
        map  - output k = sum_c map[k*nfunc + c] functional c, length nout*nfunc; copied

Return: success or error

Notes:
1) potsynth_eval() then writes nout outputs per point, out[k*n + i]
*/

int
potsynth_set_outputs(potsynth_workspace * w, const size_t nout, const double * map)
{
  if (nout < w->nfunc || map == NULL)
    {
      POTSYNTH_ERROR ("invalid output map", POTSYNTH_EINVAL);
    }

  w->omap = malloc(nout * w->nfunc * sizeof(double));
  if (w->omap == NULL)
    {
      POTSYNTH_ERROR ("unable to allocate output map", POTSYNTH_ENOMEM);
    }

  memcpy(w->omap, map, nout * w->nfunc * sizeof(double));
  w->nout = nout;

  return POTSYNTH_SUCCESS;
}

/*
potsynth_tensor_functionals()
  Define the five functionals of the gradient tensor T = grad grad V of
V = sum (a/r)^(l+1) Pbar_l^m [C cos(m phi) + S sin(m phi)], scaled by f:

T_rr = f (l+1)(l+2) u^(l+3) P,      T_rt = -f (l+2) u^(l+3) dP/dtheta,
T_rp = -f (l+2) u^(l+3) P/sin d/dphi,
D    = T_tt - T_pp = f u^(l+3) S2A,  T_tp = f u^(l+3) S2B d/dphi

Inputs: lmax - maximum SH degree
        f    - scale factor (1/a^2 for the unit potential above)
        func - (output) 5 functionals
        rfac - (output) per-degree factors referenced by func, length 5*(lmax+1)

Notes:
1) See the documentation, "The gradient tensor"; potsynth_tensor_map()
gives the six components from the five functionals
*/

void
potsynth_tensor_functionals(const size_t lmax, const double f, potsynth_functional * func, double * rfac)
{
  const size_t L1 = lmax + 1;
  const int ops[5] = { POTSYNTH_THETA_P, POTSYNTH_THETA_DP, POTSYNTH_THETA_PSIN,
                       POTSYNTH_THETA_S2A, POTSYNTH_THETA_S2B };
  const int dphi[5] = { 0, 0, 1, 0, 1 };
  size_t c, l;

  for (l = 0; l <= lmax; ++l)
    {
      const double dl = (double) l;

      rfac[0 * L1 + l] = f * (dl + 1.0) * (dl + 2.0);
      rfac[1 * L1 + l] = -f * (dl + 2.0);
      rfac[2 * L1 + l] = -f * (dl + 2.0);
      rfac[3 * L1 + l] = f;
      rfac[4 * L1 + l] = f;
    }

  for (c = 0; c < 5; ++c)
    {
      func[c].rpow = 3;
      func[c].rfac = rfac + c * L1;
      func[c].theta = ops[c];
      func[c].dphi = dphi[c];
    }
}

/*
potsynth_tensor_map()
  Rows of an output map giving the six tensor components
(rr, rt, rp, tt, tp, pp) from the five functionals of potsynth_tensor_functionals():
T_tt = (-T_rr + D)/2, T_pp = (-T_rr - D)/2 (zero trace)

Inputs: nfunc - number of functionals of the model
        c0    - index of the first tensor functional
        map   - (output) 6 rows of length nfunc, zeroed by the caller
*/

void
potsynth_tensor_map(const size_t nfunc, const size_t c0, double * map)
{
  map[0 * nfunc + c0] = 1.0;
  map[1 * nfunc + c0 + 1] = 1.0;
  map[2 * nfunc + c0 + 2] = 1.0;
  map[3 * nfunc + c0] = -0.5;
  map[3 * nfunc + c0 + 3] = 0.5;
  map[4 * nfunc + c0 + 4] = 1.0;
  map[5 * nfunc + c0] = -0.5;
  map[5 * nfunc + c0 + 3] = -0.5;
}
