/*
 * radial.c
 *
 * Optimal radial basis (step 1 of the fast method)
 *
 * The radial factors u^(l+rpow), u = a/r, l = lmin..lmax, are sampled at K+1
 * Chebyshev-Lobatto points of [u0, u1], weighted by the size of their
 * contribution to the field, and compressed with an SVD into rho functions
 * b_i(u), stored as Chebyshev series (the rank rho is the smallest whose pointwise
 * error bound, measured between the nodes too, is within the tolerance):
 *
 * u^(l+rpow_c) ~ \sum_{i=0}^{rho-1} beta_{c,i,l} b_i(u)
 *
 * The coefficients beta follow from the discrete (Chebyshev-weighted)
 * projection. The SVD is a one-sided Jacobi method, which computes small
 * singular values to high relative accuracy (needed for machine precision
 * operation).
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "potsynth_internal.h"

/* temporary arrays of the radial set-up */
typedef struct
{
  size_t nfun;          /* number of weighted functions */
  int * fpow;           /* power of u of each function, length nfun */
  double * fw;          /* weight of each function, length nfun */
  size_t K;             /* Chebyshev degree */
  double * A;           /* samples A[f*(K+1) + j] */
  double * ca;          /* Chebyshev coefficients, same layout as A */
  double * ct;          /* cos(pi i / K), length 2K */
  double * uj;          /* nodes u_j, length K+1 */
  double * sq;          /* sqrt of Lobatto weights, length K+1 */
  double * X;           /* SVD matrix, column-major nfun-by-(K+1) */
  double * V;           /* right singular vectors, (K+1)-by-(K+1) */
  double * sv;          /* singular values, length K+1 */
  double floor;         /* rounding level of the tail measure of radial_sample() */
} radial_state_t;

static void radial_state_free(radial_state_t * state);
static int radial_state_resize(const size_t K, radial_state_t * state);
static void cheb_coeffs(const size_t K, const double * f, const double * ct, double * a);
static double radial_sample(const double uc, const double ud, radial_state_t * state);
static int radial_build(const size_t lmin, const size_t lmax, const size_t nfunc, const int * rpow,
                        const double * wn, const double tol, const int maxpow, const size_t kmin,
                        potsynth_radial_workspace * w, radial_state_t * state);
static int radial_rank_bounds(const size_t lmin, const size_t lmax, const size_t nfunc, const int * rpow,
                              const double * wn, const potsynth_radial_workspace * w, double * B);

/* first functional with the radial power of functional c (c itself if none before) */
static size_t
radial_first(const size_t c, const int * rpow)
{
  size_t c2;

  for (c2 = 0; c2 < c && rpow[c2] != rpow[c]; ++c2)
    ;

  return c2;
}

/* 1 if an earlier functional has the radial power of functional c and, when wn is given,
   the same weights wn[c*L1 + l], l >= lmin: its rows need not be computed again */
static int
radial_repeat(const size_t c, const size_t lmin, const size_t L1, const int * rpow, const double * wn)
{
  size_t c2;

  for (c2 = 0; c2 < c; ++c2)
    {
      if (rpow[c2] == rpow[c] &&
          (wn == NULL || memcmp(wn + c2 * L1 + lmin, wn + c * L1 + lmin, (L1 - lmin) * sizeof(double)) == 0))
        return 1;
    }

  return 0;
}

/*
potsynth_radial_alloc()
  Compute the optimal radial basis for a set of functionals

Inputs: u0    - a / rmax
        u1    - a / rmin
        lmin  - minimum SH degree
        lmax  - maximum SH degree
        nfunc - number of functionals
        rpow  - radial power offset of each functional, length nfunc
        wn    - weights wn[c*(lmax+1) + l]: size of the contribution of degree l
                to functional c (0 = absent)
        tol   - relative tolerance of the rank truncation

Return: pointer to workspace

Notes:
1) If rmin = rmax (to rounding), or all functions are constant, a single
constant basis function b_0 = 1 is used
*/

potsynth_radial_workspace *
potsynth_radial_alloc(const double u0, const double u1, const size_t lmin, const size_t lmax,
                      const size_t nfunc, const int * rpow, const double * wn, const double tol)
{
  const size_t L1 = lmax + 1;
  int maxpow = rpow[0], minpow = rpow[0];
  potsynth_radial_workspace * w;
  radial_state_t state;
  size_t c, l;
  int status;

  for (c = 1; c < nfunc; ++c)
    {
      maxpow = (rpow[c] > maxpow) ? rpow[c] : maxpow;
      minpow = (rpow[c] < minpow) ? rpow[c] : minpow;
    }

  if ((long) lmin + minpow < 0)
    {
      POTSYNTH_ERROR_NULL ("negative powers of a/r are not supported", POTSYNTH_EINVAL);
    }
  else if (!isfinite(pow(u1, (double) lmax + maxpow)))
    {
      POTSYNTH_ERROR_NULL ("(a/r)^lmax overflows: rmin is too small", POTSYNTH_EINVAL);
    }

  w = calloc(1, sizeof(potsynth_radial_workspace));
  if (w == NULL)
    {
      POTSYNTH_ERROR_NULL ("unable to allocate radial workspace", POTSYNTH_ENOMEM);
    }

  w->uc = 0.5 * (u0 + u1);
  w->ud = 0.5 * (u1 - u0);

  /* thin shell, or constant functions: a single radial function b_0 = 1 */
  if (!(w->ud > 1.0e-15 * w->uc) || (long) lmax + maxpow == 0)
    {
      w->rho = 1;
      w->ncheb = 1;
      w->ud = 0.0;
      w->err = 0.0;
      w->basis = potsynth_aligned_alloc(sizeof(double));
      w->beta = potsynth_aligned_calloc(nfunc * L1 * sizeof(double));
      if (!w->basis || !w->beta)
        {
          potsynth_radial_free(w);
          POTSYNTH_ERROR_NULL ("unable to allocate radial basis", POTSYNTH_ENOMEM);
        }

      w->basis[0] = 1.0;
      for (c = 0; c < nfunc; ++c)
        {
          for (l = lmin; l <= lmax; ++l)
            w->beta[c * L1 + l] = pow(w->uc, (double) l + rpow[c]);
        }

      return w;
    }

  {
    /* rank selection: build a basis of full rank (all K+1 singular vectors) on K+1
     * Chebyshev nodes, measure the pointwise bound B_r of every rank r on a fine set of
     * radii (radial_rank_bounds, between the nodes too), and take the smallest rank with
     * B_r <= tol. Then raise K by about a quarter and repeat while that lowers the rank:
     * too few nodes give a basis that interpolates rather than compresses, and near the
     * rounding floor more nodes only add rounding, so the smallest rank comes from the
     * smallest K that reaches it.
     *
     * If no rank reaches tol (the bound's rounding floor is above it), the basis is taken
     * at the knee of B_r: truncation errors fall by more than 4 times per rank and
     * rounding does not, so the knee is the smallest rank within 4 times the smallest
     * bound, followed by up to 4 more ranks when they halve the bound. The passes stop
     * when a pass does not halve the smallest bound. w->err is the bound of the basis
     * returned, which is above tol only in that case. */
    const size_t Kmax = (size_t) ((long) lmax + maxpow);
    size_t kmin = 0, iter, r;
    double prevmin = HUGE_VAL;
    /* candidates over the passes: the smallest rank within tol (t), and the floor
     * candidate (f) */
    potsynth_radial_workspace tc = { 0 }, fc = { 0 };
    size_t trho = 0, prev_trho = 0, frho = 0;
    double terr = 0.0, ferr = 0.0, fmin_ = HUGE_VAL;

    for (iter = 0; iter < 16; ++iter)
      {
        double bmin = HUGE_VAL, * B;
        size_t Kused, R, rt = 0, rf = 0;
        int keep_t, keep_f;

        memset(&state, 0, sizeof(radial_state_t));
        status = radial_build(lmin, lmax, nfunc, rpow, wn, tol, maxpow, kmin, w, &state);
        Kused = state.K;
        radial_state_free(&state);

        R = w->rho;
        B = (status == POTSYNTH_SUCCESS) ? malloc(R * sizeof(double)) : NULL;
        if (status == POTSYNTH_SUCCESS && (B == NULL || radial_rank_bounds(lmin, lmax, nfunc, rpow, wn, w, B)))
          status = POTSYNTH_ENOMEM;

        if (status)
          {
            free(B);
            potsynth_aligned_free(tc.basis);
            potsynth_aligned_free(tc.beta);
            potsynth_aligned_free(fc.basis);
            potsynth_aligned_free(fc.beta);
            potsynth_radial_free(w);
            if (status == POTSYNTH_ENOMEM)
              {
                POTSYNTH_ERROR_NULL ("unable to allocate the radial basis", POTSYNTH_ENOMEM);
              }
            return NULL;
          }

        for (r = 0; r < R; ++r)
          bmin = fmin(bmin, B[r]);
        for (r = 0; r < R && !rt; ++r)
          if (B[r] <= tol)
            rt = r + 1;
        {
          /* knee, then up to 4 more ranks that halve the bound */
          size_t rk = 0;
          double m2 = HUGE_VAL;
          for (r = 0; r < R && !rk; ++r)
            if (B[r] <= 4.0 * bmin)
              rk = r + 1;
          for (r = 0; r < R && r < rk + 4; ++r)
            m2 = fmin(m2, B[r]);
          for (r = 0; r < R && !rf; ++r)
            if (B[r] <= 2.0 * m2)
              rf = r + 1;
        }

        keep_t = (rt > 0 && (trho == 0 || rt < trho));
        keep_f = (frho == 0 || bmin < 0.5 * fmin_ || (rf < frho && B[rf - 1] <= 2.0 * fmin_));

        /* keep the candidates of this pass (each owns a copy of the full basis) */
        if (keep_t || keep_f)
          {
            potsynth_radial_workspace * dst[2];
            size_t nd = 0, d;
            if (keep_t)
              dst[nd++] = &tc;
            if (keep_f)
              dst[nd++] = &fc;
            for (d = 0; d < nd; ++d)
              {
                potsynth_aligned_free(dst[d]->basis);
                potsynth_aligned_free(dst[d]->beta);
                dst[d]->basis = potsynth_aligned_alloc(R * w->ncheb * sizeof(double));
                dst[d]->beta = potsynth_aligned_alloc(nfunc * R * L1 * sizeof(double));
                if (!dst[d]->basis || !dst[d]->beta)
                  {
                    status = POTSYNTH_ENOMEM;
                    break;
                  }
                memcpy(dst[d]->basis, w->basis, R * w->ncheb * sizeof(double));
                memcpy(dst[d]->beta, w->beta, nfunc * R * L1 * sizeof(double));
                dst[d]->rho = R;
                dst[d]->ncheb = w->ncheb;
              }
            if (keep_t)
              {
                trho = rt;
                terr = B[rt - 1];
              }
            if (keep_f)
              {
                frho = rf;
                ferr = B[rf - 1];
              }
          }
        fmin_ = fmin(fmin_, bmin);
        free(B);
        potsynth_aligned_free(w->basis);
        potsynth_aligned_free(w->beta);
        w->basis = NULL;
        w->beta = NULL;

        if (status)
          {
            potsynth_aligned_free(tc.basis);
            potsynth_aligned_free(tc.beta);
            potsynth_aligned_free(fc.basis);
            potsynth_aligned_free(fc.beta);
            potsynth_radial_free(w);
            POTSYNTH_ERROR_NULL ("unable to allocate the radial basis", POTSYNTH_ENOMEM);
          }

        /* stop: no more nodes to add; a rank within tol that more nodes did not lower;
         * or no rank within tol and a pass that did not halve the smallest bound */
        if (Kused >= Kmax)
          break;
        if (trho && prev_trho && !(trho < prev_trho))
          break;
        if (!trho && iter > 0 && !(bmin < 0.5 * prevmin))
          break;

        prev_trho = trho;
        prevmin = fmin(prevmin, bmin);
        kmin = (Kused + 2 > (5 * Kused + 3) / 4) ? Kused + 2 : (5 * Kused + 3) / 4;
      }

    /* the chosen rank: its basis functions and coefficients */
    {
      potsynth_radial_workspace * cand = trho ? &tc : &fc;
      const size_t rho = trho ? trho : frho, R0 = cand->rho, nc = cand->ncheb;

      w->basis = potsynth_aligned_alloc(rho * nc * sizeof(double));
      w->beta = potsynth_aligned_alloc(nfunc * rho * L1 * sizeof(double));
      if (w->basis && w->beta)
        {
          memcpy(w->basis, cand->basis, rho * nc * sizeof(double));
          for (c = 0; c < nfunc; ++c)
            {
              size_t i;
              for (i = 0; i < rho; ++i)
                memcpy(w->beta + (c * rho + i) * L1, cand->beta + (c * R0 + i) * L1, L1 * sizeof(double));
            }
          w->rho = rho;
          w->ncheb = nc;
          w->err = trho ? terr : ferr;
        }
      potsynth_aligned_free(tc.basis);
      potsynth_aligned_free(tc.beta);
      potsynth_aligned_free(fc.basis);
      potsynth_aligned_free(fc.beta);
      if (!w->basis || !w->beta)
        {
          potsynth_radial_free(w);
          POTSYNTH_ERROR_NULL ("unable to allocate the radial basis", POTSYNTH_ENOMEM);
        }
    }
  }

  return w;
}

/*
potsynth_radial_alloc_caller()
  Radial workspace of a model with radial functions of the caller: nrad
functions that are evaluated by the caller's function at each point, so that
there is nothing to expand: all coefficients beta are 1

Inputs: nrad  - number of radial functions
        lmax  - maximum degree
        nfunc - number of functionals

Return: pointer to workspace, or NULL on error

Notes:
1) basis is 1 for each function (one Chebyshev term), which is only used to
rank the nodes of the window check
*/

potsynth_radial_workspace *
potsynth_radial_alloc_caller(const size_t nrad, const size_t lmax, const size_t nfunc)
{
  const size_t nbeta = nfunc * nrad * (lmax + 1);
  potsynth_radial_workspace * w;
  size_t i;

  w = calloc(1, sizeof(potsynth_radial_workspace));
  if (w == NULL)
    {
      POTSYNTH_ERROR_NULL ("unable to allocate radial workspace", POTSYNTH_ENOMEM);
    }

  w->err = -1.0; /* the caller's functions are used as given */
  w->uc = 1.0;
  w->ud = 0.0;
  w->rho = nrad;
  w->ncheb = 1;
  w->basis = potsynth_aligned_calloc(nrad * sizeof(double));
  w->beta = potsynth_aligned_calloc(nbeta * sizeof(double));
  if (w->basis == NULL || w->beta == NULL)
    {
      potsynth_radial_free(w);
      POTSYNTH_ERROR_NULL ("unable to allocate radial workspace", POTSYNTH_ENOMEM);
    }

  for (i = 0; i < nrad; ++i)
    w->basis[i] = 1.0;

  for (i = 0; i < nbeta; ++i)
    w->beta[i] = 1.0;

  return w;
}

void
potsynth_radial_free(potsynth_radial_workspace * w)
{
  if (w == NULL)
    return; /* as free(NULL) */

  if (w->basis)
    potsynth_aligned_free(w->basis);

  if (w->beta)
    potsynth_aligned_free(w->beta);

  free(w);
}

/*
potsynth_svd_jacobi()
  One-sided Jacobi SVD of an n-by-p matrix

Inputs: n  - number of rows
        p  - number of columns
        X  - on input, matrix (column-major); on output, U*S (orthogonal columns)
        V  - (output) right singular vectors, p-by-p column-major
        sv - (output) singular values, length p, in decreasing order

Return: success or error

Notes:
1) Columns below (eps ||X||_F)^2 are rounding noise and are not rotated
against each other (they would never satisfy the relative convergence test)
*/

int
potsynth_svd_jacobi(const size_t n, const size_t p, double * X, double * V, double * sv)
{
  const double tol = DBL_EPSILON * sqrt((double) n);
  double fro2 = 0.0, tiny;
  size_t i, j, k, sweep;

  for (i = 0; i < n * p; ++i)
    fro2 += X[i] * X[i];
  tiny = fro2 * DBL_EPSILON * DBL_EPSILON;

  for (i = 0; i < p * p; ++i)
    V[i] = 0.0;
  for (i = 0; i < p; ++i)
    V[i * p + i] = 1.0;

  for (sweep = 0; sweep < 100; ++sweep)
    {
      int rotated = 0;

      for (i = 0; i + 1 < p; ++i)
        {
          for (j = i + 1; j < p; ++j)
            {
              double * xi = X + i * n;
              double * xj = X + j * n;
              double * vi = V + i * p;
              double * vj = V + j * p;
              double alpha = 0.0, beta = 0.0, gamma = 0.0;
              double zeta, t, c, s;

              for (k = 0; k < n; ++k)
                {
                  alpha += xi[k] * xi[k];
                  beta += xj[k] * xj[k];
                  gamma += xi[k] * xj[k];
                }

              if (alpha <= tiny || beta <= tiny || fabs(gamma) <= tol * sqrt(alpha * beta))
                continue;

              rotated = 1;
              zeta = (beta - alpha) / (2.0 * gamma);
              t = ((zeta >= 0.0) ? 1.0 : -1.0) / (fabs(zeta) + sqrt(1.0 + zeta * zeta));
              c = 1.0 / sqrt(1.0 + t * t);
              s = c * t;

              for (k = 0; k < n; ++k)
                {
                  const double xik = xi[k], xjk = xj[k];

                  xi[k] = c * xik - s * xjk;
                  xj[k] = s * xik + c * xjk;
                }

              for (k = 0; k < p; ++k)
                {
                  const double vik = vi[k], vjk = vj[k];

                  vi[k] = c * vik - s * vjk;
                  vj[k] = s * vik + c * vjk;
                }
            }
        }

      if (!rotated)
        break;
    }

  for (i = 0; i < p; ++i)
    {
      double sum = 0.0;

      for (k = 0; k < n; ++k)
        sum += X[i * n + k] * X[i * n + k];

      sv[i] = sqrt(sum);
    }

  /* selection sort into decreasing order (p is small) */
  for (i = 0; i < p; ++i)
    {
      size_t imax = i;

      for (j = i + 1; j < p; ++j)
        {
          if (sv[j] > sv[imax])
            imax = j;
        }

      if (imax != i)
        {
          double tmp = sv[i];

          sv[i] = sv[imax];
          sv[imax] = tmp;

          for (k = 0; k < n; ++k)
            {
              tmp = X[i * n + k];
              X[i * n + k] = X[imax * n + k];
              X[imax * n + k] = tmp;
            }

          for (k = 0; k < p; ++k)
            {
              tmp = V[i * p + k];
              V[i * p + k] = V[imax * p + k];
              V[imax * p + k] = tmp;
            }
        }
    }

  if (sweep == 100)
    {
      POTSYNTH_ERROR ("Jacobi SVD did not converge", POTSYNTH_EINVAL);
    }

  return POTSYNTH_SUCCESS;
}

static void
radial_state_free(radial_state_t * state)
{
  if (state->fpow)
    free(state->fpow);

  if (state->fw)
    free(state->fw);

  if (state->A)
    free(state->A);

  if (state->ca)
    free(state->ca);

  if (state->ct)
    free(state->ct);

  if (state->uj)
    free(state->uj);

  if (state->sq)
    free(state->sq);

  if (state->X)
    free(state->X);

  if (state->V)
    free(state->V);

  if (state->sv)
    free(state->sv);
}

/* (re)allocate the arrays which depend on the Chebyshev degree K */
static int
radial_state_resize(const size_t K, radial_state_t * state)
{
  const size_t nfun = state->nfun;

  free(state->A);
  free(state->ca);
  free(state->ct);
  free(state->uj);
  free(state->sq);

  state->K = K;
  state->A = malloc(nfun * (K + 1) * sizeof(double));
  state->ca = malloc(nfun * (K + 1) * sizeof(double));
  state->ct = malloc(2 * K * sizeof(double));
  state->uj = malloc((K + 1) * sizeof(double));
  state->sq = malloc((K + 1) * sizeof(double));
  if (!state->A || !state->ca || !state->ct || !state->uj || !state->sq)
    {
      POTSYNTH_ERROR ("unable to allocate radial arrays", POTSYNTH_ENOMEM);
    }

  return POTSYNTH_SUCCESS;
}

/*
cheb_coeffs()
  Chebyshev coefficients a[0..K] of the interpolant through f[j] at the
Lobatto points t_j = cos(pi j / K), j = 0..K (K >= 1); ct[i] = cos(pi i / K)
*/

static void
cheb_coeffs(const size_t K, const double * f, const double * ct, double * a)
{
  size_t j, k;

  for (k = 0; k <= K; ++k)
    {
      double sum = 0.5 * (f[0] + ((k & 1) ? -f[K] : f[K]));

      for (j = 1; j < K; ++j)
        sum += f[j] * ct[(j * k) % (2 * K)];

      a[k] = (2.0 / K) * sum;
    }

  a[0] *= 0.5;
  a[K] *= 0.5;
}

/*
radial_sample()
  Sample all weighted functions at the K+1 Lobatto points and compute their
Chebyshev coefficients

Inputs: uc    - center of u interval
        ud    - half width of u interval
        state - state

Return: tail measure sum_f sqrt(2 l_f + 1) (|a_{K-1}| + |a_K|) relative to
        max_j (sum_f f(u_j)^2)^(1/2), with l_f the power of function f (the
        sqrt(2l+1) factor bounds the maximum over the sphere, see
        radial_rank_bounds())
*/

static double
radial_sample(const double uc, const double ud, radial_state_t * state)
{
  const size_t K = state->K;
  double tail = 0.0, scale = 0.0;
  size_t f, j;

  for (j = 0; j < 2 * K; ++j)
    state->ct[j] = cos(M_PI * (double) j / (double) K);

  for (j = 0; j <= K; ++j)
    {
      state->uj[j] = uc + ud * state->ct[j];
      state->sq[j] = (j == 0 || j == K) ? sqrt(0.5) : 1.0;
    }

  for (f = 0; f < state->nfun; ++f)
    {
      double * A = state->A + f * (K + 1);
      double * ca = state->ca + f * (K + 1);

      for (j = 0; j <= K; ++j)
        A[j] = state->fw[f] * pow(state->uj[j], (double) state->fpow[f]);

      cheb_coeffs(K, A, state->ct, ca);
      tail += sqrt(2.0 * state->fpow[f] + 1.0) * (fabs(ca[K]) + fabs(ca[K - 1]));
    }

  for (j = 0; j <= K; ++j)
    {
      double s2 = 0.0;
      for (f = 0; f < state->nfun; ++f)
        s2 += state->A[f * (K + 1) + j] * state->A[f * (K + 1) + j];
      scale = fmax(scale, sqrt(s2));
    }

  /* rounding level of the tail: each Chebyshev coefficient of f is computed with an
   * absolute error of about DBL_EPSILON max_j |f(u_j)| */
  {
    double fl = 0.0;
    for (f = 0; f < state->nfun; ++f)
      {
        double fm = 0.0;
        for (j = 0; j <= K; ++j)
          fm = fmax(fm, fabs(state->A[f * (K + 1) + j]));
        fl += sqrt(2.0 * state->fpow[f] + 1.0) * fm;
      }
    state->floor = (scale > 0.0) ? 2.0 * DBL_EPSILON * fl / scale : 0.0;
  }

  return (scale > 0.0) ? tail / scale : 0.0;
}

/*
radial_build()
  Build the radial basis

Inputs: lmin   - minimum SH degree
        lmax   - maximum SH degree
        nfunc  - number of functionals
        rpow   - radial power offsets, length nfunc
        wn     - weights, wn[c*(lmax+1) + l]
        tol    - truncation tolerance (Chebyshev degree and trailing terms)
        maxpow - largest rpow
        kmin   - smallest Chebyshev degree K
        w      - (output) radial workspace: all K+1 singular vectors (rho = K+1);
                 the caller chooses the rank
        state  - temporary arrays
*/

static int
radial_build(const size_t lmin, const size_t lmax, const size_t nfunc, const int * rpow,
             const double * wn, const double tol, const int maxpow, const size_t kmin,
             potsynth_radial_workspace * w, radial_state_t * state)
{
  const size_t L1 = lmax + 1;
  const size_t Kmax = (size_t) ((long) lmax + maxpow); /* exact polynomial degree */
  size_t K, Kopt, f, c, i, j, k, l, ncheb, rho, nfun = 0;
  double bmax = 0.0;
  int status;

  /* weighted functions w * u^(l+rpow): one per (distinct rpow, l), with the
   * largest weight over the functionals sharing that rpow */
  state->fpow = malloc(nfunc * L1 * sizeof(int));
  state->fw = malloc(nfunc * L1 * sizeof(double));
  if (!state->fpow || !state->fw)
    {
      POTSYNTH_ERROR ("unable to allocate radial functions", POTSYNTH_ENOMEM);
    }

  for (c = 0; c < nfunc; ++c)
    {
      size_t c2;
      if (radial_repeat(c, lmin, L1, rpow, NULL))
        continue;

      for (l = lmin; l <= lmax; ++l)
        {
          double wmax = 0.0;

          for (c2 = c; c2 < nfunc; ++c2)
            {
              if (rpow[c2] == rpow[c])
                wmax = fmax(wmax, wn[c2 * L1 + l]);
            }

          if (wmax > 0.0)
            {
              state->fpow[nfun] = (int) l + rpow[c];
              state->fw[nfun] = wmax;
              ++nfun;
            }
        }
    }

  if (nfun == 0)
    {
      /* model is identically zero: keep a single (unused) constant function */
      state->fpow[0] = 0;
      state->fw[0] = 1.0;
      nfun = 1;
    }

  state->nfun = nfun;

  /* Chebyshev degree: double K until the series are resolved to tol/2 */
  size_t prevK = 0;
  double prevtail = 0.0;

  K = (Kmax < 8) ? Kmax : 8;
  if (K < kmin)
    K = (kmin < Kmax) ? kmin : Kmax;
  if (K < 1)
    K = 1;

  while (1)
    {
      double tail;

      status = radial_state_resize(K, state);
      if (status)
        return status;

      tail = radial_sample(w->uc, w->ud, state);
      if (tail <= 0.5 * tol || K >= Kmax)
        break;

      /* rounding: a resolved function's Chebyshev tail falls much faster than 4 times
       * per doubling; if it did not and the tail is near its rounding level, the tail is
       * rounding noise: keep the previous K. Far above that level a slow fall only means
       * the functions are not resolved yet (thick shells at high degree) */
      if (prevK > 0 && tail > 0.25 * prevtail && tail <= 100.0 * state->floor)
        {
          K = prevK;
          status = radial_state_resize(K, state);
          if (status)
            return status;
          radial_sample(w->uc, w->ud, state);
          break;
        }

      prevK = K;
      prevtail = tail;
      K = (2 * K < Kmax) ? 2 * K : Kmax;
    }

  /* smallest K whose truncation tail is below tol/2 */
  {
    double scale = 0.0, tail = 0.0;

    for (j = 0; j <= K; ++j)
      {
        double s2 = 0.0;
        for (f = 0; f < nfun; ++f)
          s2 += state->A[f * (K + 1) + j] * state->A[f * (K + 1) + j];
        scale = fmax(scale, sqrt(s2));
      }

    Kopt = K;
    for (k = K; k >= 2; --k)
      {
        for (f = 0; f < nfun; ++f)
          tail += sqrt(2.0 * state->fpow[f] + 1.0) * fabs(state->ca[f * (K + 1) + k]);

        if (tail > 0.5 * tol * scale || k <= kmin)
          break;

        Kopt = k - 1;
      }
  }

  if (Kopt < K)
    {
      K = Kopt;
      status = radial_state_resize(K, state);
      if (status)
        return status;

      radial_sample(w->uc, w->ud, state);
    }

  /* weighted SVD: column j of X = sqrt(c_j) * samples at node j */
  state->X = malloc(nfun * (K + 1) * sizeof(double));
  state->V = malloc((K + 1) * (K + 1) * sizeof(double));
  state->sv = malloc((K + 1) * sizeof(double));
  if (!state->X || !state->V || !state->sv)
    {
      POTSYNTH_ERROR ("unable to allocate SVD arrays", POTSYNTH_ENOMEM);
    }

  for (j = 0; j <= K; ++j)
    {
      for (f = 0; f < nfun; ++f)
        state->X[j * nfun + f] = state->sq[j] * state->A[f * (K + 1) + j];
    }

  status = potsynth_svd_jacobi(nfun, K + 1, state->X, state->V, state->sv);
  if (status)
    return status;

  rho = K + 1; /* all singular vectors: the caller chooses the rank */

  /* basis functions b_i(u_j) = U_ji / sqrt(c_j), U_ji = V[i*(K+1) + j] */
  w->basis = potsynth_aligned_calloc(rho * (K + 1) * sizeof(double));
  w->beta = potsynth_aligned_calloc(nfunc * rho * L1 * sizeof(double));
  if (!w->basis || !w->beta)
    {
      POTSYNTH_ERROR ("unable to allocate radial basis", POTSYNTH_ENOMEM);
    }

  for (i = 0; i < rho; ++i)
    {
      double * bval = state->A; /* reuse as scratch, length K+1 */

      for (j = 0; j <= K; ++j)
        bval[j] = state->V[i * (K + 1) + j] / state->sq[j];

      cheb_coeffs(K, bval, state->ct, w->basis + i * (K + 1));

      for (k = 0; k <= K; ++k)
        bmax = fmax(bmax, fabs(w->basis[i * (K + 1) + k]));
    }

  /* drop trailing Chebyshev terms below the tolerance (rounding noise) */
  ncheb = 1;
  for (k = 0; k <= K; ++k)
    {
      for (i = 0; i < rho; ++i)
        {
          if (fabs(w->basis[i * (K + 1) + k]) > 0.1 * tol * bmax)
            ncheb = k + 1;
        }
    }

  if (ncheb < K + 1)
    {
      for (i = 0; i < rho; ++i)
        memmove(w->basis + i * ncheb, w->basis + i * (K + 1), ncheb * sizeof(double));
    }

  w->rho = rho;
  w->ncheb = ncheb;

  /* the powers are computed once per degree (not once per radial function), and the
   * rows of a functional with the radial power of an earlier one are copied from it */
  {
    double * pw = malloc((K + 1) * sizeof(double)); /* u_j^(l+rpow_c) */

    if (pw == NULL)
      {
        POTSYNTH_ERROR ("unable to allocate radial powers", POTSYNTH_ENOMEM);
      }

    for (c = 0; c < nfunc; ++c)
      {
        const size_t c2 = radial_first(c, rpow);

        if (c2 < c)
          {
            memcpy(w->beta + c * rho * L1, w->beta + c2 * rho * L1, rho * L1 * sizeof(double));
            continue;
          }

        for (l = lmin; l <= lmax; ++l)
          {
            for (j = 0; j <= K; ++j)
              pw[j] = pow(state->uj[j], (double) l + rpow[c]);

            for (i = 0; i < rho; ++i)
              {
                double sum = 0.0;

                for (j = 0; j <= K; ++j)
                  sum += state->V[i * (K + 1) + j] * state->sq[j] * pw[j];

                w->beta[(c * rho + i) * L1 + l] = sum;
              }
          }
      }

    free(pw);
  }

  return POTSYNTH_SUCCESS;
}

/*
radial_rank_bounds()
  Pointwise bound of the radial error of every rank, relative to the size of each
component, on a fine set of radii

Inputs: lmin  - minimum SH degree
        lmax  - maximum SH degree
        nfunc - number of functionals
        rpow  - radial power offsets, length nfunc
        wn    - weights, wn[c*(lmax+1) + l]
        w     - radial basis of rank w->rho
        B     - (output) B[r-1], r = 1..w->rho: largest bound over the components with
                the first r basis functions

Return: 0, or -1 if memory could not be allocated

Notes:
1) An irreducible component is a set of functionals with the same radial power and
weights. By the addition theorem (also for spin-weighted harmonics), a degree-l part
with rms w_l over the sphere is at most sqrt(2l+1) w_l anywhere on it, so the radial
error of a component at u is at most
B(u) = sum_l sqrt(2l+1) w_l |u^(l+p) - sum_{i<r} beta_il b_i(u)|, at every point of the
sphere and in any coordinate system. B is relative to the size
F = max_u (sum_l w_l^2 u^(2(l+p)))^(1/2), at most the largest magnitude of the component
in the shell.
2) B is evaluated with the basis as stored (Chebyshev series, double precision, as the
fast method uses it) at 8 ncheb + 1 Chebyshev points of the shell, so it includes the
error between the SVD nodes and the rounding of the expansion; its rounding floor is
the floor of the radial step
3) The partial sums over i are accumulated rank by rank, so all ranks cost about as
much as one; the loops over the points are innermost (vectorized)
*/

static int
radial_rank_bounds(const size_t lmin, const size_t lmax, const size_t nfunc, const int * rpow,
                   const double * wn, const potsynth_radial_workspace * w, double * B)
{
  const size_t L1 = lmax + 1, R = w->rho, nc = w->ncheb;
  const size_t M = 8 * nc + 1;
  double * bu = malloc(R * M * sizeof(double));  /* b_r at point m: bu[r*M + m] */
  double * err = malloc(R * M * sizeof(double)); /* B of rank r+1 at point m */
  double * uu = malloc(M * sizeof(double));
  double * f = malloc(M * sizeof(double));
  double * a = malloc(M * sizeof(double));
  double * f2 = malloc(M * sizeof(double));
  double * T = malloc(nc * sizeof(double));
  size_t c, m, i, k, l, r;

  if (!bu || !err || !uu || !f || !a || !f2 || !T)
    {
      free(bu);
      free(err);
      free(uu);
      free(f);
      free(a);
      free(f2);
      free(T);
      return -1;
    }

  for (r = 0; r < R; ++r)
    B[r] = 0.0;

  for (m = 0; m < M; ++m)
    {
      const double t = cos(M_PI * (double) m / (double) (M - 1));

      uu[m] = w->uc + w->ud * t;
      T[0] = 1.0;
      if (nc > 1)
        T[1] = t;
      for (k = 2; k < nc; ++k)
        T[k] = 2.0 * t * T[k - 1] - T[k - 2];

      for (i = 0; i < R; ++i)
        {
          double sum = 0.0;
          for (k = 0; k < nc; ++k)
            sum += w->basis[i * nc + k] * T[k];
          bu[i * M + m] = sum;
        }
    }

  for (c = 0; c < nfunc; ++c)
    {
      const double * wc = wn + c * L1;
      double fmax2 = 0.0;

      if (radial_repeat(c, lmin, L1, rpow, wn))
        continue; /* same component as an earlier functional */

      memset(err, 0, R * M * sizeof(double));
      memset(f2, 0, M * sizeof(double));

      for (l = lmin; l <= lmax; ++l)
        {
          const double sw = sqrt(2.0 * l + 1.0) * fabs(wc[l]);

          if (wc[l] == 0.0)
            continue;

          for (m = 0; m < M; ++m)
            {
              f[m] = pow(uu[m], (double) l + rpow[c]);
              a[m] = 0.0;
              f2[m] += wc[l] * wc[l] * f[m] * f[m];
            }

          for (r = 0; r < R; ++r)
            {
              const double be = w->beta[(c * R + r) * L1 + l];
              const double * b = bu + r * M;
              double * e = err + r * M;

              for (m = 0; m < M; ++m)
                {
                  a[m] += be * b[m];
                  e[m] += sw * fabs(f[m] - a[m]);
                }
            }
        }

      for (m = 0; m < M; ++m)
        fmax2 = fmax(fmax2, f2[m]);

      if (fmax2 > 0.0)
        {
          for (r = 0; r < R; ++r)
            {
              double bm = 0.0;
              for (m = 0; m < M; ++m)
                bm = fmax(bm, err[r * M + m]);
              B[r] = fmax(B[r], bm / sqrt(fmax2));
            }
        }
    }

  free(bu);
  free(err);
  free(uu);
  free(f);
  free(a);
  free(f2);
  free(T);

  return 0;
}
