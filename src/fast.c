/*
 * fast.c
 *
 * Construction of the fast-transform grid (steps 1-4 of the fast method)
 *
 * For each order m (in parallel):
 *
 * 1) Column recursion in l at the coarse colatitude nodes of one quadrant,
 *    theta_j in (0, pi/2), fused with the sums over l of every functional, split by
 *    the parity of l - m (potsynth_alf_accum). The other three quadrants of the
 *    doubled sphere follow from parity and sign rules, so the recursion runs on a
 *    quarter of the nodes.
 *
 * 2) Exact Fourier coefficients in theta by one FFT of length Nt per functional.
 *
 * 3) Pre-sharpening (division by the window transform), zero padding, and an
 *    inverse FFT to the fine theta grid.
 *
 * 4) The values of the stored rows are written into the grid in half-complex
 *    order along phi.
 *
 * Functionals are expanded into terms as in the direct method (column Q, E or
 * G with a weight row and a power of sin(theta)); the sums of the terms of a
 * functional are combined at each node before step 2. The second-derivative
 * operators S2A (like P) and S2B (like dP/dtheta) keep the parity and sign
 * rules of a single operator, so the quadrant symmetries apply unchanged.
 *
 * Then, for each stored row (in parallel), an inverse real FFT along phi turns
 * the half-complex spectra into grid values, and the periodic padding columns
 * are filled. Peak memory is the grid itself plus small per-thread buffers.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

#include "potsynth_internal.h"
#include "potsynth_fft.h"

#define FAST_TBLK 8 /* pairs transformed together in the theta pass (fast_order()) */

/* per-thread buffers */
typedef struct
{
  double * W;           /* weights, 2*nterm rows of length lmax+1 */
  size_t * n0;          /* first index (l - m) of each term pair, length nterm */
  unsigned char * src;  /* column of each term pair (POTSYNTH_COL_xxx), length nterm */
  double * Q;           /* column values, (lmax+1)*POTSYNTH_NB */
  double * E;           /* column values, (lmax+1)*POTSYNTH_NB */
  double * G;           /* column values, (lmax+1)*POTSYNTH_NB, or NULL */
  double * acc;         /* sums of one node block, 4*nterm*POTSYNTH_NB */
  double * S;           /* quadrant sums, 4*nterm*nqp */
  double * H;           /* doubled-sphere functions of a group of pairs, FAST_TBLK*Nt complex */
  double * Z;           /* fine theta grid of each pair, npair*nf1 complex */
  double * pm;          /* s_j^(m-d) of the order being computed, d = 0, 1, 2: mantissa, 3*nqp */
  double * p2;          /* and its power of two (with the column scaling), or 0 if not a double */
  int * pe;             /* and its exponent (with the column scaling), 3*nqp */
} fast_thread;

/* data shared by all threads during set-up */
typedef struct
{
  size_t Nt;            /* coarse theta nodes on the doubled sphere */
  size_t nq;            /* nodes per quadrant, Nt/4 */
  size_t nqp;           /* nq padded to a multiple of POTSYNTH_NB */
  double * xq;          /* cos(theta_j) of quadrant nodes, length nqp */
  double * sq;          /* sin(theta_j) of quadrant nodes, length nqp */
  double * c1;          /* pre-sharpening factors in theta, length lmax+1 */
  double * c2;          /* pre-sharpening factors in phi, length lmax+1 */
  double * shift;       /* e^{-i p pi/Nt}/Nt, p = -lmax..lmax (complex) */
  size_t nterm;         /* term pairs: (functional c, term t, radial function i) */
  size_t * tbase;       /* first term pair of functional c is tbase[c] + i; term t at + t*rho */
  int * tcol;           /* column of each term pair (POTSYNTH_COL_xxx) */
  int * tpost;          /* factor after the sum (POTSYNTH_POST_xxx) */
  size_t * tsub;        /* term number within the functional */
  size_t * tfunc;       /* functional of each term pair */
  size_t * tpart;       /* part (operator and coefficient set) of each term pair */
  int needG;            /* 1 if any term uses column G */
  int needE;            /* 1 if any term uses column E (none for potentials and radial derivatives) */
  int nthreads;         /* number of threads */
  fast_thread * threads;
} fast_state_t;

static int fast_state_alloc(potsynth_workspace * w, fast_state_t * state);
static void fast_state_free(fast_state_t * state);
static int fast_order(const size_t m, const fast_state_t * state, fast_thread * thread,
                      potsynth_workspace * w);
static int fast_row(const size_t row, const fast_state_t * state, potsynth_workspace * w);
static int func_same_sets(const potsynth_workspace * w, const size_t c, const size_t c2);
static size_t check_select(const potsynth_workspace * w, const size_t * grp, const size_t ngrp, size_t * sel);
static int check_direct(const potsynth_workspace * w, const size_t nsel, const size_t * sel, const size_t n,
                        const double * r, const double * th, const double * ph, double * out);
static long floor_mod(const long a, const long n);

/*
potsynth_fast_init()
  Construct the radial basis and the grid of the fast method

Inputs: w - workspace, with coefficients and functionals set

Return: success or error
*/

int
potsynth_fast_init(potsynth_workspace * w)
{
  const size_t L = w->lmax;
  const size_t L1 = L + 1;
  const size_t nband = 2 * L + 1;
  fast_state_t state;
  double * wn;
  double tol;
  long mm, rr;
  int status = POTSYNTH_SUCCESS;

  if (w->C == NULL || w->S == NULL)
    {
      /* freed after a fast set-up (setup.c): such a workspace cannot be set up again */
      POTSYNTH_ERROR ("the coefficients of the model were freed after its fast set-up", POTSYNTH_EINVAL);
    }

  /* window */
  w->sigma = (w->params.sigma > 0.0) ? w->params.sigma : 2.0;
  w->w = (w->params.w > 0) ? w->params.w : potsynth_window_choose(w->params.eps, w->sigma);
  if (w->w < 2 || w->w > 32)
    {
      POTSYNTH_ERROR ("window width must be between 2 and 32", POTSYNTH_EINVAL);
    }

  w->beta = M_PI * (1.0 - 1.0 / (2.0 * w->sigma)) * w->w;
  w->half = 0.5 * w->w;

  /* step 1: radial basis */
  wn = malloc(w->nfunc * L1 * sizeof(double));
  if (wn == NULL)
    {
      POTSYNTH_ERROR ("unable to allocate radial weights", POTSYNTH_ENOMEM);
    }

  status = potsynth_radial_weights(w, 1, wn);
  if (status)
    {
      free(wn);
      return status; /* reported by potsynth_radial_weights() */
    }

  tol = potsynth_radial_tol(&(w->params));

  /* window width: widen it where a component has more of its power near the band
     edge than the table assumes, or the shell reaches inside r = a
     (potsynth_window_choose_shell) */
  if (w->params.w <= 0)
    {
      const size_t nmin = (size_t) ceil(w->sigma * nband);
      const size_t nf = potsynth_next_smooth((nmin > 2 * (size_t) w->w) ? nmin : 2 * (size_t) w->w, 2);

      /* radial functions of the caller: no powers of a/r (u = 1) */
      w->w = potsynth_window_choose_shell(w->params.eps, w->sigma, L, nf, w->lmin, w->nfunc, w->rpow,
                                          wn, (w->nrad > 0) ? 1.0 : w->a / w->rmax,
                                          (w->nrad > 0) ? 1.0 : w->a / w->rmin);
      w->beta = M_PI * (1.0 - 1.0 / (2.0 * w->sigma)) * w->w;
      w->half = 0.5 * w->w;
    }

  /* the radial basis does not depend on the window: it is kept when the set-up is
     repeated with a wider one */
  if (w->radial_workspace_p == NULL)
    {
      if (w->nrad > 0)
        w->radial_workspace_p = potsynth_radial_alloc_caller(w->nrad, L, w->nfunc);
      else
        w->radial_workspace_p = potsynth_radial_alloc(w->a / w->rmax, w->a / w->rmin, w->lmin, L,
                                                      w->nfunc, w->rpow, wn, tol);
    }
  free(wn);
  if (w->radial_workspace_p == NULL)
    return POTSYNTH_EINVAL; /* reported by potsynth_radial_alloc() */

  /* grid geometry: nf >= sigma (2 lmax + 1) and >= 2 w, even and 2,3,5-smooth */
  {
    const size_t nmin = (size_t) ceil(w->sigma * nband);
    w->nf1 = potsynth_next_smooth((nmin > 2 * (size_t) w->w) ? nmin : 2 * (size_t) w->w, 2);
  }
  w->nf2 = w->nf1;
  w->h1 = 2.0 * M_PI / w->nf1;
  w->h2 = 2.0 * M_PI / w->nf2;
  w->inv_h1 = w->nf1 / (2.0 * M_PI);
  w->inv_h2 = w->nf2 / (2.0 * M_PI);

  /* rows that windows of points with theta in [th_lo,th_hi] can touch, plus one
     ([0,pi] unless a colatitude range is given) */
  w->row_lo = (long) ceil(w->th_lo * w->inv_h1 - w->half) - 1;
  w->row_hi = (long) ceil(w->th_hi * w->inv_h1 - w->half) + w->w;

  /* columns for phi in [0,2pi), with periodic padding so windows never wrap */
  w->col_lo = (long) ceil(-w->half) - 1;
  w->col_hi = (long) w->nf2 + w->w;

  w->nrows = (size_t) (w->row_hi - w->row_lo + 1);
  w->ncols = (size_t) (w->col_hi - w->col_lo + 1);
  w->nch = w->nfunc * w->radial_workspace_p->rho;
  w->nchp = (w->nch + POTSYNTH_CHAN_ALIGN - 1) / POTSYNTH_CHAN_ALIGN * POTSYNTH_CHAN_ALIGN;

  if ((double) w->nrows * w->ncols * w->nchp * sizeof(double) > (double) SIZE_MAX / 2)
    {
      POTSYNTH_ERROR ("grid too large", POTSYNTH_ENOMEM);
    }

  /* not zeroed here (a serial pass over the whole grid): the orders write the Fourier
   * slots 0..2 mmax of every row and functional, and fast_row() supplies the zeros of the
   * slots above them and of the padding functionals */
  w->grid = potsynth_aligned_alloc(w->nrows * w->ncols * w->nchp * sizeof(double));
  if (w->grid == NULL)
    {
      POTSYNTH_ERROR ("unable to allocate grid", POTSYNTH_ENOMEM);
    }

  status = fast_state_alloc(w, &state);
  if (status)
    {
      fast_state_free(&state);
      return status;
    }

  /* steps 1-4 for each order m (work decreases with m: dynamic schedule) */
#pragma omp parallel for private(mm) schedule(dynamic, 1) num_threads(state.nthreads)
  for (mm = 0; mm <= (long) w->mmax; ++mm)
    {
#ifdef _OPENMP
      fast_thread * thread = &(state.threads[omp_get_thread_num()]);
#else
      fast_thread * thread = &(state.threads[0]);
#endif
      int s = fast_order((size_t) mm, &state, thread, w);

      if (s)
        {
#pragma omp critical (potsynth_fast_status)
          status = s;
        }
    }

  /* inverse real FFT along phi of each stored row, and periodic padding */
  if (status == POTSYNTH_SUCCESS)
    {
#pragma omp parallel for private(rr) schedule(static) num_threads(state.nthreads)
      for (rr = 0; rr < (long) w->nrows; ++rr)
        {
          int s = fast_row((size_t) rr, &state, w);

          if (s)
            {
#pragma omp critical (potsynth_fast_status)
              status = s;
            }
        }
    }

  fast_state_free(&state);

  if (status)
    {
      POTSYNTH_ERROR ("FFT failed (out of memory)", status);
    }

  return POTSYNTH_SUCCESS;
}

/*
fast_state_alloc()
  Allocate and initialize the shared set-up data and per-thread buffers

Inputs: w     - workspace
        state - (output) state

Notes:
1) The coarse theta nodes on the doubled sphere are theta_j = (j + 1/2) 2 pi / Nt,
with Nt a multiple of 4 and >= 2 lmax + 2 (exact for degree lmax); the quadrant
nodes are j < Nt/4
*/

static int
fast_state_alloc(potsynth_workspace * w, fast_state_t * state)
{
  const size_t L = w->lmax;
  const size_t L1 = L + 1;
  const size_t nband = 2 * L + 1;
  size_t j;
  int i, status;

  memset(state, 0, sizeof(fast_state_t));

  /* term pairs, functional-major, then term, then radial function */
  state->tbase = malloc(w->nfunc * sizeof(size_t));
  if (state->tbase == NULL)
    {
      POTSYNTH_ERROR ("unable to allocate set-up data", POTSYNTH_ENOMEM);
    }

  /* the terms of a functional are those of its parts, one after the other */
  state->nterm = 0;
  for (j = 0; j < w->nfunc; ++j)
    {
      size_t p;

      state->tbase[j] = state->nterm;
      for (p = w->pbase[j]; p < w->pbase[j + 1]; ++p)
        state->nterm += potsynth_op_nterm(w->ptheta[p]) * w->radial_workspace_p->rho;
    }

  state->tcol = malloc(state->nterm * sizeof(int));
  state->tpost = malloc(state->nterm * sizeof(int));
  state->tsub = malloc(state->nterm * sizeof(size_t));
  state->tfunc = malloc(state->nterm * sizeof(size_t));
  state->tpart = malloc(state->nterm * sizeof(size_t));
  if (!state->tcol || !state->tpost || !state->tsub || !state->tfunc || !state->tpart)
    {
      POTSYNTH_ERROR ("unable to allocate set-up data", POTSYNTH_ENOMEM);
    }

  state->needG = 0;
  state->needE = 0;
  for (j = 0; j < w->nfunc; ++j)
    {
      const size_t rho = w->radial_workspace_p->rho;
      size_t p, t, t0 = 0, i2;

      for (p = w->pbase[j]; p < w->pbase[j + 1]; t0 += potsynth_op_nterm(w->ptheta[p]), ++p)
      for (t = 0; t < potsynth_op_nterm(w->ptheta[p]); ++t)
        {
          for (i2 = 0; i2 < rho; ++i2)
            {
              const size_t kt = state->tbase[j] + (t0 + t) * rho + i2;

              potsynth_op_term(w->ptheta[p], t, &(state->tcol[kt]), &(state->tpost[kt]));
              state->tsub[kt] = t;
              state->tfunc[kt] = j;
              state->tpart[kt] = p;
              if (state->tcol[kt] == POTSYNTH_COL_G)
                state->needG = 1;
              if (state->tcol[kt] == POTSYNTH_COL_E)
                state->needE = 1;
            }
        }
    }

  state->Nt = potsynth_next_smooth(2 * L + 2, 4);
  state->nq = state->Nt / 4;
  state->nqp = (state->nq + POTSYNTH_NB - 1) / POTSYNTH_NB * POTSYNTH_NB;
  state->nthreads = w->nthreads;

  state->xq = malloc(state->nqp * sizeof(double));
  state->sq = malloc(state->nqp * sizeof(double));
  state->c1 = malloc(L1 * sizeof(double));
  state->c2 = malloc(L1 * sizeof(double));
  state->shift = malloc(2 * nband * sizeof(double));
  state->threads = calloc((size_t) state->nthreads, sizeof(fast_thread));
  if (!state->xq || !state->sq || !state->c1 || !state->c2 || !state->shift ||
      !state->threads)
    {
      POTSYNTH_ERROR ("unable to allocate set-up data", POTSYNTH_ENOMEM);
    }

  for (j = 0; j < state->nqp; ++j)
    {
      const size_t jj = (j < state->nq) ? j : state->nq - 1; /* padding repeats the last node */
      const double theta = (jj + 0.5) * 2.0 * M_PI / state->Nt;

      state->xq[j] = cos(theta);
      state->sq[j] = sin(theta);
    }

  for (j = 0; j < nband; ++j)
    {
      const double p = (double) j - (double) L;

      state->shift[2 * j] = cos(p * M_PI / state->Nt) / state->Nt;
      state->shift[2 * j + 1] = -sin(p * M_PI / state->Nt) / state->Nt;
    }

  /* pre-sharpening factors h / phi_hat(k) */
  status = potsynth_window_ft(L1, w->w, w->beta, w->h1, state->c1);
  if (status)
    return status;

  status = potsynth_window_ft(L1, w->w, w->beta, w->h2, state->c2);
  if (status)
    return status;

  for (j = 0; j < L1; ++j)
    {
      state->c1[j] = w->h1 / state->c1[j];
      state->c2[j] = w->h2 / state->c2[j];
    }

  for (i = 0; i < state->nthreads; ++i)
    {
      fast_thread * thread = &(state->threads[i]);

      thread->W = potsynth_aligned_alloc(2 * state->nterm * L1 * sizeof(double));
      thread->n0 = malloc(state->nterm * sizeof(size_t));
      thread->src = malloc(state->nterm);
      thread->Q = potsynth_aligned_alloc(L1 * POTSYNTH_NB * sizeof(double));
      thread->E = state->needE ? potsynth_aligned_alloc(L1 * POTSYNTH_NB * sizeof(double)) : NULL;
      thread->G = state->needG ? potsynth_aligned_alloc(L1 * POTSYNTH_NB * sizeof(double)) : NULL;
      thread->acc = potsynth_aligned_alloc(4 * state->nterm * POTSYNTH_NB * sizeof(double));
      thread->S = potsynth_aligned_alloc(4 * state->nterm * state->nqp * sizeof(double));
      thread->H = potsynth_aligned_alloc(2 * FAST_TBLK * state->Nt * sizeof(double));
      thread->Z = potsynth_aligned_alloc(2 * w->nch * w->nf1 * sizeof(double));
      thread->pm = malloc(3 * state->nqp * sizeof(double));
      thread->p2 = malloc(3 * state->nqp * sizeof(double));
      thread->pe = malloc(3 * state->nqp * sizeof(int));
      if (!thread->W || !thread->n0 || !thread->src || !thread->Q || (state->needE && !thread->E) ||
          (state->needG && !thread->G) || !thread->acc || !thread->S || !thread->H ||
          !thread->Z || !thread->pm || !thread->p2 || !thread->pe)
        {
          POTSYNTH_ERROR ("unable to allocate thread buffers", POTSYNTH_ENOMEM);
        }
    }

  return POTSYNTH_SUCCESS;
}

static void
fast_state_free(fast_state_t * state)
{
  if (state->xq)
    free(state->xq);

  if (state->sq)
    free(state->sq);

  if (state->c1)
    free(state->c1);

  if (state->c2)
    free(state->c2);

  if (state->shift)
    free(state->shift);

  if (state->tbase)
    free(state->tbase);

  if (state->tcol)
    free(state->tcol);

  if (state->tpost)
    free(state->tpost);

  if (state->tsub)
    free(state->tsub);

  if (state->tfunc)
    free(state->tfunc);

  if (state->tpart)
    free(state->tpart);

  if (state->threads)
    {
      int i;

      for (i = 0; i < state->nthreads; ++i)
        {
          fast_thread * thread = &(state->threads[i]);

          if (thread->W)
            potsynth_aligned_free(thread->W);

          free(thread->pm);
          free(thread->p2);
          free(thread->pe);

          if (thread->n0)
            free(thread->n0);

          if (thread->src)
            free(thread->src);

          if (thread->Q)
            potsynth_aligned_free(thread->Q);

          if (thread->E)
            potsynth_aligned_free(thread->E);

          if (thread->G)
            potsynth_aligned_free(thread->G);

          if (thread->acc)
            potsynth_aligned_free(thread->acc);

          if (thread->S)
            potsynth_aligned_free(thread->S);

          if (thread->H)
            potsynth_aligned_free(thread->H);

          if (thread->Z)
            potsynth_aligned_free(thread->Z);

        }

      free(state->threads);
    }
}

/*
fast_order()
  Steps 1-4 for one order m: fill the half-complex slots of all stored rows

Inputs: m      - SH order
        state  - shared set-up data
        thread - per-thread buffers
        w      - workspace

Notes:
1) Pair k = (c, i) of functional c and radial function i has weight rows
A (cos m phi) and B (sin m phi); the angular function is
G = Re[ (S_A - i S_B) e^{i m phi} ]. Each pair is the sum of its terms (one
for P, P/sin, dP/dtheta; three for S2A, two for S2B), each accumulated
separately and multiplied by its power of sin(theta) (and x) at the node

2) Values on the doubled sphere from the quadrant sums E (l-m even) and
O (l-m odd) of the combined terms, with sigma = (-1)^m for P, S2A, S2B and
(-1)^(m-1) for P/sin, dP/dtheta:

theta_j       : E + O
pi - theta_j  : E - O  (P, P/sin, S2A),  O - E  (dP/dtheta, S2B)
2pi - theta_j : sigma * value at theta_j
pi + theta_j  : sigma * value at pi - theta_j

Under theta -> pi - theta, Q has the parity (-1)^(l-m) in x, E and G the
opposite parity, and a factor x changes sign; every term of S2A then
follows E - O and every term of S2B O - E. Under theta -> 2pi - theta only
the power of sin(theta) changes sign, and the powers of all terms of a
functional have the parity of m (S2A, S2B).
*/

static int
fast_order(const size_t m, const fast_state_t * state, fast_thread * thread, potsynth_workspace * w)
{
  const size_t L = w->lmax;
  const size_t L1 = L + 1;
  const size_t nn = L - m + 1;
  const size_t rho = w->radial_workspace_p->rho;
  const size_t npair = w->nch;
  const size_t Nt = state->Nt;
  const size_t nqp = state->nqp;
  const size_t nf1 = w->nf1;
  const int scale = w->alf_workspace_p->scale[m];
  const double dm = (double) m;
  size_t k, k0, b, j, l;
  long r;

  /* weight rows of each term pair */
  for (k = 0; k < state->nterm; ++k)
    {
      const size_t c = state->tfunc[k];
      const size_t i = (k - state->tbase[c]) % rho;
      const size_t p = state->tpart[k];
      const double * beta = w->radial_workspace_p->beta + (c * rho + i) * L1;
      const double * C = w->C + potsynth_part_set(w, p, i) * w->nalf;
      const double * S = w->S + potsynth_part_set(w, p, i) * w->nalf;
      double * wa = thread->W + 2 * k * L1;
      double * wb = wa + L1;
      /* factors per degree in wa, then the coefficients */
      const size_t lstart = potsynth_term_factors(w, p, state->tsub[k], m, beta, wa);

      thread->src[k] = (unsigned char) state->tcol[k];
      thread->n0[k] = lstart - m;

      for (l = m; l <= L; ++l)
        {
          const size_t idx = potsynth_alf_idx(l, m, L);
          const double f = wa[l - m];
          double A, B;

          if (w->pdphi[p] == 0)
            {
              A = C[idx];
              B = S[idx];
            }
          else
            {
              A = dm * S[idx];
              B = -dm * C[idx];
            }

          wa[l - m] = f * A;
          wb[l - m] = f * B;
        }
    }

  /* recursion + accumulation over node blocks of the first quadrant */
  for (b = 0; b < nqp; b += POTSYNTH_NB)
    {
      if (state->needE)
        potsynth_alf_column(m, state->xq + b, state->sq + b, thread->Q, thread->E, thread->G,
                            w->alf_workspace_p);
      else
        potsynth_alf_column_q(m, state->xq + b, thread->Q, thread->G, w->alf_workspace_p);
      potsynth_alf_accum(nn, state->nterm, thread->W, L1, thread->src, thread->n0, thread->Q,
                         thread->E, thread->G, thread->acc);

      for (k = 0; k < 4 * state->nterm; ++k)
        {
          for (j = 0; j < POTSYNTH_NB; ++j)
            thread->S[k * nqp + b + j] = thread->acc[k * POTSYNTH_NB + j];
        }
    }

  /* s_j^(m-d), d = 0, 1, 2, at each node, as mantissa and exponent (with the column
   * scaling of order m); they depend only on the node and the order, not on the term.
   * The power of two is applied by multiplying with 2^e, which gives the value of
   * ldexp() bit for bit (one rounding of the exact product) when 2^e is a double;
   * otherwise ldexp() is used */
  if (m > 0)
    {
      for (j = 0; j < state->nq; ++j)
        {
          size_t d;

          for (d = 0; d <= 2 && d <= m; ++d)
            {
              int expo;

              potsynth_pow_scaled(state->sq[j], m - d, &thread->pm[3 * j + d], &expo);
              expo += scale;
              thread->pe[3 * j + d] = expo;
              thread->p2[3 * j + d] = (expo >= -1074 && expo <= 1023) ? ldexp(1.0, expo) : 0.0;
            }
        }
    }

  /* the pairs in groups of FAST_TBLK: combine terms and build the doubled-sphere
   * functions, forward FFTs of the group in one call (vectorized across its pairs),
   * pre-sharpen, backward FFTs; a group's functions stay in cache */
  for (k0 = 0; k0 < npair; k0 += FAST_TBLK)
    {
      const size_t k1 = (k0 + FAST_TBLK < npair) ? k0 + FAST_TBLK : npair;

      for (k = k0; k < k1; ++k)
        {
          const size_t c = k / rho;
          const size_t i = k % rho;
          /* the parts of a functional have the same sign rule under theta -> 2pi - theta
             (potsynth_alloc_sets()); the rule under theta -> pi - theta is that of each part */
          const int op0 = w->ptheta[w->pbase[c]];
          const int even = (op0 == POTSYNTH_THETA_P || op0 == POTSYNTH_THETA_S2A || op0 == POTSYNTH_THETA_S2B);
          const double sgn = even ? ((m & 1) ? -1.0 : 1.0) : ((m & 1) ? 1.0 : -1.0);
          double * H = thread->H + 2 * (k - k0) * Nt;

          for (j = 0; j < state->nq; ++j)
            {
              double vA = 0.0, vB = 0.0, uA = 0.0, uB = 0.0;
              size_t part, t0 = 0;

              for (part = w->pbase[c]; part < w->pbase[c + 1]; ++part)
              {
              const int op = w->ptheta[part];
              const size_t nt = potsynth_op_nterm(op);
              const int typeE = (op == POTSYNTH_THETA_DP || op == POTSYNTH_THETA_S2B);
              double ae = 0.0, ao = 0.0, be = 0.0, bo = 0.0;
              size_t t;

              for (t = 0; t < nt; ++t)
                {
                  const size_t kt = state->tbase[c] + (t0 + t) * rho + i;
                  const int post = state->tpost[kt];
                  const double * S = thread->S + 4 * kt * nqp;
                  const int xf = (post == POTSYNTH_POST_XSM || post == POTSYNTH_POST_XSM2);
                  double f = xf ? state->xq[j] : 1.0;

                  /* s^(m-2) terms have zero weights for m < 2 */
                  if (m < 2 && (post == POTSYNTH_POST_SM2 || post == POTSYNTH_POST_XSM2))
                    continue;

                  if (m > 0)
                    {
                      const size_t d = (post == POTSYNTH_POST_SM1) ? 1 :
                                       (post == POTSYNTH_POST_SM2 || post == POTSYNTH_POST_XSM2) ? 2 : 0;
                      const double p2 = thread->p2[3 * j + d];

                      f *= thread->pm[3 * j + d];
                      if (p2 != 0.0)
                        {
                          ae += (f * S[j]) * p2;
                          ao += (f * S[nqp + j]) * p2;
                          be += (f * S[2 * nqp + j]) * p2;
                          bo += (f * S[3 * nqp + j]) * p2;
                        }
                      else
                        {
                          const int expo = thread->pe[3 * j + d];

                          ae += ldexp(f * S[j], expo);
                          ao += ldexp(f * S[nqp + j], expo);
                          be += ldexp(f * S[2 * nqp + j], expo);
                          bo += ldexp(f * S[3 * nqp + j], expo);
                        }
                    }
                  else
                    {
                      ae += f * S[j];
                      ao += f * S[nqp + j];
                      be += f * S[2 * nqp + j];
                      bo += f * S[3 * nqp + j];
                    }
                }

              /* theta_j and pi - theta_j */
              if (part == w->pbase[c])
                {
                  vA = ae + ao;
                  vB = be + bo;
                  uA = typeE ? (ao - ae) : (ae - ao);
                  uB = typeE ? (bo - be) : (be - bo);
                }
              else
                {
                  vA += ae + ao;
                  vB += be + bo;
                  uA += typeE ? (ao - ae) : (ae - ao);
                  uB += typeE ? (bo - be) : (be - bo);
                }

              t0 += nt;
              }

              /* H = vA - i vB */
              H[2 * j] = vA;
              H[2 * j + 1] = -vB;
              H[2 * (Nt / 2 - 1 - j)] = uA;
              H[2 * (Nt / 2 - 1 - j) + 1] = -uB;

              /* 2pi - theta_j and pi + theta_j: sin(theta) changes sign */
              H[2 * (Nt - 1 - j)] = sgn * vA;
              H[2 * (Nt - 1 - j) + 1] = -sgn * vB;
              H[2 * (Nt / 2 + j)] = sgn * uA;
              H[2 * (Nt / 2 + j) + 1] = -sgn * uB;
            }

        }

      if (potsynth_fft_c2c(Nt, k1 - k0, Nt, 1, thread->H, 1.0) != 0)
        return POTSYNTH_ENOMEM;

      /* pre-sharpen each pair of the group onto the fine theta grid */
      {
        const double pfac = (m == 0) ? 1.0 : 0.5;

        for (k = k0; k < k1; ++k)
          {
            const double * H = thread->H + 2 * (k - k0) * Nt;
            double * Z = thread->Z + 2 * k * nf1;
            long p;

            memset(Z, 0, 2 * nf1 * sizeof(double));
            for (p = -(long) L; p <= (long) L; ++p)
              {
                const size_t ip = (size_t) floor_mod(p, (long) Nt);
                const size_t iz = (size_t) floor_mod(p, (long) nf1);
                const double * sh = state->shift + 2 * (size_t) (p + (long) L);
                const double f = state->c1[labs(p)] * state->c2[m] * pfac;
                const double hr = H[2 * ip], hi = H[2 * ip + 1];

                Z[2 * iz] = f * (sh[0] * hr - sh[1] * hi);
                Z[2 * iz + 1] = f * (sh[0] * hi + sh[1] * hr);
              }
          }
      }

      if (potsynth_fft_c2c(nf1, k1 - k0, nf1, 0, thread->Z + 2 * k0 * nf1, 1.0) != 0)
        return POTSYNTH_ENOMEM;
    }

  /* half-complex slots: m = 0 -> column 0 (real part), m >= 1 -> columns 2m-1 (real), 2m (imag) */
  for (r = w->row_lo; r <= w->row_hi; ++r)
    {
      const size_t lrow = (size_t) floor_mod(r, (long) nf1);
      double * row = w->grid + (size_t) (r - w->row_lo) * w->ncols * w->nchp;

      if (m == 0)
        {
          double * g = row + (size_t) (0 - w->col_lo) * w->nchp;

          for (k = 0; k < npair; ++k)
            g[potsynth_layer(w, k / rho, k % rho, rho)] = thread->Z[2 * k * nf1 + 2 * lrow];
        }
      else
        {
          double * gr = row + (size_t) ((long) (2 * m - 1) - w->col_lo) * w->nchp;
          double * gi = gr + w->nchp;

          for (k = 0; k < npair; ++k)
            {
              const size_t q = potsynth_layer(w, k / rho, k % rho, rho);

              gr[q] = thread->Z[2 * k * nf1 + 2 * lrow];
              gi[q] = thread->Z[2 * k * nf1 + 2 * lrow + 1];
            }
        }
    }

  return POTSYNTH_SUCCESS;
}

/*
fast_row()
  Inverse real FFT along phi of each functional of one stored row, and
periodic padding columns

Notes:
1) The grid is not zeroed at allocation: the Fourier slots above 2 mmax, which no
order writes, are set to zero here, and so are the padding functionals nch..nchp-1
2) The row is transformed in place, all its functionals in one call: element k of
functional q is at column k, q (stride nchp between columns), and the transforms
are vectorized across the functionals, which are next to each other

Inputs: row   - stored row index
        state - shared set-up data
        w     - workspace
*/

static int
fast_row(const size_t row, const fast_state_t * state, potsynth_workspace * w)
{
  const size_t nf2 = w->nf2;
  const size_t nchp = w->nchp;
  const size_t nch = w->nch;
  const size_t kmax = 2 * w->mmax; /* last Fourier slot written by the orders */
  double * g = w->grid + row * w->ncols * nchp;
  double * g0 = g + (size_t) (0 - w->col_lo) * nchp; /* column 0 */
  size_t k;
  long col;

  (void) state;

  for (k = kmax + 1; k < nf2; ++k)
    memset(g0 + k * nchp, 0, nch * sizeof(double));

  if (potsynth_fft_r_backward(nf2, nch, (ptrdiff_t) nchp, 1, g0, 1.0) != 0)
    return POTSYNTH_ENOMEM;

  /* periodic padding columns */
  for (col = w->col_lo; col <= w->col_hi; ++col)
    {
      if (col < 0 || col >= (long) nf2)
        memcpy(g + (size_t) (col - w->col_lo) * nchp, g0 + (size_t) floor_mod(col, (long) nf2) * nchp,
               nch * sizeof(double));
    }

  /* padding functionals */
  if (nchp > nch)
    {
      for (k = 0; k < w->ncols; ++k)
        memset(g + k * nchp + nch, 0, (nchp - nch) * sizeof(double));
    }

  return POTSYNTH_SUCCESS;
}

/*
sets_msq()
  Mean square over the sphere of the degree-l, order-m part of a component at
r = a, for a model with several coefficient sets

Inputs: w   - workspace
        g   - component
        grp - component of each functional
        l   - degree
        m   - order
        idx - index of (l, m) in a coefficient set

Notes:
1) In each functional, the parts with the same theta and phi operator are summed
with their coefficients before squaring, so that sets which cancel (e.g. the
difference of two models) give the power of their sum; parts with different
operators, and different functionals, add in power. With radial functions of the
caller, each function adds in power.
*/

static double
sets_msq(const potsynth_workspace * w, const size_t g, const size_t * grp, const size_t l, const size_t m,
         const size_t idx)
{
  const size_t L1 = w->lmax + 1;
  const size_t nsetr = potsynth_nsetr(w);
  double msq = 0.0;
  size_t c, p, p2, j;

  for (c = 0; c < w->nfunc; ++c)
    {
      double wt;

      if (grp[c] != g)
        continue;

      wt = potsynth_func_weight(w, c);

      for (p = w->pbase[c]; p < w->pbase[c + 1]; ++p)
        {
          double q;
          int first = 1;

          /* the first part of the functional with this operator takes all of them */
          for (p2 = w->pbase[c]; p2 < p && first; ++p2)
            first = !(w->ptheta[p2] == w->ptheta[p] && w->pdphi[p2] == w->pdphi[p]);

          if (!first)
            continue;

          q = potsynth_op_msq(w->ptheta[p], w->pdphi[p], l, m);
          if (q == 0.0)
            continue;

          for (j = 0; j < nsetr; ++j)
            {
              double ac = 0.0, as = 0.0;

              for (p2 = p; p2 < w->pbase[c + 1]; ++p2)
                {
                  if (w->ptheta[p2] == w->ptheta[p] && w->pdphi[p2] == w->pdphi[p])
                    {
                      const double rf = w->prfac[p2 * L1 + l];
                      const size_t set = w->pset[p2] + j;

                      ac += rf * w->C[set * w->nalf + idx];
                      as += rf * w->S[set * w->nalf + idx];
                    }
                }

              msq += wt * wt * q * (ac * ac + as * as);
            }
        }
    }

  return msq;
}

/*
potsynth_radial_weights()
  Weights of the radial basis (step 1). Functionals with the same radial power and
field family (scalar, surface gradient, spin-2 tensor) form one irreducible component; its
weight for degree l is the rms over the sphere of its degree-l part at r = a

Inputs: w         - workspace, with coefficients and functionals set
        normalize - 1: normalize the weights of each component over l to unit
                    sum of squares, so that the rank truncation is relative to
                    the size of each component (e.g. potential and acceleration
                    differ by a factor ~a); 0: rms values
        wn        - (output) weights wn[c*(lmax+1) + l], length nfunc*(lmax+1);
                    all functionals of a component get the same weights

Notes:
1) The radial basis takes the largest weight over the components sharing a
radial power (potsynth_radial_alloc)
*/

int
potsynth_radial_weights(const potsynth_workspace * w, const int normalize, double * wn)
{
  const size_t L = w->lmax;
  const size_t L1 = L + 1;
  const size_t nfunc = w->nfunc;
  const size_t nsetr = potsynth_nsetr(w); /* coefficient sets of a part */
  size_t * grp = malloc(nfunc * sizeof(size_t));
  size_t * plist = malloc(w->npart * sizeof(size_t));
  size_t * pfunc = malloc(w->npart * sizeof(size_t));
  size_t ngrp, g, c, l, m;

  if (grp == NULL || plist == NULL || pfunc == NULL)
    {
      free(grp);
      free(plist);
      free(pfunc);
      POTSYNTH_ERROR ("unable to allocate components", POTSYNTH_ENOMEM);
    }

  if (w->C == NULL || w->S == NULL)
    {
      /* freed after a fast set-up (setup.c) */
      free(grp);
      free(plist);
      free(pfunc);
      POTSYNTH_ERROR ("the coefficients of the model were freed after its fast set-up", POTSYNTH_EINVAL);
    }

  for (c = 0; c < nfunc; ++c)
    {
      size_t p;

      for (p = w->pbase[c]; p < w->pbase[c + 1]; ++p)
        pfunc[p] = c;
    }

  memset(wn, 0, nfunc * L1 * sizeof(double));

  ngrp = potsynth_components(w, grp);

  for (g = 0; g < ngrp; ++g)
    {
      size_t c0 = nfunc, np = 0, s, i;
      double norm = 0.0;

      /* parts of the functionals of the component, ordered by coefficient set and
         within a set by functional and part (insertion sort: the lists are short) */
      for (c = 0; c < nfunc; ++c)
        {
          size_t p;

          if (grp[c] != g)
            continue;

          if (c0 == nfunc)
            c0 = c;

          for (p = w->pbase[c]; p < w->pbase[c + 1]; ++p)
            {
              for (i = np; i > 0 && w->pset[plist[i - 1]] > w->pset[p]; --i)
                plist[i] = plist[i - 1];
              plist[i] = p;
              ++np;
            }
        }

      for (l = w->lmin; l <= L; ++l)
        {
          const size_t M = (l < w->mmax) ? l : w->mmax;
          double msq = 0.0;

          for (m = 0; m <= M; ++m)
            {
              const size_t idx = potsynth_alf_idx(l, m, L);

              if (w->nset > nsetr)
                {
                  /* several coefficient sets: parts of a functional with the same operator
                     add coherently (sets that nearly cancel, such as the difference of
                     two epochs, give a small weight); different operators and functionals
                     add in power */
                  msq += sets_msq(w, g, grp, l, m, idx);
                  continue;
                }

              /* sum over the sets of (power of the set) x (sum over its parts) */
              for (i = 0; i < np; i = s)
                {
                  const size_t set = w->pset[plist[i]];
                  double cs, q = 0.0;

                  if (nsetr == 1)
                    {
                      const double C = w->C[set * w->nalf + idx], S = w->S[set * w->nalf + idx];
                      cs = C * C + S * S;
                    }
                  else
                    {
                      /* radial functions of the caller: the power of all their sets */
                      size_t j;

                      cs = 0.0;
                      for (j = 0; j < nsetr; ++j)
                        {
                          const double C = w->C[(set + j) * w->nalf + idx], S = w->S[(set + j) * w->nalf + idx];
                          cs += C * C + S * S;
                        }
                    }

                  for (s = i; s < np && w->pset[plist[s]] == set; ++s)
                    {
                      const size_t p = plist[s];
                      const double rf = w->prfac[p * L1 + l];
                      /* weight of the functional in its component; spin 2: the component is
                         (S2A, 2 S2B), with an m-independent norm */
                      const double wt = potsynth_func_weight(w, pfunc[p]);

                      q += wt * wt * rf * rf * potsynth_op_msq(w->ptheta[p], w->pdphi[p], l, m);
                    }

                  msq += cs * q;
                }
            }

          wn[c0 * L1 + l] = sqrt(msq);
          norm += msq;
        }

      if (normalize && norm > 0.0)
        {
          norm = 1.0 / sqrt(norm);
          for (l = w->lmin; l <= L; ++l)
            wn[c0 * L1 + l] *= norm;
        }

      /* the other functionals of the component share its weights */
      for (c = c0 + 1; c < nfunc; ++c)
        {
          if (grp[c] == g)
            memcpy(wn + c * L1, wn + c0 * L1, L1 * sizeof(double));
        }
    }

  free(grp);
  free(plist);
  free(pfunc);

  return POTSYNTH_SUCCESS;
}

/*
potsynth_func_weight()
  Weight of a functional in the magnitude of its component: given by the caller, or
1, and 2 for the off-diagonal part S2B of a spin-2 tensor (which counts twice)
*/

double
potsynth_func_weight(const potsynth_workspace * w, const size_t c)
{
  if (w->cwt[c] > 0.0)
    return w->cwt[c];

  return (potsynth_func_family(w, c) == 2 && w->ptheta[w->pbase[c]] == POTSYNTH_THETA_S2B) ? 2.0 : 1.0;
}

/*
potsynth_func_family()
  Field family of a functional (potsynth_op_family()): that of its parts, which
all have the same operators or are all of family 1 (potsynth_alloc_sets())
*/

int
potsynth_func_family(const potsynth_workspace * w, const size_t c)
{
  const size_t p = w->pbase[c];
  return potsynth_op_family(w->ptheta[p], w->pdphi[p]);
}

/*
potsynth_components()
  Irreducible components of the functionals: functionals with the same
component given by the caller; otherwise functionals with the same radial power,
field family (scalar, surface gradient, spin-2 tensor) and coefficient sets
form one component, and a functional without a family is a component of its own

Inputs: w   - workspace, with functionals set
        grp - (output) component of each functional, 0..ngrp-1, numbered in the order
              of their first functionals; length nfunc

Return: number of components
*/

size_t
potsynth_components(const potsynth_workspace * w, size_t * grp)
{
  size_t ngrp = 0;
  size_t c, c2;

  for (c = 0; c < w->nfunc; ++c)
    {
      const int fam = potsynth_func_family(w, c);

      grp[c] = ngrp;

      if (w->ccomp[c] > 0)
        {
          /* component given by the caller */
          for (c2 = 0; c2 < c; ++c2)
            {
              if (w->ccomp[c2] == w->ccomp[c])
                {
                  grp[c] = grp[c2];
                  break;
                }
            }
        }
      else
        {
          for (c2 = 0; c2 < c && fam >= 0; ++c2)
            {
              if (w->ccomp[c2] == 0 && w->rpow[c2] == w->rpow[c] && potsynth_func_family(w, c2) == fam &&
                  func_same_sets(w, c, c2))
                {
                  grp[c] = grp[c2];
                  break;
                }
            }
        }

      if (grp[c] == ngrp)
        ++ngrp;
    }

  return ngrp;
}

/* 1 if the parts of the functionals c and c2 use the same coefficient sets */
static int
func_same_sets(const potsynth_workspace * w, const size_t c, const size_t c2)
{
  size_t p, q;

  if (w->nset == 1)
    return 1;

  for (p = w->pbase[c]; p < w->pbase[c + 1]; ++p)
    {
      for (q = w->pbase[c2]; q < w->pbase[c2 + 1] && w->pset[q] != w->pset[p]; ++q)
        ;
      if (q == w->pbase[c2 + 1])
        return 0;
    }

  for (q = w->pbase[c2]; q < w->pbase[c2 + 1]; ++q)
    {
      for (p = w->pbase[c]; p < w->pbase[c + 1] && w->pset[p] != w->pset[q]; ++p)
        ;
      if (p == w->pbase[c + 1])
        return 0;
    }

  return 1;
}

/*
potsynth_radial_tol()
  Relative tolerance of the radial step: rank_tol, or eps/2 by default (the
other half of the error budget goes to the window), and at least 1e-16. The
radial basis bounds the maximum of each component's radial error over the shell
by it (see radial_rank_check() in radial.c)
*/

double
potsynth_radial_tol(const potsynth_parameters * params)
{
  const double tol = (params->rank_tol > 0.0) ? params->rank_tol : 0.5 * params->eps;
  return (tol < 1.0e-16) ? 1.0e-16 : tol;
}

/*
potsynth_op_msq()
  Mean square over the sphere of a functional operator applied to one
4-pi normalized harmonic Pbar_l^m(cos theta) cos(m phi) (or sin(m phi))

Inputs: op   - operator in theta (POTSYNTH_THETA_xxx)
        dphi - number of phi derivatives (0 or 1)
        l    - degree
        m    - order

Return: mean square (1 for op = P, dphi = 0)

Notes:
1) Exact, from the orthogonality of the spherical harmonics. Spin 1:
dP and Psin d/dphi (the theta and phi components of the surface gradient)
give l(l+1) - m(2l+1)/2 and m(2l+1)/2, with sum l(l+1). Spin 2: S2A and
S2B d/dphi give (l-1)l(l+1)(l+2) - 4 b and b, b = m(2l+1)(l(l+1) - 1 - m^2)/4,
so that S2A^2 + 4 S2B^2 has mean square (l-1)l(l+1)(l+2). A phi derivative
multiplies by m^2.
2) Psin and S2B vanish for m = 0.
*/

double
potsynth_op_msq(const int op, const int dphi, const size_t l, const size_t m)
{
  const double dl = (double) l, dm = (double) m;
  const double ll = dl * (dl + 1.0);
  const double m2 = (dphi > 0) ? dm * dm : 1.0;
  double q;

  if (op == POTSYNTH_THETA_P)
    q = 1.0;
  else if (op == POTSYNTH_THETA_DP)
    q = ll - 0.5 * dm * (2.0 * dl + 1.0);
  else if (op == POTSYNTH_THETA_PSIN)
    q = (m == 0) ? 0.0 : 0.5 * (2.0 * dl + 1.0) / dm;
  else if (op == POTSYNTH_THETA_S2A)
    q = (dl - 1.0) * ll * (dl + 2.0) - dm * (2.0 * dl + 1.0) * (ll - 1.0 - dm * dm);
  else if (op == POTSYNTH_THETA_S2B)
    q = (m == 0) ? 0.0 : 0.25 * (2.0 * dl + 1.0) * (ll - 1.0 - dm * dm) / dm;
  else
    q = 0.0;

  return (l < m) ? 0.0 : fmax(q, 0.0) * m2;
}

/*
potsynth_op_family()
  Field family of a functional operator: 0 = scalar (P), 1 = surface gradient
(dP, Psin d/dphi), 2 = spin-2 horizontal tensor (S2A, S2B d/dphi), -1 = none
(the functional forms a family of its own)
*/

int
potsynth_op_family(const int op, const int dphi)
{
  if (op == POTSYNTH_THETA_P && dphi == 0)
    return 0;
  else if ((op == POTSYNTH_THETA_DP && dphi == 0) || (op == POTSYNTH_THETA_PSIN && dphi == 1))
    return 1;
  else if ((op == POTSYNTH_THETA_S2A && dphi == 0) || (op == POTSYNTH_THETA_S2B && dphi == 1))
    return 2;

  return -1;
}

/* a mod n in [0, n) */
static long
floor_mod(const long a, const long n)
{
  const long r = a % n;
  return (r < 0) ? r + n : r;
}

/*
potsynth_fast_check()
  Largest error of any irreducible component of the fast method, against the
direct sum, at points chosen where the window error is likely to be largest,
relative to the component's largest magnitude at the points

Inputs: w - workspace after potsynth_fast_init()

Return: worst error ratio; 0 if there is nothing to check; minus the error
        code if the check failed (out of memory, or the radial functions of the
        caller failed at a check point)

Notes:
1) The window error is largest where the field has most of its power near the
band edge, which for a real model is local (mountains, a buried mass). For each
functional, the CHECK_TOP grid nodes with the largest five-point Laplacian of the
functional at the bottom of the shell (where the highest degrees are largest) are
taken, and around each the four cell midpoints (half a grid step from the nodes
in both angles, where the window error of a band-edge mode peaks), at the bottom,
middle and top of the shell. CHECK_NRAND points uniform in the shell are added
for the sizes of the components.

2) The candidate scan is split over the threads, each keeping its own
CHECK_TOP list per functional; the lists are merged afterwards, so the result does
not depend on the number of threads (ties apart).

3) An irreducible component is a set of functionals (potsynth_components(); eq. 10 of the
paper): its error and size at a point are the lengths of the vectors its
functionals form, with the weights of the functionals (potsynth_func_weight()).

4) A model with more than CHECK_MAXCOMP components (many fields on one grid) is
checked for the CHECK_MAXCOMP components with the largest share of their power
in the highest eighth of the degrees, where the window error is largest. A model
with more than CHECK_BUDGET / CHECK_TOP functionals to check takes CHECK_TOP candidate
nodes per irreducible component instead of per functional, ranked by the weighted
root sum of squares of the Laplacians of its functionals; the direct sum of a
component is computed at its own candidate points and the random points only. The direct sum is computed for the functionals of
the checked components only, so that the cost of the check does not grow with
the number of fields.
*/

#define CHECK_TOP      16   /* candidate nodes per functional */
#define CHECK_NRAND    256  /* random points */
#define CHECK_MAXCOMP  64   /* components checked of a model with more components */
#define CHECK_BUDGET   256  /* candidate nodes of all functionals; with more functionals the
                               candidates are taken per component (CHECK_TOP each) */
#define CHECK_CHUNK    256  /* points per evaluation of the fast method */

/* offer a candidate node to a list of ctop entries whose smallest entry is at *jmin;
 * the list is searched only when it changes (the first smallest, as a search before
 * each candidate would find) */
static void
check_offer(const double lap, const long row, const long col, const size_t ctop, double * tops,
            long * topr, long * topc, size_t * jmin)
{
  size_t q;

  if (!(lap > tops[*jmin]))
    return;

  tops[*jmin] = lap;
  topr[*jmin] = row;
  topc[*jmin] = col;

  *jmin = 0;
  for (q = 1; q < ctop; ++q)
    if (tops[q] < tops[*jmin])
      *jmin = q;
}

/*
check_many()
  Errors of the components of a model with many checked functionals: each component
against the direct sum of its own functionals at its own candidate points and the
random points (pown = -1)

Inputs: of - fast values of the checked functionals, of[ck*npt + i]
        worst  - (input/output) largest relative error
        failed - (output) error code if a direct sum failed

Return: 1 on success, 0 on failure
*/

static int
check_many(potsynth_workspace * w, const size_t * grp, const size_t nsel, const size_t * sel,
           const size_t * cg, const size_t ncg, const size_t npt, const double * r, const double * th,
           const double * ph, const long * pown, const double * of, double * worst, int * failed)
{
  size_t * idx = malloc(npt * sizeof(size_t));
  size_t * fs = malloc(nsel * sizeof(size_t));   /* functionals of the component */
  size_t * fk = malloc(nsel * sizeof(size_t));   /* their positions in sel */
  double * rr = malloc(npt * sizeof(double));
  double * tt = malloc(npt * sizeof(double));
  double * pp = malloc(npt * sizeof(double));
  double * od = malloc(nsel * npt * sizeof(double));
  int ok = (idx && fs && fk && rr && tt && pp && od);
  size_t kg;

  (void) grp;
  if (!ok)
    *failed = POTSYNTH_ENOMEM;

  for (kg = 0; ok && kg < ncg; ++kg)
    {
      size_t n = 0, nf = 0, i, j;
      double emax = 0.0, vmax = 0.0;
      int st;

      for (j = 0; j < nsel; ++j)
        {
          if (cg[j] == kg)
            {
              fk[nf] = j;
              fs[nf++] = sel[j];
            }
        }

      for (i = 0; i < npt; ++i)
        {
          if (pown[i] < 0 || (size_t) pown[i] == kg)
            {
              idx[n] = i;
              rr[n] = r[i];
              tt[n] = th[i];
              pp[n] = ph[i];
              ++n;
            }
        }

      st = check_direct(w, nf, fs, n, rr, tt, pp, od);
      if (st != POTSYNTH_SUCCESS)
        {
          *failed = st;
          ok = 0;
          break;
        }

      for (i = 0; i < n; ++i)
        {
          double e2 = 0.0, v2 = 0.0;

          for (j = 0; j < nf; ++j)
            {
              const double wt = potsynth_func_weight(w, fs[j]);
              const double dv = od[j * n + i];
              const double d = of[fk[j] * npt + idx[i]] - dv;

              e2 += wt * wt * d * d;
              v2 += wt * wt * dv * dv;
            }

          emax = fmax(emax, sqrt(e2));
          vmax = fmax(vmax, sqrt(v2));
        }

      if (vmax > 0.0)
        *worst = fmax(*worst, emax / vmax);
    }

  free(idx); free(fs); free(fk); free(rr); free(tt); free(pp); free(od);

  return ok;
}

double
potsynth_fast_check(potsynth_workspace * w)
{
  const size_t nfunc = w->nfunc, rho = w->radial_workspace_p->rho, nchp = w->nchp;
  const size_t nc = w->radial_workspace_p->ncheb;
  /* grid rows inside the colatitude range of the points, without the poles */
  const long rpole = (long) floor(M_PI * w->inv_h1);
  const long rlo = (long) ceil(w->th_lo * w->inv_h1), rhi = (long) floor(w->th_hi * w->inv_h1) + 1;
  const long rfirst = (rlo > 1) ? rlo : 1;
  const long rlast = (rhi < rpole) ? rhi : rpole;
  const int nth = (w->nthreads > 0) ? w->nthreads : 1;
  size_t * grp = malloc(nfunc * sizeof(size_t));
  size_t * sel = malloc(nfunc * sizeof(size_t));
  double * b1 = malloc(rho * sizeof(double));
  size_t * tmin = NULL;
  double * tops = NULL, * r = NULL, * th = NULL, * ph = NULL, * of = NULL, * od = NULL;
  double * buf = NULL, * acc = NULL;
  long * topr = NULL, * topc = NULL;
  double * omap = w->omap;
  size_t ngrp = 0, nsel = 0, ctop = CHECK_TOP, ntop = 0, npt = 0, i, c, k;
  size_t * cg = NULL;   /* many: compact component number of each checked functional */
  double * lsqall = NULL; /* many: per thread, sums of squares of the components' Laplacians */
  size_t ncg = 0, nlist = 0;
  long * pown = NULL;   /* many: list (component) whose candidate a point is, -1: random */
  int many = 0;
  double worst = 0.0;
  int ok = (grp != NULL && sel != NULL && b1 != NULL);
  int failed = ok ? POTSYNTH_SUCCESS : POTSYNTH_ENOMEM;
  long row;
  int th0;

  /* components, and the functionals that are checked */
  if (ok)
    {
      ngrp = potsynth_components(w, grp);
      nsel = check_select(w, grp, ngrp, sel);
      ok = (nsel > 0);
      if (!ok && nfunc > 0)
        failed = POTSYNTH_ENOMEM; /* check_select() returns 0 only when out of memory */
      /* candidates: CHECK_TOP per functional; with many functionals CHECK_TOP per
         component, so that the number of check points stays bounded without fewer
         candidates for each component */
      many = (nsel > CHECK_BUDGET / CHECK_TOP);
      nlist = nsel;
      if (ok && many)
        {
          cg = malloc(nsel * sizeof(size_t));
          ok = (cg != NULL);
          if (!ok)
            failed = POTSYNTH_ENOMEM;
          for (c = 0; ok && c < nsel; ++c)
            {
              size_t c2;

              cg[c] = ncg;
              for (c2 = 0; c2 < c; ++c2)
                {
                  if (grp[sel[c2]] == grp[sel[c]])
                    {
                      cg[c] = cg[c2];
                      break;
                    }
                }
              if (cg[c] == ncg)
                ++ncg;
            }
          nlist = ncg;
        }
      ntop = nlist * ctop;
    }

  if (ok)
    {
      npt = ntop * 12 + CHECK_NRAND;
      tops = calloc(ntop * (size_t) nth, sizeof(double));
      tmin = calloc(nlist * (size_t) nth, sizeof(size_t));
      lsqall = many ? calloc(ncg * (size_t) nth, sizeof(double)) : NULL;
      pown = malloc(npt * sizeof(long));
      topr = calloc(ntop * (size_t) nth, sizeof(long));
      topc = calloc(ntop * (size_t) nth, sizeof(long));
      r = malloc(npt * sizeof(double));
      th = malloc(npt * sizeof(double));
      ph = malloc(npt * sizeof(double));
      of = malloc(nsel * npt * sizeof(double));
      od = malloc(nsel * npt * sizeof(double));
      buf = malloc(nfunc * CHECK_CHUNK * sizeof(double));
      acc = calloc(4 * ngrp, sizeof(double));
      ok = (tops && tmin && (!many || lsqall) && pown && topr && topc && r && th && ph && of && od && buf && acc);
      if (!ok)
        failed = POTSYNTH_ENOMEM;
    }

  if (ok)
    {
      /* b_i at the bottom of the shell (t = +1: all T_k = 1) */
      for (i = 0; i < rho; ++i)
        {
          double s = 0.0;
          for (k = 0; k < nc; ++k)
            s += w->radial_workspace_p->basis[i * nc + k];
          b1[i] = (w->radial_workspace_p->ud > 0.0) ? s : w->radial_workspace_p->basis[i * nc];
        }

      /* candidates: largest |Laplacian| of each checked functional at r = rmin, one list per thread.
       * Each node's radial sum is computed once: a thread keeps the sums of three
       * consecutive rows (columns -1..nf2) and adds one row per step; its rows are
       * consecutive under the static schedule (otherwise the three rows are recomputed) */
#pragma omp parallel num_threads(nth)
      {
#ifdef _OPENMP
        const size_t base = (size_t) omp_get_thread_num() * ntop;
#else
        const size_t base = 0;
#endif
        const size_t nrc = w->nf2 + 2;              /* columns -1 .. nf2 */
        double * ring = malloc(3 * nrc * nsel * sizeof(double));
        double * lsq = many ? lsqall + (base / ntop) * ncg : NULL; /* per component, many */
        long have = LONG_MIN;                       /* centre row of the sums in the ring */

#pragma omp for private(row) schedule(static)
        for (row = rfirst; row < rlast; ++row)
          {
            long col, rr;
            size_t ck;

            if (ring != NULL)
              {
                /* rows row-1, row, row+1 in slots r mod 3 */
                for (rr = (have == row - 1) ? row + 1 : row - 1; rr <= row + 1; ++rr)
                  {
                    double * slotp = ring + (size_t) (((rr % 3) + 3) % 3) * nrc * nsel;

                    for (col = -1; col <= (long) w->nf2; ++col)
                      {
                        const double * g = w->grid + ((size_t) (rr - w->row_lo) * w->ncols +
                                                      (size_t) (col - w->col_lo)) * nchp;

                        for (ck = 0; ck < nsel; ++ck)
                          {
                            const size_t cc = sel[ck];
                            double sum = 0.0;
                            size_t ii;

                            for (ii = 0; ii < rho; ++ii)
                              sum += b1[ii] * g[potsynth_layer(w, cc, ii, rho)];
                            slotp[(size_t) (col + 1) * nsel + ck] = sum;
                          }
                      }
                  }
                have = row;
              }

            for (col = 0; col < (long) w->nf2; ++col)
              {
                for (ck = 0; ck < nsel; ++ck)
                  {
                    const size_t cc = sel[ck];
                    const long dr[5] = { 0, -1, 1, 0, 0 }, dc[5] = { 0, 0, 0, -1, 1 };
                    double v[5];
                    size_t q, ii, slot = base + ck * ctop;
                    size_t * jmin = tmin + base / ctop + ck; /* position of the list's smallest entry */
                    double lap;

                    for (q = 0; q < 5; ++q)
                      {
                        const long rq = row + dr[q], cq = col + dc[q];

                        if (ring != NULL)
                          {
                            v[q] = ring[((size_t) (((rq % 3) + 3) % 3) * nrc + (size_t) (cq + 1)) * nsel + ck];
                          }
                        else
                          {
                            const double * g = w->grid + ((size_t) (rq - w->row_lo) * w->ncols +
                                                          (size_t) (cq - w->col_lo)) * nchp;
                            double sum = 0.0;
                            for (ii = 0; ii < rho; ++ii)
                              sum += b1[ii] * g[potsynth_layer(w, cc, ii, rho)];
                            v[q] = sum;
                          }
                      }

                    lap = fabs(v[1] + v[2] + v[3] + v[4] - 4.0 * v[0]);

                    if (many)
                      {
                        /* weighted sum of squares over the functionals of the component */
                        const double wt = potsynth_func_weight(w, cc);
                        lsq[cg[ck]] += wt * wt * lap * lap;
                        continue;
                      }

                    check_offer(lap, row, col, ctop, tops + slot, topr + slot, topc + slot, jmin);
                  }

                if (many)
                  {
                    size_t kg;

                    for (kg = 0; kg < ncg; ++kg)
                      {
                        const size_t slot = base + kg * ctop;

                        check_offer(sqrt(lsq[kg]), row, col, ctop, tops + slot, topr + slot, topc + slot,
                                    tmin + base / ctop + kg);
                        lsq[kg] = 0.0;
                      }
                  }
              }
          }

        free(ring);
      }

      /* merge the threads' lists into the first */
      for (th0 = 1; th0 < nth; ++th0)
        {
          const size_t base = (size_t) th0 * ntop;
          for (c = 0; c < nlist; ++c)
            {
              size_t q, q2;
              for (q = 0; q < ctop; ++q)
                {
                  const size_t src = base + c * ctop + q;
                  size_t jmin = c * ctop;
                  for (q2 = 1; q2 < ctop; ++q2)
                    if (tops[c * ctop + q2] < tops[jmin])
                      jmin = c * ctop + q2;
                  if (tops[src] > tops[jmin])
                    {
                      tops[jmin] = tops[src];
                      topr[jmin] = topr[src];
                      topc[jmin] = topc[src];
                    }
                }
            }
        }

      /* points */
      {
        const double rad[3] = { w->rmin, 0.5 * (w->rmin + w->rmax), w->rmax };
        unsigned long long s = 88172645463325252ULL;
        size_t n = 0, t, d, q;

        for (t = 0; t < ntop; ++t)
          {
            if (tops[t] <= 0.0)
              continue;
            for (d = 0; d < 4; ++d)
              {
                const double tt = (topr[t] + ((d & 1) ? 0.5 : -0.5)) * w->h1;
                const double pp = (topc[t] + ((d & 2) ? 0.5 : -0.5)) * w->h2;
                if (tt < w->th_lo || tt > w->th_hi)
                  continue;
                for (q = 0; q < 3; ++q)
                  {
                    r[n] = rad[q];
                    th[n] = tt;
                    ph[n] = pp;
                    pown[n] = many ? (long) (t / ctop) : -1;
                    ++n;
                  }
              }
          }

        for (t = 0; t < CHECK_NRAND; ++t)
          {
            double x[3];
            for (q = 0; q < 3; ++q)
              {
                s ^= s << 13; s ^= s >> 7; s ^= s << 17;
                x[q] = (double) (s >> 11) / 9007199254740992.0;
              }
            r[n] = cbrt(pow(w->rmin, 3) + x[0] * (pow(w->rmax, 3) - pow(w->rmin, 3)));
            th[n] = (w->th_lo > 0.0 || w->th_hi < M_PI) ?
                    fmin(fmax(acos(cos(w->th_lo) - x[1] * (cos(w->th_lo) - cos(w->th_hi))), w->th_lo), w->th_hi) :
                    acos(1.0 - 2.0 * x[1]);
            ph[n] = 2.0 * M_PI * x[2];
            pown[n] = -1;
            ++n;
          }
        npt = n;
      }

      /* fast method: functional values (no output map) of the checked functionals */
      {
        potsynth_eval_workspace * ew;
        size_t i0;

        w->omap = NULL;
        ew = potsynth_eval_alloc(0, w);
        ok = (ew != NULL);
        if (!ok)
          failed = POTSYNTH_ENOMEM;

        for (i0 = 0; i0 < npt && ok; i0 += CHECK_CHUNK)
          {
            const size_t m = (npt - i0 < CHECK_CHUNK) ? npt - i0 : CHECK_CHUNK;

            const int st = potsynth_eval(m, r + i0, th + i0, ph + i0, buf, ew);

            ok = (st == POTSYNTH_SUCCESS);
            if (!ok)
              failed = st;

            for (c = 0; c < nsel && ok; ++c)
              memcpy(of + c * npt + i0, buf + sel[c] * m, m * sizeof(double));
          }

        if (ew)
          potsynth_eval_free(ew);
        w->omap = omap;
      }
    }

  /* many functionals: each component against the direct sum of its own functionals at
   * its own candidate points and the random points */
  if (ok && many)
    {
      ok = check_many(w, grp, nsel, sel, cg, ncg, npt, r, th, ph, pown, of, &worst, &failed);
    }

  /* direct sum of the checked functionals */
  if (ok && !many)
    {
      const int st = check_direct(w, nsel, sel, npt, r, th, ph, od);

      ok = (st == POTSYNTH_SUCCESS);
      if (!ok)
        failed = st;
    }

  /* components: largest error and largest value of each over the points */
  if (ok && !many)
    {
      double * e2 = acc, * v2 = acc + ngrp, * emax = acc + 2 * ngrp, * vmax = acc + 3 * ngrp;
      size_t g;

      for (i = 0; i < npt; ++i)
        {
          for (c = 0; c < nsel; ++c)
            e2[grp[sel[c]]] = v2[grp[sel[c]]] = 0.0;

          for (c = 0; c < nsel; ++c)
            {
              const double wt = potsynth_func_weight(w, sel[c]);
              const double d = of[c * npt + i] - od[c * npt + i];

              e2[grp[sel[c]]] += wt * wt * d * d;
              v2[grp[sel[c]]] += wt * wt * od[c * npt + i] * od[c * npt + i];
            }

          for (c = 0; c < nsel; ++c)
            {
              g = grp[sel[c]];
              emax[g] = fmax(emax[g], sqrt(e2[g]));
              vmax[g] = fmax(vmax[g], sqrt(v2[g]));
            }
        }

      for (g = 0; g < ngrp; ++g)
        {
          if (vmax[g] > 0.0)
            worst = fmax(worst, emax[g] / vmax[g]);
        }
    }

  free(grp); free(sel); free(b1); free(tops); free(tmin); free(topr); free(topc); free(cg); free(pown); free(lsqall);
  free(r); free(th); free(ph); free(of); free(od); free(buf); free(acc);

  return failed ? -(double) failed : worst;
}

/*
check_select()
  Functionals checked by potsynth_fast_check(): all functionals, or for a model with
more than CHECK_MAXCOMP components the functionals of the CHECK_MAXCOMP components
with the largest share of their power in the highest eighth of the degrees

Inputs: w    - workspace
        grp  - component of each functional (potsynth_components())
        ngrp - number of components
        sel  - (output) checked functionals, ascending; length nfunc

Return: number of checked functionals (0 if memory could not be allocated)
*/

static size_t
check_select(const potsynth_workspace * w, const size_t * grp, const size_t ngrp, size_t * sel)
{
  const size_t nfunc = w->nfunc, L1 = w->lmax + 1;
  const size_t ltail = w->lmax - w->lmax / 8;
  size_t nsel = 0, c;

  if (ngrp <= CHECK_MAXCOMP)
    {
      for (c = 0; c < nfunc; ++c)
        sel[c] = c;

      return nfunc;
    }
  else
    {
      double * wn = malloc(nfunc * L1 * sizeof(double));
      double * tail = calloc(ngrp, sizeof(double));
      unsigned char * pick = calloc(ngrp, 1);
      size_t g, k, l;

      if (wn != NULL && tail != NULL && pick != NULL && potsynth_radial_weights(w, 1, wn) == POTSYNTH_SUCCESS)
        {
          /* the weights of a component (the same for its functionals) have unit sum of squares */
          for (c = 0; c < nfunc; ++c)
            {
              double s = 0.0;

              for (l = ltail; l < L1; ++l)
                s += wn[c * L1 + l] * wn[c * L1 + l];

              tail[grp[c]] = s;
            }

          for (k = 0; k < CHECK_MAXCOMP; ++k)
            {
              size_t gmax = ngrp;

              for (g = 0; g < ngrp; ++g)
                {
                  if (!pick[g] && (gmax == ngrp || tail[g] > tail[gmax]))
                    gmax = g;
                }

              pick[gmax] = 1;
            }

          for (c = 0; c < nfunc; ++c)
            {
              if (pick[grp[c]])
                sel[nsel++] = c;
            }
        }

      free(wn);
      free(tail);
      free(pick);

      return nsel;
    }
}

/*
check_direct()
  Direct sum of some functionals of a model at points, for potsynth_fast_check():
a model of the direct method with these functionals only, which shares the
coefficients of w

Inputs: w    - workspace
        nsel - number of functionals
        sel  - functionals, length nsel
        n    - number of points
        r    - radii, length n
        th   - colatitudes, length n
        ph   - longitudes, length n
        out  - (output) functional values out[c*n + i], length nsel*n

Return: success or error
*/

static int
check_direct(const potsynth_workspace * w, const size_t nsel, const size_t * sel, const size_t n,
             const double * r, const double * th, const double * ph, double * out)
{
  const size_t L1 = w->lmax + 1;
  potsynth_workspace ws = *w; /* scalars, coefficients and Legendre workspace of w */
  potsynth_eval_workspace * ew = NULL;
  size_t np = 0, c, p, k;
  int status = POTSYNTH_ENOMEM;

  for (c = 0; c < nsel; ++c)
    np += w->pbase[sel[c] + 1] - w->pbase[sel[c]];

  ws.method = POTSYNTH_METHOD_DIRECT;
  ws.nfunc = nsel;
  ws.nout = nsel;
  ws.omap = NULL;
  ws.npart = np;
  ws.grid = NULL;
  ws.radial_workspace_p = NULL;
  ws.doff = NULL;
  ws.tcol = NULL;
  ws.tpost = NULL;
  ws.tfunc = NULL;
  ws.dW = NULL;
  ws.rpow = malloc(nsel * sizeof(int));
  ws.ccomp = malloc(nsel * sizeof(int));
  ws.cwt = malloc(nsel * sizeof(double));
  ws.pbase = malloc((nsel + 1) * sizeof(size_t));
  ws.pset = malloc(np * sizeof(size_t));
  ws.ptheta = malloc(np * sizeof(int));
  ws.pdphi = malloc(np * sizeof(int));
  ws.prfac = potsynth_aligned_calloc(np * L1 * sizeof(double));

  if (ws.rpow && ws.ccomp && ws.cwt && ws.pbase && ws.pset && ws.ptheta && ws.pdphi && ws.prfac)
    {
      for (c = 0, k = 0; c < nsel; ++c)
        {
          ws.rpow[c] = w->rpow[sel[c]];
          ws.ccomp[c] = w->ccomp[sel[c]];
          ws.cwt[c] = w->cwt[sel[c]];
          ws.pbase[c] = k;

          for (p = w->pbase[sel[c]]; p < w->pbase[sel[c] + 1]; ++p, ++k)
            {
              ws.pset[k] = w->pset[p];
              ws.ptheta[k] = w->ptheta[p];
              ws.pdphi[k] = w->pdphi[p];
              memcpy(ws.prfac + k * L1, w->prfac + p * L1, L1 * sizeof(double));
            }
        }
      ws.pbase[nsel] = k;

      status = potsynth_direct_init(&ws);
      if (status == POTSYNTH_SUCCESS)
        {
          ew = potsynth_eval_alloc(0, &ws);
          status = (ew != NULL) ? potsynth_eval(n, r, th, ph, out, ew) : POTSYNTH_ENOMEM;
        }
    }

  if (ew)
    potsynth_eval_free(ew);

  potsynth_direct_free(&ws);
  free(ws.rpow);
  free(ws.ccomp);
  free(ws.cwt);
  free(ws.pbase);
  free(ws.pset);
  free(ws.ptheta);
  free(ws.pdphi);
  if (ws.prfac)
    potsynth_aligned_free(ws.prfac);

  return status;
}
