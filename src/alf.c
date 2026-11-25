/*
 * alf.c
 *
 * 4-pi fully normalized associated Legendre functions
 *
 * The recursions are those of gsl_sf_alf_vsh_array() (GSL specfunc/alf_P.c),
 * which compute P(l,m)/sin(theta) and dP(l,m)/dtheta without dividing by
 * sin(theta), following
 *
 * [1] S. A. Holmes and W. E. Featherstone, A unified approach to the Clenshaw
 *     summation and the recursive computation of very high degree and order
 *     normalised associated Legendre functions, J. Geodesy 76, 279-299, 2002.
 *
 * Differences from GSL:
 *
 * 1) for m >= 1 the columns are computed for P/sin^m and (dP/dtheta)/sin^(m-1),
 *    which are polynomials in x = cos(theta); the powers of sin(theta) are
 *    applied afterwards by the caller, so the recursion itself never underflows
 *
 * 2) instead of a fixed 1e-280 factor, each column is scaled by an exact power
 *    of 2 derived from the analytic maximum of P/sin^m (attained at the poles),
 *    so most columns are not scaled at all and none can overflow
 *
 * 3) sin(theta) may be negative (doubled sphere); only s^2 enters for m >= 1
 *
 * 4) the recursion runs over a block of POTSYNTH_NB nodes at once (vectorized),
 *    and the values are consumed immediately by a register-blocked accumulation
 *    kernel
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "potsynth_internal.h"

/*
potsynth_alf_alloc()
  Allocate ALF workspace and precompute recurrence factors

Inputs: lmax - maximum SH degree

Return: pointer to workspace

Notes:
1) The diagonal Q(m,m) = 2^-scale[m] prod_{k=1}^m d_k, with d_1 = sqrt(3),
d_k = sqrt(1 + 1/(2k)). The scaling is chosen from the maximum of Q(L,m)
at x = 1,

Q(l,m)(1) = sqrt(2(2l+1)) sqrt((l+m)!/(l-m)!) / (2^m m!)

E(l,m) is bounded by about (l+1) times that, hence an extra margin.

2) Returns NULL with POTSYNTH_EUNSUP if lmax is too large for the scaling
*/

potsynth_alf_workspace *
potsynth_alf_alloc(const size_t lmax)
{
  const size_t nlm = (lmax + 1) * (lmax + 2) / 2;
  const double L = (double) lmax;
  potsynth_alf_workspace * w;
  double prod = 1.0;
  size_t l, m;

  w = calloc(1, sizeof(potsynth_alf_workspace));
  if (w == NULL)
    {
      POTSYNTH_ERROR_NULL ("unable to allocate ALF workspace", POTSYNTH_ENOMEM);
    }

  w->lmax = lmax;
  w->ab = potsynth_aligned_calloc(2 * nlm * sizeof(double));
  w->cl = potsynth_aligned_calloc((lmax + 1) * sizeof(double));
  w->qmm = potsynth_aligned_calloc((lmax + 1) * sizeof(double));
  w->scale = calloc(lmax + 1, sizeof(int));
  if (!w->ab || !w->cl || !w->qmm || !w->scale)
    {
      potsynth_alf_free(w);
      POTSYNTH_ERROR_NULL ("unable to allocate ALF arrays", POTSYNTH_ENOMEM);
    }

  /* recurrence factors (GSL_SF_ALF_FOURPI) */
  for (m = 0; m <= lmax; ++m)
    {
      const double dm = (double) m;

      for (l = m + 2; l <= lmax; ++l)
        {
          const size_t idx = potsynth_alf_idx(l, m, lmax);
          const double dl = (double) l;

          w->ab[2 * idx] = (sqrt(2.0 * dl + 1.0) / sqrt(dl + dm)) *
                           (sqrt(2.0 * dl - 1.0) / sqrt(dl - dm));
          w->ab[2 * idx + 1] = -(sqrt(dl + dm - 1.0) / sqrt(dl + dm)) *
                                (sqrt(dl - dm - 1.0) / sqrt(dl - dm)) *
                                (sqrt(2.0 * dl + 1.0) / sqrt(2.0 * dl - 3.0));
        }

      w->cl[m] = sqrt(2.0 * dm + 3.0);
    }

  /* diagonal terms and column scaling */
  w->qmm[0] = 1.0;
  w->scale[0] = 0;
  for (m = 1; m <= lmax; ++m)
    {
      const double dm = (double) m;
      double log2q, need;

      prod *= (m == 1) ? sqrt(3.0) : sqrt(1.0 + 0.5 / dm);

      log2q = 0.5 * log2(2.0 * (2.0 * L + 1.0)) +
              (0.5 * (lgamma(L + dm + 1.0) - lgamma(L - dm + 1.0)) - lgamma(dm + 1.0)) / M_LN2 - dm;
      need = log2q + log2(L + 1.0) + 1.0 - POTSYNTH_ALF_EXPMAX;

      w->scale[m] = (need > 0.0) ? (int) ceil(need) : 0;
      if (w->scale[m] > POTSYNTH_ALF_SCALEMAX)
        {
          potsynth_alf_free(w);
          POTSYNTH_ERROR_NULL ("lmax too large for double precision ALFs", POTSYNTH_EUNSUP);
        }

      w->qmm[m] = ldexp(prod, -w->scale[m]);
    }

  return w;
}

void
potsynth_alf_free(potsynth_alf_workspace * w)
{
  if (w == NULL)
    return; /* as free(NULL) */

  if (w->ab)
    potsynth_aligned_free(w->ab);

  if (w->cl)
    potsynth_aligned_free(w->cl);

  if (w->qmm)
    potsynth_aligned_free(w->qmm);

  if (w->scale)
    free(w->scale);

  free(w);
}

/* potsynth_alf_column() without E: the Q (and G) recursion alone, for models whose
 * terms do not use E (a separate function, so that potsynth_alf_column() compiles as
 * before and its values do not change in the last bits) */
void
potsynth_alf_column_q(const size_t m, const double * restrict x, double * restrict Q, double * restrict G,
             const potsynth_alf_workspace * w)
{
  const size_t lmax = w->lmax;
  const size_t nn = lmax - m + 1;
  const double * ab = w->ab + 2 * potsynth_alf_idx(m, m, lmax);
  size_t k, j;

  for (j = 0; j < POTSYNTH_NB; ++j)
    {
      Q[j] = w->qmm[m];
      if (G)
        G[j] = 0.0;
    }

  if (nn == 1)
    return;

  for (j = 0; j < POTSYNTH_NB; ++j)
    {
      Q[POTSYNTH_NB + j] = w->cl[m] * x[j] * Q[j];
      if (G)
        G[POTSYNTH_NB + j] = -w->cl[m] * Q[j];
    }

  for (k = 2; k < nn; ++k)
    {
      const double a = ab[2 * k], b = ab[2 * k + 1];
      const double * restrict q1 = Q + (k - 1) * POTSYNTH_NB;
      const double * restrict q2 = Q + (k - 2) * POTSYNTH_NB;
      double * restrict q0 = Q + k * POTSYNTH_NB;

      for (j = 0; j < POTSYNTH_NB; ++j)
        q0[j] = a * x[j] * q1[j] + b * q2[j];

      if (G)
        {
          const double * restrict g1 = G + (k - 1) * POTSYNTH_NB;
          const double * restrict g2 = G + (k - 2) * POTSYNTH_NB;
          double * restrict g0 = G + k * POTSYNTH_NB;

          for (j = 0; j < POTSYNTH_NB; ++j)
            g0[j] = a * (x[j] * g1[j] - q1[j]) + b * g2[j];
        }
    }
}

/*
potsynth_alf_column()
  Compute one column (fixed m, l = m..lmax) of the ALF recursion at
POTSYNTH_NB nodes

Inputs: m - SH order
        x - cos(theta) at the nodes, length POTSYNTH_NB
        s - sin(theta) at the nodes, length POTSYNTH_NB (may be negative)
        Q - (output) Q[(l-m)*POTSYNTH_NB + j] = P(l,0) (m = 0) or scaled Q(l,m) (m > 0)
        E - (output) E[(l-m)*POTSYNTH_NB + j] = dP(l,0)/dtheta (m = 0) or scaled E(l,m) (m > 0)
            (potsynth_alf_column_q() computes Q and G alone)
        G - (output) G[(l-m)*POTSYNTH_NB + j] = -dQ(l,m)/dx (scaled as Q), or NULL if
            not needed
        w - workspace

Notes:
1) Recurrences (m >= 1, s2 = s^2):

Q(m,m)   = qmm[m],               E(m,m)   = m x qmm[m]
Q(m+1,m) = c_m x Q(m,m),         E(m+1,m) = c_m (x E(m,m) - s2 Q(m,m))
Q(l,m)   = a x Q(l-1,m) + b Q(l-2,m)
E(l,m)   = a (x E(l-1,m) - s2 Q(l-1,m)) + b E(l-2,m)

For m = 0 the same recursions apply to P and dP/dtheta with s in place of s2.

2) G = -dQ/dx has the same recurrence factors (differentiate the Q recursion):

G(m,m) = 0, G(m+1,m) = c_m (x G(m,m) - Q(m,m)), G(l,m) = a (x G(l-1,m) - Q(l-1,m)) + b G(l-2,m)
*/

void
potsynth_alf_column(const size_t m, const double * restrict x, const double * restrict s,
                    double * restrict Q, double * restrict E, double * restrict G,
                    const potsynth_alf_workspace * w)
{
  const size_t lmax = w->lmax;
  const size_t nn = lmax - m + 1;
  const double * ab = w->ab + 2 * potsynth_alf_idx(m, m, lmax);
  double t[POTSYNTH_NB];   /* s (m = 0) or s^2 (m > 0) */
  size_t k, j;

  for (j = 0; j < POTSYNTH_NB; ++j)
    {
      t[j] = (m == 0) ? s[j] : s[j] * s[j];
      Q[j] = w->qmm[m];
      E[j] = (double) m * x[j] * w->qmm[m];
      if (G)
        G[j] = 0.0;
    }

  if (nn == 1)
    return;

  /* l = m + 1 */
  for (j = 0; j < POTSYNTH_NB; ++j)
    {
      Q[POTSYNTH_NB + j] = w->cl[m] * x[j] * Q[j];
      E[POTSYNTH_NB + j] = w->cl[m] * (x[j] * E[j] - t[j] * Q[j]);
      if (G)
        G[POTSYNTH_NB + j] = -w->cl[m] * Q[j];
    }

  /* l = m + 2 ... lmax */
  for (k = 2; k < nn; ++k)
    {
      const double a = ab[2 * k], b = ab[2 * k + 1];
      const double * restrict q1 = Q + (k - 1) * POTSYNTH_NB;
      const double * restrict q2 = Q + (k - 2) * POTSYNTH_NB;
      const double * restrict e1 = E + (k - 1) * POTSYNTH_NB;
      const double * restrict e2 = E + (k - 2) * POTSYNTH_NB;
      double * restrict q0 = Q + k * POTSYNTH_NB;
      double * restrict e0 = E + k * POTSYNTH_NB;

      for (j = 0; j < POTSYNTH_NB; ++j)
        {
          q0[j] = a * x[j] * q1[j] + b * q2[j];
          e0[j] = a * (x[j] * e1[j] - t[j] * q1[j]) + b * e2[j];
        }

      if (G)
        {
          const double * restrict g1 = G + (k - 1) * POTSYNTH_NB;
          const double * restrict g2 = G + (k - 2) * POTSYNTH_NB;
          double * restrict g0 = G + k * POTSYNTH_NB;

          for (j = 0; j < POTSYNTH_NB; ++j)
            g0[j] = a * (x[j] * g1[j] - q1[j]) + b * g2[j];
        }
    }
}

/* a*b + c with one rounding where the hardware has fused multiply-add (fma() is
 * then one instruction), so that no sum depends on which products the compiler
 * contracts; without it nothing is contracted and a*b + c is used */
#ifdef FP_FAST_FMA
#define ALF_FMA(a, b, c) fma((a), (b), (c))
#else
#define ALF_FMA(a, b, c) ((a) * (b) + (c))
#endif

/*
alf_accum_blocks()
  Kernel of potsynth_alf_accum() for nblk (1 or 2) node blocks; the weights of
each degree are loaded once for all nblk*POTSYNTH_NB nodes. Each node's sums are
formed in the same order whatever nblk is (ALF_FMA)
*/

static inline void
alf_accum_blocks(const size_t nblk, const size_t nn, const size_t npair, const double * restrict W,
                 const size_t ldw, const unsigned char * src, const size_t * n0,
                 const double * restrict Q, const double * restrict E, const double * restrict G,
                 const size_t bs, double * restrict acc)
{
  enum { NB2 = 2 * POTSYNTH_NB };
  const size_t nj = nblk * POTSYNTH_NB;
  size_t k;

  for (k = 0; k < npair; ++k)
    {
      const double * restrict V = (src[k] == POTSYNTH_COL_Q) ? Q : (src[k] == POTSYNTH_COL_E) ? E : G;
      const double * restrict wa = W + 2 * k * ldw;
      const double * restrict wb = wa + ldw;
      double ae[NB2], ao[NB2], be[NB2], bo[NB2];
      size_t l = n0[k], j, h;

      for (j = 0; j < nj; ++j)
        ae[j] = ao[j] = be[j] = bo[j] = 0.0;

      /* make l even so that (l, l+1) = (even, odd) */
      if ((l & 1) && l < nn)
        {
          for (h = 0; h < nblk; ++h)
            {
              const double * restrict v = V + h * bs + l * POTSYNTH_NB;
              double * restrict po = ao + h * POTSYNTH_NB;
              double * restrict qo = bo + h * POTSYNTH_NB;

              for (j = 0; j < POTSYNTH_NB; ++j)
                {
                  po[j] = ALF_FMA(wa[l], v[j], po[j]);
                  qo[j] = ALF_FMA(wb[l], v[j], qo[j]);
                }
            }

          ++l;
        }

      for (; l + 1 < nn; l += 2)
        {
          const double wa0 = wa[l], wa1 = wa[l + 1];
          const double wb0 = wb[l], wb1 = wb[l + 1];

          for (h = 0; h < nblk; ++h)
            {
              const double * restrict v0 = V + h * bs + l * POTSYNTH_NB;
              const double * restrict v1 = v0 + POTSYNTH_NB;
              double * restrict pe = ae + h * POTSYNTH_NB;
              double * restrict po = ao + h * POTSYNTH_NB;
              double * restrict qe = be + h * POTSYNTH_NB;
              double * restrict qo = bo + h * POTSYNTH_NB;

              for (j = 0; j < POTSYNTH_NB; ++j)
                {
                  pe[j] = ALF_FMA(wa0, v0[j], pe[j]);
                  po[j] = ALF_FMA(wa1, v1[j], po[j]);
                  qe[j] = ALF_FMA(wb0, v0[j], qe[j]);
                  qo[j] = ALF_FMA(wb1, v1[j], qo[j]);
                }
            }
        }

      if (l < nn)
        {
          for (h = 0; h < nblk; ++h)
            {
              const double * restrict v = V + h * bs + l * POTSYNTH_NB;
              double * restrict pe = ae + h * POTSYNTH_NB;
              double * restrict qe = be + h * POTSYNTH_NB;

              for (j = 0; j < POTSYNTH_NB; ++j)
                {
                  pe[j] = ALF_FMA(wa[l], v[j], pe[j]);
                  qe[j] = ALF_FMA(wb[l], v[j], qe[j]);
                }
            }
        }

      for (j = 0; j < nj; ++j)
        {
          acc[(4 * k + 0) * nj + j] = ae[j];
          acc[(4 * k + 1) * nj + j] = ao[j];
          acc[(4 * k + 2) * nj + j] = be[j];
          acc[(4 * k + 3) * nj + j] = bo[j];
        }
    }
}

/*
potsynth_alf_accum()
  Accumulate weighted sums over l of column values, split by the parity
of l - m, for npair pairs of weight rows (A and B) sharing a source array,
over one or two blocks of POTSYNTH_NB nodes

Inputs: nblk  - number of node blocks, 1 or 2
        nn    - column length, lmax - m + 1
        npair - number of pairs
        W     - weights, 2*npair rows of length ldw (row 2k = A, 2k+1 = B)
        ldw   - leading dimension of W
        src   - source of pair k: POTSYNTH_COL_Q, POTSYNTH_COL_E or POTSYNTH_COL_G
        n0    - first index (l - m) of pair k
        Q     - column values from potsynth_alf_column(), block h at Q + h*bs
        E     - same for E, or NULL if no pair uses E
        G     - same for G, or NULL if no pair uses G
        bs    - offset between the blocks of Q, E and G
        acc   - (output) sums, length 4*npair*nj with nj = nblk*POTSYNTH_NB, node
                j = h*POTSYNTH_NB + i (node i of block h):
                acc[(4k+0)*nj + j] = sum_{l-m even} WA_k[l-m] src_k[l-m][j]
                acc[(4k+1)*nj + j] = sum_{l-m odd}  WA_k[l-m] src_k[l-m][j]
                acc[(4k+2)*nj + j], acc[(4k+3)*nj + j] = same for WB_k

Notes:
1) The weight rows are the largest operand (2*npair*(lmax+1) doubles, beyond L1
and often L2); with two blocks they are read once per 2*POTSYNTH_NB nodes.
The sums are bit for bit the same for one or two blocks
*/

void
potsynth_alf_accum(const size_t nblk, const size_t nn, const size_t npair, const double * restrict W,
                   const size_t ldw, const unsigned char * src, const size_t * n0,
                   const double * restrict Q, const double * restrict E, const double * restrict G,
                   const size_t bs, double * restrict acc)
{
  if (nblk == 2)
    alf_accum_blocks(2, nn, npair, W, ldw, src, n0, Q, E, G, bs, acc);
  else
    alf_accum_blocks(1, nn, npair, W, ldw, src, n0, Q, E, G, bs, acc);
}

/*
potsynth_op_nterm()
  Number of terms of a functional with a given theta operator
*/

size_t
potsynth_op_nterm(const int op)
{
  if (op == POTSYNTH_THETA_S2A)
    return 3;
  else if (op == POTSYNTH_THETA_S2B)
    return 2;
  else
    return 1;
}

/*
potsynth_op_term()
  Column and post-factor of term sub of a functional with a given theta operator

Inputs: op   - theta operator
        sub  - term number, < potsynth_op_nterm(op)
        col  - (output) POTSYNTH_COL_xxx
        post - (output) POTSYNTH_POST_xxx
*/

void
potsynth_op_term(const int op, const size_t sub, int * col, int * post)
{
  switch (op)
    {
      case POTSYNTH_THETA_P:
        *col = POTSYNTH_COL_Q;
        *post = POTSYNTH_POST_SM;
        break;

      case POTSYNTH_THETA_PSIN:
        *col = POTSYNTH_COL_Q;
        *post = POTSYNTH_POST_SM1;
        break;

      case POTSYNTH_THETA_DP:
        *col = POTSYNTH_COL_E;
        *post = POTSYNTH_POST_SM1;
        break;

      case POTSYNTH_THETA_S2A:
        /* 2m(m-1) s^(m-2) Q + (2m - l(l+1)) s^m Q - 2 x s^m G */
        *col = (sub == 2) ? POTSYNTH_COL_G : POTSYNTH_COL_Q;
        *post = (sub == 0) ? POTSYNTH_POST_SM2 : (sub == 1) ? POTSYNTH_POST_SM : POTSYNTH_POST_XSM;
        break;

      default: /* POTSYNTH_THETA_S2B */
        /* (m-1) x s^(m-2) Q + s^m G */
        *col = (sub == 0) ? POTSYNTH_COL_Q : POTSYNTH_COL_G;
        *post = (sub == 0) ? POTSYNTH_POST_XSM2 : POTSYNTH_POST_SM;
        break;
    }
}

/*
potsynth_op_factor()
  Factor multiplying the weights of term sub at degree l and order m
*/

double
potsynth_op_factor(const int op, const size_t sub, const size_t l, const size_t m)
{
  const double dl = (double) l;
  const double dm = (double) m;

  if (op == POTSYNTH_THETA_S2A)
    {
      if (sub == 0)
        return 2.0 * dm * (dm - 1.0);
      else if (sub == 1)
        return 2.0 * dm - dl * (dl + 1.0);
      else
        return -2.0;
    }
  else if (op == POTSYNTH_THETA_S2B)
    {
      return (sub == 0) ? dm - 1.0 : 1.0;
    }

  return 1.0;
}
