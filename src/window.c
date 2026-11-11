/*
 * window.c
 *
 * Exponential of semicircle window, its Fourier transform, and the choice of
 * window width:
 *
 * phi(z) = exp(beta (sqrt(1 - z^2) - 1)),  |z| <= 1,  z = x / (w h / 2)
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "potsynth_internal.h"

/*
 * Measured maximum error of the complete transform, relative to max |B|,
 * against a 19-digit reference (tests/gen_reference.py), worst case over
 * L = 13, L = 185 (250-800 km shell) and L = 185 (r in [a, 2a]); indexed by w.
 * The error is essentially independent of L and of the radial range; ~5e-15 is
 * the rounding floor of double precision.
 */
static const double err_sigma2[POTSYNTH_WMAX + 1] =
{
  0, 0, 0, 0, 8e-4, 2e-5, 9e-6, 7e-7, 3e-8, 5e-9, 6e-10, 3e-11, 6e-12, 4e-13, 3e-14, 1e-14, 7e-15
};

static const double err_sigma3[POTSYNTH_WMAX + 1] =
{
  0, 0, 0, 0, 2e-4, 3e-5, 1e-6, 2e-7, 7e-9, 7e-10, 6e-11, 3e-12, 4e-13, 1e-14, 1e-14, 7e-15, 9e-15
};

/*
potsynth_window_choose()
  Choose the smallest window width whose measured error is below eps

Inputs: eps   - requested accuracy
        sigma - oversampling factor

Return: window width w
*/

int
potsynth_window_choose(const double eps, const double sigma)
{
  const double * table = (sigma >= 3.0) ? err_sigma3 : err_sigma2;
  int w;

  for (w = 4; w <= POTSYNTH_WMAX; ++w)
    {
      if (table[w] <= eps)
        return w;
    }

  return POTSYNTH_WMAX;
}

/*
potsynth_window_choose_shell()
  Window width for a given model and shell: the smallest width, at least that of
potsynth_window_choose(), at which the estimated window error of every quantity
is within eps of its largest size in the shell

Inputs: eps    - requested accuracy
        sigma  - oversampling factor
        L      - maximum SH degree
        nf     - grid points per period (both directions)
        lmin   - minimum SH degree
        nchan  - number of channels
        rpow   - radial power offsets, length nchan
        wn     - weights wn[c*(L+1) + l] (potsynth_radial_weights), the size
                 of degree l of the quantity of channel c at r = a
        u0, u1 - a / rmax, a / rmin

Return: window width w

Notes:
1) The pre-sharpened window sum returns a Fourier mode k with relative aliasing
error E(k) = sum_{j != 0} |phi_hat(k + j nf)| / phi_hat(k), nearly flat over
most of the band and rising steeply near the band edge |k| = L (by ~30 at
sigma = 2). The table of potsynth_window_choose() was measured for the
magnetic field at satellite altitude, whose power lies well inside the band.
A quantity with most of its power at the highest degrees (the horizontal
gravity gradients), or a shell reaching inside r = a (u > 1, where the highest
degrees grow as u^l), sees more of the band-edge error.
2) A degree-l harmonic has torus modes |k| <= l and |m| <= l, so its relative
error is at most Ebar(l) = max_{k<=l} E(k) + max_{m<=l} E(m). The rms over the
sphere of the error of a quantity at u is then at most
(sum_l wn_l^2 u^(2(l+p)) Ebar(l)^2)^(1/2), against its rms size
(sum_l wn_l^2 u^(2(l+p)))^(1/2); both are taken at the Chebyshev points of the
shell, and the ratio of their maxima R is the quantity's window error relative
to its largest size.
3) R is a bound: Ebar(l) charges every mode of degree l with the band-edge
error in both directions. In tests (synthetic gravity model, L = 360 and 720, full
tensor, shells from r = a and from 21 km inside it; tests/test_window.c) the largest measured error
at a point was 0.20-0.27 R; WINDOW_PEAK = 0.5 keeps a factor 2 on that. The
width is raised until WINDOW_PEAK * R <= eps for every quantity, or up to
POTSYNTH_WMAX. Quantities whose power lies well inside the band, such as the
magnetic field at satellite altitude for which the table was measured, keep
the table's width.
*/

#define WINDOW_PEAK 0.5
#define WINDOW_NALIAS 2       /* aliases j = +-1, +-2: phi_hat decays fast beyond */
#define WINDOW_NU 9           /* Chebyshev points of the shell */

int
potsynth_window_choose_shell(const double eps, const double sigma, const size_t L, const size_t nf,
                             const size_t lmin, const size_t nchan, const int * rpow, const double * wn,
                             const double u0, const double u1)
{
  const size_t L1 = L + 1;
  const size_t nk = L + WINDOW_NALIAS * nf + 1;
  const double h = 2.0 * M_PI / nf;
  double * ph = malloc(nk * sizeof(double));
  double * eb = malloc(L1 * sizeof(double));
  double u[WINDOW_NU];
  int w = potsynth_window_choose(eps, sigma);
  size_t i;

  if (!ph || !eb)
    {
      free(ph);
      free(eb);
      return w; /* fall back on the table */
    }

  for (i = 0; i < WINDOW_NU; ++i)
    u[i] = 0.5 * (u0 + u1) + 0.5 * (u1 - u0) * cos(M_PI * i / (WINDOW_NU - 1));

  for (; w <= POTSYNTH_WMAX; ++w)
    {
      const double beta = M_PI * (1.0 - 1.0 / (2.0 * sigma)) * w;
      double emax = 0.0, worst = 0.0;
      size_t c, k, j;

      if (potsynth_window_ft(nk, w, beta, h, ph))
        break;

      /* Ebar(l): running maximum of E(k), in both directions */
      for (k = 0; k <= L; ++k)
        {
          double e = 0.0;

          for (j = 1; j <= WINDOW_NALIAS; ++j)
            e += fabs(ph[k + j * nf]) + fabs(ph[j * nf - k]);

          e /= ph[k];
          emax = (e > emax) ? e : emax;
          eb[k] = 2.0 * emax;
        }

      for (c = 0; c < nchan; ++c)
        {
          double err = 0.0, size = 0.0;

          for (i = 0; i < WINDOW_NU; ++i)
            {
              double e2 = 0.0, s2 = 0.0;
              size_t l;

              for (l = lmin; l <= L; ++l)
                {
                  const double t = wn[c * L1 + l] * pow(u[i], (double) l + rpow[c]);
                  s2 += t * t;
                  e2 += t * t * eb[l] * eb[l];
                }

              err = fmax(err, sqrt(e2));
              size = fmax(size, sqrt(s2));
            }

          if (size > 0.0)
            worst = fmax(worst, err / size);
        }

      if (getenv("POTSYNTH_DEBUG_WINDOW"))
        fprintf(stderr, "potsynth: window w = %d: largest estimated error %.2e (x%.1f = %.2e), eps %.1e\n",
                w, worst, WINDOW_PEAK, WINDOW_PEAK * worst, eps);

      if (WINDOW_PEAK * worst <= eps)
        break;
    }

  free(ph);
  free(eb);

  return (w <= POTSYNTH_WMAX) ? w : POTSYNTH_WMAX;
}

/*
potsynth_gauss_legendre()
  Gauss-Legendre nodes and weights on [-1,1] (Newton iteration on P_n)

Inputs: n - number of nodes
        x - (output) nodes, length n, increasing
        w - (output) weights, length n
*/

void
potsynth_gauss_legendre(const size_t n, double * x, double * w)
{
  size_t i, k, iter;

  for (i = 0; i < (n + 1) / 2; ++i)
    {
      double z = cos(M_PI * (i + 0.75) / (n + 0.5));
      double p0, p1, dp;

      for (iter = 0; iter < 100; ++iter)
        {
          double dz;

          p0 = 1.0;
          p1 = z;
          for (k = 2; k <= n; ++k)
            {
              const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;

              p0 = p1;
              p1 = p2;
            }

          /* p1 = P_n(z), p0 = P_{n-1}(z) */
          dp = n * (z * p1 - p0) / (z * z - 1.0);
          dz = p1 / dp;
          z -= dz;

          if (fabs(dz) <= 1.0e-16 * fabs(z))
            break;
        }

      p0 = 1.0;
      p1 = z;
      for (k = 2; k <= n; ++k)
        {
          const double p2 = ((2.0 * k - 1.0) * z * p1 - (k - 1.0) * p0) / k;

          p0 = p1;
          p1 = p2;
        }

      dp = n * (z * p1 - p0) / (z * z - 1.0);

      x[i] = -z;
      x[n - 1 - i] = z;
      w[i] = w[n - 1 - i] = 2.0 / ((1.0 - z * z) * dp * dp);
    }
}

/*
potsynth_window_ft()
  Fourier transform of the window

phi_hat(k) = alpha \int_{-1}^{1} phi(z) cos(k alpha z) dz,  alpha = w h / 2

Inputs: nk   - number of frequencies, k = 0..nk-1
        w    - window width (grid steps)
        beta - window shape parameter
        h    - grid spacing
        out  - (output) phi_hat(k), length nk

Return: success or error

Notes:
1) phi has a square-root singularity at z = +-1 (of size e^-beta, not
negligible for small w), which limits Gauss-Legendre in z to algebraic
convergence. With z = sin(t), sqrt(1 - z^2) = cos(t) and the integrand

exp(beta (cos t - 1)) cos(k alpha sin t) cos t,  t in [-pi/2, pi/2]

is entire, so Gauss-Legendre in t converges exponentially
*/

int
potsynth_window_ft(const size_t nk, const int w, const double beta, const double h, double * out)
{
  const size_t nq = 8 * (size_t) w + 100;
  const double alpha = 0.5 * w * h;
  double * zq = malloc(nq * sizeof(double));
  double * wq = malloc(nq * sizeof(double));
  size_t k, q;

  if (!zq || !wq)
    {
      free(zq);
      free(wq);
      POTSYNTH_ERROR ("unable to allocate quadrature arrays", POTSYNTH_ENOMEM);
    }

  potsynth_gauss_legendre(nq, zq, wq);

  for (q = 0; q < nq; ++q)
    {
      const double t = 0.5 * M_PI * zq[q];
      const double ct = cos(t);

      zq[q] = sin(t);
      wq[q] *= 0.5 * M_PI * ct * exp(beta * (ct - 1.0));
    }

  for (k = 0; k < nk; ++k)
    {
      double sum = 0.0;

      for (q = 0; q < nq; ++q)
        sum += wq[q] * cos((double) k * alpha * zq[q]);

      out[k] = alpha * sum;
    }

  free(zq);
  free(wq);

  return POTSYNTH_SUCCESS;
}
