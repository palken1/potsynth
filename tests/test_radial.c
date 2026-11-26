/*
 * test_radial.c
 */

/* evaluate radial basis function i at u */
static double
test_radial_eval(const size_t i, const double u, const potsynth_radial_workspace * w)
{
  const double t = (w->ud > 0.0) ? (u - w->uc) / w->ud : 0.0;
  double T0 = 1.0, T1 = t, sum = w->basis[i * w->ncheb];
  size_t k;

  if (w->ncheb > 1)
    sum += w->basis[i * w->ncheb + 1] * T1;

  for (k = 2; k < w->ncheb; ++k)
    {
      const double T2 = 2.0 * t * T1 - T0;

      sum += w->basis[i * w->ncheb + k] * T2;
      T0 = T1;
      T1 = T2;
    }

  return sum;
}

/*
test_radial_error()
  Weighted approximation error of u^(l+rpow) over a fine grid, relative to
the largest weighted function value
*/

static double
test_radial_error(const double u0, const double u1, const size_t lmin, const size_t lmax,
                  const size_t nfunc, const int * rpow, const double * wn,
                  const potsynth_radial_workspace * w)
{
  const size_t L1 = lmax + 1;
  const size_t ng = 500;
  double emax = 0.0, fmaxv = 0.0;
  size_t c, l, g, i;

  for (c = 0; c < nfunc; ++c)
    {
      for (l = lmin; l <= lmax; ++l)
        {
          for (g = 0; g <= ng; ++g)
            {
              const double u = u0 + (u1 - u0) * g / (double) ng;
              const double f = pow(u, (double) l + rpow[c]);
              double approx = 0.0;

              for (i = 0; i < w->rho; ++i)
                approx += w->beta[(c * w->rho + i) * L1 + l] * test_radial_eval(i, u, w);

              emax = fmax(emax, wn[c * L1 + l] * fabs(approx - f));
              fmaxv = fmax(fmaxv, wn[c * L1 + l] * fabs(f));
            }
        }
    }

  return emax / fmaxv;
}

static void test_radial_quantities(void);
static void test_radial_surface(void);
static void test_radial_bound(void);

static void
test_radial(void)
{
  const double a = 6371.2;
  const double rmin[] = { a + 250.0, a, a, a + 450.0 };
  const double rmax[] = { a + 800.0, 2.0 * a, a + 1000.0, a + 450.0 };
  const size_t lmin[] = { 1, 1, 0, 1 };
  const size_t lmax[] = { 185, 185, 13, 60 };
  const double tols[] = { 1.0e-6, 1.0e-10, 1.0e-14 };
  const int rpow[2] = { 2, 1 };
  size_t ic, it, l;

  for (ic = 0; ic < sizeof(lmax) / sizeof(lmax[0]); ++ic)
    {
      const size_t L = lmax[ic];
      const size_t L1 = L + 1;
      const double u0 = a / rmax[ic];
      const double u1 = a / rmin[ic];
      double * wn = malloc(2 * L1 * sizeof(double));
      size_t lastrho = 0;

      for (l = 0; l <= L; ++l)
        {
          wn[l] = (l + 1.0) * pow(0.8, (double) l); /* functional 0 */
          wn[L1 + l] = pow(0.8, (double) l);        /* functional 1 */
        }

      for (it = 0; it < sizeof(tols) / sizeof(tols[0]); ++it)
        {
          potsynth_radial_workspace * w = potsynth_radial_alloc(u0, u1, lmin[ic], L, 2, rpow, wn, tols[it]);
          double err;

          ps_test(w == NULL, "radial: case %zu tol=%g allocation", ic, tols[it]);
          if (w == NULL)
            continue;

          err = test_radial_error(u0, u1, lmin[ic], L, 2, rpow, wn, w);

          ps_test(err > 20.0 * tols[it] + 1.0e-15,
                  "radial: case %zu tol=%g rho=%zu ncheb=%zu (error %.2e)",
                  ic, tols[it], w->rho, w->ncheb, err);
          ps_test(w->rho < lastrho, "radial: case %zu rank increases as tol decreases", ic);

          lastrho = w->rho;
          potsynth_radial_free(w);
        }

      free(wn);
    }

  test_radial_quantities();
  test_radial_surface();
  test_radial_bound();
}

/*
test_radial_quantities()
  Each quantity within tol on its own. Two quantities share a radial power, as
the gravity gradient parts T_rr and the spin-2 part (T_tt - T_pp, 2 T_tp) do:
T_rr is dominated by degree 0, the spin-2 part (which has no degrees 0 and 1)
by degree 2, and both have a Kaula tail up to L = 2190 on a GOCE shell. The
singular value rule alone gives rho = 4 at tol = 1e-7, with the error of the
spin-2 part at 2.9 tol (with EGM2008 itself, the maximum error of that part
reached 2.6 eps at eps = 1e-6).
*/

static void
test_radial_quantities(void)
{
  const double a = 6378.1363, R = 6371.2;
  const double u0 = a / (R + 281.0), u1 = a / (R + 220.0);
  const size_t L = 2190, L1 = L + 1, ng = 400;
  const double tols[] = { 1.0e-7, 1.0e-11 };
  const int rpow[2] = { 3, 3 };
  double * wn = malloc(2 * L1 * sizeof(double));
  size_t it, l, c, g, i;

  /* degree rms of the two parts at r = a: sigma_l = 1 (l = 0), 4.8e-4 (C20), 1e-5/l^2 */
  for (l = 0; l <= L; ++l)
    {
      const double dl = (double) l;
      const double sig = (l == 0) ? 1.0 : (l == 1) ? 0.0 : (l == 2) ? 4.8e-4 : 1.0e-5 / (dl * dl);
      const double s2l = sqrt(2.0 * dl + 1.0) * sig;

      wn[l] = (dl + 1.0) * (dl + 2.0) * s2l;                                            /* T_rr */
      wn[L1 + l] = (l < 2) ? 0.0 : sqrt((dl - 1.0) * dl * (dl + 1.0) * (dl + 2.0)) * s2l; /* spin 2 */
    }

  for (c = 0; c < 2; ++c)
    {
      double norm = 0.0;

      for (l = 0; l <= L; ++l)
        norm += wn[c * L1 + l] * wn[c * L1 + l];
      for (l = 0; l <= L; ++l)
        wn[c * L1 + l] /= sqrt(norm);
    }

  for (it = 0; it < sizeof(tols) / sizeof(tols[0]); ++it)
    {
      potsynth_radial_workspace * w = potsynth_radial_alloc(u0, u1, 0, L, 2, rpow, wn, tols[it]);

      ps_test(w == NULL, "radial quantities: tol=%g allocation", tols[it]);
      if (w == NULL)
        continue;

      /* per quantity: max over the shell of the rms error over the sphere, relative to
       * the largest rms in the shell */
      for (c = 0; c < 2; ++c)
        {
          double emax2 = 0.0, fmax2 = 0.0;

          for (g = 0; g <= ng; ++g)
            {
              const double u = u0 + (u1 - u0) * g / (double) ng;
              double b[64], e2 = 0.0, f2 = 0.0;

              for (i = 0; i < w->rho && i < 64; ++i)
                b[i] = test_radial_eval(i, u, w);

              for (l = 0; l <= L; ++l)
                {
                  const double f = pow(u, (double) l + rpow[c]);
                  double approx = 0.0;

                  for (i = 0; i < w->rho && i < 64; ++i)
                    approx += w->beta[(c * w->rho + i) * L1 + l] * b[i];

                  e2 += wn[c * L1 + l] * wn[c * L1 + l] * (approx - f) * (approx - f);
                  f2 += wn[c * L1 + l] * wn[c * L1 + l] * f * f;
                }

              emax2 = fmax(emax2, e2);
              fmax2 = fmax(fmax2, f2);
            }

          ps_test(sqrt(emax2 / fmax2) > 1.5 * tols[it],
                  "radial quantities: tol=%g quantity %zu rho=%zu (error %.2e)",
                  tols[it], c, w->rho, sqrt(emax2 / fmax2));
        }

      potsynth_radial_free(w);
    }

  free(wn);
}

/*
test_radial_surface()
  Pointwise radial error on a thin shell reaching inside r = a, at high degree:
g_r and g_h of a Kaula model to L = 2190 (with C20), from 21 km inside r = a to
10 km above it, as for Earth's topographic surface. Over this shell u^(l+2)
varies by a factor of about 1800 at l = 2190. The largest error of a degree-l
part anywhere on the sphere is sqrt(2l+1) times its rms (addition theorem), so
the radial error of a quantity at u is at most
B(u) = sum_l sqrt(2l+1) wn_l |u^(l+p) - sum_i beta_il b_i(u)|, here evaluated at
2001 radii (between the SVD nodes too) and required to be within tol of the
quantity's largest rms in the shell. With the rms rank check and the Chebyshev
degree chosen from the summed tails (potsynth 9558b56), B reached about 10 tol
for g_h at tol = 5e-7 (EGM2008 itself: 5.1 eps for g_h at eps = 1e-6 on the
ETOPO 2022 surface).
*/

static void
test_radial_surface(void)
{
  const double a = 6378.1363;
  const double u0 = a / (a + 10.0), u1 = a / (a - 21.0);
  const size_t L = 2190, L1 = L + 1, ng = 2000;
  const double tol = 5.0e-7;
  const int rpow[2] = { 2, 2 };
  double * wn = malloc(2 * L1 * sizeof(double));
  potsynth_radial_workspace * w;
  size_t l, c, g, i;

  for (l = 0; l <= L; ++l)
    {
      const double dl = (double) l;
      const double sig = (l == 0) ? 1.0 : (l == 1) ? 0.0 : (l == 2) ? 4.8e-4 : 1.0e-5 / (dl * dl);
      const double s2l = sqrt(2.0 * dl + 1.0) * sig;

      wn[l] = (dl + 1.0) * s2l;                          /* g_r */
      wn[L1 + l] = sqrt(dl * (dl + 1.0)) * s2l;          /* g_h */
    }

  for (c = 0; c < 2; ++c)
    {
      double norm = 0.0;

      for (l = 0; l <= L; ++l)
        norm += wn[c * L1 + l] * wn[c * L1 + l];
      for (l = 0; l <= L; ++l)
        wn[c * L1 + l] /= sqrt(norm);
    }

  w = potsynth_radial_alloc(u0, u1, 0, L, 2, rpow, wn, tol);
  ps_test(w == NULL, "radial surface: allocation");
  if (w != NULL)
    {
      for (c = 0; c < 2; ++c)
        {
          double bmax = 0.0, fmax2 = 0.0;

          for (g = 0; g <= ng; ++g)
            {
              const double u = u0 + (u1 - u0) * g / (double) ng;
              double b[64], bsum = 0.0, f2 = 0.0;

              for (i = 0; i < w->rho && i < 64; ++i)
                b[i] = test_radial_eval(i, u, w);

              for (l = 0; l <= L; ++l)
                {
                  const double f = pow(u, (double) l + rpow[c]);
                  double approx = 0.0;

                  for (i = 0; i < w->rho && i < 64; ++i)
                    approx += w->beta[(c * w->rho + i) * L1 + l] * b[i];

                  bsum += sqrt(2.0 * l + 1.0) * wn[c * L1 + l] * fabs(approx - f);
                  f2 += wn[c * L1 + l] * wn[c * L1 + l] * f * f;
                }

              bmax = fmax(bmax, bsum);
              fmax2 = fmax(fmax2, f2);
            }

          ps_test(bmax / sqrt(fmax2) > tol,
                  "radial surface: quantity %zu rho=%zu ncheb=%zu (pointwise bound %.2e, tol %.0e)",
                  c, w->rho, w->ncheb, bmax / sqrt(fmax2), tol);
        }

      potsynth_radial_free(w);
    }

  free(wn);
}

/*
test_radial_bound_ld()
  Independent check of err_radial: the pointwise bound of the radial error,
max over the components and over 801 equally spaced radii of
B(u) = sum_l sqrt(2l+1) w_l |u^(l+p) - sum_i beta_il b_i(u)| relative to
max_u (sum_l w_l^2 u^(2(l+p)))^(1/2), in long double (potsynth: double, on
Chebyshev points)
*/

static double
test_radial_bound_ld(const double u0, const double u1, const size_t lmin, const size_t lmax,
                     const size_t nfunc, const int * rpow, const double * wn,
                     const potsynth_radial_workspace * w)
{
  const size_t L1 = lmax + 1, ng = 800;
  double worst = 0.0;
  size_t c, c2, g, i, l;

  for (c = 0; c < nfunc; ++c)
    {
      double bmax = 0.0, fmax2 = 0.0;
      int rep = 0;

      for (c2 = 0; c2 < c && !rep; ++c2)
        rep = (rpow[c2] == rpow[c] &&
               memcmp(wn + c2 * L1 + lmin, wn + c * L1 + lmin, (L1 - lmin) * sizeof(double)) == 0);
      if (rep)
        continue;

      for (g = 0; g <= ng; ++g)
        {
          const double u = u0 + (u1 - u0) * g / (double) ng;
          long double b[128], f, bsum = 0.0L, f2 = 0.0L;

          for (i = 0; i < w->rho && i < 128; ++i)
            b[i] = test_radial_eval(i, u, w);

          f = powl((long double) u, (long double) ((long) lmin + rpow[c]));
          for (l = lmin; l <= lmax; ++l, f *= u)
            {
              long double approx = 0.0L;

              if (wn[c * L1 + l] == 0.0)
                continue;
              for (i = 0; i < w->rho && i < 128; ++i)
                approx += (long double) w->beta[(c * w->rho + i) * L1 + l] * b[i];

              bsum += sqrtl(2.0L * l + 1.0L) * fabsl((long double) wn[c * L1 + l] * (approx - f));
              f2 += (long double) wn[c * L1 + l] * wn[c * L1 + l] * f * f;
            }

          bmax = fmax(bmax, (double) bsum);
          fmax2 = fmax(fmax2, (double) f2);
        }

      if (fmax2 > 0.0)
        worst = fmax(worst, bmax / sqrt(fmax2));
    }

  return worst;
}

/*
test_radial_bound()
  The documented bound of the radial step: err_radial <= rank_tol wherever a
basis reaches it, for the functionals of the paper's GOCE benchmark, g and T of
a Kaula model to L = 2190 (random coefficients of rms 1e-5/l^2, weights from
potsynth_radial_weights() as in the set-up), on the GOCE shell 220-281 km.
potsynth 6255bf0 returned err_radial at 10 to 46 times rank_tol here at
eps = 1e-6 to 1e-11 (its Chebyshev degree was too small for the rank, and its
refinement stopped after one pass that did not halve the bound); with EGM2008
the server rerun of 2026-11-26 found 1.9 eps at eps = 1e-10. Also: the rank does not fall
as eps is tightened; err_radial agrees with an independent long double bound
(between the Chebyshev points too); a thick shell (0-2000 km, g); and two
shells whose bound has a rounding floor above rank_tol (inside r = a, and
eps = 1e-14), where the rank must stay near the knee of the bound rather than
grow towards the number of nodes.
*/

static void
test_radial_bound(void)
{
  const size_t L = 2190, L1 = L + 1;
  const double GM = 3.986004415e14, a = 6378136.3;
  const double epss[] = { 1e-6, 1e-8, 1e-9, 1e-10, 1e-11 };
  double * coef = malloc(L1 * L1 * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  uint64_t state = 2026;
  size_t l, it, lastrho = 0;
  long m;
  int k;

  if (coef == NULL)
    {
      ps_test(1, "radial bound: allocation");
      return;
    }

  /* Kaula model: normal random coefficients (Box-Muller), rms 1e-5/l^2, C00 = 1 */
  coef[0] = 1.0;
  for (l = 1; l <= L; ++l)
    {
      const double sig = (l == 1) ? 0.0 : 1.0e-5 / ((double) l * (double) l);

      for (m = -(long) l; m <= (long) l; ++m)
        {
          double u1, u2;

          state = state * 6364136223846793005ULL + 1442695040888963407ULL;
          u1 = ((double) (state >> 11) + 0.5) / 9007199254740992.0;
          state = state * 6364136223846793005ULL + 1442695040888963407ULL;
          u2 = ((double) (state >> 11) + 0.5) / 9007199254740992.0;
          coef[l * (l + 1) + m] = sig * sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
        }
    }

  params.method = POTSYNTH_METHOD_DIRECT;
  params.nthreads = 1;

  for (k = 0; k < 2; ++k)
    {
      /* k = 0: g and T; k = 1: g */
      const int q = (k == 0) ? (POTSYNTH_GRAV_G | POTSYNTH_GRAV_T) : POTSYNTH_GRAV_G;
      potsynth_workspace * w = potsynth_grav_alloc(0, L, L, coef, GM, a, a, a, q, &params);
      double * wn;

      ps_test(w == NULL, "radial bound: model %d", k);
      if (w == NULL)
        continue;

      wn = malloc(w->nfunc * L1 * sizeof(double));
      potsynth_radial_weights(w, 1, wn);

      if (k == 0)
        {
          const double u0 = a / (a + 281.0e3), u1 = a / (a + 220.0e3);

          /* GOCE shell, eps = 1e-6 .. 1e-11 */
          for (it = 0; it < sizeof(epss) / sizeof(epss[0]); ++it)
            {
              const double tol = 0.5 * epss[it];
              potsynth_radial_workspace * rw = potsynth_radial_alloc(u0, u1, 0, L, w->nfunc, w->rpow, wn, tol);

              ps_test(rw == NULL, "radial bound: GOCE g,T eps=%g allocation", epss[it]);
              if (rw == NULL)
                continue;

              ps_test(!(rw->err > 0.0 && rw->err <= tol),
                      "radial bound: GOCE g,T eps=%g rho=%zu err_radial %.2e <= rank_tol %.0e",
                      epss[it], rw->rho, rw->err, tol);
              ps_test(rw->rho < lastrho, "radial bound: GOCE g,T rank %zu does not fall as eps tightens",
                      rw->rho);
              lastrho = rw->rho;

              if (epss[it] == 1e-10)
                {
                  const double ind = test_radial_bound_ld(u0, u1, 0, L, w->nfunc, w->rpow, wn, rw);
                  ps_test(!(ind <= 1.1 * rw->err && ind >= 0.5 * rw->err),
                          "radial bound: GOCE g,T independent bound %.2e, err_radial %.2e", ind, rw->err);
                }

              potsynth_radial_free(rw);
            }

          /* rounding floor above rank_tol: inside r = a (the topographic surface, floor
           * about 2e-10), and eps = 1e-14 on the GOCE shell (floor about 3e-13); the
           * rank stays near the knee of the bound instead of growing with the number of
           * nodes (at eps = 1e-14 it is 31, against 13 at the knee of eps = 1e-13: the
           * first pass starts at a larger Chebyshev degree) */
          {
            const double s0[2] = { a / (a + 8.9e3), u0 }, s1[2] = { a / (a - 21.8e3), u1 };
            const size_t rmax[2] = { 20, 40 };
            int j;

            for (j = 0; j < 2; ++j)
              {
                const double tol = (j == 0) ? 5e-11 : 5e-15;
                potsynth_radial_workspace * rw =
                  potsynth_radial_alloc(s0[j], s1[j], 0, L, w->nfunc, w->rpow, wn, tol);

                ps_test(rw == NULL, "radial bound: floor case %d allocation", j);
                if (rw == NULL)
                  continue;
                ps_test(!(rw->rho <= rmax[j] && rw->err > 0.0 && rw->err < 1e-8),
                        "radial bound: floor case %d rho=%zu (at most %zu) err_radial %.2e", j, rw->rho,
                        rmax[j], rw->err);
                potsynth_radial_free(rw);
              }
          }
        }
      else
        {
          /* thick shell 0-2000 km, g (the rounding floor of this model's bound there is
           * 5.3e-11, so eps = 1e-9 is the tightest well above it) */
          const double u0 = a / (a + 2000.0e3), u1 = 1.0;
          const double e2[2] = { 1e-6, 1e-9 };

          for (it = 0; it < 2; ++it)
            {
              const double tol = 0.5 * e2[it];
              potsynth_radial_workspace * rw = potsynth_radial_alloc(u0, u1, 0, L, w->nfunc, w->rpow, wn, tol);

              ps_test(rw == NULL, "radial bound: thick shell eps=%g allocation", e2[it]);
              if (rw == NULL)
                continue;
              ps_test(!(rw->err > 0.0 && rw->err <= tol),
                      "radial bound: thick shell g eps=%g rho=%zu err_radial %.2e <= rank_tol %.0e", e2[it],
                      rw->rho, rw->err, tol);
              if (it == 1)
                {
                  const double ind = test_radial_bound_ld(u0, u1, 0, L, w->nfunc, w->rpow, wn, rw);
                  ps_test(!(ind <= 1.1 * rw->err && ind >= 0.5 * rw->err),
                          "radial bound: thick shell independent bound %.2e, err_radial %.2e", ind, rw->err);
                }
              potsynth_radial_free(rw);
            }
        }

      free(wn);
      potsynth_free(w);
    }

  free(coef);
}
