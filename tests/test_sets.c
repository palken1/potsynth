/*
 * test_sets.c
 *
 * Models with several coefficient sets (potsynth_alloc_sets()): outputs that
 * are sums of terms of different sets against sums of single-set models
 */

/* coefficients of nset sets, set-major */
static double *
test_sets_coef(const size_t nset, const size_t lmin, const size_t lmax, const size_t mmax)
{
  const size_t nlm = potsynth_nlm(lmin, lmax, mmax);
  double * coef = malloc(nset * nlm * sizeof(double));
  size_t s;

  for (s = 0; s < nset; ++s)
    {
      double * c = test_random_coef(lmin, lmax, mmax);

      memcpy(coef + s * nlm, c, nlm * sizeof(double));
      free(c);
    }

  return coef;
}

/* outputs of a single-set model with the direct method, out[c*n + i] */
static void
test_sets_single(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef,
                 const double a, const size_t nchan, const potsynth_channel * chan,
                 const size_t n, const double * r, const double * theta, const double * phi, double * out)
{
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;

  params.method = POTSYNTH_METHOD_DIRECT;
  w = potsynth_alloc(lmin, lmax, mmax, coef, a, 0.0, 0.0, nchan, chan, &params);

  if (w)
    {
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);

      potsynth_eval(n, r, theta, phi, out, ew);
      potsynth_eval_free(ew);
      potsynth_free(w);
    }
}

/*
test_sets_vector()
  A scalar and a horizontal vector of three sets (the components of a vector
spherical harmonic expansion, with factors):

out 0 = [set 0] P
out 1 = [set 1] dP - f [set 2] (P/sin) d/dphi
out 2 = [set 1] (P/sin) d/dphi + g rf(l) [set 2] dP

Inputs: desc - description
        lmin - minimum degree
        lmax - maximum degree
        mmax - maximum order
        rmin - minimum radius
        rmax - maximum radius
*/

static void
test_sets_vector(const char * desc, const size_t lmin, const size_t lmax, const size_t mmax,
                 const double rmin, const double rmax)
{
  const double a = 6371.2, f = -1.0, g = 2.5;
  const size_t n = 3000, nlm = potsynth_nlm(lmin, lmax, mmax);
  const double epsv[3] = { 1.0e-6, 1.0e-10, 1.0e-13 };
  double * coef = test_sets_coef(3, lmin, lmax, mmax);
  double * rf = malloc((lmax + 1) * sizeof(double));
  double * grf = malloc((lmax + 1) * sizeof(double));
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * ref = calloc(3 * n, sizeof(double));
  double * out = malloc(3 * n * sizeof(double));
  double * tmp = malloc(2 * n * sizeof(double));
  potsynth_term t0[1] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_term t1[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_DP, 0 }, { 2, f, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_term t2[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 }, { 2, g, rf, POTSYNTH_THETA_DP, 0 } };
  potsynth_output o[3] = { { 2, 1, t0, 0, 0.0 }, { 1, 2, t1, 0, 0.0 }, { 1, 2, t2, 0, 0.0 } };
  potsynth_channel ch_r[1] = { { 2, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_channel ch_h[2] = { { 1, NULL, POTSYNTH_THETA_DP, 0 }, { 1, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_channel ch_g[1] = { { 1, grf, POTSYNTH_THETA_DP, 0 } };
  potsynth_parameters params = potsynth_default_parameters();
  double mx[3] = { 0.0, 0.0, 0.0 };
  size_t i, c, e;

  for (i = 0; i <= lmax; ++i)
    {
      rf[i] = 1.0 + 0.01 * i;
      grf[i] = g * rf[i];
    }

  test_random_points(n, rmin, rmax, r, theta, phi);

  /* reference from single-set models */
  test_sets_single(lmin, lmax, mmax, coef, a, 1, ch_r, n, r, theta, phi, ref);
  test_sets_single(lmin, lmax, mmax, coef + nlm, a, 2, ch_h, n, r, theta, phi, ref + n);
  test_sets_single(lmin, lmax, mmax, coef + 2 * nlm, a, 2, ch_h, n, r, theta, phi, tmp);
  for (i = 0; i < n; ++i)
    ref[n + i] += f * tmp[n + i];
  test_sets_single(lmin, lmax, mmax, coef + 2 * nlm, a, 1, ch_g, n, r, theta, phi, tmp);
  for (i = 0; i < n; ++i)
    ref[2 * n + i] += tmp[i];

  for (c = 0; c < 3; ++c)
    {
      for (i = 0; i < n; ++i)
        mx[c] = fmax(mx[c], fabs(ref[c * n + i]));
    }

  for (e = 0; e < 4; ++e)
    {
      const int fast = (e < 3);
      const double tol = fast ? 3.0 * epsv[e] : 1.0e-13;
      potsynth_workspace * w;

      params.method = fast ? POTSYNTH_METHOD_FAST : POTSYNTH_METHOD_DIRECT;
      params.eps = fast ? epsv[e] : 1.0e-10;

      w = potsynth_alloc_sets(lmin, lmax, mmax, 3, coef, a, rmin, rmax, NULL, 3, o, &params);
      ps_test(w == NULL, "sets: %s %s allocation", desc, fast ? "fast" : "direct");

      if (w)
        {
          potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
          potsynth_info info;

          potsynth_get_info(&info, w);
          ps_test(info.nout != 3 || info.nchan != 3, "sets: %s number of outputs", desc);
          ps_test_int(potsynth_eval(n, r, theta, phi, out, ew), POTSYNTH_SUCCESS, "sets: %s %s evaluation",
                      desc, fast ? "fast" : "direct");

          for (c = 0; c < 3; ++c)
            {
              double emax = 0.0;

              for (i = 0; i < n; ++i)
                emax = fmax(emax, fabs(out[c * n + i] - ref[c * n + i]));

              ps_test(emax / mx[c] > tol, "sets: %s %s eps %.0e output %zu against single-set models (%.2e)",
                      desc, fast ? "fast" : "direct", params.eps, c, emax / mx[c]);
            }

          potsynth_eval_free(ew);
          potsynth_free(w);
        }
    }

  free(coef);
  free(rf);
  free(grf);
  free(r);
  free(theta);
  free(phi);
  free(ref);
  free(out);
  free(tmp);
}

/*
test_sets_same_op()
  Terms with the same operator: the output is that of the model with the
combined coefficients, for a scalar (P) and for the tensor operators
*/

static void
test_sets_same_op(void)
{
  const size_t lmin = 2, lmax = 36, mmax = 36;
  const double a = 6371.2, rmin = 6700.0, rmax = 6900.0;
  const size_t n = 2000, nlm = potsynth_nlm(lmin, lmax, mmax);
  double * coef = test_sets_coef(2, lmin, lmax, mmax);
  double * comb = malloc(nlm * sizeof(double));
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * ref = malloc(3 * n * sizeof(double));
  double * out = malloc(3 * n * sizeof(double));
  potsynth_term tp[2] = { { 0, 2.0, NULL, POTSYNTH_THETA_P, 0 }, { 1, -3.0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_term ta[2] = { { 0, 2.0, NULL, POTSYNTH_THETA_S2A, 0 }, { 1, -3.0, NULL, POTSYNTH_THETA_S2A, 0 } };
  potsynth_term tb[2] = { { 0, 2.0, NULL, POTSYNTH_THETA_S2B, 1 }, { 1, -3.0, NULL, POTSYNTH_THETA_S2B, 1 } };
  potsynth_output o[3] = { { 1, 2, tp, 0, 0.0 }, { 3, 2, ta, 0, 0.0 }, { 3, 2, tb, 0, 0.0 } };
  potsynth_channel ch[3] = { { 1, NULL, POTSYNTH_THETA_P, 0 }, { 3, NULL, POTSYNTH_THETA_S2A, 0 },
                             { 3, NULL, POTSYNTH_THETA_S2B, 1 } };
  potsynth_parameters params = potsynth_default_parameters();
  size_t i, c;
  int method;

  for (i = 0; i < nlm; ++i)
    comb[i] = 2.0 * coef[i] - 3.0 * coef[nlm + i];

  test_random_points(n, rmin, rmax, r, theta, phi);
  test_sets_single(lmin, lmax, mmax, comb, a, 3, ch, n, r, theta, phi, ref);

  for (method = POTSYNTH_METHOD_FAST; method <= POTSYNTH_METHOD_DIRECT; ++method)
    {
      const int fast = (method == POTSYNTH_METHOD_FAST);
      potsynth_workspace * w;

      params.method = method;
      params.eps = 1.0e-12;
      w = potsynth_alloc_sets(lmin, lmax, mmax, 2, coef, a, rmin, rmax, NULL, 3, o, &params);
      ps_test(w == NULL, "sets: same operator %s allocation", fast ? "fast" : "direct");

      if (w)
        {
          potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);

          potsynth_eval(n, r, theta, phi, out, ew);

          for (c = 0; c < 3; ++c)
            {
              double emax = 0.0, mx = 0.0;

              for (i = 0; i < n; ++i)
                {
                  emax = fmax(emax, fabs(out[c * n + i] - ref[c * n + i]));
                  mx = fmax(mx, fabs(ref[c * n + i]));
                }

              ps_test(emax / mx > (fast ? 1.0e-11 : 1.0e-13),
                      "sets: same operator %s output %zu against combined coefficients (%.2e)",
                      fast ? "fast" : "direct", c, emax / mx);
            }

          potsynth_eval_free(ew);
          potsynth_free(w);
        }
    }

  free(coef);
  free(comb);
  free(r);
  free(theta);
  free(phi);
  free(ref);
  free(out);
}

/*
test_sets_wrapper()
  potsynth_alloc() is potsynth_alloc_sets() with one set and one term of
factor 1 per output: the results are equal bit for bit. Quantities of a model
of two vector expansions (white box)
*/

static void
test_sets_wrapper(void)
{
  const size_t lmin = 1, lmax = 30, mmax = 20;
  const double a = 6371.2, rmin = 6650.0, rmax = 6950.0;
  const size_t n = 1500, nlm = potsynth_nlm(lmin, lmax, mmax);
  double * coef = test_sets_coef(6, lmin, lmax, mmax);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * o1 = malloc(3 * n * sizeof(double));
  double * o2 = malloc(3 * n * sizeof(double));
  potsynth_channel ch[3] = { { 2, NULL, POTSYNTH_THETA_P, 0 }, { 2, NULL, POTSYNTH_THETA_DP, 0 },
                             { 2, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_term t[3] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 }, { 0, 1.0, NULL, POTSYNTH_THETA_DP, 0 },
                         { 0, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_output o[3] = { { 2, 1, &t[0], 0, 0.0 }, { 2, 1, &t[1], 0, 0.0 }, { 2, 1, &t[2], 0, 0.0 } };
  potsynth_parameters params = potsynth_default_parameters();
  int method;
  size_t i;

  (void) nlm;

  test_random_points(n, rmin, rmax, r, theta, phi);

  for (method = POTSYNTH_METHOD_FAST; method <= POTSYNTH_METHOD_DIRECT; ++method)
    {
      potsynth_workspace * w1, * w2;

      params.method = method;
      w1 = potsynth_alloc(lmin, lmax, mmax, coef, a, rmin, rmax, 3, ch, &params);
      w2 = potsynth_alloc_sets(lmin, lmax, mmax, 1, coef, a, rmin, rmax, NULL, 3, o, &params);
      ps_test(w1 == NULL || w2 == NULL, "sets: wrapper allocation");

      if (w1 && w2)
        {
          potsynth_eval_workspace * e1 = potsynth_eval_alloc(0, w1);
          potsynth_eval_workspace * e2 = potsynth_eval_alloc(0, w2);
          size_t ndiff = 0;

          potsynth_eval(n, r, theta, phi, o1, e1);
          potsynth_eval(n, r, theta, phi, o2, e2);

          for (i = 0; i < 3 * n; ++i)
            ndiff += (o1[i] != o2[i]);

          ps_test(ndiff != 0, "sets: potsynth_alloc() equal to one set and one term (%s, %zu values differ)",
                  (method == POTSYNTH_METHOD_FAST) ? "fast" : "direct", ndiff);

          potsynth_eval_free(e1);
          potsynth_eval_free(e2);
        }

      if (w1)
        potsynth_free(w1);
      if (w2)
        potsynth_free(w2);
    }

  /* two vector expansions (sets a, b, c each): V_r is a quantity of its own, and
     (V_theta, V_phi) of each expansion one quantity */
  {
    potsynth_term tr0[1] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
    potsynth_term tt0[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_DP, 0 }, { 2, -1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
    potsynth_term tp0[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 }, { 2, 1.0, NULL, POTSYNTH_THETA_DP, 0 } };
    potsynth_term tr1[1] = { { 3, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
    potsynth_term tt1[2] = { { 4, 1.0, NULL, POTSYNTH_THETA_DP, 0 }, { 5, -1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
    potsynth_term tp1[2] = { { 4, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 }, { 5, 1.0, NULL, POTSYNTH_THETA_DP, 0 } };
    potsynth_output ov[6] = { { 0, 1, tr0, 0, 0.0 }, { 0, 2, tt0, 0, 0.0 }, { 0, 2, tp0, 0, 0.0 }, { 0, 1, tr1, 0, 0.0 }, { 0, 2, tt1, 0, 0.0 }, { 0, 2, tp1, 0, 0.0 } };
    const size_t expected[6] = { 0, 1, 1, 2, 3, 3 };
    size_t grp[6];
    potsynth_workspace * w;

    params.method = POTSYNTH_METHOD_DIRECT;
    w = potsynth_alloc_sets(lmin, lmax, mmax, 6, coef, a, rmin, rmax, NULL, 6, ov, &params);
    ps_test(w == NULL, "sets: two vector expansions allocation");

    if (w)
      {
        const size_t ngrp = potsynth_quantities(w, grp);
        size_t nbad = 0;

        for (i = 0; i < 6; ++i)
          nbad += (grp[i] != expected[i]);

        ps_test(ngrp != 4 || nbad != 0, "sets: quantities of two vector expansions (%zu quantities)", ngrp);
        potsynth_free(w);
      }
  }

  free(coef);
  free(r);
  free(theta);
  free(phi);
  free(o1);
  free(o2);
}

/* invalid arguments of potsynth_alloc_sets() */
static void
test_sets_errors(void)
{
  const size_t L = 8;
  double * coef = test_sets_coef(2, 1, L, L);
  potsynth_term t[2] = { { 0, 1.0, NULL, POTSYNTH_THETA_DP, 0 }, { 1, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_output o[1] = { { 1, 2, t, 0, 0.0 } };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;

  params.method = POTSYNTH_METHOD_DIRECT;
  test_last_error();

  w = potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params);
  ps_test(w == NULL, "sets: valid model accepted");
  if (w)
    potsynth_free(w);

  ps_test(potsynth_alloc_sets(1, L, L, 0, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: nset = 0 rejected");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "sets: nset = 0 error code");
  ps_test(potsynth_alloc_sets(1, L, L, 2, NULL, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: coef = NULL rejected");
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 0, o, &params) != NULL, "sets: nout = 0 rejected");
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, NULL, &params) != NULL, "sets: out = NULL rejected");

  o[0].nterm = 0;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: output without terms rejected");
  o[0].nterm = 2;

  t[1].set = 2;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: set out of range rejected");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "sets: set out of range error code");
  t[1].set = 1;

  t[1].fac = NAN;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: factor NaN rejected");
  t[1].fac = 1.0;

  t[1].theta = 99;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: invalid operator rejected");
  t[1].theta = POTSYNTH_THETA_PSIN;

  o[0].rpow = -2;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: lmin + rpow < 0 rejected");
  o[0].rpow = 1;

  /* terms that do not form one field component */
  t[0].theta = POTSYNTH_THETA_P;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: scalar + vector term rejected");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "sets: scalar + vector term error code");

  t[0].theta = POTSYNTH_THETA_DP;
  t[1].dphi = 0;
  ps_test(potsynth_alloc_sets(1, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: dP + P/sin without d/dphi rejected");
  t[1].dphi = 1;

  t[0].theta = POTSYNTH_THETA_S2A;
  t[1].theta = POTSYNTH_THETA_S2B;
  ps_test(potsynth_alloc_sets(2, L, L, 2, coef, 1.0, 1.0, 1.1, NULL, 1, o, &params) != NULL, "sets: S2A + S2B in one output rejected");

  test_last_error();
  free(coef);
}


/*
test_sets_rows()
  potsynth_eval_rows() against potsynth_eval() for the front-ends (with output
maps), fast and direct
*/

static void
test_sets_rows(void)
{
  const size_t L = 40, n = 1003;
  const double a = 6371.2, rmin = a + 300.0, rmax = a + 700.0;
  double * gh = test_random_coef(1, L, L);
  double * cf = test_random_coef(0, L, L);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * o1 = malloc(16 * n * sizeof(double));
  double * o2 = malloc(16 * n * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  int method, model;
  size_t i, k;

  test_random_points(n, rmin, rmax, r, theta, phi);
  r[20] = NAN;

  for (method = POTSYNTH_METHOD_FAST; method <= POTSYNTH_METHOD_DIRECT; ++method)
    {
      for (model = 0; model < 2; ++model)
        {
          potsynth_workspace * w;

          params.method = method;
          w = (model == 0) ?
              potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params) :
              potsynth_grav_alloc(0, L, L, cf, 3.986e5, a, rmin, rmax,
                                  POTSYNTH_GRAV_V | POTSYNTH_GRAV_G | POTSYNTH_GRAV_T, &params);
          ps_test(w == NULL, "rows: allocation");

          if (w)
            {
              potsynth_eval_workspace * ew = potsynth_eval_alloc(100, w);
              potsynth_info info;
              size_t ndiff = 0, nnan = 0;

              potsynth_get_info(&info, w);
              ps_test_int(potsynth_eval(n, r, theta, phi, o1, ew), POTSYNTH_SUCCESS, "rows: output-major");
              ps_test_int(potsynth_eval_rows(n, r, theta, phi, o2, ew), POTSYNTH_SUCCESS, "rows: point-major");

              for (i = 0; i < n; ++i)
                {
                  for (k = 0; k < info.nout; ++k)
                    {
                      const double x = o1[k * n + i], y = o2[i * info.nout + k];

                      if (isnan(x) && isnan(y))
                        ++nnan;
                      else
                        ndiff += (x != y);
                    }
                }

              ps_test(ndiff != 0 || nnan != info.nout,
                      "rows: %s %s, %zu outputs: equal to the output-major values (%zu differ, %zu NaN)",
                      (model == 0) ? "magnetic" : "gravity", (method == POTSYNTH_METHOD_FAST) ? "fast" : "direct",
                      info.nout, ndiff, nnan);

              potsynth_eval_free(ew);
              potsynth_free(w);
            }
        }
    }

  free(gh);
  free(cf);
  free(r);
  free(theta);
  free(phi);
  free(o1);
  free(o2);
}

/*
test_sets_components()
  Components given by the caller, and a horizontal tensor of two coefficient
sets (mixed operators, which need a component and weights)
*/

static void
test_sets_components(void)
{
  const size_t lmin = 2, lmax = 24, mmax = 24, n = 1500, nlm = potsynth_nlm(lmin, lmax, mmax);
  const double a = 6371.2, rmin = 6700.0, rmax = 6900.0;
  double * coef = test_sets_coef(2, lmin, lmax, mmax);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * ref = calloc(2 * n, sizeof(double));
  double * tmp = malloc(2 * n * sizeof(double));
  double * out = malloc(4 * n * sizeof(double));
  /* tensor of sets e (0) and b (1): D = e A - b B d/dphi, T = e B d/dphi + b A */
  potsynth_term td[2] = { { 0, 1.0, NULL, POTSYNTH_THETA_S2A, 0 }, { 1, -1.0, NULL, POTSYNTH_THETA_S2B, 1 } };
  potsynth_term tt[2] = { { 0, 1.0, NULL, POTSYNTH_THETA_S2B, 1 }, { 1, 1.0, NULL, POTSYNTH_THETA_S2A, 0 } };
  potsynth_term ts[1] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_output o[4] = { { 3, 2, td, 7, 1.0 }, { 3, 2, tt, 7, 2.0 }, { 3, 1, ts, 9, 1.0 }, { 3, 1, ts, 0, 0.0 } };
  potsynth_channel ch[2] = { { 3, NULL, POTSYNTH_THETA_S2A, 0 }, { 3, NULL, POTSYNTH_THETA_S2B, 1 } };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;
  int method;
  size_t i, c;

  test_random_points(n, rmin, rmax, r, theta, phi);

  test_sets_single(lmin, lmax, mmax, coef, a, 2, ch, n, r, theta, phi, ref);
  test_sets_single(lmin, lmax, mmax, coef + nlm, a, 2, ch, n, r, theta, phi, tmp);
  for (i = 0; i < n; ++i)
    {
      ref[i] -= tmp[n + i];
      ref[n + i] += tmp[i];
    }

  for (method = POTSYNTH_METHOD_FAST; method <= POTSYNTH_METHOD_DIRECT; ++method)
    {
      params.method = method;
      params.eps = 1.0e-11;
      w = potsynth_alloc_sets(lmin, lmax, mmax, 2, coef, a, rmin, rmax, NULL, 4, o, &params);
      ps_test(w == NULL, "components: tensor of two sets allocation");

      if (w)
        {
          potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
          potsynth_info info;
          size_t grp[4];

          potsynth_get_info(&info, w);
          potsynth_quantities(w, grp);
          ps_test(info.ncomp != 3 || grp[0] != grp[1] || grp[2] == grp[0] || grp[3] == grp[2],
                  "components: three components (%zu)", info.ncomp);

          potsynth_eval(n, r, theta, phi, out, ew);

          for (c = 0; c < 2; ++c)
            {
              double emax = 0.0, mx = 0.0;

              for (i = 0; i < n; ++i)
                {
                  emax = fmax(emax, fabs(out[c * n + i] - ref[c * n + i]));
                  mx = fmax(mx, fabs(ref[c * n + i]));
                }

              ps_test(emax / mx > ((method == POTSYNTH_METHOD_FAST) ? 3.0e-11 : 1.0e-13),
                      "components: tensor of two sets %s output %zu (%.2e)",
                      (method == POTSYNTH_METHOD_FAST) ? "fast" : "direct", c, emax / mx);
            }

          potsynth_eval_free(ew);
          potsynth_free(w);
        }
    }

  /* mixed tensor operators without a component and weight; a component with two radial powers */
  params.method = POTSYNTH_METHOD_DIRECT;
  test_last_error();
  o[0].component = 0;
  ps_test(potsynth_alloc_sets(lmin, lmax, mmax, 2, coef, a, rmin, rmax, NULL, 1, o, &params) != NULL,
          "components: mixed tensor operators need a component");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "components: error code");
  o[0].component = 7;
  o[1].rpow = 2;
  ps_test(potsynth_alloc_sets(lmin, lmax, mmax, 2, coef, a, rmin, rmax, NULL, 2, o, &params) != NULL,
          "components: one component with two radial powers rejected");
  o[1].rpow = 3;
  o[2].component = -1;
  ps_test(potsynth_alloc_sets(lmin, lmax, mmax, 2, coef, a, rmin, rmax, NULL, 3, o, &params) != NULL,
          "components: negative component rejected");
  test_last_error();

  free(coef);
  free(r);
  free(theta);
  free(phi);
  free(ref);
  free(tmp);
  free(out);
}

static void
test_sets(void)
{
  test_sets_vector("L = 40, shell", 1, 40, 40, 6621.2, 7171.2);
  test_sets_vector("L = 30, mmax = 12, lmin = 3", 3, 30, 12, 6800.0, 6900.0);
  test_sets_vector("L = 24, one radius", 1, 24, 24, 6900.0, 6900.0);
  test_sets_same_op();
  test_sets_wrapper();
  test_sets_errors();
  test_sets_rows();
  test_sets_components();
}
