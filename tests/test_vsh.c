/*
 * test_vsh.c
 *
 * Vector spherical harmonic front-end: fields against single-set models,
 * sign conventions against derivatives of scalar fields, radial functions
 */

/* outputs of a single-set model on the unit sphere with the direct method, out[c*n + i] */
static void
test_vsh_single(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef,
                const size_t nfunc, const potsynth_functional * func,
                const size_t n, const double * theta, const double * phi, double * out)
{
  double * ones = malloc(n * sizeof(double));
  size_t i;

  for (i = 0; i < n; ++i)
    ones[i] = 1.0;

  test_sets_single(lmin, lmax, mmax, coef, 1.0, nfunc, func, n, ones, theta, phi, out);
  free(ones);
}

/* coefficients of full normalization from those of normalization norm */
static double *
test_vsh_full(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef, const int norm)
{
  const size_t nlm = potsynth_nlm(lmin, lmax, mmax);
  double * c = malloc(nlm * sizeof(double));
  size_t l;

  for (l = lmin; l <= lmax; ++l)
    {
      const long M = (long) ((l < mmax) ? l : mmax);
      const double f = (norm == POTSYNTH_VSH_NORM_SCHMIDT) ? 1.0 / sqrt(2.0 * l + 1.0) : 1.0;
      long m;

      for (m = -M; m <= M; ++m)
        {
          const size_t idx = potsynth_idx((long) l, m, (long) lmin, (long) mmax);
          c[idx] = f * coef[idx];
        }
    }

  return c;
}

/*
test_vsh_field()
  Two fields on the sphere with all quantities against single-set models

Inputs: desc   - description
        lmin   - minimum degree
        lmax   - maximum degree
        mmax   - maximum order
        norm   - normalization
        have   - bit s set: array s (a, b, c) is given
        method - method
        eps    - accuracy
        tol    - tolerance
*/

static void
test_vsh_field(const char * desc, const size_t lmin, const size_t lmax, const size_t mmax,
               const int norm, const int have, const int method, const double eps, const double tol)
{
  const size_t nexp = 2, n = 2000, nlm = potsynth_nlm(lmin, lmax, mmax);
  const int quantities = POTSYNTH_VSH_V | POTSYNTH_VSH_DIV | POTSYNTH_VSH_CURL;
  const size_t nq = potsynth_vsh_nout(quantities);
  double * abc[3];
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * ref = calloc(nexp * nq * n, sizeof(double));
  double * out = malloc(nexp * nq * n * sizeof(double));
  double * tmp = malloc(2 * n * sizeof(double));
  double * rl = malloc((lmax + 1) * sizeof(double));
  potsynth_functional ch_p[1] = { { 0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_functional ch_h[2] = { { 0, NULL, POTSYNTH_THETA_DP, 0 }, { 0, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_functional ch_l[1] = { { 0, rl, POTSYNTH_THETA_P, 0 } };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;
  size_t e, i, q, s;

  for (i = 0; i <= lmax; ++i)
    rl[i] = -(double) i * (double) (i + 1);

  for (s = 0; s < 3; ++s)
    abc[s] = test_sets_coef(nexp, lmin, lmax, mmax);

  test_random_points(n, 1.0, 1.0, r, theta, phi);

  /* reference: V_r, V_theta, V_phi, div, curl of each field */
  for (e = 0; e < nexp; ++e)
    {
      double * R = ref + e * nq * n;
      double * fa = test_vsh_full(lmin, lmax, mmax, abc[0] + e * nlm, norm);
      double * fb = test_vsh_full(lmin, lmax, mmax, abc[1] + e * nlm, norm);
      double * fc = test_vsh_full(lmin, lmax, mmax, abc[2] + e * nlm, norm);

      if (have & 1)
        test_vsh_single(lmin, lmax, mmax, fa, 1, ch_p, n, theta, phi, R);

      if (have & 2)
        {
          test_vsh_single(lmin, lmax, mmax, fb, 2, ch_h, n, theta, phi, tmp);
          for (i = 0; i < n; ++i)
            {
              R[n + i] += tmp[i];
              R[2 * n + i] += tmp[n + i];
            }
          test_vsh_single(lmin, lmax, mmax, fb, 1, ch_l, n, theta, phi, R + 3 * n);
        }

      if (have & 4)
        {
          test_vsh_single(lmin, lmax, mmax, fc, 2, ch_h, n, theta, phi, tmp);
          for (i = 0; i < n; ++i)
            {
              R[n + i] -= tmp[n + i];
              R[2 * n + i] += tmp[i];
            }
          test_vsh_single(lmin, lmax, mmax, fc, 1, ch_l, n, theta, phi, R + 4 * n);
        }

      free(fa);
      free(fb);
      free(fc);
    }

  params.method = method;
  params.eps = eps;
  w = potsynth_vsh_alloc(lmin, lmax, mmax, nexp, (have & 1) ? abc[0] : NULL, (have & 2) ? abc[1] : NULL,
                         (have & 4) ? abc[2] : NULL, norm, quantities, 0.0, 0.0, NULL, &params);
  ps_test(w == NULL, "vsh: %s allocation", desc);

  if (w)
    {
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
      potsynth_info info;

      potsynth_get_info(&info, w);
      ps_test(info.nout != nexp * nq || info.ncomp != nexp * 4, "vsh: %s outputs and components (%zu, %zu)",
              desc, info.nout, info.ncomp);

      /* fields on the sphere: no radii */
      ps_test_int(potsynth_eval(n, NULL, theta, phi, out, ew), POTSYNTH_SUCCESS, "vsh: %s evaluation", desc);

      for (q = 0; q < nq; ++q)
        {
          double emax = 0.0, mx = 0.0;

          for (e = 0; e < nexp; ++e)
            {
              for (i = 0; i < n; ++i)
                {
                  emax = fmax(emax, fabs(out[(e * nq + q) * n + i] - ref[(e * nq + q) * n + i]));
                  mx = fmax(mx, fabs(ref[(e * nq + q) * n + i]));
                }
            }

          if (mx > 0.0)
            ps_test(emax / mx > tol, "vsh: %s quantity %zu against single-set models (%.2e)", desc, q, emax / mx);
          else
            ps_test(emax != 0.0, "vsh: %s quantity %zu is zero", desc, q);
        }

      potsynth_eval_free(ew);
      potsynth_free(w);
    }

  for (s = 0; s < 3; ++s)
    free(abc[s]);
  free(r);
  free(theta);
  free(phi);
  free(ref);
  free(out);
  free(tmp);
  free(rl);
}

/*
test_vsh_signs()
  Sign conventions against numerical derivatives: with the scalar fields
S_b = sum b Y, S_c = sum c Y,

V_theta = dS_b/dtheta - (1/sin) dS_c/dphi,  V_phi = (1/sin) dS_b/dphi + dS_c/dtheta
div  = (1/sin) [ d(sin V_theta)/dtheta + dV_phi/dphi ]
curl = (1/sin) [ d(sin V_phi)/dtheta - dV_theta/dphi ]
*/

static void
test_vsh_signs(void)
{
  const size_t L = 10, n = 300, nlm = potsynth_nlm(1, L, L);
  const double h = 1.0e-5;
  const int quantities = POTSYNTH_VSH_V | POTSYNTH_VSH_DIV | POTSYNTH_VSH_CURL;
  double * b = test_random_coef(1, L, L);
  double * c = test_random_coef(1, L, L);
  double * theta = malloc(5 * n * sizeof(double));
  double * phi = malloc(5 * n * sizeof(double));
  double * sb = malloc(5 * n * sizeof(double));
  double * sc = malloc(5 * n * sizeof(double));
  double * out = malloc(5 * 5 * n * sizeof(double));
  potsynth_functional ch_p[1] = { { 0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;
  double e[4] = { 0.0, 0.0, 0.0, 0.0 }, mx[4] = { 0.0, 0.0, 0.0, 0.0 };
  size_t i;

  (void) nlm;

  /* points i, and shifted by +-h in theta (n + i, 2n + i) and in phi (3n + i, 4n + i) */
  for (i = 0; i < n; ++i)
    {
      const double t = 0.3 + 2.5 * test_urand(), p = 2.0 * M_PI * test_urand();

      theta[i] = t; phi[i] = p;
      theta[n + i] = t + h; phi[n + i] = p;
      theta[2 * n + i] = t - h; phi[2 * n + i] = p;
      theta[3 * n + i] = t; phi[3 * n + i] = p + h;
      theta[4 * n + i] = t; phi[4 * n + i] = p - h;
    }

  test_vsh_single(1, L, L, b, 1, ch_p, 5 * n, theta, phi, sb);
  test_vsh_single(1, L, L, c, 1, ch_p, 5 * n, theta, phi, sc);

  params.method = POTSYNTH_METHOD_DIRECT;
  w = potsynth_vsh_alloc(1, L, L, 1, NULL, b, c, POTSYNTH_VSH_NORM_FULL, quantities, 0.0, 0.0, NULL, &params);
  ps_test(w == NULL, "vsh: signs allocation");

  if (w)
    {
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
      const size_t N = 5 * n;
      const double * Vt = out + N, * Vp = out + 2 * N, * dv = out + 3 * N, * cu = out + 4 * N;

      potsynth_eval(N, NULL, theta, phi, out, ew);

      for (i = 0; i < n; ++i)
        {
          const double s = sin(theta[i]), sp = sin(theta[i] + h), sm = sin(theta[i] - h);
          const double dbt = (sb[n + i] - sb[2 * n + i]) / (2.0 * h), dbp = (sb[3 * n + i] - sb[4 * n + i]) / (2.0 * h);
          const double dct = (sc[n + i] - sc[2 * n + i]) / (2.0 * h), dcp = (sc[3 * n + i] - sc[4 * n + i]) / (2.0 * h);
          const double vt = dbt - dcp / s, vp = dbp / s + dct;
          const double div = ((sp * Vt[n + i] - sm * Vt[2 * n + i]) + (Vp[3 * n + i] - Vp[4 * n + i])) / (2.0 * h * s);
          const double curl = ((sp * Vp[n + i] - sm * Vp[2 * n + i]) - (Vt[3 * n + i] - Vt[4 * n + i])) / (2.0 * h * s);

          e[0] = fmax(e[0], fabs(Vt[i] - vt)); mx[0] = fmax(mx[0], fabs(vt));
          e[1] = fmax(e[1], fabs(Vp[i] - vp)); mx[1] = fmax(mx[1], fabs(vp));
          e[2] = fmax(e[2], fabs(dv[i] - div)); mx[2] = fmax(mx[2], fabs(div));
          e[3] = fmax(e[3], fabs(cu[i] - curl)); mx[3] = fmax(mx[3], fabs(curl));
        }

      ps_test(e[0] / mx[0] > 1.0e-6, "vsh: V_theta = dS_b/dtheta - (1/sin) dS_c/dphi (%.2e)", e[0] / mx[0]);
      ps_test(e[1] / mx[1] > 1.0e-6, "vsh: V_phi = (1/sin) dS_b/dphi + dS_c/dtheta (%.2e)", e[1] / mx[1]);
      ps_test(e[2] / mx[2] > 1.0e-6, "vsh: divergence against derivatives of the field (%.2e)", e[2] / mx[2]);
      ps_test(e[3] / mx[3] > 1.0e-6, "vsh: curl against derivatives of the field (%.2e)", e[3] / mx[3]);

      potsynth_eval_free(ew);
      potsynth_free(w);
    }

  free(b);
  free(c);
  free(theta);
  free(phi);
  free(sb);
  free(sc);
  free(out);
}

/* a cubic polynomial and its derivative */
static double
test_vsh_poly(const double r, const int deriv)
{
  const double x = r - 6800.0;
  return deriv ? (0.3 - 2.0e-3 * x + 1.5e-5 * x * x) : (2.0 + 0.3 * x - 1.0e-3 * x * x + 5.0e-6 * x * x * x);
}

/* radial functions: reproduction of polynomials and of their derivatives, range */
static void
test_vsh_radial(void)
{
  const double rn[4] = { 6700.0, 6780.0, 6850.0, 6900.0 };
  double wb[4], phi[8], dphi[8];
  test_radial_nodes nodes = { 4, rn, wb };
  size_t first = 0, cnt = 0, i, j;
  double ev[3] = { 0.0, 0.0, 0.0 }, ed[3] = { 0.0, 0.0, 0.0 }, es = 0.0;
  int k;

  ps_test_int(test_bary_weights(4, rn, wb), POTSYNTH_SUCCESS, "vsh: barycentric weights");

  for (i = 0; i <= 40; ++i)
    {
      /* includes the nodes 6700 and 6900 and passes 6780, 6850 closely */
      const double r = (i == 13) ? rn[1] : (i == 27) ? rn[2] : 6700.0 + 5.0 * i;
      double v, d, s;

      /* barycentric: cubic polynomials are reproduced */
      test_radial_bary(r, &first, &cnt, phi, dphi, &nodes);
      v = d = s = 0.0;
      for (j = 0; j < cnt; ++j)
        {
          v += phi[j] * test_vsh_poly(rn[first + j], 0);
          d += dphi[j] * test_vsh_poly(rn[first + j], 0);
          s += phi[j];
        }
      ev[0] = fmax(ev[0], fabs(v - test_vsh_poly(r, 0)));
      ed[0] = fmax(ed[0], fabs(d - test_vsh_poly(r, 1)));
      es = fmax(es, fabs(s - 1.0));

      /* linear: linear functions are reproduced */
      test_radial_linear(r, &first, &cnt, phi, dphi, &nodes);
      v = d = 0.0;
      for (j = 0; j < cnt; ++j)
        {
          v += phi[j] * (3.0 + 0.5 * rn[first + j]);
          d += dphi[j] * (3.0 + 0.5 * rn[first + j]);
        }
      ev[1] = fmax(ev[1], fabs(v - (3.0 + 0.5 * r)));
      ed[1] = fmax(ed[1], fabs(d - 0.5));
      ps_test(cnt != 2, "vsh: linear radial functions: two are nonzero");

      /* Hermite: cubic polynomials are reproduced from values and derivatives */
      test_radial_hermite(r, &first, &cnt, phi, dphi, &nodes);
      v = d = 0.0;
      for (j = 0; j < cnt; ++j)
        {
          const size_t node = (first + j) / 2;
          const double cf = test_vsh_poly(rn[node], (int) ((first + j) % 2));
          v += phi[j] * cf;
          d += dphi[j] * cf;
        }
      ev[2] = fmax(ev[2], fabs(v - test_vsh_poly(r, 0)));
      ed[2] = fmax(ed[2], fabs(d - test_vsh_poly(r, 1)));
      ps_test(cnt != 4 || first % 2 != 0, "vsh: Hermite radial functions: four are nonzero");
    }

  ps_test(ev[0] > 1.0e-9, "vsh: barycentric radial functions reproduce a cubic (%.2e)", ev[0]);
  ps_test(ed[0] > 1.0e-10, "vsh: barycentric radial derivatives (%.2e)", ed[0]);
  ps_test(es > 1.0e-13, "vsh: barycentric radial functions sum to 1 (%.2e)", es);
  ps_test(ev[1] > 1.0e-9, "vsh: linear radial functions reproduce a linear function (%.2e)", ev[1]);
  ps_test(ed[1] > 1.0e-12, "vsh: linear radial derivatives (%.2e)", ed[1]);
  ps_test(ev[2] > 1.0e-9, "vsh: Hermite radial functions reproduce a cubic (%.2e)", ev[2]);
  ps_test(ed[2] > 1.0e-10, "vsh: Hermite radial derivatives (%.2e)", ed[2]);

  for (k = 0; k < 2; ++k)
    {
      const double r = k ? 6900.001 : 6699.999;

      ps_test_int(test_radial_bary(r, &first, &cnt, phi, NULL, &nodes), POTSYNTH_EDOM,
                  "vsh: barycentric radius outside of the nodes");
      ps_test_int(test_radial_linear(r, &first, &cnt, phi, NULL, &nodes), POTSYNTH_EDOM,
                  "vsh: linear radius outside of the nodes");
      ps_test_int(test_radial_hermite(r, &first, &cnt, phi, NULL, &nodes), POTSYNTH_EDOM,
                  "vsh: Hermite radius outside of the nodes");
    }
}

/*
test_vsh_shell()
  Fields with radial functions: coefficient sets p(r_j) x (a, b, c) (and
p'(r_j) x (a, b, c) for the Hermite functions) of a cubic p give
V(r, theta, phi) = p(r) V_0(theta, phi)

Inputs: desc   - description
        radial - radial functions
        nper   - coefficient sets per node (2 for Hermite)
        method - method
*/

static void
test_vsh_shell(const char * desc, potsynth_radial_fn * radial, const size_t nper, const int method)
{
  const size_t lmin = 1, lmax = 20, mmax = 20, nn = 4, nrad = nn * nper, n = 1000;
  const size_t nlm = potsynth_nlm(lmin, lmax, mmax);
  const double rn[4] = { 6700.0, 6780.0, 6850.0, 6900.0 };
  const int linear = (radial == test_radial_linear);
  double wb[4];
  test_radial_nodes nodes = { 4, rn, wb };
  double * abc0[3], * abc[3];
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * v0 = malloc(3 * n * sizeof(double));
  double * out = malloc(3 * n * sizeof(double));
  double * out2 = malloc(3 * n * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  const potsynth_radial rad = { nrad, radial, &nodes };
  potsynth_workspace * w0, * w;
  size_t i, j, k, s, q;

  test_bary_weights(4, rn, wb);

  for (s = 0; s < 3; ++s)
    {
      abc0[s] = test_random_coef(lmin, lmax, mmax);
      abc[s] = malloc(nrad * nlm * sizeof(double));

      for (j = 0; j < nrad; ++j)
        {
          const double x = rn[j / nper];
          const double f = linear ? ((j % nper) ? 0.0 : 3.0 + 0.5 * x) : test_vsh_poly(x, (int) (j % nper));

          for (k = 0; k < nlm; ++k)
            abc[s][j * nlm + k] = f * abc0[s][k];
        }
    }

  test_random_points(n, rn[0], rn[3], r, theta, phi);
  for (i = 0; i < n; ++i)
    r[i] = (i % 97 == 0) ? rn[(i / 97) % 4] : rn[0] + (rn[3] - rn[0]) * test_urand();

  params.method = method;
  params.eps = 1.0e-11;
  w0 = potsynth_vsh_alloc(lmin, lmax, mmax, 1, abc0[0], abc0[1], abc0[2], POTSYNTH_VSH_NORM_SCHMIDT,
                          POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params);
  w = potsynth_vsh_alloc(lmin, lmax, mmax, 1, abc[0], abc[1], abc[2], POTSYNTH_VSH_NORM_SCHMIDT,
                         POTSYNTH_VSH_V, rn[0], rn[3], &rad, &params);
  ps_test(w0 == NULL || w == NULL, "vsh: %s allocation", desc);

  if (w0 && w)
    {
      potsynth_eval_workspace * e0 = potsynth_eval_alloc(0, w0);
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
      potsynth_eval_workspace * e7 = potsynth_eval_alloc(7, w);
      potsynth_info info;
      double emax = 0.0, mx = 0.0, sav;
      size_t ndiff = 0;

      potsynth_get_info(&info, w);
      if (method == POTSYNTH_METHOD_FAST)
        ps_test(info.rho != nrad, "vsh: %s radial functions (%zu)", desc, info.rho);
      ps_test(info.nout != 3 || info.ncomp != 2, "vsh: %s outputs and components", desc);

      potsynth_eval(n, NULL, theta, phi, v0, e0);
      ps_test_int(potsynth_eval(n, r, theta, phi, out, ew), POTSYNTH_SUCCESS, "vsh: %s evaluation", desc);
      ps_test_int(potsynth_eval(n, r, theta, phi, out2, e7), POTSYNTH_SUCCESS, "vsh: %s chunks of 7 points", desc);

      for (i = 0; i < n; ++i)
        {
          const double p = linear ? 3.0 + 0.5 * r[i] : test_vsh_poly(r[i], 0);

          for (q = 0; q < 3; ++q)
            {
              emax = fmax(emax, fabs(out[q * n + i] - p * v0[q * n + i]));
              mx = fmax(mx, fabs(p * v0[q * n + i]));
              ndiff += (out[q * n + i] != out2[q * n + i]);
            }
        }

      ps_test(emax / mx > 1.0e-10, "vsh: %s V(r) = p(r) V_0 (%.2e)", desc, emax / mx);
      ps_test(ndiff != 0, "vsh: %s independent of the chunk size (%zu values differ)", desc, ndiff);

      /* point-major output */
      ndiff = 0;
      ps_test_int(potsynth_eval_rows(n, r, theta, phi, out2, e7), POTSYNTH_SUCCESS, "vsh: %s rows", desc);
      for (i = 0; i < n; ++i)
        {
          for (q = 0; q < 3; ++q)
            ndiff += (out[q * n + i] != out2[i * 3 + q]);
        }
      ps_test(ndiff != 0, "vsh: %s rows equal to the output-major values (%zu values differ)", desc, ndiff);

      /* the radii are needed with radial functions */
      test_last_error();
      ps_test_int(potsynth_eval(n, NULL, theta, phi, out2, ew), POTSYNTH_EINVAL, "vsh: %s radii not given", desc);
      test_last_error();

      /* radius outside of the nodes: error, output untouched; NaN radius: NaN */
      test_last_error();
      sav = r[5];
      r[5] = rn[3] + 1.0;
      out[0] = -12345.0;
      ps_test_int(potsynth_eval(n, r, theta, phi, out, ew), POTSYNTH_EDOM, "vsh: %s radius outside", desc);
      ps_test(out[0] != -12345.0, "vsh: %s output untouched on error", desc);
      ps_test_int(test_last_error(), POTSYNTH_EDOM, "vsh: %s error reported", desc);
      r[5] = NAN;
      ps_test_int(potsynth_eval(n, r, theta, phi, out, ew), POTSYNTH_SUCCESS, "vsh: %s NaN radius", desc);
      ps_test(!isnan(out[5]) || !isnan(out[n + 5]) || isnan(out[6]), "vsh: %s NaN radius gives NaN", desc);
      r[5] = sav;

      potsynth_eval_free(e0);
      potsynth_eval_free(ew);
      potsynth_eval_free(e7);
    }

  if (w0)
    potsynth_free(w0);
  if (w)
    potsynth_free(w);

  for (s = 0; s < 3; ++s)
    {
      free(abc0[s]);
      free(abc[s]);
    }
  free(r);
  free(theta);
  free(phi);
  free(v0);
  free(out);
  free(out2);
}

/* invalid arguments, colatitude range */
static void
test_vsh_errors(void)
{
  const size_t L = 8;
  double * a = test_sets_coef(2, 1, L, L);
  const potsynth_radial norad = { 2, NULL, NULL };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;

  params.method = POTSYNTH_METHOD_FAST;
  test_last_error();

  ps_test_int((int) potsynth_vsh_nout(POTSYNTH_VSH_V | POTSYNTH_VSH_CURL), 4, "vsh: number of outputs");

  ps_test(potsynth_vsh_alloc(1, L, L, 0, a, a, a, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params) != NULL,
          "vsh: nexp = 0 rejected");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "vsh: nexp = 0 error code");
  ps_test(potsynth_vsh_alloc(1, L, L, 1, NULL, NULL, NULL, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params) != NULL,
          "vsh: no coefficients rejected");
  ps_test(potsynth_vsh_alloc(1, L, L, 1, a, a, a, POTSYNTH_VSH_NORM_FULL, 0, 0.0, 0.0, NULL, &params) != NULL,
          "vsh: no quantities rejected");
  ps_test(potsynth_vsh_alloc(1, L, L, 1, a, a, a, 7, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params) != NULL,
          "vsh: invalid normalization rejected");
  ps_test(potsynth_vsh_alloc(1, L, L, 1, a, a, a, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V, 1.0, 1.0, &norad, &params) != NULL,
          "vsh: two radial functions without a function rejected");
  ps_test(potsynth_vsh_alloc(1, L, L + 1, 1, a, a, a, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params) != NULL,
          "vsh: mmax > lmax rejected");
  test_last_error();

  /* polar cap */
  params.theta_max = 0.5;
  w = potsynth_vsh_alloc(1, L, L, 2, a, NULL, a, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params);
  ps_test(w == NULL, "vsh: polar cap allocation");

  if (w)
    {
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, w);
      double th[2] = { 0.5, 0.2 }, ph[2] = { 1.0, 2.0 }, out[12];
      potsynth_info info;

      potsynth_get_info(&info, w);
      ps_test_int((int) info.nout, 6, "vsh: outputs of the model");
      ps_test_int(potsynth_eval(2, NULL, th, ph, out, ew), POTSYNTH_SUCCESS, "vsh: points on the cap");
      th[1] = 0.6;
      ps_test_int(potsynth_eval(2, NULL, th, ph, out, ew), POTSYNTH_EDOM, "vsh: point outside of the cap");
      ps_test_int(potsynth_eval_rows(2, NULL, th, ph, out, ew), POTSYNTH_EDOM, "vsh: rows outside of the cap");
      test_last_error();

      potsynth_eval_free(ew);
      potsynth_free(w);
    }

  free(a);
}

/* more components than the window check tests: nexp fields of degree L against the
   direct method (40 fields; 350 fields have more than 1024 grid layers, for which
   the window sum adds node after node) */
static void
test_vsh_many(const size_t L, const size_t nexp)
{
  const size_t n = 600, nq = 3;
  double * abc[3];
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * r = malloc(n * sizeof(double));
  double * of = malloc(nexp * nq * n * sizeof(double));
  double * od = malloc(nexp * nq * n * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * wf, * wd;
  size_t e, i, q, s;

  for (s = 0; s < 3; ++s)
    abc[s] = test_sets_coef(nexp, 1, L, L);

  /* fields of very different sizes: the accuracy is relative to each */
  for (e = 0; e < nexp; ++e)
    {
      const double f = pow(10.0, -8.0 * (double) e / (double) nexp);

      for (s = 0; s < 3; ++s)
        {
          for (i = 0; i < potsynth_nlm(1, L, L); ++i)
            abc[s][e * potsynth_nlm(1, L, L) + i] *= f;
        }
    }

  test_random_points(n, 1.0, 1.0, r, theta, phi);

  params.eps = 1.0e-10;
  params.method = POTSYNTH_METHOD_FAST;
  wf = potsynth_vsh_alloc(1, L, L, nexp, abc[0], abc[1], abc[2], POTSYNTH_VSH_NORM_SCHMIDT, POTSYNTH_VSH_V,
                          0.0, 0.0, NULL, &params);
  params.method = POTSYNTH_METHOD_DIRECT;
  wd = potsynth_vsh_alloc(1, L, L, nexp, abc[0], abc[1], abc[2], POTSYNTH_VSH_NORM_SCHMIDT, POTSYNTH_VSH_V,
                          0.0, 0.0, NULL, &params);
  ps_test(wf == NULL || wd == NULL, "vsh: %zu fields allocation", nexp);

  if (wf && wd)
    {
      potsynth_eval_workspace * ef = potsynth_eval_alloc(0, wf);
      potsynth_eval_workspace * ed = potsynth_eval_alloc(0, wd);
      potsynth_info info;
      double worst = 0.0;

      potsynth_get_info(&info, wf);
      ps_test(info.ncomp != 2 * nexp || info.nout != nq * nexp, "vsh: %zu fields components (%zu)", nexp, info.ncomp);

      potsynth_eval(n, NULL, theta, phi, of, ef);
      potsynth_eval(n, NULL, theta, phi, od, ed);

      /* each field: V_r, and the horizontal vector, relative to its own largest value */
      for (e = 0; e < nexp; ++e)
        {
          double er = 0.0, mr = 0.0, eh = 0.0, mh = 0.0;

          for (i = 0; i < n; ++i)
            {
              for (q = 0; q < nq; ++q)
                {
                  const double d = fabs(of[(e * nq + q) * n + i] - od[(e * nq + q) * n + i]);
                  const double v = fabs(od[(e * nq + q) * n + i]);

                  if (q == 0)
                    {
                      er = fmax(er, d);
                      mr = fmax(mr, v);
                    }
                  else
                    {
                      eh = fmax(eh, d);
                      mh = fmax(mh, v);
                    }
                }
            }

          worst = fmax(worst, fmax(er / mr, eh / mh));
        }

      ps_test(worst > 3.0e-10, "vsh: %zu fields over 8 decades, each against the direct method (%.2e)", nexp, worst);

      potsynth_eval_free(ef);
      potsynth_eval_free(ed);
    }

  if (wf)
    potsynth_free(wf);
  if (wd)
    potsynth_free(wd);

  for (s = 0; s < 3; ++s)
    free(abc[s]);
  free(theta);
  free(phi);
  free(r);
  free(of);
  free(od);
}

/* radial and horizontal parts requested on their own: the values of the whole field */
static void
test_vsh_parts(void)
{
  const size_t L = 18, nexp = 3, n = 500;
  double * b = test_sets_coef(nexp, 1, L, L);
  double * c = test_sets_coef(nexp, 1, L, L);
  double * a = test_sets_coef(nexp, 1, L, L);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * ov = malloc(nexp * 3 * n * sizeof(double));
  double * op = malloc(nexp * 3 * n * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  const int quant[3] = { POTSYNTH_VSH_VH, POTSYNTH_VSH_VR, POTSYNTH_VSH_VR | POTSYNTH_VSH_CURL };
  const size_t nq[3] = { 2, 1, 2 };
  /* output q of each request in the whole field (V_r, V_theta, V_phi); 3: the curl, not compared */
  const size_t qmap[3][2] = { { 1, 2 }, { 0, 0 }, { 0, 3 } };
  potsynth_workspace * wv;
  size_t e, i, k, q;

  test_random_points(n, 1.0, 1.0, r, theta, phi);
  params.method = POTSYNTH_METHOD_DIRECT;
  wv = potsynth_vsh_alloc(1, L, L, nexp, a, b, c, POTSYNTH_VSH_NORM_SCHMIDT, POTSYNTH_VSH_V, 0.0, 0.0, NULL, &params);
  ps_test(wv == NULL, "vsh: parts allocation");

  if (wv)
    {
      potsynth_eval_workspace * ev = potsynth_eval_alloc(0, wv);

      potsynth_eval(n, NULL, theta, phi, ov, ev);

      for (k = 0; k < 3; ++k)
        {
          potsynth_workspace * wp = potsynth_vsh_alloc(1, L, L, nexp, a, b, c, POTSYNTH_VSH_NORM_SCHMIDT, quant[k],
                                                       0.0, 0.0, NULL, &params);

          ps_test(wp == NULL || potsynth_vsh_nout(quant[k]) != nq[k], "vsh: part %zu allocation", k);

          if (wp)
            {
              potsynth_eval_workspace * ep = potsynth_eval_alloc(0, wp);
              potsynth_info info;
              size_t ndiff = 0;

              potsynth_get_info(&info, wp);
              potsynth_eval(n, NULL, theta, phi, op, ep);

              for (e = 0; e < nexp; ++e)
                {
                  for (q = 0; q < nq[k]; ++q)
                    {
                      if (qmap[k][q] == 3)
                        continue;

                      for (i = 0; i < n; ++i)
                        ndiff += (op[(e * nq[k] + q) * n + i] != ov[(e * 3 + qmap[k][q]) * n + i]);
                    }
                }

              ps_test(info.nout != nexp * nq[k] || ndiff != 0,
                      "vsh: part %zu, %zu outputs, equal to the whole field (%zu differ)", k, info.nout, ndiff);

              potsynth_eval_free(ep);
              potsynth_free(wp);
            }
        }

      potsynth_eval_free(ev);
      potsynth_free(wv);
    }

  free(a);
  free(b);
  free(c);
  free(r);
  free(theta);
  free(phi);
  free(ov);
  free(op);
}

static void
test_vsh(void)
{
  test_vsh_field("L = 30 Schmidt fast", 1, 30, 30, POTSYNTH_VSH_NORM_SCHMIDT, 7, POTSYNTH_METHOD_FAST, 1.0e-10, 3.0e-10);
  test_vsh_field("L = 30 Schmidt direct", 1, 30, 30, POTSYNTH_VSH_NORM_SCHMIDT, 7, POTSYNTH_METHOD_DIRECT, 1.0e-10, 1.0e-13);
  test_vsh_field("L = 24 full, lmin = 0, mmax = 10, fast", 0, 24, 10, POTSYNTH_VSH_NORM_FULL, 7, POTSYNTH_METHOD_FAST, 1.0e-13, 3.0e-13);
  test_vsh_field("L = 20 no radial part", 1, 20, 20, POTSYNTH_VSH_NORM_SCHMIDT, 6, POTSYNTH_METHOD_FAST, 1.0e-10, 3.0e-10);
  test_vsh_field("L = 20 toroidal only", 1, 20, 20, POTSYNTH_VSH_NORM_SCHMIDT, 4, POTSYNTH_METHOD_FAST, 1.0e-10, 3.0e-10);
  test_vsh_field("L = 20 radial only", 1, 20, 20, POTSYNTH_VSH_NORM_FULL, 1, POTSYNTH_METHOD_DIRECT, 1.0e-10, 1.0e-13);
  test_vsh_signs();
  test_vsh_radial();
  test_vsh_shell("barycentric shell fast", test_radial_bary, 1, POTSYNTH_METHOD_FAST);
  test_vsh_shell("barycentric shell direct", test_radial_bary, 1, POTSYNTH_METHOD_DIRECT);
  test_vsh_shell("linear shell", test_radial_linear, 1, POTSYNTH_METHOD_FAST);
  test_vsh_shell("Hermite shell", test_radial_hermite, 2, POTSYNTH_METHOD_FAST);
  test_vsh_many(16, 40);
  test_vsh_many(8, 350);
  test_vsh_parts();
  test_vsh_errors();
}
