/*
 * test_region.c
 *
 * Colatitude range of the evaluation points (theta_min, theta_max of the
 * parameters): the fast method stores only the grid rows of the range
 */

/* random points with colatitudes in [t0,t1], the first four at the ends of the range */
static void
test_region_points(const size_t n, const double rmin, const double rmax, const double t0, const double t1,
                   double * r, double * theta, double * phi)
{
  size_t i;

  for (i = 0; i < n; ++i)
    {
      r[i] = rmin + (rmax - rmin) * test_urand();
      theta[i] = acos(cos(t0) - (cos(t0) - cos(t1)) * test_urand());
      theta[i] = fmin(fmax(theta[i], t0), t1);
      phi[i] = 2.0 * M_PI * test_urand();
    }

  theta[0] = t0;
  theta[1] = t1;
  theta[2] = t0;
  phi[2] = 0.0;
  theta[3] = t1;
  phi[3] = 2.0 * M_PI - 1.0e-14;
}

/*
test_region_range()
  A model with a colatitude range against the model of the whole sphere and
against the direct method, at points inside the range

Inputs: desc - description
        t0   - smallest colatitude
        t1   - largest colatitude
*/

static void
test_region_range(const char * desc, const double t0, const double t1)
{
  const size_t L = 60, n = 3000;
  const size_t nout = potsynth_mag_nout(POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB);
  const double a = 6371.2, rmin = a + 300.0, rmax = a + 800.0;
  const double epsv[3] = { 1.0e-6, 1.0e-10, 1.0e-13 };
  double * gh = test_random_coef(1, L, L);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * od = malloc(nout * n * sizeof(double));
  double * ow = malloc(nout * n * sizeof(double));
  double * oc = malloc(nout * n * sizeof(double));
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * wd;
  size_t i, e;

  test_region_points(n, rmin, rmax, t0, t1, r, theta, phi);

  /* direct method: the range is ignored */
  params.method = POTSYNTH_METHOD_DIRECT;
  params.theta_min = t0;
  params.theta_max = t1;
  wd = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params);
  ps_test(wd == NULL, "region: %s direct allocation", desc);

  if (wd)
    {
      potsynth_eval_workspace * ew = potsynth_eval_alloc(0, wd);
      double rr = rmin, tt = (t0 > 0.1) ? 0.5 * t0 : fmin(t1 + 0.1, M_PI), pp = 1.0, B[12];

      potsynth_eval(n, r, theta, phi, od, ew);
      ps_test_int(potsynth_eval(1, &rr, &tt, &pp, B, ew), POTSYNTH_SUCCESS,
                  "region: %s direct method outside of the range", desc);
      potsynth_eval_free(ew);
      potsynth_free(wd);
    }

  for (e = 0; e < 3; ++e)
    {
      potsynth_workspace * ww, * wc;
      potsynth_info iw, ic;

      params.method = POTSYNTH_METHOD_FAST;
      params.eps = epsv[e];

      params.theta_min = 0.0;
      params.theta_max = 0.0;
      ww = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params);

      params.theta_min = t0;
      params.theta_max = t1;
      wc = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params);

      ps_test(ww == NULL || wc == NULL, "region: %s eps %.0e allocation", desc, epsv[e]);

      if (ww && wc)
        {
          potsynth_eval_workspace * ew = potsynth_eval_alloc(0, ww);
          potsynth_eval_workspace * ec = potsynth_eval_alloc(0, wc);
          double emax = 0.0, dmax = 0.0, bmax = 0.0;
          size_t ndiff = 0;

          potsynth_get_info(&iw, ww);
          potsynth_get_info(&ic, wc);

          ps_test_int(potsynth_eval(n, r, theta, phi, ow, ew), POTSYNTH_SUCCESS, "region: %s whole sphere", desc);
          ps_test_int(potsynth_eval(n, r, theta, phi, oc, ec), POTSYNTH_SUCCESS,
                      "region: %s points inside and at the ends of the range", desc);

          for (i = 0; i < nout * n; ++i)
            {
              ndiff += (ow[i] != oc[i]);
              dmax = fmax(dmax, fabs(ow[i] - oc[i]));
              emax = fmax(emax, fabs(oc[i] - od[i]));
              bmax = fmax(bmax, fabs(od[i]));
            }

          ps_test(emax / bmax > 3.0 * epsv[e], "region: %s eps %.0e against direct method (%.2e)",
                  desc, epsv[e], emax / bmax);

          /* same window: the stored rows hold the same values */
          if (iw.w == ic.w)
            ps_test(ndiff != 0, "region: %s eps %.0e equal to the whole sphere (%zu values differ, %.2e)",
                    desc, epsv[e], ndiff, dmax / bmax);

          if (t0 > 0.0 || t1 < M_PI)
            ps_test(ic.nrows >= iw.nrows || ic.grid_bytes >= iw.grid_bytes,
                    "region: %s eps %.0e fewer rows (%zu of %zu)", desc, epsv[e], ic.nrows, iw.nrows);
          else
            ps_test(ic.nrows != iw.nrows, "region: %s eps %.0e all rows (%zu of %zu)", desc, epsv[e],
                    ic.nrows, iw.nrows);
          ps_test((double) ic.nrows > (t1 - t0) / M_PI * (double) iw.nrows + 2.0 * ic.w + 4.0,
                  "region: %s eps %.0e number of rows (%zu)", desc, epsv[e], ic.nrows);

          /* points outside of the range */
          {
            double rr = rmin, pp = 1.0, B[12], tt;

            test_last_error();

            if (t0 > 0.0)
              {
                tt = t0 * (1.0 - 1.0e-12);
                ps_test_int(potsynth_eval(1, &rr, &tt, &pp, B, ec), POTSYNTH_EDOM,
                            "region: %s point below the range", desc);
                ps_test_int(test_last_error(), POTSYNTH_EDOM, "region: %s error reported", desc);
              }

            if (t1 < M_PI)
              {
                tt = t1 * (1.0 + 1.0e-12);
                ps_test_int(potsynth_eval(1, &rr, &tt, &pp, B, ec), POTSYNTH_EDOM,
                            "region: %s point above the range", desc);
              }

            test_last_error();
          }

          potsynth_eval_free(ew);
          potsynth_eval_free(ec);
        }

      if (ww)
        potsynth_free(ww);
      if (wc)
        potsynth_free(wc);
    }

  free(gh);
  free(r);
  free(theta);
  free(phi);
  free(od);
  free(ow);
  free(oc);
}

/* a model with several coefficient sets on a polar cap, fast against direct */
static void
test_region_sets(void)
{
  const size_t L = 36, n = 2000, nlm = potsynth_nlm(1, L, L);
  const double t1 = 50.0 * M_PI / 180.0;
  double * coef = test_sets_coef(3, 1, L, L);
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  double * of = malloc(3 * n * sizeof(double));
  double * od = malloc(3 * n * sizeof(double));
  potsynth_term t0[1] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
  potsynth_term tt[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_DP, 0 }, { 2, -1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
  potsynth_term tp[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_PSIN, 1 }, { 2, 1.0, NULL, POTSYNTH_THETA_DP, 0 } };
  potsynth_output o[3] = { { 0, 1, t0, 0, 0.0 }, { 0, 2, tt, 0, 0.0 }, { 0, 2, tp, 0, 0.0 } };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * wf, * wd;
  size_t i;

  (void) nlm;

  test_region_points(n, 1.0, 1.0, 0.0, t1, r, theta, phi);

  params.eps = 1.0e-11;
  params.theta_max = t1;
  params.method = POTSYNTH_METHOD_FAST;
  wf = potsynth_alloc_sets(1, L, L, 3, coef, 1.0, 1.0, 1.0, NULL, 3, o, &params);
  params.method = POTSYNTH_METHOD_DIRECT;
  wd = potsynth_alloc_sets(1, L, L, 3, coef, 1.0, 1.0, 1.0, NULL, 3, o, &params);
  ps_test(wf == NULL || wd == NULL, "region: sets on a polar cap allocation");

  if (wf && wd)
    {
      potsynth_eval_workspace * ef = potsynth_eval_alloc(0, wf);
      potsynth_eval_workspace * ed = potsynth_eval_alloc(0, wd);
      double emax = 0.0, bmax = 0.0;

      potsynth_eval(n, r, theta, phi, of, ef);
      potsynth_eval(n, r, theta, phi, od, ed);

      for (i = 0; i < 3 * n; ++i)
        {
          emax = fmax(emax, fabs(of[i] - od[i]));
          bmax = fmax(bmax, fabs(od[i]));
        }

      ps_test(emax / bmax > 3.0e-11, "region: sets on a polar cap fast against direct (%.2e)", emax / bmax);

      potsynth_eval_free(ef);
      potsynth_eval_free(ed);
    }

  if (wf)
    potsynth_free(wf);
  if (wd)
    potsynth_free(wd);

  free(coef);
  free(r);
  free(theta);
  free(phi);
  free(of);
  free(od);
}

/* invalid colatitude ranges */
static void
test_region_errors(void)
{
  const size_t L = 8;
  const double a = 6371.2;
  double * gh = test_random_coef(1, L, L);
  potsynth_parameters params = potsynth_default_parameters();
  const double bad[5][2] = { { 1.0, 0.5 }, { 0.5, 0.0 }, { -0.1, 1.0 }, { 0.0, 3.2 }, { NAN, 1.0 } };
  size_t k;

  params.method = POTSYNTH_METHOD_DIRECT;
  test_last_error();

  for (k = 0; k < 5; ++k)
    {
      params.theta_min = bad[k][0];
      params.theta_max = bad[k][1];
      ps_test(potsynth_mag_alloc(1, L, L, gh, a, a, a + 100.0, POTSYNTH_MAG_B, &params) != NULL,
              "region: range [%g, %g] rejected", bad[k][0], bad[k][1]);
      ps_test_int(test_last_error(), POTSYNTH_EINVAL, "region: range [%g, %g] error code", bad[k][0], bad[k][1]);
    }

  free(gh);
}

static void
test_region(void)
{
  test_region_range("mid-latitude band", 0.6, 1.3);
  test_region_range("north polar cap", 0.0, 50.0 * M_PI / 180.0);
  test_region_range("south polar cap", M_PI - 0.4, M_PI);
  test_region_range("narrow band at the equator", 1.55, 1.59);
  test_region_range("whole sphere given explicitly", 0.0, M_PI);
  test_region_sets();
  test_region_errors();
}
