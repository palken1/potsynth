/*
 * test_robust.c
 *
 * Parameter checks, damaged model files, radial functions of the caller that
 * fail after validation, the achieved errors in potsynth_info, and coefficient
 * sets that nearly cancel
 */

/* radial functions that fail after a given number of calls (barycentric otherwise) */
typedef struct
{
  const test_radial_nodes * nodes;
  long calls;           /* calls left before failing; < 0: never fail */
} test_robust_failing;

static int
test_robust_failing_fn(const double r, size_t * first, size_t * n, double * f, double * df, void * params)
{
  test_robust_failing * p = (test_robust_failing *) params;

  if (p->calls == 0)
    return POTSYNTH_EDOM;
  if (p->calls > 0)
    --(p->calls);

  return test_radial_bary(r, first, n, f, df, (void *) p->nodes);
}

/* overwrite the 8 bytes at offset off of a file */
static void
test_robust_patch(const char * name, const long off, const void * val)
{
  FILE * fp = fopen(name, "r+b");

  if (fp != NULL)
    {
      fseek(fp, off, SEEK_SET);
      fwrite(val, 8, 1, fp);
      fclose(fp);
    }
}

static void
test_robust(void)
{
  const size_t L = 24;
  const double a = 6371.2, rmin = a + 300.0, rmax = a + 700.0;
  double * gh = test_random_coef(1, L, L);
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;
  potsynth_info info;

  params.method = POTSYNTH_METHOD_FAST;

  /* the former name of potsynth_functional still compiles */
  {
    potsynth_channel old = { 0, NULL, POTSYNTH_THETA_P, 0 };
    potsynth_functional * f = &old;
    ps_test(f->theta != POTSYNTH_THETA_P, "robust: potsynth_channel is potsynth_functional");
  }

  /* the free functions accept NULL, as free() does */
  potsynth_free(NULL);
  potsynth_eval_free(NULL);
  potsynth_shc_free(NULL);
  potsynth_gfc_free(NULL);
  ps_test(0, "robust: free functions accept NULL");

  /* sigma below 1 is rejected (it overflowed the grid); 1 is accepted */
  test_last_error();
  params.sigma = 0.6;
  ps_test(potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params) != NULL, "robust: sigma 0.6 rejected");
  ps_test_int(test_last_error(), POTSYNTH_EINVAL, "robust: sigma 0.6 error code");
  params.sigma = 0.9;
  ps_test(potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params) != NULL, "robust: sigma 0.9 rejected");
  test_last_error();
  params.sigma = 1.0;
  w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
  ps_test(w == NULL, "robust: sigma 1 accepted");
  potsynth_free(w);
  params.sigma = 0.0;

  /* just past the poles (within the colatitude tolerance): fast and direct agree, and
   * the values continue smoothly from those at the pole */
  {
    potsynth_workspace * wf, * wdir;
    const double tq[6] = { 0.0, -5.0e-7, -9.9e-7, M_PI, M_PI + 5.0e-7, (double) (float) M_PI };
    double rq[6], pq[6], of[18], od[18], big = 0.0, dfd = 0.0, dpole = 0.0;
    size_t jj;

    for (jj = 0; jj < 6; ++jj)
      {
        rq[jj] = rmin + 100.0;
        pq[jj] = 1.1;
      }
    params.eps = 1e-10;
    wf = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
    params.method = POTSYNTH_METHOD_DIRECT;
    wdir = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
    params.method = POTSYNTH_METHOD_FAST;
    if (wf && wdir)
      {
        potsynth_eval_workspace * ef = potsynth_eval_alloc(0, wf), * ed = potsynth_eval_alloc(0, wdir);

        ps_test_int(potsynth_eval(6, rq, tq, pq, of, ef), POTSYNTH_SUCCESS, "robust: fast just past the poles");
        ps_test_int(potsynth_eval(6, rq, tq, pq, od, ed), POTSYNTH_SUCCESS, "robust: direct just past the poles");
        for (jj = 0; jj < 18; ++jj)
          big = fmax(big, fabs(od[jj]));
        for (jj = 0; jj < 18; ++jj)
          dfd = fmax(dfd, fabs(of[jj] - od[jj]));
        for (jj = 0; jj < 3; ++jj)
          {
            /* component c of point k is out[c*6 + k]: compare k = 1, 2 with 0 and 4, 5 with 3 */
            size_t k;
            for (k = 1; k < 6; ++k)
              if (k != 3)
                dpole = fmax(dpole, fabs(od[jj * 6 + k] - od[jj * 6 + ((k < 3) ? 0 : 3)]));
          }
        ps_test(!(dfd <= 3.0 * params.eps * big), "robust: fast and direct agree past the poles (%.1e)", dfd / big);
        ps_test(!(dpole <= 1.0e-4 * big), "robust: values past the poles close to those at the poles (%.1e)",
                dpole / big);
        potsynth_eval_free(ef);
        potsynth_eval_free(ed);
      }
    potsynth_free(wf);
    potsynth_free(wdir);
  }

  /* direct method: a point so far inside a that (a/r)^lmax overflows is refused
   * (it gave NaN); one where it does not is evaluated */
  {
    const size_t L2 = 300;
    double * c3 = test_random_coef(1, L2, L2);
    potsynth_workspace * wd;

    params.method = POTSYNTH_METHOD_DIRECT;
    wd = potsynth_mag_alloc(1, L2, L2, c3, a, 0.0, 0.0, POTSYNTH_MAG_B, &params);
    params.method = POTSYNTH_METHOD_FAST;
    ps_test(wd == NULL, "robust: direct model at L = 300");
    if (wd)
      {
        potsynth_eval_workspace * ed = potsynth_eval_alloc(0, wd);
        double rr = a * exp(-3.0), tt = 1.0, pp = 0.5, o3[3];

        test_last_error();
        ps_test_int(potsynth_eval(1, &rr, &tt, &pp, o3, ed), POTSYNTH_EDOM,
                    "robust: direct refuses (a/r)^300 = e^900");
        test_last_error();
        rr = 0.5 * a;
        ps_test_int(potsynth_eval(1, &rr, &tt, &pp, o3, ed), POTSYNTH_SUCCESS, "robust: direct at r = a/2");
        ps_test(!(isfinite(o3[0]) && isfinite(o3[1]) && isfinite(o3[2])), "robust: direct at r = a/2 finite");
        potsynth_eval_free(ed);
        potsynth_free(wd);
      }
    free(c3);
  }

  /* achieved errors: within the tolerances for the automatic choice, -1 without a check */
  params.eps = 1e-10;
  w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params);
  ps_test(w == NULL, "robust: allocation");
  if (w)
    {
      potsynth_get_info(&info, w);
      ps_test(!(info.err_radial > 0.0 && info.err_radial <= 0.5 * params.eps), "robust: err_radial %g within eps/2",
              info.err_radial);
      ps_test(!(info.err_check > 0.0 && info.err_check <= params.eps), "robust: err_check %g within eps",
              info.err_check);
      potsynth_free(w);
    }

  params.w = 9;
  w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
  if (w)
    {
      potsynth_get_info(&info, w);
      ps_test(info.err_check != -1.0, "robust: err_check -1 for a given width");
      potsynth_free(w);
    }
  params.w = 0;

  params.method = POTSYNTH_METHOD_DIRECT;
  w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
  if (w)
    {
      potsynth_get_info(&info, w);
      ps_test(info.err_radial != 0.0 || info.err_check != 0.0, "robust: error fields 0 for the direct method");
      potsynth_free(w);
    }
  params.method = POTSYNTH_METHOD_FAST;

  /* a model file with damaged sizes is refused: window wider than the stored rows, and
   * a grid size that overflows. Header: 8 bytes magic, 8 byte order mark, 17 unsigned
   * (nrows is the 11th, ncols the 12th), then 4 signed (row_lo, row_hi, col_lo, col_hi) */
  {
    char name[1024];
    const long hu = 16, hi = 16 + 17 * 8;

    test_io_name(name, sizeof(name));
    w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
    if (w)
      {
        uint64_t u;
        int64_t lo, hi1;
        FILE * fp;
        potsynth_workspace * w2;

        potsynth_get_info(&info, w);
        ps_test_int(potsynth_write(name, w), POTSYNTH_SUCCESS, "robust: write");

        fp = fopen(name, "rb");
        fseek(fp, hi, SEEK_SET);
        if (fread(&lo, 8, 1, fp) != 1)
          lo = 0;
        fclose(fp);

        /* fewer rows than the window width (row_hi kept consistent) */
        u = (uint64_t) info.w - 1;
        hi1 = lo + (int64_t) u - 1;
        test_robust_patch(name, hu + 10 * 8, &u);
        test_robust_patch(name, hi + 8, &hi1);
        test_last_error();
        w2 = potsynth_read(name, NULL, 1);
        ps_test(w2 != NULL, "robust: file with fewer rows than the window refused");
        ps_test_int(test_last_error(), POTSYNTH_EINVAL, "robust: damaged file error code");
        if (w2)
          potsynth_free(w2);

        /* nrows = ncols = 2^32: the grid size overflows */
        ps_test_int(potsynth_write(name, w), POTSYNTH_SUCCESS, "robust: write again");
        u = (uint64_t) 1 << 32;
        hi1 = lo + (int64_t) u - 1;
        test_robust_patch(name, hu + 10 * 8, &u);
        test_robust_patch(name, hu + 11 * 8, &u);
        test_robust_patch(name, hi + 8, &hi1);
        fp = fopen(name, "rb");
        fseek(fp, hi + 16, SEEK_SET);
        if (fread(&lo, 8, 1, fp) != 1)
          lo = 0;
        fclose(fp);
        hi1 = lo + (int64_t) u - 1;
        test_robust_patch(name, hi + 24, &hi1);
        w2 = potsynth_read(name, NULL, 1);
        ps_test(w2 != NULL, "robust: file with an overflowing grid size refused");
        test_last_error();

        remove(name);
        potsynth_free(w);
      }
  }

  /* radial functions of the caller that fail after the points were validated: NaN for
   * those points and POTSYNTH_EDOM, both methods; points before the failure keep their
   * values */
  {
    const double rn[4] = { 6671.2, 6800.0, 6950.0, 7071.2 };
    double wb[4];
    test_radial_nodes nodes = { 4, rn, wb };
    test_robust_failing fail = { &nodes, -1 };
    const potsynth_radial rad = { 4, test_robust_failing_fn, &fail };
    double * abc = test_sets_coef(4, 1, L, L);
    const size_t n = 40;
    double r[40], th[40], ph[40], out[40 * 3];
    int meth;

    test_bary_weights(4, rn, wb);
    test_random_points(n, rn[0], rn[3], r, th, ph);

    for (meth = 0; meth < 2; ++meth)
      {
        potsynth_eval_workspace * ew;
        const char * desc = meth ? "direct" : "fast";
        size_t i, nnan = 0;

        params.method = meth ? POTSYNTH_METHOD_DIRECT : POTSYNTH_METHOD_FAST;
        params.nthreads = 1;
        fail.calls = -1;
        w = potsynth_vsh_alloc(1, L, L, 1, abc, NULL, NULL, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V,
                               rn[0], rn[3], &rad, &params);
        ps_test(w == NULL, "robust: %s model with radial functions", desc);
        if (w == NULL)
          continue;

        ew = potsynth_eval_alloc(0, w);

        /* validation calls the functions once per point; fail on the second pass */
        fail.calls = (long) n + 5;
        test_last_error();
        ps_test_int(potsynth_eval(n, r, th, ph, out, ew), POTSYNTH_EDOM,
                    "robust: %s radial function failing after validation", desc);
        test_last_error();
        for (i = 0; i < 3 * n; ++i)
          nnan += isnan(out[i]);
        ps_test(nnan == 0 || nnan == 3 * n, "robust: %s NaN at the failed points only (%zu NaN)", desc, nnan);

        /* and a following call that does not fail succeeds */
        fail.calls = -1;
        ps_test_int(potsynth_eval(n, r, th, ph, out, ew), POTSYNTH_SUCCESS, "robust: %s next call", desc);

        potsynth_eval_free(ew);
        potsynth_free(w);
      }

    params.method = POTSYNTH_METHOD_FAST;
    params.nthreads = 0;
    free(abc);
  }

  /* sets that nearly cancel: the output set1 - set0 with set1 = set0 + d gets the
   * radial basis of d alone (a model of the one set d), not of the sets */
  {
    const size_t nlm = potsynth_nlm(1, L, L);
    double * c2 = malloc(2 * nlm * sizeof(double));
    double * d = malloc(nlm * sizeof(double));
    potsynth_term t2[2] = { { 1, 1.0, NULL, POTSYNTH_THETA_P, 0 }, { 0, -1.0, NULL, POTSYNTH_THETA_P, 0 } };
    potsynth_term t1[1] = { { 0, 1.0, NULL, POTSYNTH_THETA_P, 0 } };
    potsynth_output o2[1] = { { 1, 2, t2, 0, 0.0 } };
    potsynth_output o1[1] = { { 1, 1, t1, 0, 0.0 } };
    potsynth_workspace * wd, * w1;
    size_t k, l;

    /* steep common part, flat difference */
    for (l = 1; l <= L; ++l)
      {
        long m;
        for (m = -(long) l; m <= (long) l; ++m)
          {
            k = potsynth_idx(l, m, 1, L);
            c2[k] = 1.0e3 * pow((double) l, -4.0) * (1.0 + 0.1 * cos((double) k));
            d[k] = 1.0e-4 * (0.5 + 0.5 * sin(3.7 * (double) k)); /* flat spectrum */
            c2[nlm + k] = c2[k] + d[k];
            d[k] = c2[nlm + k] - c2[k]; /* the difference as the sets hold it */
          }
      }

    params.eps = 1e-8;
    wd = potsynth_alloc_sets(1, L, L, 2, c2, a, rmin, rmax, NULL, 1, o2, &params);
    w1 = potsynth_alloc_sets(1, L, L, 1, d, a, rmin, rmax, NULL, 1, o1, &params);
    ps_test(wd == NULL || w1 == NULL, "robust: cancelling sets allocation");
    if (wd && w1)
      {
        potsynth_info i1;

        potsynth_get_info(&info, wd);
        potsynth_get_info(&i1, w1);
        ps_test(info.rho != i1.rho, "robust: cancelling sets get the rank of their difference (%zu, %zu)",
                info.rho, i1.rho);
      }
    potsynth_free(wd);
    potsynth_free(w1);
    free(c2);
    free(d);
  }

  free(gh);
}
