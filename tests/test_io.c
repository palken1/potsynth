/*
 * test_io.c
 *
 * Model files: a model read from a file gives the same values as the model
 * that was written
 */

/* name of a scratch file */
static void
test_io_name(char * name, const size_t len)
{
  const char * dir = getenv("TMPDIR");

  snprintf(name, len, "%s/potsynth_test_model.bin", (dir != NULL && dir[0] != '\0') ? dir : ".");
}

/* write w, read it and compare the outputs at n points; frees w */
static void
test_io_roundtrip(const char * desc, potsynth_workspace * w, const potsynth_radial * radial,
                  const double rmin, const double rmax, const double t0, const double t1)
{
  const size_t n = 1200;
  char name[1024];
  double * r = malloc(n * sizeof(double));
  double * theta = malloc(n * sizeof(double));
  double * phi = malloc(n * sizeof(double));
  potsynth_info i1, i2;
  potsynth_workspace * w2;
  size_t i;

  ps_test(w == NULL, "io: %s allocation", desc);
  if (w == NULL)
    {
      free(r);
      free(theta);
      free(phi);
      return;
    }

  test_io_name(name, sizeof(name));
  if (t0 > 0.0 || t1 < M_PI)
    test_region_points(n, rmin, rmax, t0, t1, r, theta, phi);
  else
    test_random_points(n, rmin, rmax, r, theta, phi);
  potsynth_get_info(&i1, w);

  ps_test_int(potsynth_write(name, w), POTSYNTH_SUCCESS, "io: %s write", desc);
  w2 = potsynth_read(name, radial, 1);
  ps_test(w2 == NULL, "io: %s read", desc);

  if (w2)
    {
      potsynth_eval_workspace * e1 = potsynth_eval_alloc(0, w);
      potsynth_eval_workspace * e2 = potsynth_eval_alloc(0, w2);
      double * o1 = malloc(i1.nout * n * sizeof(double));
      double * o2 = malloc(i1.nout * n * sizeof(double));
      size_t ndiff = 0;

      potsynth_get_info(&i2, w2);
      ps_test(i1.lmax != i2.lmax || i1.nfunc != i2.nfunc || i1.nout != i2.nout || i1.w != i2.w ||
              i1.rho != i2.rho || i1.nrows != i2.nrows || i1.grid_bytes != i2.grid_bytes ||
              i1.ncomp != i2.ncomp || i2.method != POTSYNTH_METHOD_FAST, "io: %s parameters of the model", desc);

      ps_test_int(potsynth_eval(n, r, theta, phi, o1, e1), POTSYNTH_SUCCESS, "io: %s evaluation", desc);
      ps_test_int(potsynth_eval(n, r, theta, phi, o2, e2), POTSYNTH_SUCCESS, "io: %s evaluation after reading", desc);

      for (i = 0; i < i1.nout * n; ++i)
        ndiff += (o1[i] != o2[i]);

      ps_test(ndiff != 0, "io: %s same values after reading (%zu differ)", desc, ndiff);

      /* a point outside of the shell is still refused */
      r[3] = 0.5 * rmin;
      test_last_error();
      ps_test_int(potsynth_eval(n, r, theta, phi, o2, e2), POTSYNTH_EDOM, "io: %s domain after reading", desc);
      test_last_error();

      potsynth_eval_free(e1);
      potsynth_eval_free(e2);
      potsynth_free(w2);
      free(o1);
      free(o2);
    }

  remove(name);
  potsynth_free(w);
  free(r);
  free(theta);
  free(phi);
}

static void
test_io(void)
{
  const size_t L = 30, nlm = potsynth_nlm(1, L, L);
  const double a = 6371.2, rmin = a + 300.0, rmax = a + 700.0;
  const double rn[4] = { 6671.2, 6800.0, 6950.0, 7071.2 };
  double wb[4];
  const test_radial_nodes nodes = { 4, rn, wb };
  const potsynth_radial rad = { 4, test_radial_bary, (void *) &nodes };
  const potsynth_radial rad3 = { 3, test_radial_bary, (void *) &nodes };
  double * gh = test_random_coef(1, L, L);
  double * cf = test_random_coef(0, L, L);
  double * abc = test_sets_coef(2 * 4, 1, L, L);
  potsynth_parameters params = potsynth_default_parameters();
  char name[1024];

  (void) nlm;
  test_bary_weights(4, rn, wb);
  params.method = POTSYNTH_METHOD_FAST;

  test_io_roundtrip("magnetic field and gradient",
                    potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params),
                    NULL, rmin, rmax, 0.0, M_PI);
  test_io_roundtrip("gravity",
                    potsynth_grav_alloc(0, L, L, cf, 3.986e5, a, rmin, rmax, POTSYNTH_GRAV_V | POTSYNTH_GRAV_G, &params),
                    NULL, rmin, rmax, 0.0, M_PI);

  params.theta_min = 0.4;
  params.theta_max = 2.9;
  test_io_roundtrip("vector fields with radial functions",
                    potsynth_vsh_alloc(1, L, L, 2, abc, abc, NULL, POTSYNTH_VSH_NORM_SCHMIDT,
                                       POTSYNTH_VSH_V | POTSYNTH_VSH_DIV, rn[0], rn[3], &rad, &params),
                    &rad, rn[0], rn[3], 0.4, 2.9);
  params.theta_min = 0.0;
  params.theta_max = 0.0;

  /* errors */
  test_io_name(name, sizeof(name));
  test_last_error();

  {
    potsynth_workspace * w;

    params.method = POTSYNTH_METHOD_DIRECT;
    w = potsynth_mag_alloc(1, L, L, gh, a, rmin, rmax, POTSYNTH_MAG_B, &params);
    ps_test_int(potsynth_write(name, w), POTSYNTH_EUNSUP, "io: model of the direct method not written");
    potsynth_free(w);

    ps_test(potsynth_read("/nonexistent/potsynth_model.bin", NULL, 0) != NULL, "io: missing file");
    ps_test_int(test_last_error(), POTSYNTH_EINVAL, "io: missing file error code");

    /* radial functions that do not match; a truncated file; a file of another kind */
    params.method = POTSYNTH_METHOD_FAST;
    w = potsynth_vsh_alloc(1, L, L, 1, abc, NULL, NULL, POTSYNTH_VSH_NORM_FULL, POTSYNTH_VSH_V,
                           rn[0], rn[3], &rad, &params);
    ps_test(w == NULL, "io: allocation");
    if (w)
      {
        FILE * fp;
        long size;

        ps_test_int(potsynth_write(name, w), POTSYNTH_SUCCESS, "io: write");
        ps_test(potsynth_read(name, NULL, 0) != NULL, "io: radial functions missing");
        ps_test(potsynth_read(name, &rad3, 0) != NULL, "io: wrong number of radial functions");

        fp = fopen(name, "rb");
        fseek(fp, 0, SEEK_END);
        size = ftell(fp);
        fclose(fp);

        fp = fopen(name, "r+b");
        if (fp != NULL && size > 4096)
          {
            /* cut the grid short by rewriting the file with its first part */
            char * buf = malloc(4096);
            size_t got = fread(buf, 1, 4096, fp);

            fclose(fp);
            fp = fopen(name, "wb");
            fwrite(buf, 1, got, fp);
            fclose(fp);
            free(buf);
            ps_test(potsynth_read(name, &rad, 0) != NULL, "io: truncated file");
          }

        fp = fopen(name, "wb");
        fwrite("not a model file of potsynth", 1, 28, fp);
        fclose(fp);
        ps_test(potsynth_read(name, &rad, 0) != NULL, "io: file of another kind");

        potsynth_free(w);
      }

    test_last_error();
    remove(name);
  }

  free(gh);
  free(cf);
  free(abc);
}
