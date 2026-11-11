/*
 * test_window.c
 *
 * Window width per quantity (potsynth_window_choose_shell): the horizontal gravity
 * gradients, whose power lies near the band edge, and a shell reaching inside r = a
 * (u > 1) must meet eps, while the magnetic field at satellite altitude keeps the
 * width of the measured table
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

/* largest error of each irreducible gravity quantity relative to its largest magnitude */
static double
test_window_grav_error(const size_t n, const double * f, const double * d)
{
  double emax[6] = { 0 }, vmax[6] = { 0 }, worst = 0.0;
  size_t i, q;

  for (i = 0; i < n; ++i)
    {
#define F(c) f[(c) * n + i]
#define D(c) d[(c) * n + i]
      const double e[6] = { fabs(F(0) - D(0)), fabs(F(1) - D(1)), hypot(F(2) - D(2), F(3) - D(3)),
                            fabs(F(4) - D(4)), hypot(F(5) - D(5), F(6) - D(6)),
                            hypot((F(7) - F(9)) - (D(7) - D(9)), 2.0 * (F(8) - D(8))) };
      const double v[6] = { fabs(D(0)), fabs(D(1)), hypot(D(2), D(3)), fabs(D(4)), hypot(D(5), D(6)),
                            hypot(D(7) - D(9), 2.0 * D(8)) };
#undef F
#undef D
      for (q = 0; q < 6; ++q)
        {
          emax[q] = fmax(emax[q], e[q]);
          vmax[q] = fmax(vmax[q], v[q]);
        }
    }

  for (q = 0; q < 6; ++q)
    worst = fmax(worst, emax[q] / vmax[q]);

  return worst;
}

static void
test_window_shell(void)
{
  const size_t L = 360, n = 600;
  const double GM = 3.986004415e14, a = 6378136.3, rmax = a + 300.0e3;
  const double rmin_list[] = { a, a - 21.0e3 };   /* from r = a, and from inside it (u > 1) */
  const double eps_list[] = { 1.0e-6, 1.0e-10 };
  const int quantities = POTSYNTH_GRAV_V | POTSYNTH_GRAV_G | POTSYNTH_GRAV_T;
  const size_t nchan = potsynth_grav_nout(quantities);
  double * coef = calloc(potsynth_nlm(0, L, L), sizeof(double));
  double * r = malloc(n * sizeof(double)), * th = malloc(n * sizeof(double)), * ph = malloc(n * sizeof(double));
  double * out = malloc(nchan * n * sizeof(double)), * ref = malloc(nchan * n * sizeof(double));
  size_t l, m, i, s, k;

  /* Kaula-like spectrum, as in examples/example_grav.c */
  for (l = 2; l <= L; ++l)
    {
      const double sd = sqrt(3.0) * 1.0e-5 / ((double) l * l);

      for (m = 0; m <= l; ++m)
        {
          coef[potsynth_idx((long) l, (long) m, 0L, (long) L)] = sd * (2.0 * test_urand() - 1.0);
          if (m > 0)
            coef[potsynth_idx((long) l, -(long) m, 0L, (long) L)] = sd * (2.0 * test_urand() - 1.0);
        }
    }
  coef[potsynth_idx(0L, 0L, 0L, (long) L)] = 1.0;
  coef[potsynth_idx(2L, 0L, 0L, (long) L)] = -0.484165143790e-03;

  for (s = 0; s < 2; ++s)
    {
      const double rmin = rmin_list[s];
      potsynth_parameters p = potsynth_default_parameters();
      potsynth_workspace * wd;
      potsynth_eval_workspace * ewd;

      /* half the points uniform in the shell, half on its bottom */
      for (i = 0; i < n; ++i)
        {
          r[i] = (i < n / 2) ? rmin + (rmax - rmin) * test_urand() : rmin;
          th[i] = acos(2.0 * test_urand() - 1.0);
          ph[i] = 2.0 * M_PI * test_urand();
        }

      p.method = POTSYNTH_METHOD_DIRECT;
      wd = potsynth_grav_alloc(0, L, L, coef, GM, a, rmin, rmax, quantities, &p);
      ewd = potsynth_eval_alloc(0, wd);
      potsynth_eval(n, r, th, ph, ref, ewd);

      for (k = 0; k < 2; ++k)
        {
          potsynth_workspace * w;
          potsynth_eval_workspace * ew;
          double worst;

          p = potsynth_default_parameters();
          p.method = POTSYNTH_METHOD_FAST;
          p.eps = eps_list[k];
          w = potsynth_grav_alloc(0, L, L, coef, GM, a, rmin, rmax, quantities, &p);
          ew = potsynth_eval_alloc(0, w);
          potsynth_eval(n, r, th, ph, out, ew);

          worst = test_window_grav_error(n, out, ref);
          ps_test(!(worst <= eps_list[k]),
                  "window shell: gravity tensor, rmin = a%+.0f km, eps = %.0e, worst error %.2e",
                  (rmin - a) / 1.0e3, eps_list[k], worst);

          potsynth_eval_free(ew);
          potsynth_free(w);
        }

      potsynth_eval_free(ewd);
      potsynth_free(wd);
    }

  /* magnetic field at satellite altitude: the width of the measured table */
  {
    const size_t Lm = 185;
    double * gh = malloc(potsynth_nlm(1, Lm, Lm) * sizeof(double));

    for (l = 1; l <= Lm; ++l)
      {
        const double sd = sqrt(3.0 * 0.5 / ((l + 1.0) * (2.0 * l + 1.0)));

        for (m = 0; m <= l; ++m)
          {
            gh[potsynth_idx((long) l, (long) m, 1L, (long) Lm)] = sd * (2.0 * test_urand() - 1.0);
            if (m > 0)
              gh[potsynth_idx((long) l, -(long) m, 1L, (long) Lm)] = sd * (2.0 * test_urand() - 1.0);
          }
      }

    for (k = 0; k < 2; ++k)
      {
        potsynth_parameters p = potsynth_default_parameters();
        potsynth_workspace * w;
        potsynth_info info;

        p.method = POTSYNTH_METHOD_FAST;
        p.eps = eps_list[k];
        w = potsynth_mag_alloc(1, Lm, Lm, gh, 6371.2, 6371.2 + 250.0, 6371.2 + 800.0,
                               POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &p);
        potsynth_get_info(&info, w);
        ps_test_int(info.w, potsynth_window_choose(eps_list[k], 2.0),
                    "window shell: magnetic field at 250-800 km keeps the table width, eps = %.0e",
                    eps_list[k]);
        potsynth_free(w);
      }

    free(gh);
  }

  free(coef); free(r); free(th); free(ph); free(out); free(ref);
}
