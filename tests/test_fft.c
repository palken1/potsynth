/*
 * test_fft.c
 */

/*
test_fft()
  Test the FFT conventions relied on by the set-up (fast.c, potsynth_fft.h):

forward:  F_k = sum_j c_j e^{-2 pi i jk/N}
backward: c_j = sum_k F_k e^{+2 pi i jk/N}
real backward, half-complex input r0, r1, i1, r2, i2, ...:
  x_k = r0 + 2 sum_m (r_m cos - i_m sin)(2 pi mk/N) + r_{N/2} (-1)^k

with several transforms per call, as the set-up does: complex transforms one after
the other, and real transforms whose elements are nb apart (the functionals of a
grid row), in place
*/

static void
test_fft(void)
{
  const size_t ns[] = { 12, 30, 64, 375, 750 };
  const size_t nb = 3; /* transforms per call */
  size_t in;

  for (in = 0; in < sizeof(ns) / sizeof(ns[0]); ++in)
    {
      const size_t N = ns[in];
      double * c = malloc(2 * N * nb * sizeof(double));
      double * c0 = malloc(2 * N * nb * sizeof(double));
      double * r = malloc(N * nb * sizeof(double));
      double * r0 = malloc(N * nb * sizeof(double));
      double err_fwd = 0.0, err_inv = 0.0, err_real = 0.0;
      size_t j, k, m, t;

      for (j = 0; j < 2 * N * nb; ++j)
        c[j] = c0[j] = 2.0 * test_urand() - 1.0;

      ps_test_int(potsynth_fft_c2c(N, nb, N, 1, c, 1.0), 0, "fft: forward N=%zu call", N);
      for (t = 0; t < nb; ++t)
        {
          const double * ct = c + 2 * t * N, * c0t = c0 + 2 * t * N;

          for (k = 0; k < N; ++k)
            {
              long double re = 0.0L, im = 0.0L;

              for (j = 0; j < N; ++j)
                {
                  const long double a = -2.0L * M_PI * (long double) ((j * k) % N) / N;

                  re += c0t[2 * j] * cosl(a) - c0t[2 * j + 1] * sinl(a);
                  im += c0t[2 * j] * sinl(a) + c0t[2 * j + 1] * cosl(a);
                }

              err_fwd = fmax(err_fwd, (double) fmaxl(fabsl(ct[2 * k] - re), fabsl(ct[2 * k + 1] - im)) / N);
            }
        }

      ps_test_int(potsynth_fft_c2c(N, nb, N, 0, c, 1.0 / N), 0, "fft: backward N=%zu call", N);
      for (j = 0; j < 2 * N * nb; ++j)
        err_inv = fmax(err_inv, fabs(c[j] - c0[j]));

      /* real: element k of transform t at r[k*nb + t] */
      for (j = 0; j < N * nb; ++j)
        r[j] = r0[j] = 2.0 * test_urand() - 1.0;

      ps_test_int(potsynth_fft_r_backward(N, nb, (ptrdiff_t) nb, 1, r, 1.0), 0,
                  "fft: real backward N=%zu call", N);
      for (t = 0; t < nb; ++t)
        {
          for (k = 0; k < N; ++k)
            {
              long double sum = r0[t];

              for (m = 1; 2 * m < N; ++m)
                {
                  const long double a = 2.0L * M_PI * (long double) ((m * k) % N) / N;
                  sum += 2.0L * (r0[(2 * m - 1) * nb + t] * cosl(a) - r0[2 * m * nb + t] * sinl(a));
                }

              if (N % 2 == 0)
                sum += r0[(N - 1) * nb + t] * ((k & 1) ? -1.0L : 1.0L);

              err_real = fmax(err_real, (double) fabsl(r[k * nb + t] - sum) / N);
            }
        }

      ps_test(err_fwd > 1.0e-14, "fft: forward N=%zu (error %.2e)", N, err_fwd);
      ps_test(err_inv > 1.0e-14, "fft: backward N=%zu (error %.2e)", N, err_inv);
      ps_test(err_real > 1.0e-14, "fft: real backward N=%zu (error %.2e)", N, err_real);

      free(c);
      free(c0);
      free(r);
      free(r0);
    }
}
