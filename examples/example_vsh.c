/*
 * example_vsh.c
 *
 * Vector fields that are not gradients of a potential: nexp fields given by
 * vector spherical harmonic coefficients at a few radii, evaluated at scattered
 * points with radial functions of the caller, and written as the rows of a
 * design matrix.
 *
 * Each field is
 *
 *   V(r,theta,phi) = sum_j f_j(r) sum_lm [ a^j_lm Y_lm rhat + b^j_lm Psi_lm + c^j_lm Phi_lm ]
 *
 * where the coefficients a^j, b^j, c^j are the values at the radial node r_j and
 * f_j is the hat function of that node (linear interpolation between the nodes).
 * Any other radial dependence that is linear in the coefficients is used in the
 * same way: only radial_linear() below changes.
 *
 * Here the coefficients are random numbers with a decaying spectrum.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <potsynth.h>
#include <potsynth_vsh.h>

typedef struct
{
  size_t n;             /* number of radial nodes */
  const double * r;     /* nodes, ascending */
} nodes_t;

/* radial functions: the hat functions of the nodes, two of which are nonzero at a radius */
static int
radial_linear(const double r, size_t * first, size_t * n, double * f, double * df, void * params)
{
  const nodes_t * nodes = (const nodes_t *) params;
  size_t k = 0;

  if (!(r >= nodes->r[0] && r <= nodes->r[nodes->n - 1]))
    return POTSYNTH_EDOM;

  while (k + 2 < nodes->n && r >= nodes->r[k + 1])
    ++k;

  {
    const double h = nodes->r[k + 1] - nodes->r[k];
    const double t = (r - nodes->r[k]) / h;

    *first = k;
    *n = 2;
    f[0] = 1.0 - t;
    f[1] = t;

    if (df != NULL)
      {
        df[0] = -1.0 / h;
        df[1] = 1.0 / h;
      }
  }

  return POTSYNTH_SUCCESS;
}

static double
urand(void)
{
  static unsigned long long s = 88172645463325252ULL;

  s ^= s << 13;
  s ^= s >> 7;
  s ^= s << 17;

  return (double) (s >> 11) / 9007199254740992.0;
}

int
main(void)
{
  const size_t lmax = 60, nexp = 20, nrad = 4, npts = 100000;
  const size_t nlm = potsynth_nlm(1, lmax, lmax);
  const double rnode[4] = { 6671.2, 6871.2, 7071.2, 7221.2 };    /* km */
  const double theta_max = 50.0 * M_PI / 180.0;                   /* points on a polar cap */
  nodes_t nodes = { 4, rnode };
  potsynth_radial radial = { 4, radial_linear, &nodes };
  potsynth_parameters params = potsynth_default_parameters();
  potsynth_workspace * w;
  potsynth_eval_workspace * ew;
  potsynth_info info;
  double * a = malloc(nexp * nrad * nlm * sizeof(double));
  double * b = malloc(nexp * nrad * nlm * sizeof(double));
  double * c = malloc(nexp * nrad * nlm * sizeof(double));
  double * r = malloc(npts * sizeof(double));
  double * theta = malloc(npts * sizeof(double));
  double * phi = malloc(npts * sizeof(double));
  double * G = malloc(npts * 3 * nexp * sizeof(double));
  size_t i, l;
  int status;

  /* coefficients of field e at radial node j: a + (e*nrad + j)*nlm, in the layout of potsynth_idx() */
  for (i = 0; i < nexp * nrad; ++i)
    {
      for (l = 1; l <= lmax; ++l)
        {
          const double sd = pow(0.9, (double) l);
          long m;

          for (m = -(long) l; m <= (long) l; ++m)
            {
              const size_t k = i * nlm + potsynth_idx((long) l, m, 1L, (long) lmax);

              a[k] = sd * (2.0 * urand() - 1.0);
              b[k] = sd * (2.0 * urand() - 1.0);
              c[k] = sd * (2.0 * urand() - 1.0);
            }
        }
    }

  /* scattered points on the cap, between the first and the last node */
  for (i = 0; i < npts; ++i)
    {
      r[i] = rnode[0] + (rnode[3] - rnode[0]) * urand();
      theta[i] = acos(1.0 - (1.0 - cos(theta_max)) * urand());
      phi[i] = 2.0 * M_PI * urand();
    }

  /* set-up: all fields on one grid, stored only for the colatitudes of the points */
  params.eps = 1.0e-8;
  params.theta_max = theta_max;

  w = potsynth_vsh_alloc(1, lmax, lmax, nexp, a, b, c, POTSYNTH_VSH_NORM_SCHMIDT, rnode[0], rnode[3], &radial,
                         POTSYNTH_VSH_V, &params);
  if (w == NULL)
    return 1;

  ew = potsynth_eval_alloc(0, w);
  potsynth_get_info(&info, w);

  /* rows of the design matrix: (V_r, V_theta, V_phi) of field e at point i in
     G[(i*nexp + e)*3 + 0..2] */
  status = potsynth_eval_rows(npts, r, theta, phi, G, ew);
  if (status)
    {
      fprintf(stderr, "potsynth_eval_rows: %s\n", potsynth_strerror(status));
      return 1;
    }

  printf("%zu fields of degree %zu with %zu radial functions: %zu outputs, %zu components\n",
         nexp, lmax, nrad, info.nout, info.ncomp);
  printf("method %s, window width %d, grid %.1f MB (%zu rows of %zu)\n",
         (info.method == POTSYNTH_METHOD_FAST) ? "fast" : "direct", info.w,
         info.grid_bytes / 1048576.0, info.nrows, info.nf1 / 2);
  printf("field 0 at the first point (r = %.1f km, theta = %.2f deg): V_r %.6f, V_theta %.6f, V_phi %.6f\n",
         r[0], theta[0] * 180.0 / M_PI, G[0], G[1], G[2]);

  /* a model of the fast method can be written and read instead of being set up again */
  status = potsynth_write("example_vsh_model.bin", w);
  if (status == POTSYNTH_SUCCESS)
    {
      potsynth_workspace * w2 = potsynth_read("example_vsh_model.bin", &radial, 0);

      if (w2 != NULL)
        {
          potsynth_eval_workspace * ew2 = potsynth_eval_alloc(0, w2);
          double V[3 * 20];

          potsynth_eval_rows(1, r, theta, phi, V, ew2);
          printf("the same from the model file: V_r %.6f, V_theta %.6f, V_phi %.6f\n", V[0], V[1], V[2]);
          potsynth_eval_free(ew2);
          potsynth_free(w2);
        }

      remove("example_vsh_model.bin");
    }

  potsynth_eval_free(ew);
  potsynth_free(w);
  free(a);
  free(b);
  free(c);
  free(r);
  free(theta);
  free(phi);
  free(G);

  return 0;
}
