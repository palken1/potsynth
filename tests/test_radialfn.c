/*
 * test_radialfn.c
 *
 * Radial functions for the tests of models with a radial basis of the caller
 * (potsynth_radial): interpolation between radial nodes, so that the coefficient
 * sets of a model are the values of its fields at the nodes (barycentric and
 * piecewise linear interpolation) or their values and radial derivatives
 * (piecewise cubic Hermite interpolation). They are also examples of how such
 * functions are written: the library itself has none, since the radial
 * dependence of a field is the caller's.
 */

/* radial nodes: the parameters of the functions below */
typedef struct
{
  size_t n;             /* number of nodes */
  const double * r;     /* nodes, length n (linear, Hermite: ascending) */
  const double * w;     /* barycentric weights, length n (test_bary_weights()); not used by the
                           linear and Hermite functions */
} test_radial_nodes;

static int test_radial_interval(const double r, const test_radial_nodes * nodes, size_t * k);

/*
test_bary_weights()
  Barycentric weights of n distinct nodes, w_j = 1 / prod_{k != j} (r_j - r_k),
scaled by the capacity (r_max - r_min)/4 so that they do not overflow

Inputs: n - number of nodes
        r - nodes, length n
        w - (output) weights, length n

Return: success or error
*/

static int
test_bary_weights(const size_t n, const double * r, double * w)
{
  double rmin, rmax, cap;
  size_t j, k;

  if (n == 0 || r == NULL || w == NULL)
    {
      return POTSYNTH_EINVAL;
    }

  rmin = rmax = r[0];
  for (j = 1; j < n; ++j)
    {
      rmin = fmin(rmin, r[j]);
      rmax = fmax(rmax, r[j]);
    }

  cap = (n > 1) ? 0.25 * (rmax - rmin) : 1.0;

  for (j = 0; j < n; ++j)
    {
      double p = 1.0;

      for (k = 0; k < n; ++k)
        {
          if (k == j)
            continue;

          if (r[j] == r[k])
            {
              return POTSYNTH_EINVAL;
            }

          p *= (r[j] - r[k]) / cap;
        }

      w[j] = 1.0 / p;
    }

  return POTSYNTH_SUCCESS;
}

/*
test_radial_bary()
  Radial functions of barycentric (polynomial) interpolation between the nodes:
function j is the Lagrange polynomial of node j, so the coefficient sets are
the values of the field at the nodes. params: test_radial_nodes with weights.
The radius must lie between the smallest and the largest node
*/

static int
test_radial_bary(const double r, size_t * first, size_t * n, double * f, double * df,
                         void * params)
{
  const test_radial_nodes * nodes = (const test_radial_nodes *) params;
  const size_t nn = nodes->n;
  double rmin = nodes->r[0], rmax = nodes->r[0], sum = 0.0, sum2 = 0.0;
  size_t j, hit = nn;

  for (j = 1; j < nn; ++j)
    {
      rmin = fmin(rmin, nodes->r[j]);
      rmax = fmax(rmax, nodes->r[j]);
    }

  if (!(r >= rmin && r <= rmax))
    return POTSYNTH_EDOM;

  *first = 0;
  *n = nn;

  for (j = 0; j < nn; ++j)
    {
      if (r == nodes->r[j])
        hit = j;
    }

  if (hit < nn)
    {
      /* at a node: l_j = delta_ij, l_j' = (w_j/w_i)/(r_i - r_j), l_i' = -sum of the others */
      double s = 0.0;

      for (j = 0; j < nn; ++j)
        f[j] = (j == hit) ? 1.0 : 0.0;

      if (df != NULL)
        {
          for (j = 0; j < nn; ++j)
            {
              if (j == hit)
                continue;

              df[j] = (nodes->w[j] / nodes->w[hit]) / (nodes->r[hit] - nodes->r[j]);
              s += df[j];
            }

          df[hit] = -s;
        }

      return POTSYNTH_SUCCESS;
    }

  for (j = 0; j < nn; ++j)
    {
      const double d = 1.0 / (r - nodes->r[j]);

      f[j] = nodes->w[j] * d;
      sum += f[j];
      sum2 += f[j] * d;
    }

  for (j = 0; j < nn; ++j)
    f[j] /= sum;

  /* l_j' = l_j [ (sum_k w_k/(r-r_k)^2) / (sum_k w_k/(r-r_k)) - 1/(r-r_j) ] */
  if (df != NULL)
    {
      for (j = 0; j < nn; ++j)
        df[j] = f[j] * (sum2 / sum - 1.0 / (r - nodes->r[j]));
    }

  return POTSYNTH_SUCCESS;
}

/*
test_radial_linear()
  Radial functions of piecewise linear interpolation between ascending nodes:
function j is the hat function of node j, and the coefficient sets are the
values of the field at the nodes. params: test_radial_nodes
*/

static int
test_radial_linear(const double r, size_t * first, size_t * n, double * f, double * df,
                           void * params)
{
  const test_radial_nodes * nodes = (const test_radial_nodes *) params;
  size_t k;

  if (nodes->n == 1)
    {
      if (r != nodes->r[0])
        return POTSYNTH_EDOM;

      *first = 0;
      *n = 1;
      f[0] = 1.0;
      if (df != NULL)
        df[0] = 0.0;

      return POTSYNTH_SUCCESS;
    }

  if (test_radial_interval(r, nodes, &k))
    return POTSYNTH_EDOM;

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

/*
test_radial_hermite()
  Radial functions of piecewise cubic Hermite interpolation between ascending
nodes: functions 2k and 2k+1 belong to the value and to the radial derivative
at node k, so there are nrad = 2 n functions and the coefficient sets are the
values and the derivatives of the field at the nodes (e.g. of a cubic spline).
params: test_radial_nodes with at least two nodes
*/

static int
test_radial_hermite(const double r, size_t * first, size_t * n, double * f, double * df,
                            void * params)
{
  const test_radial_nodes * nodes = (const test_radial_nodes *) params;
  size_t k;

  if (nodes->n < 2 || test_radial_interval(r, nodes, &k))
    return POTSYNTH_EDOM;

  {
    const double h = nodes->r[k + 1] - nodes->r[k];
    const double t = (r - nodes->r[k]) / h;
    const double t2 = t * t, t3 = t2 * t;

    *first = 2 * k;
    *n = 4;
    f[0] = 2.0 * t3 - 3.0 * t2 + 1.0;
    f[1] = h * (t3 - 2.0 * t2 + t);
    f[2] = -2.0 * t3 + 3.0 * t2;
    f[3] = h * (t3 - t2);

    if (df != NULL)
      {
        df[0] = (6.0 * t2 - 6.0 * t) / h;
        df[1] = 3.0 * t2 - 4.0 * t + 1.0;
        df[2] = (-6.0 * t2 + 6.0 * t) / h;
        df[3] = 3.0 * t2 - 2.0 * t;
      }
  }

  return POTSYNTH_SUCCESS;
}

/******************************************
 * INTERNAL ROUTINES                      *
 ******************************************/

/* interval k of ascending nodes with r in [r_k, r_k+1] (the last one for r = r_n-1);
   nonzero if r is outside of the nodes */
static int
test_radial_interval(const double r, const test_radial_nodes * nodes, size_t * k)
{
  size_t lo = 0, hi = nodes->n - 1;

  if (!(r >= nodes->r[0] && r <= nodes->r[hi]))
    return 1;

  while (hi - lo > 1)
    {
      const size_t mid = (lo + hi) / 2;

      if (r >= nodes->r[mid])
        lo = mid;
      else
        hi = mid;
    }

  *k = lo;

  return 0;
}
