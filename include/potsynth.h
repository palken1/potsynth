/*
 * potsynth.h
 *
 * Fast synthesis of internal potential-field models at scattered points:
 * generic (field-agnostic) interface. Physical front-ends (potsynth_mag.h,
 * potsynth_grav.h, potsynth_vsh.h) translate their quantities into the
 * functionals defined here.
 *
 * Each functional c is
 *
 * f_c(r,theta,phi) = \sum_{l=lmin}^{lmax} \sum_{m=0}^{min(l,mmax)}
 *                      rfac_c(l) (a/r)^(l + rpow_c) T_c[ Pbar_l^m(cos theta) ]
 *                      D_c[ C_lm cos(m phi) + S_lm sin(m phi) ]
 *
 * where Pbar_l^m are 4-pi fully normalized associated Legendre functions
 * (no Condon-Shortley phase), T_c an operator in theta and D_c = (d/dphi)^dphi_c.
 * Constant factors (a, GM/a, ...) are folded into rfac.
 *
 * A model may also hold several coefficient sets, with outputs that are sums
 * of such terms of different sets (potsynth_alloc_sets()), and take its radial
 * dependence from radial functions of the caller in place of the powers of
 * a/r (potsynth_radial): fields that are not gradients of one potential.
 *
 * Copyright (C) 2026 Patrick Alken and Johnathan Burchill, BSD 3-Clause License
 */

#ifndef INCLUDED_potsynth_h
#define INCLUDED_potsynth_h

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

/* symbol export: only the public API is exported from the shared library */
#if defined(_WIN32) && defined(POTSYNTH_SHARED)
#  ifdef POTSYNTH_BUILDING
#    define POTSYNTH_API __declspec(dllexport)
#  else
#    define POTSYNTH_API __declspec(dllimport)
#  endif
#elif defined(__GNUC__) || defined(__clang__)
#  define POTSYNTH_API __attribute__((visibility("default")))
#else
#  define POTSYNTH_API
#endif

#define POTSYNTH_VERSION_MAJOR       0
#define POTSYNTH_VERSION_MINOR       1
#define POTSYNTH_VERSION_PATCH       0
/* layout of the public structs; raised when a struct changes, so that a program (or the
   Python package) can refuse a library built from other headers (potsynth_abi()) */
#define POTSYNTH_ABI                 1

/* largest supported SH degree (double precision Holmes-Featherstone scaling) */
#define POTSYNTH_LMAX                2800

/* status codes */
#define POTSYNTH_SUCCESS             0
#define POTSYNTH_EINVAL              1 /* invalid argument */
#define POTSYNTH_ENOMEM              2 /* memory allocation failed */
#define POTSYNTH_EDOM                3 /* evaluation point outside valid domain */
#define POTSYNTH_EUNSUP              4 /* unsupported option (e.g. lmax too large) */

/* operators in theta */
#define POTSYNTH_THETA_P             0 /* Pbar_l^m(cos theta) */
#define POTSYNTH_THETA_DP            1 /* d/dtheta Pbar_l^m */
#define POTSYNTH_THETA_PSIN          2 /* Pbar_l^m / sin(theta), m >= 1, computed without dividing by sin(theta) */
#define POTSYNTH_THETA_S2A           3 /* d2Pbar/dtheta2 - cot(theta) dPbar/dtheta + m^2 Pbar / sin^2(theta)
                                          (Hessian T_thth - T_phph), pole-safe */
#define POTSYNTH_THETA_S2B           4 /* (dPbar/dtheta - cot(theta) Pbar) / sin(theta), m >= 1
                                          (Hessian T_thph), pole-safe */

/* evaluation methods */
#define POTSYNTH_METHOD_AUTO         0 /* fast or direct, chosen from lmax (and npts_hint) */
#define POTSYNTH_METHOD_FAST         1 /* gridded transform */
#define POTSYNTH_METHOD_DIRECT       2 /* point-by-point sum */

/*
potsynth_idx()
  Index of SH coefficient (l,m) in a coefficient vector spanning degrees
[lmin,lmax] and orders up to mmax:

[ (lmin,-M) ... (lmin,M) (lmin+1,-M') ... (lmin+1,M') ... (lmax,M'') ]

with M = min(l,mmax). Entry (l,m), m >= 0 holds the cosine coefficient
(C_lm or g_l^m), entry (l,-m) the sine coefficient (S_lm or h_l^m).

Input conditions: 0 <= lmin <= l, |m| <= min(l,mmax); use signed arguments
*/
#define potsynth_idx(l, m, lmin, mmax) ( ((lmin) <= (mmax)) ? \
                                           ( ((l) <= (mmax)) ? \
                                               (l)*((l)+1) + (m) - (lmin)*(lmin) : \
                                               (l)*(2*(mmax)+1) + (m) - (mmax)*(mmax) - (lmin)*(lmin) ) \
                                         : ( ((l)-(lmin))*(2*(mmax)+1) + (m) + (mmax) ) )

typedef struct
{
  int rpow;             /* radial power offset: term carries (a/r)^(l + rpow) */
  const double * rfac;  /* per-degree factors rfac[l], l = 0..lmax, or NULL for 1; copied at allocation */
  int theta;            /* operator in theta (POTSYNTH_THETA_xxx) */
  int dphi;             /* number of phi derivatives: 0 or 1 */
} potsynth_functional;


/* one term of an output of a model with several coefficient sets (potsynth_alloc_sets()):
   fac * rfac[l] * T[Pbar_l^m] D[C_lm cos(m phi) + S_lm sin(m phi)] with the coefficients of one set */
typedef struct
{
  size_t set;           /* coefficient set, 0..nset-1 */
  double fac;           /* constant factor */
  const double * rfac;  /* per-degree factors rfac[l], l = 0..lmax, or NULL for 1; copied at allocation */
  int theta;            /* operator in theta (POTSYNTH_THETA_xxx) */
  int dphi;             /* number of phi derivatives: 0 or 1 */
} potsynth_term;

/* one output of such a model: (a/r)^(l + rpow) times the sum of its terms */
typedef struct
{
  int rpow;                     /* radial power offset of all terms (not used with radial functions
                                   of the caller) */
  size_t nterm;                 /* number of terms */
  const potsynth_term * term;   /* terms, length nterm; copied at allocation */
  int component;                /* field component for the accuracy: outputs with the same value > 0
                                   form one component, and eps is relative to the largest magnitude
                                   of each component; 0: chosen by the library (outputs with the same
                                   radial power, kind of operators and coefficient sets) */
  double weight;                /* weight of the output in the magnitude of its component; 0: chosen
                                   by the library (1, and 2 for POTSYNTH_THETA_S2B outputs) */
} potsynth_output;

/*
 * Radial functions of the caller, in place of the powers of a/r: on output
 * f[0..*n-1] = f_first(r), ..., f_{first+n-1}(r), the functions that are nonzero
 * at r (all others are zero), *n <= number of functions; if df is not NULL, also
 * their derivatives (reserved for radial derivatives; not used at present).
 * Returns POTSYNTH_SUCCESS, or POTSYNTH_EDOM if r is outside of the range of the
 * functions. Must be safe to call from several threads.
 */
typedef int potsynth_radial_fn(const double r, size_t * first, size_t * n,
                               double * f, double * df, void * params);

/*
 * A radial basis of the caller for potsynth_alloc_sets(): each output is
 * \sum_j f_j(r) x (sum of its terms with the coefficient sets set + j), i.e. the
 * set of a term is the first of n consecutive sets, one for each radial function
 */
typedef struct
{
  size_t n;                     /* number of radial functions */
  potsynth_radial_fn * f;       /* the functions; NULL with n = 1: no radial dependence
                                   (a field on the sphere; r is then not used and may be NULL) */
  void * params;                /* passed to f; must remain valid while the model is used */
} potsynth_radial;

typedef struct
{
  double eps;           /* requested accuracy relative to max of each irreducible component, default 1e-10 */
  int method;           /* POTSYNTH_METHOD_xxx, default POTSYNTH_METHOD_AUTO */
  double sigma;         /* grid oversampling factor, at least 1; 0 = default (2) */
  int w;                /* window width; 0 = chosen from eps */
  double rank_tol;      /* radial SVD tolerance; 0 = 0.5 * eps */
  int nthreads;         /* OpenMP threads; 0 = all available */
  size_t npts_hint;     /* expected number of evaluation points, used by METHOD_AUTO to
                           choose the faster method including set-up; 0 = unknown */
  double theta_min;     /* smallest colatitude of the evaluation points (radians), default 0 */
  double theta_max;     /* largest colatitude of the evaluation points (radians); 0 = pi (default).
                           The fast method stores only the grid rows of this range, and points
                           outside of it are a domain error; ignored by the direct method */
  double reserved[8];   /* zero (potsynth_default_parameters()); room for later parameters */
} potsynth_parameters;

typedef struct
{
  size_t lmin;          /* minimum SH degree */
  size_t lmax;          /* maximum SH degree */
  size_t mmax;          /* maximum SH order */
  size_t nfunc;         /* number of functionals */
  size_t nout;         /* outputs per point written by potsynth_eval() (= nfunc, or more with
                           a front-end output map, e.g. 6 tensor components from 5 functionals) */
  int method;           /* method selected (POTSYNTH_METHOD_FAST or POTSYNTH_METHOD_DIRECT) */
  int w;                /* window width */
  double sigma;         /* oversampling factor */
  size_t rho;           /* number of radial functions */
  size_t ncheb;         /* Chebyshev terms per radial function */
  size_t nf1;           /* fine grid size in theta (doubled sphere) */
  size_t nf2;           /* fine grid size in phi */
  size_t nrows;         /* stored grid rows */
  size_t grid_bytes;    /* memory used by the grid */
  size_t ncomp;         /* number of irreducible components (eps is relative to each) */
  double err_radial;    /* bound of the error of the radial basis over the shell, relative to
                           the size of each component (at most rank_tol where it is reachable);
                           -1 if not known (radial functions of the caller, model read from a file) */
  double err_check;     /* largest error of the window check after set-up, against the direct
                           sum, relative to the size of each component (at most eps when the
                           check passed); -1 if there was no check (window width given, or a
                           model read from a file) */
  double reserved[8];   /* zero; room for later fields */
} potsynth_info;

/*
 * potsynth_workspace: everything that depends on the model (coefficients,
 * radial basis, grid). Read-only after allocation, so it may be shared by any
 * number of threads.
 *
 * potsynth_eval_workspace: scratch space for evaluation; each thread calling
 * potsynth_eval() concurrently needs its own.
 */
typedef struct potsynth_workspace potsynth_workspace;
typedef struct potsynth_eval_workspace potsynth_eval_workspace;

/* error handler, called by all functions on error before returning */
typedef void potsynth_error_handler_t(const char * reason, const char * file, int line,
                                      int potsynth_errno);

/*
 * Prototypes
 */

/* error.c */
POTSYNTH_API potsynth_error_handler_t * potsynth_set_error_handler(potsynth_error_handler_t * new_handler);
POTSYNTH_API potsynth_error_handler_t * potsynth_set_error_handler_off(void);
POTSYNTH_API const char * potsynth_strerror(const int potsynth_errno);
POTSYNTH_API const char * potsynth_version(void);
POTSYNTH_API int potsynth_abi(void);

/* setup.c */
POTSYNTH_API potsynth_parameters potsynth_default_parameters(void);
POTSYNTH_API size_t potsynth_nlm(const size_t lmin, const size_t lmax, const size_t mmax);
POTSYNTH_API potsynth_workspace * potsynth_alloc(const size_t lmin, const size_t lmax, const size_t mmax,
                                                 const double * coef, const double a,
                                                 const double rmin, const double rmax,
                                                 const size_t nfunc, const potsynth_functional * func,
                                                 const potsynth_parameters * params);
POTSYNTH_API potsynth_workspace * potsynth_alloc_sets(const size_t lmin, const size_t lmax, const size_t mmax,
                                                      const size_t nset, const double * coef,
                                                      const double a, const double rmin, const double rmax,
                                                      const potsynth_radial * radial,
                                                      const size_t nout, const potsynth_output * out,
                                                      const potsynth_parameters * params);
POTSYNTH_API void potsynth_free(potsynth_workspace * w);
POTSYNTH_API int potsynth_get_info(potsynth_info * info, const potsynth_workspace * w);

/* io.c: files of models of the fast method */
POTSYNTH_API int potsynth_write(const char * filename, const potsynth_workspace * w);
POTSYNTH_API potsynth_workspace * potsynth_read(const char * filename, const potsynth_radial * radial,
                                                const int nthreads);

/* frame.c */
POTSYNTH_API void potsynth_tensor_rtp2nec(const size_t n, double * T);
POTSYNTH_API void potsynth_tensor_rtp2xyz(const size_t n, const double * theta, const double * phi, double * T);

/* util.c */
POTSYNTH_API int potsynth_cpu_avx2_fma(void);

/* eval.c */
POTSYNTH_API potsynth_eval_workspace * potsynth_eval_alloc(const size_t chunk, const potsynth_workspace * w);
POTSYNTH_API void potsynth_eval_free(potsynth_eval_workspace * ew);
POTSYNTH_API int potsynth_eval_set_threads(potsynth_eval_workspace * ew, const int nthreads);
POTSYNTH_API int potsynth_eval(const size_t n, const double * r, const double * theta, const double * phi,
                               double * out, potsynth_eval_workspace * ew);
POTSYNTH_API int potsynth_eval_rows(const size_t n, const double * r, const double * theta, const double * phi,
                                    double * out, potsynth_eval_workspace * ew);

#ifdef __cplusplus
}
#endif

#endif /* INCLUDED_potsynth_h */
