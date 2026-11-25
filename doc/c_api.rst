.. _c_api:

*****
C API
*****

There are four headers:

- ``potsynth.h``: the generic core, with functionals of an internal potential.
- ``potsynth_mag.h``: the magnetic front-end. It includes ``potsynth.h``.
- ``potsynth_grav.h``: the gravity front-end. It includes ``potsynth.h``.
- ``potsynth_vsh.h``: the vector spherical harmonic front-end, for vector fields that
  need not be gradients of a potential. It includes ``potsynth.h``.

Typical use
===========

.. code-block:: c

   #include <potsynth_mag.h>

   potsynth_parameters params = potsynth_default_parameters();
   potsynth_workspace * w;
   potsynth_eval_workspace * ew;
   int status;

   params.eps = 1.0e-10;                        /* requested accuracy */

   w = potsynth_mag_alloc(1, lmax, lmax, gh, 6371.2, rmin, rmax, POTSYNTH_MAG_B, &params);
   ew = potsynth_eval_alloc(0, w);

   /* B[0..n) = B_r, B[n..2n) = B_theta, B[2n..3n) = B_phi */
   status = potsynth_eval(n, r, theta, phi, B, ew);

   potsynth_eval_free(ew);
   potsynth_free(w);

Objects and thread safety
=========================

The library uses two opaque objects:

.. c:type:: potsynth_workspace

   Holds everything that depends on the model: radial basis, grids, and for the direct
   method the coefficients. It is built by one of the allocation functions
   (:c:func:`potsynth_alloc`, :c:func:`potsynth_alloc_sets`,
   :c:func:`potsynth_mag_alloc`, :c:func:`potsynth_grav_alloc`,
   :c:func:`potsynth_vsh_alloc`) or read from a file (:c:func:`potsynth_read`), and is
   read-only afterwards, so any number of threads can share it. Radial functions of the
   caller are called from all of them.

.. c:type:: potsynth_eval_workspace

   Scratch space for evaluation: point sorting and per-thread buffers. Each thread
   that calls :c:func:`potsynth_eval` or :c:func:`potsynth_eval_rows` at the same time
   as others needs its own.

Evaluation never allocates memory. All allocation happens in the ``*_alloc``
functions. Within one call, :c:func:`potsynth_eval` uses the number of OpenMP
threads given in the parameters.

Errors
======

Errors are reported in the same way as in GSL. Every function that detects an
error calls the current error handler, then returns an error code or ``NULL``.
Unlike GSL, the default handler does not abort. It prints a message such as::

   potsynth: src/eval.c:181: ERROR: theta outside [0,pi], ... (evaluation point outside the valid domain)

to ``stderr`` and lets the function return.

.. c:type:: potsynth_error_handler_t

   ``void potsynth_error_handler_t(const char * reason, const char * file, int line, int potsynth_errno)``

.. c:function:: potsynth_error_handler_t * potsynth_set_error_handler(potsynth_error_handler_t * new_handler)

   Installs a new error handler and returns the previous one.

.. c:function:: potsynth_error_handler_t * potsynth_set_error_handler_off(void)

   Disables error messages. Functions still return error codes or ``NULL``.

.. c:function:: const char * potsynth_strerror(const int potsynth_errno)

   Returns a short description of an error code.

The error codes are:

.. c:macro:: POTSYNTH_SUCCESS

   0: success.

.. c:macro:: POTSYNTH_EINVAL

   An argument is invalid.

.. c:macro:: POTSYNTH_ENOMEM

   A memory allocation failed.

.. c:macro:: POTSYNTH_EDOM

   An evaluation point is outside the valid domain.

.. c:macro:: POTSYNTH_EUNSUP

   The option is not supported, for example :math:`L >` :c:macro:`POTSYNTH_LMAX`.

Coefficient layout
==================

Coefficients are passed as a single vector. For degrees :math:`l_{min} \le l \le l_{max}` and orders
:math:`|m| \le \min(l, m_{max})` the order is

.. code-block:: text

   [ (lmin,-M) ... (lmin,M)  (lmin+1,-M') ... (lmin+1,M')  ...  (lmax,M'') ]

Entry :math:`(l, m)` with :math:`m \ge 0` holds the cosine coefficient
(:math:`C_{lm}` or :math:`g_l^m`). Entry :math:`(l, -m)` holds the sine coefficient
(:math:`S_{lm}` or :math:`h_l^m`). The layout is read only once, during set-up.
Internally the library reorders everything to match its recursions.

.. c:macro:: potsynth_idx(l, m, lmin, mmax)

   Index of coefficient :math:`(l, m)` in the vector, with
   :math:`-\min(l,m_{max}) \le m \le \min(l,m_{max})`. Pass signed integers.

.. c:function:: size_t potsynth_nlm(const size_t lmin, const size_t lmax, const size_t mmax)

   Length of the coefficient vector.

.. c:macro:: POTSYNTH_LMAX

   Largest supported degree (2800).

Parameters
==========

.. c:struct:: potsynth_parameters

   Get the defaults from :c:func:`potsynth_default_parameters` and change fields as
   needed. A ``NULL`` pointer in place of the parameters means all defaults.

   .. c:member:: double eps

      Requested accuracy, relative to the largest magnitude of each irreducible component: scalar
      functionals (:math:`V`, :math:`B_r`, :math:`g_r`, :math:`T_{rr}`) relative to their own
      maximum, horizontal components relative to the maximum of the horizontal vector
      (e.g. :math:`B_\theta` and :math:`B_\phi` relative to :math:`\max|\mathbf{B}_h|`), see
      :ref:`accuracy`. Default :math:`10^{-10}`. Any value from about :math:`10^{-4}` down to
      :math:`5 \times 10^{-15}` (the double-precision floor) is meaningful.

   .. c:member:: int method

      :c:macro:`POTSYNTH_METHOD_AUTO` (default),
      :c:macro:`POTSYNTH_METHOD_FAST` or :c:macro:`POTSYNTH_METHOD_DIRECT`.

   .. c:member:: int nthreads

      Number of OpenMP threads; 0 (the default) means all available.

   .. c:member:: double sigma

      Oversampling factor of the grid; 0 selects the default, 2. Values below 1 are
      rejected.

   .. c:member:: int w

      Window width; 0 means choose it from :c:member:`eps`.

   .. c:member:: double rank_tol

      Tolerance of the radial SVD; 0 means :math:`0.5\,\varepsilon`.

   .. c:member:: size_t npts_hint

      Expected number of points that will be evaluated with this model; 0 (the
      default) means unknown. :c:macro:`POTSYNTH_METHOD_AUTO` uses it to choose the
      method with the smaller total time, set-up plus evaluation. The fast method's
      set-up only pays off for enough points; see :ref:`choosing-method`.

   .. c:member:: double theta_min

      Smallest colatitude of the evaluation points (radians). The default, 0, with
      :c:member:`theta_max` = 0 means the whole sphere.

   .. c:member:: double theta_max

      Largest colatitude of the evaluation points (radians); 0 means :math:`\pi`.

      The fast method stores only the grid rows that points in
      :math:`[\theta_{min}, \theta_{max}]` can reach, so the memory of the grid shrinks
      roughly in proportion to the range. The values are the same as those of the model
      of the whole sphere. A point outside the range is a domain error
      (:c:macro:`POTSYNTH_EDOM`), as a radius outside :math:`[r_{min}, r_{max}]` is. The
      direct method ignores the range. Use it for data in a region, such as a polar cap
      or a regional survey.

.. c:function:: potsynth_parameters potsynth_default_parameters(void)

   Returns the default parameters.

.. c:macro:: POTSYNTH_METHOD_AUTO

   With :c:member:`potsynth_parameters.npts_hint` set: use the method with the smaller estimated total
   time for that many points (see :ref:`choosing-method`). Without it: use the fast
   method for :math:`l_{max} \ge 24` and the direct method otherwise. If no valid
   radial range :math:`0 < r_{min} \le r_{max}` is given, the direct method is used.

.. c:macro:: POTSYNTH_METHOD_FAST

   Gridded transform (see :ref:`algorithm`). Needs the radial range
   :math:`[r_{min}, r_{max}]` of the points.

.. c:macro:: POTSYNTH_METHOD_DIRECT

   Point-by-point sum. Accepts any :math:`r > 0`.

Magnetic front-end
==================

.. c:macro:: POTSYNTH_MAG_B

   The field :math:`(B_r, B_\theta, B_\phi) = -\nabla V` (nT).

.. c:macro:: POTSYNTH_MAG_GRADB

   The gradient tensor :math:`\nabla\mathbf{B} = -\nabla\nabla V` (nT/km), as six
   components :math:`\partial B_{rr}, \partial B_{r\theta}, \partial B_{r\phi},
   \partial B_{\theta\theta}, \partial B_{\theta\phi}, \partial B_{\phi\phi}` in the
   spherical basis (see :ref:`gradient-tensor`). The tensor is symmetric and traceless;
   it is computed from five functionals.

The flags are combined with ``|``. The outputs are in the order :math:`B_r, B_\theta,
B_\phi`, then the six tensor components, restricted to those requested. At the poles,
the :math:`\theta` and :math:`\phi` components refer to the local frame defined by
the longitude of the point; :c:func:`potsynth_tensor_rtp2xyz` gives components that are
unique there.

.. c:function:: potsynth_workspace * potsynth_mag_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const double * gh, const double a, const double rmin, const double rmax, const int quantities, const potsynth_parameters * params)

   Builds a magnetic field model for the potential :eq:`mag_potential`: the field
   :math:`\mathbf{B} = -\nabla V` and/or its gradient tensor.

   :param lmin: minimum degree, :math:`\ge 1`
   :param lmax: maximum degree, :math:`\le` :c:macro:`POTSYNTH_LMAX`
   :param mmax: maximum order, :math:`\le` ``lmax``
   :param gh: Schmidt semi-normalized Gauss coefficients in nT, indexed by
      :c:macro:`potsynth_idx`. The array is copied.
   :param a: reference radius (km), usually 6371.2
   :param rmin: smallest radius (km) at which the model will be evaluated. Required
      by the fast method; ignored by the direct method.
   :param rmax: largest radius (km) at which the model will be evaluated, as for
      ``rmin``
   :param quantities: ``POTSYNTH_MAG_xxx`` flags
   :param params: parameters, or ``NULL`` for the defaults
   :returns: the workspace, or ``NULL`` on error (reported to the error handler)

   Choose :math:`[r_{min}, r_{max}]` as tight as the data allow. A narrower shell
   needs fewer radial functions, which makes evaluation faster.

   For example, the field and its gradient tensor, 9 outputs per point:

   .. code-block:: c

      w = potsynth_mag_alloc(1, lmax, lmax, gh, 6371.2, rmin, rmax,
                             POTSYNTH_MAG_B | POTSYNTH_MAG_GRADB, &params);
      ew = potsynth_eval_alloc(0, w);
      status = potsynth_eval(n, r, theta, phi, out, ew);  /* out: 9*n doubles */

      potsynth_mag_rtp2nec(n, out);              /* B: rows 0-2 */
      potsynth_tensor_rtp2nec(n, out + 3 * n);   /* grad B: rows 3-8 */

.. c:function:: size_t potsynth_mag_nout(const int quantities)

   The number of outputs per point for a set of flags, or 0 if the flags are invalid.

.. c:function:: void potsynth_mag_rtp2nec(const size_t n, double * B)

   Converts output in place, from :math:`(B_r, B_\theta, B_\phi)` to
   :math:`(B_N, B_E, B_C) = (-B_\theta, B_\phi, -B_r)`.

SHC files
---------

Static models in the Swarm SHC format can be read directly, for example the
lithospheric model LCS-1. The file contains ``#`` comment lines, a header line
``lmin lmax ntimes [spline_order nstep]``, a line with the epochs, and one line
``l m value`` per coefficient (:math:`m \ge 0` for :math:`g_l^m`, :math:`m < 0` for
:math:`h_l^{|m|}`). Only static models (one epoch) are supported. Time-dependent
core field models are of low degree, where a direct sum is already fast.

.. c:struct:: potsynth_shc

   .. c:member:: size_t lmin

      minimum degree

   .. c:member:: size_t lmax

      maximum degree

   .. c:member:: size_t mmax

      maximum order (equal to ``lmax``)

   .. c:member:: double epoch

      model epoch (decimal year)

   .. c:member:: double * gh

      Schmidt semi-normalized Gauss coefficients (nT), indexed by
      :c:macro:`potsynth_idx` ``(l, m, lmin, mmax)``, ready for
      :c:func:`potsynth_mag_alloc`

.. c:function:: potsynth_shc * potsynth_shc_read(const char * filename)

   Reads a static model from an SHC file. All coefficients for
   :math:`l_{min} \le l \le l_{max}`, :math:`|m| \le l` must be present exactly once.
   Malformed files are reported to the error handler with the line number. Problems
   detected include a missing header or epoch, out-of-range or duplicate entries,
   missing coefficients, extra columns, non-finite values and over-long lines.
   Time-dependent files give :c:macro:`POTSYNTH_EUNSUP`.

   :returns: the model, or ``NULL`` on error

.. c:function:: void potsynth_shc_free(potsynth_shc * shc)

   Frees a model read by :c:func:`potsynth_shc_read`. A NULL pointer is ignored.

.. c:function:: potsynth_workspace * potsynth_mag_alloc_shc(const char * filename, const double a, const double rmin, const double rmax, const int quantities, const potsynth_parameters * params)

   Reads an SHC file and builds the magnetic field model, as
   :c:func:`potsynth_mag_alloc`.

Gravity front-end
=================

The gravity front-end evaluates the gravitational potential :eq:`grav_potential` and the
acceleration :math:`\mathbf{g} = \nabla V` of a model given by 4-:math:`\pi` fully
normalized coefficients :math:`\bar{C}_{lm}, \bar{S}_{lm}`. It computes the synthesis
only:

- no centrifugal term and no normal field are included;
- coordinates are geocentric spherical, with :math:`r` in metres;
- coefficients are used as given (no tide-system conversion).

For disturbing quantities, subtract the normal field's even zonal coefficients before
set-up.

.. c:macro:: POTSYNTH_GRAV_V

   The potential :math:`V` (m\ :sup:`2`/s\ :sup:`2`).

.. c:macro:: POTSYNTH_GRAV_GR

   The radial acceleration :math:`g_r = \partial V/\partial r` (m/s\ :sup:`2`), which is
   negative (towards the centre).

.. c:macro:: POTSYNTH_GRAV_G

   The acceleration :math:`(g_r, g_\theta, g_\phi) = \nabla V` (m/s\ :sup:`2`), with
   :math:`g_\theta` positive southwards. It includes :c:macro:`POTSYNTH_GRAV_GR`.

.. c:macro:: POTSYNTH_GRAV_T

   The gravity gradient tensor :math:`T = \nabla\nabla V` (s\ :sup:`-2`;
   :math:`1\ \mathrm{E} = 10^{-9}` s\ :sup:`-2`), as six components :math:`T_{rr},
   T_{r\theta}, T_{r\phi}, T_{\theta\theta}, T_{\theta\phi}, T_{\phi\phi}` in the
   spherical basis (see :ref:`gradient-tensor`), computed from five functionals.

The flags are combined with ``|``. The outputs are in the order :math:`V, g_r,
g_\theta, g_\phi`, then the six tensor components, restricted to those requested. Each
functional adds to the set-up time and grid memory of the fast method, so request only
what is needed: for example :c:macro:`POTSYNTH_GRAV_GR` alone for radial gravity.

.. c:function:: potsynth_workspace * potsynth_grav_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef, const double GM, const double a, const double rmin, const double rmax, const int quantities, const potsynth_parameters * params)

   Builds a gravity field model.

   :param lmin: minimum degree, :math:`\ge 0`
   :param lmax: maximum degree, :math:`\le` :c:macro:`POTSYNTH_LMAX`
   :param mmax: maximum order, :math:`\le` ``lmax``
   :param coef: fully normalized coefficients, indexed by :c:macro:`potsynth_idx`
      (:math:`\bar{C}_{lm}` at :math:`m \ge 0`, :math:`\bar{S}_{lm}` at :math:`-m`). The
      array is copied.
   :param GM: gravitational constant times mass (m\ :sup:`3`/s\ :sup:`2`)
   :param a: reference radius (m)
   :param rmin: smallest radius (m) at which the model will be evaluated. Required by
      the fast method; ignored by the direct method.
   :param rmax: largest radius (m) at which the model will be evaluated, as for ``rmin``
   :param quantities: ``POTSYNTH_GRAV_xxx`` flags
   :param params: parameters, or ``NULL`` for the defaults
   :returns: the workspace, or ``NULL`` on error (reported to the error handler)

   The accuracy :c:member:`potsynth_parameters.eps` applies to each irreducible component
   separately (:math:`V`, :math:`g_r`, the horizontal acceleration
   :math:`(g_\theta, g_\phi)`, and the three parts of the tensor; see :ref:`accuracy`).
   For a full model, each of them is dominated by its central
   (:math:`l = 0`) term; for :math:`10^{-10}` the potential is accurate to about
   :math:`10^{-3}` m\ :sup:`2`/s\ :sup:`2` (0.1 mm of geoid height) and
   :math:`g_r` to about :math:`2 \times 10^{-5}` mGal, measured with EGM2008 at satellite
   altitudes.

.. c:function:: size_t potsynth_grav_nout(const int quantities)

   The number of outputs per point for a set of flags, or 0 if the flags are invalid.

ICGEM files
-----------

Static models in the ICGEM format (``.gfc`` files, as distributed by the International
Centre for Global Earth Models, http://icgem.gfz-potsdam.de) can be read directly. The
file has a header ending with a line ``end_of_head``, and one line
``gfc l m C S [sigma_C sigma_S]`` per coefficient. The reader uses these header keywords:

- ``earth_gravity_constant``, ``radius`` and ``max_degree``, which are required;
- ``norm``: ``fully_normalized`` (the default) or ``unnormalized``, which is converted;
- ``modelname``, ``tide_system`` and ``product_type`` (which must be ``gravity_field``
  if present).

Other header lines are ignored. Numbers may use Fortran exponents (``1.0D-09``).
Coefficients absent from the file are zero: many models omit degrees 0 and 1. Files with
time-variable terms (``gfct``, ``trnd``, ``dot``, ``acos``, ``asin``) give
:c:macro:`POTSYNTH_EUNSUP`.

.. c:struct:: potsynth_gfc

   .. c:member:: size_t lmin

      minimum degree (0)

   .. c:member:: size_t lmax

      maximum degree (``max_degree``, or the truncation degree)

   .. c:member:: size_t mmax

      maximum order (equal to ``lmax``)

   .. c:member:: double GM

      ``earth_gravity_constant`` (m\ :sup:`3`/s\ :sup:`2`)

   .. c:member:: double a

      ``radius`` (m)

   .. c:member:: char modelname[128]

      ``modelname`` from the header

   .. c:member:: char tide_system[64]

      ``tide_system`` from the header, or ``unknown``

   .. c:member:: double * coef

      fully normalized coefficients, indexed by :c:macro:`potsynth_idx`
      ``(l, m, 0, lmax)``, ready for :c:func:`potsynth_grav_alloc`

.. c:function:: potsynth_gfc * potsynth_gfc_read(const char * filename, const size_t lmax)

   Reads a static model from a ``.gfc`` file, truncated at degree ``lmax`` (0 for all
   degrees of the file). Malformed files are reported to the error handler with the line
   number. Problems detected include a missing ``end_of_head`` or required keyword,
   unknown keywords, invalid numbers, :math:`m > l`, :math:`l >` ``max_degree``,
   duplicate entries, a nonzero :math:`S_{l0}` and over-long lines.

   :returns: the model, or ``NULL`` on error

.. c:function:: void potsynth_gfc_free(potsynth_gfc * gfc)

   Frees a model read by :c:func:`potsynth_gfc_read`. A NULL pointer is ignored.

.. c:function:: potsynth_workspace * potsynth_grav_alloc_gfc(const char * filename, const size_t lmax, const double rmin, const double rmax, const int quantities, const potsynth_parameters * params)

   Reads a ``.gfc`` file, truncated at degree ``lmax`` (0 for all degrees), and builds
   the gravity field model with the file's :math:`GM` and reference radius, as
   :c:func:`potsynth_grav_alloc`. For example, radial gravity from EGM2008 to degree 720
   at satellite altitudes:

   .. code-block:: c

      #include <potsynth_grav.h>

      w = potsynth_grav_alloc_gfc("EGM2008.gfc", 720, 6.55e6, 6.9e6, POTSYNTH_GRAV_GR, NULL);
      ew = potsynth_eval_alloc(0, w);
      status = potsynth_eval(n, r, theta, phi, g_r, ew);   /* r in m, g_r in m/s^2 */

.. _c-geodesy:

Geodetic helpers
----------------

Small helpers for the recipes of :ref:`gravity-recipes`. The constants derived from an
ellipsoid (:math:`J_2`, normal gravity at the equator and the poles) are computed from
its four defining parameters with the closed formulas of Moritz (2000); they reproduce
the published WGS84 and GRS80 values.

.. c:struct:: potsynth_ellipsoid

   .. c:member:: double a

      semi-major axis (m)

   .. c:member:: double f

      flattening

   .. c:member:: double GM

      geocentric gravitational constant (m\ :sup:`3`/s\ :sup:`2`)

   .. c:member:: double omega

      angular velocity (rad/s)

.. c:var:: const potsynth_ellipsoid potsynth_wgs84

   WGS84: :math:`a = 6378137` m, :math:`1/f = 298.257223563`,
   :math:`GM = 3.986004418 \times 10^{14}`, :math:`\omega = 7.292115 \times 10^{-5}`.

.. c:var:: const potsynth_ellipsoid potsynth_grs80

   GRS80: :math:`a = 6378137` m, :math:`1/f = 298.257222101`,
   :math:`GM = 3.986005 \times 10^{14}`, :math:`\omega = 7.292115 \times 10^{-5}`.

.. c:function:: int potsynth_geodetic2geocentric(const size_t n, const double * lat, const double * h, const potsynth_ellipsoid * ell, double * r, double * theta)

   Converts geodetic latitude ``lat`` (radians) and ellipsoidal height ``h`` (m) to
   geocentric radius ``r`` (m) and colatitude ``theta`` (radians). The longitude is the
   same in both systems. The outputs may be the input arrays.

.. c:function:: int potsynth_normal_zonals(const potsynth_ellipsoid * ell, const size_t lmax, double * C)

   Fills ``C[l]``, :math:`l = 0 \ldots l_{max}`, with the 4-:math:`\pi` normalized zonal
   coefficients of the ellipsoid's normal gravitational potential, in its own scaling:
   :math:`C_0 = 1`, :math:`C_{2n} = -J_{2n}/\sqrt{4n+1}` for :math:`2n \le 20`, and 0
   otherwise.

.. c:function:: int potsynth_grav_subtract_normal(const size_t lmin, const size_t lmax, const size_t mmax, double * coef, const double GM, const double a, const potsynth_ellipsoid * ell)

   Subtracts the normal gravitational potential from a model in place, giving the
   coefficients of the disturbing potential :math:`T = V - U`. The normal zonals are
   rescaled to the model, :math:`C_l\, (GM_e/GM)(a_e/a)^l`. Degree 0 is included when
   :math:`l_{min} = 0`, so :math:`T` contains the term :math:`(GM - GM_e)/r`.

.. c:function:: int potsynth_normal_gravity(const size_t n, const double * lat, const potsynth_ellipsoid * ell, double * gamma)

   Normal gravity (m/s\ :sup:`2`) on the surface of the ellipsoid at geodetic latitude
   ``lat``, by Somigliana's formula.

Generic core
============

.. c:struct:: potsynth_functional

   One functional, as in :eq:`functional`.
   ``potsynth_channel`` is its former name, kept as a deprecated alias.

   .. c:member:: int rpow

      Radial power offset :math:`p`: the term carries :math:`(a/r)^{l+p}`.
      :math:`l_{min} + p \ge 0` is required.

   .. c:member:: const double * rfac

      Per-degree factors :math:`\rho(l)` for :math:`l = 0 \ldots l_{max}` (length
      :math:`l_{max}+1`), or ``NULL`` for 1. The array is copied.

   .. c:member:: int theta

      The operator in colatitude: one of the ``POTSYNTH_THETA_xxx`` constants below.

   .. c:member:: int dphi

      0, or 1 to apply :math:`\partial/\partial\phi` to
      :math:`C\cos m\phi + S\sin m\phi`.

.. c:macro:: POTSYNTH_THETA_P

   :math:`\bar{P}_l^m(\cos\theta)`

.. c:macro:: POTSYNTH_THETA_DP

   :math:`d\bar{P}_l^m/d\theta`

.. c:macro:: POTSYNTH_THETA_PSIN

   :math:`\bar{P}_l^m/\sin\theta`, computed without dividing by :math:`\sin\theta`.
   It is defined for :math:`m \ge 1`; the :math:`m = 0` terms are omitted.

.. c:macro:: POTSYNTH_THETA_S2A

   :math:`\mathcal{A}_l^m = d^2\bar{P}_l^m/d\theta^2 - \cot\theta\, d\bar{P}_l^m/d\theta
   + m^2 \bar{P}_l^m/\sin^2\theta`, the angular part of :math:`T_{\theta\theta} -
   T_{\phi\phi}` in the gradient tensor (see :ref:`gradient-tensor`). It is computed
   without dividing by :math:`\sin\theta` and is exact at the poles.

.. c:macro:: POTSYNTH_THETA_S2B

   :math:`\mathcal{B}_l^m = (d\bar{P}_l^m/d\theta - \cot\theta\, \bar{P}_l^m)/\sin\theta`,
   the angular part of :math:`T_{\theta\phi}` (used with ``dphi = 1``). It is defined
   for :math:`m \ge 1`; the :math:`m = 0` terms are omitted.

Both methods support all operators.

The front-ends provide the gradient tensor directly (:c:macro:`POTSYNTH_MAG_GRADB`,
:c:macro:`POTSYNTH_GRAV_T`). For other potentials, the gradient tensor
:math:`T = \nabla\nabla V` of
:math:`V = \sum (a/r)^{l+1} \bar{P}_l^m [\ldots]` is given by five functionals, all with
``rpow = 3``:

.. list-table::
   :header-rows: 1

   * - component
     - ``rfac``
     - ``theta``
     - ``dphi``
   * - :math:`T_{rr}`
     - :math:`(l+1)(l+2)/a^2`
     - ``P``
     - 0
   * - :math:`T_{r\theta}`
     - :math:`-(l+2)/a^2`
     - ``DP``
     - 0
   * - :math:`T_{r\phi}`
     - :math:`-(l+2)/a^2`
     - ``PSIN``
     - 1
   * - :math:`D = T_{\theta\theta} - T_{\phi\phi}`
     - :math:`1/a^2`
     - ``S2A``
     - 0
   * - :math:`T_{\theta\phi}`
     - :math:`1/a^2`
     - ``S2B``
     - 1

with :math:`T_{\theta\theta} = (-T_{rr} + D)/2` and :math:`T_{\phi\phi} = (-T_{rr} - D)/2`.

.. c:function:: potsynth_workspace * potsynth_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const double * coef, const double a, const double rmin, const double rmax, const size_t nfunc, const potsynth_functional * chan, const potsynth_parameters * params)

   Builds a model with user-defined functionals.

   ``coef`` holds the 4-:math:`\pi` fully normalized coefficients, indexed by
   :c:macro:`potsynth_idx`. ``lmin`` may be 0. The remaining arguments are as for
   :c:func:`potsynth_mag_alloc`. All arrays are copied.

   For example, the magnetic front-end is equivalent to the three functionals
   ``{2, rfac=l+1, THETA_P, 0}``, ``{2, rfac=-1, THETA_DP, 0}`` and
   ``{2, rfac=-1, THETA_PSIN, 1}``, with ``coef = gh / sqrt(2l+1)``.

Several coefficient sets
------------------------

A model may hold several coefficient sets. Each output is then a sum of terms, and each
term applies one operator to one set. All outputs share one grid of the fast method, so
one window sum per point gives all of them. This is how fields that are not gradients of
one potential are evaluated, for example the components of a vector spherical harmonic
expansion (see :ref:`vsh-frontend`), and how many fields are evaluated at once.

The terms of an output are evaluated separately and added, by both methods. When the
terms nearly cancel, as for the difference of two models, rounding errors of the size of
the separate terms remain, about :math:`10^{-16}` of them; relative to the output this
is larger by the ratio of the terms to the output, and ``eps`` cannot be met below that
(``err_check`` in :c:struct:`potsynth_info` shows it). Form the difference of the
coefficients first where possible.

.. c:struct:: potsynth_term

   One term of an output:
   :math:`f\,\rho(l)\; T\!\left[\bar{P}_l^m\right] D\!\left[C_{lm}\cos m\phi + S_{lm}\sin m\phi\right]`
   with the coefficients :math:`C_{lm}, S_{lm}` of one set.

   .. c:member:: size_t set

      Coefficient set of the term, :math:`0 \ldots n_{set}-1`.

   .. c:member:: double fac

      Constant factor :math:`f`.

   .. c:member:: const double * rfac

      Per-degree factors :math:`\rho(l)` for :math:`l = 0 \ldots l_{max}`, or ``NULL``
      for 1. The array is copied.

   .. c:member:: int theta

      The operator in colatitude, one of the ``POTSYNTH_THETA_xxx`` constants.

   .. c:member:: int dphi

      0, or 1 to apply :math:`\partial/\partial\phi`.

.. c:struct:: potsynth_output

   One output: :math:`(a/r)^{l+p}` times the sum of its terms.

   .. c:member:: int rpow

      Radial power offset :math:`p` of all terms of the output.

   .. c:member:: size_t nterm

      Number of terms.

   .. c:member:: const potsynth_term * term

      The terms, length ``nterm``. They are copied.

   .. c:member:: int component

      Field component of the output for the accuracy. Outputs with the same value
      :math:`> 0` form one component, and ``eps`` is relative to the largest magnitude
      of each component. 0 (the default of a zero-initialized structure) lets the
      library choose: outputs with the same radial power, the same kind of operators
      and the same coefficient sets form a component.

   .. c:member:: double weight

      Weight of the output in the magnitude of its component; 0 lets the library choose
      (1, and 2 for ``POTSYNTH_THETA_S2B`` outputs, which count twice in a tensor).

.. c:function:: potsynth_workspace * potsynth_alloc_sets(const size_t lmin, const size_t lmax, const size_t mmax, const size_t nset, const double * coef, const double a, const double rmin, const double rmax, const potsynth_radial * radial, const size_t nout, const potsynth_output * out, const potsynth_parameters * params)

   Builds a model with ``nset`` coefficient sets and ``nout`` outputs. ``radial`` is
   ``NULL`` for the powers of :math:`a/r`, or the radial functions of the caller (see
   :ref:`radial-functions`).

   ``coef`` holds the sets one after the other: set :math:`s` starts at
   ``coef + s * potsynth_nlm(lmin, lmax, mmax)`` and is indexed by
   :c:macro:`potsynth_idx`. The other arguments are as for :c:func:`potsynth_alloc`,
   which is the case of one set and one term with factor 1 per output.

   The terms of an output must form one field component, because the accuracy is
   controlled per component (see :ref:`accuracy`). Either all terms have the same
   operators, or each term is one of ``(POTSYNTH_THETA_DP, dphi = 0)`` and
   ``(POTSYNTH_THETA_PSIN, dphi = 1)``, the two operators of a horizontal vector; or,
   for an output with a ``component`` and a ``weight`` given, each term is one of
   ``(POTSYNTH_THETA_S2A, dphi = 0)`` and ``(POTSYNTH_THETA_S2B, dphi = 1)``, the two
   operators of a horizontal tensor. Other combinations give
   :c:macro:`POTSYNTH_EINVAL`. The outputs of a component must have the same radial
   power. :c:member:`potsynth_info.ncomp` is the number of components of a model.

   For example, the horizontal components of
   :math:`\sum_{lm} b_{lm}\boldsymbol{\Psi}_{lm} + c_{lm}\boldsymbol{\Phi}_{lm}` with the
   sets :math:`b` (0) and :math:`c` (1) are the two outputs

   .. code-block:: c

      potsynth_term tt[2] = { { 0,  1.0, NULL, POTSYNTH_THETA_DP,   0 },
                              { 1, -1.0, NULL, POTSYNTH_THETA_PSIN, 1 } };
      potsynth_term tp[2] = { { 0,  1.0, NULL, POTSYNTH_THETA_PSIN, 1 },
                              { 1,  1.0, NULL, POTSYNTH_THETA_DP,   0 } };
      potsynth_output out[2] = { { 0, 2, tt, 0, 0.0 }, { 0, 2, tp, 0, 0.0 } };

      w = potsynth_alloc_sets(1, lmax, lmax, 2, coef, a, rmin, rmax, NULL, 2, out, &params);

.. _radial-functions:

Radial functions of the caller
------------------------------

By default the radial dependence of an output is :math:`(a/r)^{l+p}`, which the fast
method expands in a few radial functions (step 1 of :ref:`algorithm`). A model may
instead take its radial functions :math:`f_j(r)` from the caller. Each output is then

.. math::

   f(r,\theta,\phi) = \sum_{j=0}^{n-1} f_j(r)\; G_j(\theta,\phi),

where :math:`G_j` is the sum of the output's terms with the coefficient sets
``set + j``: the set of a term names the first of :math:`n` consecutive sets, one for
each radial function. The coefficients may be values at radial nodes with any linear
interpolation rule, spline coefficients, or anything else that is linear in them. This
is the separated form of the fast method with the caller's functions in place of its
own, so the grid, the window sum and the accuracy control are unchanged; ``a`` and
``rpow`` are not used. At a point only the functions that are nonzero there are
evaluated and summed, so the cost per point does not grow with the number of
functions of a local basis.

.. c:type:: int potsynth_radial_fn(const double r, size_t * first, size_t * n, double * f, double * df, void * params)

   Evaluates the radial functions that are nonzero at ``r``: on output ``f[0..*n-1]``
   are :math:`f_{first}(r), \ldots, f_{first+n-1}(r)`, and all other functions are
   zero. ``df``, if not ``NULL``, receives their derivatives (reserved for radial
   derivatives; the library does not use it at present). Returns
   :c:macro:`POTSYNTH_SUCCESS`, or :c:macro:`POTSYNTH_EDOM` if ``r`` is outside the
   range of the functions. The function must be safe to call from several threads.

.. c:struct:: potsynth_radial

   A radial basis of the caller.

   .. c:member:: size_t n

      Number of radial functions.

   .. c:member:: potsynth_radial_fn * f

      The functions. ``NULL`` with ``n = 1`` means no radial dependence, a field on
      the sphere: the radii of the evaluation points are then not used and may be
      ``NULL``.

   .. c:member:: void * params

      Passed to ``f``. Must remain valid while the model is used.

``rmin`` and ``rmax`` of the allocation function must lie in the range of the
functions. A point whose radius the functions refuse is a domain error.

The library has no radial functions of its own: how a field depends on radius, and
how it is interpolated between the radii at which it is known, is the caller's. A
radial function is short. This one interpolates linearly between ascending nodes, so
that the coefficient sets of the model are the values of its fields at the nodes:

.. code-block:: c

   typedef struct { size_t n; const double * r; } nodes_t;

   /* hat functions of the nodes: two are nonzero at a radius */
   int
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

   nodes_t nodes = { nrad, rnode };
   potsynth_radial radial = { nrad, radial_linear, &nodes };

Other schemes differ only in the functions returned: the Lagrange polynomials of the
nodes for polynomial (barycentric) interpolation, with all functions nonzero; the four
cubic Hermite functions of an interval for coefficient sets that hold values and radial
derivatives, as of a cubic spline; B-splines; or functions that involve no nodes at
all, such as a few analytic profiles. ``tests/test_radialfn.c`` has barycentric, linear
and Hermite functions as further examples.

.. c:function:: void potsynth_free(potsynth_workspace * w)

   Frees a workspace. A NULL pointer is ignored.

.. _vsh-frontend:

Vector spherical harmonic front-end
===================================

Header ``potsynth_vsh.h``. A model holds ``nexp`` vector fields, each

.. math::

   \mathbf{V}(r,\theta,\phi) = \sum_{j=0}^{n_{rad}-1} f_j(r) \sum_{l,m}
   \left[ a^j_{lm} Y_{lm}\,\hat{\mathbf{r}} + b^j_{lm} \boldsymbol{\Psi}_{lm}
        + c^j_{lm} \boldsymbol{\Phi}_{lm} \right],

with the real harmonics :math:`Y_{lm} = P_l^m(\cos\theta)\cos m\phi` for
:math:`m \ge 0` and :math:`P_l^{|m|}(\cos\theta)\sin |m|\phi` for :math:`m < 0`, and the
vector harmonics of Barrera et al. (1985),

.. math::

   \boldsymbol{\Psi}_{lm} = \frac{\partial Y_{lm}}{\partial\theta}\,\hat{\boldsymbol\theta}
     + \frac{1}{\sin\theta}\frac{\partial Y_{lm}}{\partial\phi}\,\hat{\boldsymbol\phi},
   \qquad
   \boldsymbol{\Phi}_{lm} = \hat{\mathbf{r}} \times \boldsymbol{\Psi}_{lm}
     = -\frac{1}{\sin\theta}\frac{\partial Y_{lm}}{\partial\phi}\,\hat{\boldsymbol\theta}
     + \frac{\partial Y_{lm}}{\partial\theta}\,\hat{\boldsymbol\phi}.

The fields need not be gradients of a potential. The radial functions :math:`f_j` are
those of the caller (:ref:`radial-functions`); without them the fields are functions on
the sphere. The model is an ordinary workspace, evaluated with :c:func:`potsynth_eval`
or :c:func:`potsynth_eval_rows`. All fields share one grid, so one window sum per point
gives every field.

.. code-block:: c

   #include <potsynth_vsh.h>

   potsynth_radial radial = { nrad, my_radial_functions, &my_nodes };   /* of the caller */
   potsynth_workspace * w;
   potsynth_eval_workspace * ew;

   w = potsynth_vsh_alloc(1, lmax, lmax, nexp, a, b, c, POTSYNTH_VSH_NORM_SCHMIDT,
                          POTSYNTH_VSH_V, rmin, rmax, &radial, NULL);
   ew = potsynth_eval_alloc(0, w);

   status = potsynth_eval(n, r, theta, phi, V, ew);   /* V: 3*nexp*n doubles */

   potsynth_eval_free(ew);
   potsynth_free(w);

.. c:macro:: POTSYNTH_VSH_V

   The field :math:`(V_r, V_\theta, V_\phi)`.

.. c:macro:: POTSYNTH_VSH_VR

   The radial component :math:`V_r` alone.

.. c:macro:: POTSYNTH_VSH_VH

   The horizontal components :math:`(V_\theta, V_\phi)` alone, for example of fields
   without a radial part. ``POTSYNTH_VSH_V`` is both.

.. c:macro:: POTSYNTH_VSH_DIV

   The surface divergence of the tangential part on the unit sphere,
   :math:`-\sum_j f_j \sum_{lm} l(l+1)\, b^j_{lm} Y_{lm}`. Divide by :math:`r` for the
   divergence on the sphere of radius :math:`r`.

.. c:macro:: POTSYNTH_VSH_CURL

   The radial component of the surface curl of the tangential part on the unit sphere,
   :math:`-\sum_j f_j \sum_{lm} l(l+1)\, c^j_{lm} Y_{lm}`. Divide by :math:`r` as for the
   divergence.

The quantities are bit flags. The outputs of each field are ordered :math:`V_r`,
:math:`V_\theta`, :math:`V_\phi`, divergence, curl: quantity ``q`` of field ``e`` is
output ``e*nout + q``. The divergence and the curl are evaluated as outputs of their
own, which is far more accurate than differencing the evaluated field (see
:ref:`accuracy`).

.. c:macro:: POTSYNTH_VSH_NORM_FULL

   The coefficients refer to 4-:math:`\pi` fully normalized :math:`P_l^m`.

.. c:macro:: POTSYNTH_VSH_NORM_SCHMIDT

   The coefficients refer to Schmidt semi-normalized :math:`P_l^m`.

.. c:function:: size_t potsynth_vsh_nout(const int quantities)

   Number of outputs per point of each field: 1 for the radial component, 2 for the
   horizontal components, 1 each for the divergence and the curl.

.. c:function:: potsynth_workspace * potsynth_vsh_alloc(const size_t lmin, const size_t lmax, const size_t mmax, const size_t nexp, const double * a, const double * b, const double * c, const int norm, const int quantities, const double rmin, const double rmax, const potsynth_radial * radial, const potsynth_parameters * params)

   Builds a model of ``nexp`` fields.

   :param a: coefficients of :math:`Y_{lm}\hat{\mathbf r}`: field ``e`` and radial
      function ``j`` at ``a + (e*nrad + j) * potsynth_nlm(lmin, lmax, mmax)``, indexed
      by :c:macro:`potsynth_idx`, with ``nrad`` the number of radial functions (1
      without); or ``NULL`` for zero
   :param b: coefficients of :math:`\boldsymbol\Psi_{lm}`, as ``a``; or ``NULL``
   :param c: coefficients of :math:`\boldsymbol\Phi_{lm}`, as ``a``; or ``NULL``
   :param norm: normalization of the coefficients
   :param quantities: quantities to compute
   :param rmin: smallest radius of the evaluation points, within the range of the
      radial functions; not used without radial functions
   :param rmax: largest radius of the evaluation points, as ``rmin``
   :param radial: radial functions, or ``NULL`` for fields on the sphere (the radii of
      the evaluation points may then be ``NULL``)
   :param params: parameters, or ``NULL`` for the defaults. ``theta_min`` and
      ``theta_max`` restrict the grid to the colatitudes of the points
   :returns: the model, or ``NULL`` on error

   The coefficients are copied once, into the model, and the arrays of the caller are
   not needed afterwards. An output of an array that is ``NULL`` is zero but is still
   stored in the grid, so request only the quantities that are needed. The accuracy ``eps`` is relative to the largest
   magnitude of each component of each field: of :math:`V_r`, of the horizontal vector
   :math:`(V_\theta, V_\phi)`, of the divergence and of the curl. A weak field in a model
   with strong ones is therefore as accurate, relative to itself, as the strong ones.

   The grid holds ``nexp * nrad * potsynth_vsh_nout(quantities)`` values per node;
   :c:func:`potsynth_get_info` returns its size.

Evaluation
==========

.. c:function:: potsynth_eval_workspace * potsynth_eval_alloc(const size_t chunk, const potsynth_workspace * w)

   Allocates scratch space for evaluating ``w``.

   Points are processed in chunks of at most ``chunk`` points; 0 selects the default,
   262144. :c:func:`potsynth_eval` accepts any number of points regardless. The result
   does not depend on the chunk size.

.. c:function:: void potsynth_eval_free(potsynth_eval_workspace * ew)

   Frees an evaluation workspace. A NULL pointer is ignored.

.. c:function:: int potsynth_eval(const size_t n, const double * r, const double * theta, const double * phi, double * out, potsynth_eval_workspace * ew)

   Evaluates the model at ``n`` points.

   :param n: number of points
   :param r: geocentric radius (the same unit as ``a``), length ``n``
   :param theta: colatitude in radians, :math:`0 \le \theta \le \pi`. Both poles are
      valid.
   :param phi: longitude in radians; any finite value
   :param out: output array of ``nout * n`` doubles, output-major: ``out[k*n + i]``.
      ``nout`` is the number of functionals, or, for the front-ends, the number of
      requested outputs (see :c:func:`potsynth_mag_nout`, :c:func:`potsynth_grav_nout`
      and ``nout`` in :c:struct:`potsynth_info`).
   :param ew: evaluation workspace
   :returns: a status code

   Error handling:

   - **Validation first.** All points are checked before anything is written.
   - **Domain errors.** :c:macro:`POTSYNTH_EDOM` is returned, and the output left
     untouched, if any point has :math:`\theta \notin [0,\pi]` or :math:`r \le 0`;
     for the fast method, if :math:`r` lies outside :math:`[r_{min}, r_{max}]` or
     :math:`\theta` outside the colatitude range of the model; or, with radial functions
     of the caller, if the functions refuse :math:`r`; for the direct method, if
     :math:`(a/r)^{l_{max}+p}` would exceed about :math:`10^{304}` (a point far inside
     the reference sphere at high degree, where the sum would overflow). Colatitudes up to
     :math:`10^{-6}` rad outside :math:`[0,\pi]` are accepted (the ends of a colatitude
     range stay strict, except at a pole),
     so that :math:`\pi` rounded to single precision is; both methods then give the
     smooth continuation of the field across the pole. A radial function of the caller
     that fails after validation gives NaN for its points and :c:macro:`POTSYNTH_EDOM`.
   - **Non-finite inputs.** A NaN or infinite coordinate is not an error; that
     point's outputs are set to NaN.

   The results are bitwise identical for any number of threads and any chunk size.

   ``r`` may be ``NULL`` for a model without radial dependence (see
   :ref:`radial-functions`).

.. c:function:: int potsynth_eval_rows(const size_t n, const double * r, const double * theta, const double * phi, double * out, potsynth_eval_workspace * ew)

   As :c:func:`potsynth_eval`, with the outputs of a point next to each other:
   ``out[i*nout + k]``. With the outputs as the basis functions of a fit, these are the
   rows of its design matrix. For models with many outputs this layout is also the
   faster one to write. The values are those of :c:func:`potsynth_eval` bit for bit, for
   both methods and all front-ends.

Model files
===========

The set-up of a large model of the fast method (a high degree, or many fields on one
grid) takes seconds to minutes and is the same on every run. A model can be written to
a file after its set-up and read instead of being set up again.

.. c:function:: int potsynth_write(const char * filename, const potsynth_workspace * w)

   Writes a model of the fast method: the grid, the radial basis, the window and the
   output map. The file is about as large as the grid
   (:c:member:`potsynth_info.grid_bytes`). Models of the direct method are not written
   (:c:macro:`POTSYNTH_EUNSUP`); their set-up is cheap.

.. c:function:: potsynth_workspace * potsynth_read(const char * filename, const potsynth_radial * radial, const int nthreads)

   Reads a model. ``radial`` must be the radial functions of the caller for a model that
   was set up with them (they are not stored), and ``NULL`` otherwise. ``nthreads`` is
   the number of threads of the evaluation, 0 for all. The model gives the same values
   as the one that was written, bit for bit. Returns ``NULL`` if the file cannot be
   read, is not a model file, is damaged, or does not match ``radial``.

The file is in the binary format of the machine and is not portable between byte
orders. It does not hold the coefficients.

Tensor frames
=============

Gradient tensors are returned as six components per point, in the order (11, 12, 13,
22, 23, 33) of the spherical basis :math:`(\hat{r}, \hat{\theta}, \hat{\phi})`, laid
out like any other output: ``T[k*n + i]``.

.. c:function:: void potsynth_tensor_rtp2nec(const size_t n, double * T)

   Converts a tensor in place to North-East-Center components (NN, NE, NC, EE, EC, CC),
   with :math:`\hat{e}_N = -\hat{\theta}`, :math:`\hat{e}_E = \hat{\phi}`,
   :math:`\hat{e}_C = -\hat{r}`:
   NN = tt, NE = -tp, NC = rt, EE = pp, EC = -rp, CC = rr.

.. c:function:: void potsynth_tensor_rtp2xyz(const size_t n, const double * theta, const double * phi, double * T)

   Converts a tensor in place to Earth-fixed Cartesian components (xx, xy, xz, yy, yz,
   zz): :math:`T^{xyz} = R^T T^{r\theta\phi} R`, where :math:`R` is the rotation from
   Cartesian to spherical components. Unlike the spherical components, these are
   unique at the poles: they do not depend on the longitude given for a polar point.

CPU features
============

.. c:function:: int potsynth_cpu_avx2_fma(void)

   Returns 1 if the processor and operating system support AVX2 and FMA, and 0
   otherwise (always 0 on other architectures). The Python package uses it to choose
   between the AVX2/FMA build and the portable build of the library in the x86-64
   wheels. It is safe to call on any processor.

Introspection
=============

.. c:struct:: potsynth_info

   The parameters chosen during set-up. Its members are:

   - ``lmin``, ``lmax``, ``mmax``, ``nfunc``;
   - ``nout``, the number of outputs per point written by :c:func:`potsynth_eval`. It
     equals ``nfunc``, except for front-end gradient tensors, whose six components are
     computed from five functionals;
   - ``method``, the method actually selected;
   - ``w`` and ``sigma``;
   - ``rho``, the radial rank;
   - ``ncheb``, the number of Chebyshev terms per radial function;
   - ``nf1``, ``nf2``, the fine grid size;
   - ``nrows``, the number of stored rows;
   - ``grid_bytes``;
   - ``ncomp``, the number of irreducible components (``eps`` is relative to each);
   - ``err_radial``, the bound of the error of the radial basis over the shell, relative
     to the size of each component. It is at most ``rank_tol`` (default ``eps/2``) unless
     the tolerance is below the double-precision floor of the bound (about
     :math:`10^{-12}` for EGM2008 at degree 2190), which it then shows; -1 if not known
     (radial functions of the caller, or a model read from a file);
   - ``err_check``, the largest error found by the window check after set-up, against the
     direct sum at the points where the window error is likely to be largest, relative to
     the size of each component. It is at most ``eps`` when the check passed; it can exceed
     ``eps`` near the rounding floor or with a non-default ``sigma``, after the largest
     width or four widenings. -1 if there was no check (window width given, or a model
     read from a file).

   The fast-method fields are zero for the direct method.

.. c:function:: int potsynth_get_info(potsynth_info * info, const potsynth_workspace * w)

   Fills ``info``.

Version
=======

The macros ``POTSYNTH_VERSION_MAJOR``, ``POTSYNTH_VERSION_MINOR`` and
``POTSYNTH_VERSION_PATCH`` give the version. Before 1.0 the interface may still
change.
