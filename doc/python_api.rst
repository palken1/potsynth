.. _python_api:

**********
Python API
**********

.. py:module:: potsynth

The ``potsynth`` package wraps the C library (see :ref:`c_api`). Arrays are passed
without copying, and the Python interpreter lock is released while points are
evaluated. To build and load the library, see :ref:`installation <python-install>`.

.. code-block:: python

   import numpy as np
   import potsynth

   model = potsynth.MagModel((g, h), rmin=6371.2 + 250, rmax=6371.2 + 800, eps=1e-10)
   Br, Btheta, Bphi = model(r, theta, phi)          # nT; theta = colatitude [rad]
   BN, BE, BC = potsynth.rtp_to_nec(Br, Btheta, Bphi)

   grav = potsynth.GravModel.from_gfc("EGM2008.gfc", lmax=720, quantities="V,g",
                                      rmin=6.55e6, rmax=6.9e6)
   field = grav(r, theta, phi)                       # r in m
   field.V, field.g_r                                # m^2/s^2, m/s^2

Magnetic models
===============

.. py:class:: MagModel(gh, lmax=None, lmin=1, mmax=None, a=6371.2, rmin=0.0, rmax=0.0, quantities="B", **options)

   An internal magnetic field model :math:`\mathbf{B} = -\nabla V`, where :math:`V` is
   given by Schmidt semi-normalized Gauss coefficients (:eq:`mag_potential`), and
   optionally its gradient tensor :math:`\nabla\mathbf{B}`.

   :param gh: Gauss coefficients in nT, in one of two forms:

      - a tuple ``(g, h)`` of 2-D arrays indexed ``[l, m]``; or
      - a 1-D vector in the layout of :py:func:`idx` (``g`` at :math:`m \ge 0`, ``h``
        at :math:`m < 0`).
   :param int lmax: maximum degree. By default it is taken from ``gh``.
   :param int lmin: minimum degree, :math:`\ge 1`
   :param int mmax: maximum order; default ``lmax``
   :param float a: reference radius in km
   :param float rmin: smallest radius (km) at which the model will be evaluated.
      Required by the fast method; ignored by the direct method.
   :param float rmax: largest radius (km) at which the model will be evaluated, as for
      ``rmin``
   :param options: see :ref:`python-options`
   :param quantities: ``"B"`` (the field, nT) and/or ``"gradB"`` (the gradient tensor
      :math:`\nabla\mathbf{B} = -\nabla\nabla V`, nT/km), as a list or a
      comma-separated string

   The set-up is done in the constructor. Choose :math:`[r_{\min}, r_{\max}]` as
   tight as the data allow.

   .. py:method:: __call__(r, theta, phi)

      Evaluates the field.

      :param r: radius in km
      :param theta: colatitude in radians, :math:`[0, \pi]`
      :param phi: longitude in radians
      :returns: a named tuple with the requested quantities, each with the broadcast
         shape of the inputs: ``B_r``, ``B_theta``, ``B_phi`` (nT), then ``dB_rr``,
         ``dB_rt``, ``dB_rp``, ``dB_tt``, ``dB_tp``, ``dB_pp`` (nT/km, spherical basis).
         With ``quantities="B"`` it unpacks as ``Br, Btheta, Bphi = model(r, theta, phi)``.

      Inputs may be scalars or arrays of any broadcastable shapes. Points with a NaN
      coordinate give NaN outputs. :py:exc:`DomainError` is raised if a point lies
      outside the valid domain.

   .. py:attribute:: quantities

      The names of the outputs, for example ``("B_r", "B_theta", "B_phi")``.

   .. py:attribute:: info

      A dictionary of the parameters chosen during set-up: ``method``, ``w``,
      ``sigma``, ``rho``, ``ncheb``, ``nf1``, ``nf2``, ``nrows``, ``grid_bytes``,
      ``lmin``, ``lmax``, ``mmax``, ``nfunc`` and ``nout``.

   .. py:classmethod:: from_shc(filename, a=6371.2, rmin=0.0, rmax=0.0, **kwargs)

      Builds a model from a static Swarm SHC file (see :py:func:`read_shc`). The other
      arguments are as for :py:class:`MagModel`, for example
      ``MagModel.from_shc("LCS-1.shc", rmin=a + 250, rmax=a + 800)``.

.. py:function:: read_shc(filename)

   Reads a static model (one epoch, e.g. a lithospheric model such as LCS-1) from a
   Swarm SHC file, using the C reader (see :c:func:`potsynth_shc_read`). Malformed
   files raise :py:exc:`PotsynthError` with the line number. Time-dependent files are
   not supported.

   :returns: :py:class:`ShcModel`

.. py:class:: ShcModel

   The result of :py:func:`read_shc`, with attributes ``lmin``, ``lmax``, ``mmax``,
   ``epoch`` (decimal year), and ``gh``: the Schmidt Gauss coefficients in nT, indexed by
   :py:func:`idx`.

.. py:function:: rtp_to_nec(Br, Btheta, Bphi)

   Returns :math:`(B_N, B_E, B_C) = (-B_\theta, B_\phi, -B_r)`.

Gravity models
==============

.. py:class:: GravModel(coef, GM, a, lmax=None, lmin=0, mmax=None, quantities="V,g", rmin=0.0, rmax=0.0, **options)

   A gravity field model: the gravitational potential :eq:`grav_potential` and the
   acceleration :math:`\mathbf{g} = \nabla V`. Only the synthesis is computed: no
   centrifugal term or normal field, and no tide-system conversion (see
   :c:func:`potsynth_grav_alloc`).

   :param coef: 4-:math:`\pi` fully normalized coefficients, in one of two forms:

      - a tuple ``(C, S)`` of 2-D arrays indexed ``[l, m]``; or
      - a 1-D vector in the layout of :py:func:`idx` (``C`` at :math:`m \ge 0`, ``S``
        at :math:`m < 0`).
   :param float GM: gravitational constant times mass in m\ :sup:`3`/s\ :sup:`2`
   :param float a: reference radius in m
   :param int lmax: maximum degree. By default it is taken from ``coef``.
   :param int lmin: minimum degree, :math:`\ge 0`
   :param int mmax: maximum order; default ``lmax``
   :param quantities: any of ``"V"`` (potential, m\ :sup:`2`/s\ :sup:`2`), ``"g_r"``
      (radial acceleration, m/s\ :sup:`2`), ``"g"`` (acceleration
      :math:`g_r, g_\theta, g_\phi`, m/s\ :sup:`2`) and ``"T"`` (gravity gradient tensor,
      s\ :sup:`-2`; :math:`1\ \mathrm{E} = 10^{-9}` s\ :sup:`-2`), as a list or a
      comma-separated string. Each quantity adds to the set-up time and memory of the
      fast method.
   :param float rmin: smallest radius (m) at which the model will be evaluated.
      Required by the fast method; ignored by the direct method.
   :param float rmax: largest radius (m), as for ``rmin``

   :param options: see :ref:`python-options`; ``eps`` applies to each irreducible component
      separately

   .. py:method:: __call__(r, theta, phi)

      Evaluates the field.

      :param r: radius in m
      :param theta: colatitude in radians, :math:`[0, \pi]`
      :param phi: longitude in radians
      :returns: a named tuple with the requested quantities, in the order ``V``,
         ``g_r``, ``g_theta``, ``g_phi``, ``T_rr``, ``T_rt``, ``T_rp``, ``T_tt``,
         ``T_tp``, ``T_pp``, each with the broadcast shape of the inputs.
         :math:`g_r` is negative (towards the centre) and :math:`g_\theta` is positive
         southwards.

   .. py:attribute:: quantities

      The names of the output quantities, for example ``("V", "g_r", "g_theta",
      "g_phi")``.

   .. py:attribute:: info

      As for :py:attr:`MagModel.info`.

   .. py:classmethod:: from_gfc(filename, lmax=0, quantities="V,g", rmin=0.0, rmax=0.0, **kwargs)

      Builds a model from an ICGEM ``.gfc`` file (see :py:func:`read_gfc`), with the
      file's :math:`GM` and reference radius, truncated at degree ``lmax`` (0 for all
      degrees). The other arguments are as for :py:class:`GravModel`.

.. py:function:: read_gfc(filename, lmax=0)

   Reads a static gravity field model from an ICGEM ``.gfc`` file, using the C reader
   (see :c:func:`potsynth_gfc_read`), truncated at degree ``lmax`` (0 for all degrees).
   Unnormalized files are converted. Malformed files raise :py:exc:`PotsynthError` with
   the line number. Time-variable models are not supported.

   :returns: :py:class:`GfcModel`

.. py:class:: GfcModel

   The result of :py:func:`read_gfc`, with attributes ``lmin`` (0), ``lmax``, ``mmax``,
   ``GM`` (m\ :sup:`3`/s\ :sup:`2`), ``a`` (m), ``modelname``, ``tide_system``, and
   ``coef``: the fully normalized coefficients, indexed by :py:func:`idx`, with absent
   coefficients set to 0.

   To compute disturbing quantities, subtract the normal field's even zonal coefficients
   from ``coef`` before building the model, for example
   ``coef[potsynth.idx(2, 0, 0, lmax)] -= C20_normal``. The normal coefficients must
   first be rescaled to the model's ``GM`` and ``a``.

Geodetic helpers
================

For the recipes of :ref:`gravity-recipes` (see also :ref:`c-geodesy`).

.. py:class:: Ellipsoid(a, f, GM, omega)

   A reference ellipsoid by its defining parameters: semi-major axis (m), flattening,
   geocentric gravitational constant (m\ :sup:`3`/s\ :sup:`2`) and angular velocity
   (rad/s). ``potsynth.WGS84`` and ``potsynth.GRS80`` are predefined.

.. py:function:: geodetic_to_geocentric(lat, h, ellipsoid=WGS84)

   Converts geodetic latitude (radians) and ellipsoidal height (m) to ``(r, theta)``:
   geocentric radius (m) and colatitude (radians).

.. py:function:: normal_zonals(lmax, ellipsoid=WGS84)

   The zonal coefficients ``C[l]`` of the normal gravitational potential, in the
   ellipsoid's own scaling.

.. py:function:: subtract_normal_field(coef, GM, a, lmin=0, lmax=None, mmax=None, ellipsoid=WGS84)

   Returns a copy of the coefficient vector with the normal field subtracted: the
   coefficients of the disturbing potential :math:`T`.

.. py:function:: normal_gravity(lat, ellipsoid=WGS84)

   Normal gravity (m/s\ :sup:`2`) on the ellipsoid at geodetic latitude ``lat``
   (radians).

Gradient tensors
================

``MagModel(..., quantities="gradB")`` and ``GravModel(..., quantities="T")`` return the
gradient tensor as six components in the spherical basis (see :ref:`gradient-tensor`).
These functions convert them; all accept arrays of any broadcastable shapes.

.. code-block:: python

   field = model(r, theta, phi)                       # MagModel, quantities="B,gradB"
   dB = field[3:]                                     # dB_rr, dB_rt, ..., dB_pp
   NN, NE, NC, EE, EC, CC = potsynth.tensor_rtp_to_nec(*dB)
   M = potsynth.tensor_matrix(*dB)                    # shape (..., 3, 3)

.. py:function:: tensor_rtp_to_nec(rr, rt, rp, tt, tp, pp)

   Returns the North-East-Center components (NN, NE, NC, EE, EC, CC), with
   :math:`\hat{e}_N = -\hat{\theta}`, :math:`\hat{e}_E = \hat{\phi}`,
   :math:`\hat{e}_C = -\hat{r}` (see :c:func:`potsynth_tensor_rtp2nec`).

.. py:function:: tensor_rtp_to_xyz(theta, phi, rr, rt, rp, tt, tp, pp)

   Returns the Earth-fixed Cartesian components (xx, xy, xz, yy, yz, zz) at colatitude
   ``theta`` and longitude ``phi`` (radians). Unlike the spherical components, these are
   unique at the poles (see :c:func:`potsynth_tensor_rtp2xyz`).

.. py:function:: tensor_matrix(t11, t12, t13, t22, t23, t33)

   Assembles the six components into symmetric matrices of shape ``(..., 3, 3)``.

Generic models
==============

.. py:class:: Functional(rpow, theta="P", dphi=0, rfac=None)

   One functional of a generic model, as in :eq:`functional`.

   :param int rpow: radial power offset :math:`p`; the term carries
      :math:`(a/r)^{l+p}`
   :param str theta: the operator in colatitude: ``"P"``, ``"dP"``
      (:math:`d\bar{P}/d\theta`), ``"Psin"`` (:math:`\bar{P}/\sin\theta`), or
      ``"S2A"``, ``"S2B"`` (the second-derivative operators of the gradient tensor, see
      :c:macro:`POTSYNTH_THETA_S2A`)
   :param int dphi: 1 to apply :math:`\partial/\partial\phi`
   :param rfac: per-degree factors :math:`\rho(l)`, length :math:`L+1`, or ``None``
      for 1

.. py:class:: Model(coef, lmin, lmax, functionals, a, rmin=0.0, rmax=0.0, mmax=None, **options)

   A model with user-defined functionals. ``coef`` holds the 4-:math:`\pi` fully
   normalized coefficients as a vector (see :py:func:`coef_vector`), and ``lmin``
   may be 0. ``a``, ``rmin`` and ``rmax`` are as for :py:class:`MagModel`, in the
   unit of the radii; the options are described in :ref:`python-options`.

   Calling the model with ``(r, theta, phi)`` returns an array of shape
   ``(len(functionals),) + shape``. :py:attr:`MagModel.info` is also available.

   Example: a potential :math:`V/a` and its radial derivative:

   .. code-block:: python

      l = np.arange(L + 1)
      funcs = [potsynth.Functional(rpow=1),                        # sum (a/r)^(l+1) ...
               potsynth.Functional(rpow=2, rfac=l + 1.0)]          # -dV/dr
      model = potsynth.Model(coef, 0, L, funcs, a, rmin, rmax)
      V, dV = model(r, theta, phi)

.. py:class:: SetsModel(coef, lmin, lmax, outputs, a=1.0, rmin=None, rmax=None, mmax=None, radial=None, **options)

   A model with several coefficient sets (see :c:func:`potsynth_alloc_sets`): each
   output is a sum of terms, and each term applies one operator to one set. All outputs
   share one grid, so many fields, or fields that are not gradients of one potential,
   are evaluated with one sum per point.

   :param coef: array of shape ``(nset, nlm)``: the 4-:math:`\pi` fully normalized
      coefficients of each set in the layout of :py:func:`idx`
   :param outputs: list of :py:class:`Output`
   :param radial: radial functions of the caller (:py:class:`Radial`) in place of the
      powers of :math:`a/r`, or ``None``. With radial functions a term's ``set`` is the
      first of ``radial.n`` consecutive sets, one for each function, and ``a`` and
      ``rpow`` are not used
   :param rmin, rmax: radius range of the points; with radial functions by default
      theirs

   Calling the model returns an array of shape ``(len(outputs),) + shape``.

.. py:class:: Output(terms, rpow=0, component=0, weight=0.0)

   One output: :math:`(a/r)^{l+p}` times the sum of its terms. The terms must form one
   field component: all with the same operator, or ``("dP", dphi=0)`` and ``("Psin",
   dphi=1)`` (a horizontal vector), or, with ``component`` and ``weight`` given,
   ``("S2A", dphi=0)`` and ``("S2B", dphi=1)`` (a horizontal tensor). Outputs with the
   same ``component`` :math:`> 0` form one component for the accuracy; 0 lets the
   library choose. ``weight`` is the weight of the output in the magnitude of its
   component (0: chosen by the library).

.. py:class:: Term(set, theta="P", dphi=0, fac=1.0, rfac=None)

   One term: ``fac * rfac[l]`` times the operator ``theta`` (as for :py:class:`Functional`)
   applied to coefficient set ``set``.

   Example: the horizontal components of a vector field
   :math:`\sum b_{lm}\boldsymbol\Psi_{lm} + c_{lm}\boldsymbol\Phi_{lm}` with ``b`` and
   ``c`` as sets 0 and 1:

   .. code-block:: python

      T, O = potsynth.Term, potsynth.Output
      outputs = [O([T(0, "dP"), T(1, "Psin", dphi=1, fac=-1.0)]),     # V_theta
                 O([T(0, "Psin", dphi=1), T(1, "dP")])]               # V_phi
      model = potsynth.SetsModel(np.stack([b, c]), 1, L, outputs, a, rmin, rmax)

Vector fields
=============

.. py:class:: VshModel(a=None, b=None, c=None, lmax=None, lmin=1, mmax=None, norm="schmidt", quantities="V", radial=None, rmin=None, rmax=None, **options)

   Vector fields given by vector spherical harmonic coefficients, which need not be
   gradients of a potential (see :ref:`vsh-frontend` for the definition of the fields
   and of :math:`\boldsymbol\Psi_{lm}`, :math:`\boldsymbol\Phi_{lm}`). A model holds
   any number of fields, which share one grid of the fast method.

   :param a: coefficients of :math:`Y_{lm}\hat{\mathbf r}`, ``None`` (zero) or an array
      of shape ``(nexp, nrad, nlm)``: ``nexp`` fields, ``nrad`` radial functions, and
      ``nlm = nlm(lmin, lmax, mmax)`` coefficients in the layout of :py:func:`idx`.
      ``(nexp, nlm)`` is accepted for ``nrad = 1`` and ``(nlm,)`` for one field
   :param b: coefficients of :math:`\boldsymbol\Psi_{lm}`, as ``a``
   :param c: coefficients of :math:`\boldsymbol\Phi_{lm}`, as ``a``
   :param str norm: ``"schmidt"`` or ``"full"`` normalization of the coefficients
   :param quantities: any of ``"V"`` (:math:`V_r, V_\theta, V_\phi`), ``"V_r"``,
      ``"V_h"`` (:math:`V_\theta, V_\phi`), ``"div"`` and ``"curl"`` (surface divergence
      and radial component of the surface curl of the tangential part on the unit
      sphere; divide by :math:`r`), as a list or a comma-separated string
   :param radial: radial functions (a :py:class:`Radial`), or ``None`` for fields on
      the sphere
   :param rmin, rmax: radius range of the points; by default that of ``radial``
   :param options: see :ref:`python-options`; ``eps`` is relative to the largest
      magnitude of each component of each field (:math:`V_r`; the horizontal vector;
      the divergence; the curl)

   Calling the model with ``(r, theta, phi)`` returns an array of shape
   ``(nexp, nout) + shape``, with the outputs of each field in the order
   :math:`V_r, V_\theta, V_\phi`, divergence, curl; :py:attr:`outputs` names them. For
   fields on the sphere ``r`` is not used and may be ``None``.

   .. py:attribute:: outputs

      Names of the outputs of each field.

   .. py:attribute:: nexp

      Number of fields.

   Example: 300 fields given at 8 Chebyshev radii, as a design matrix:

   .. code-block:: python

      radial = potsynth.Radial(8, weights, rmin=rlo, rmax=rhi)   # the caller's 8 functions
      model = potsynth.VshModel(a, b, c, radial=radial,     # arrays (300, 8, nlm)
                                theta_max=np.radians(50.0)) # points on a polar cap
      G = model.rows(r, theta, phi)                         # (npoints, 900)

Radial functions
----------------

By default the radial dependence of a model is the power of :math:`a/r` of each
degree. A :py:class:`VshModel` takes its radial functions from the caller (see
:ref:`radial-functions`).

.. py:class:: Radial(n, func, params=None, rmin=None, rmax=None)

   ``n`` radial functions :math:`f_j(r)`: a model that uses them has ``n`` coefficient
   sets for each field.

   :param func: the functions, as a Python function or a compiled function (below)
   :param params: for a compiled function, what it receives as ``params``: ``None``, an
      integer address, a ctypes pointer or structure, or a NumPy array (its data)
   :param rmin, rmax: range of radii in which the functions are defined; the default
      radius range of a model that uses them

   ``func=None`` with ``n = 1`` means no radial dependence: fields on the sphere, for
   which the radii of the points are not used and may be ``None``.

   **A Python function** ``func(r)`` returns ``(first, values)``: the values of the
   functions ``first, first + 1, ...`` that are nonzero at the radius ``r``; all others
   are zero. It returns ``None`` if ``r`` is outside the range of the functions. An
   exception it raises is raised by the evaluation.

   .. code-block:: python

      def hats(r):                                   # linear interpolation between nodes
          if not nodes[0] <= r <= nodes[-1]:
              return None
          k = min(np.searchsorted(nodes, r, side="right") - 1, len(nodes) - 2)
          t = (r - nodes[k]) / (nodes[k + 1] - nodes[k])
          return k, (1.0 - t, t)

      radial = potsynth.Radial(len(nodes), hats, rmin=nodes[0], rmax=nodes[-1])
      model = potsynth.VshModel(a, b, c, radial=radial)       # arrays (nexp, len(nodes), nlm)

   A Python function is called for every evaluation point, one call at a time whatever
   the number of threads, and costs 10 to 20 microseconds per point, against well
   under a microsecond for the rest of the evaluation. It is suited to trying a scheme
   out and to thousands of points.

   **A compiled function** has the C signature of :c:type:`potsynth_radial_fn` and runs
   at the speed of C, from several threads at once. It is given as a ctypes function
   pointer, as an object with an ``address`` attribute, or as an integer address. With
   numba the Python function above becomes a compiled one:

   .. code-block:: python

      from numba import cfunc, types, carray

      sig = types.intc(types.double, types.CPointer(types.uintp), types.CPointer(types.uintp),
                       types.CPointer(types.double), types.CPointer(types.double), types.voidptr)
      x = tuple(nodes)

      @cfunc(sig)
      def hats(r, first, n, f, df, params):
          if not (r >= x[0] and r <= x[-1]):
              return 3                               # POTSYNTH_EDOM
          k = 0
          while k + 2 < len(x) and r >= x[k + 1]:
              k += 1
          t = (r - x[k]) / (x[k + 1] - x[k])
          out = carray(f, 2)
          out[0] = 1.0 - t
          out[1] = t
          first[0] = k
          n[0] = 2
          return 0

      radial = potsynth.Radial(len(nodes), hats, rmin=nodes[0], rmax=nodes[-1])

   A function in a C or Cython extension is passed in the same way, by its address,
   with its node positions or other data as ``params``.

.. _python-regions:

Rows, regions and model files
=============================

These apply to every model (:py:class:`MagModel`, :py:class:`GravModel`,
:py:class:`Model`, :py:class:`VshModel`).

.. py:method:: Model.rows(r, theta, phi)

   Evaluates the model with the outputs of a point next to each other: an array of
   shape ``shape + (nout,)``. With the outputs as the basis functions of a fit, these
   are the rows of its design matrix. The values are those of calling the model. For
   models with many outputs this is also the faster layout.

**Colatitude range.** The arguments ``theta_min`` and ``theta_max`` (radians;
``theta_max = 0`` means :math:`\pi`) give the colatitude range of the evaluation
points. The fast method then stores only the part of its grid that the range needs,
with the same values, and a point outside the range raises
:py:class:`DomainError`. Use it for data on a polar cap or in a region.

.. py:method:: Model.write(filename)

   Writes a model of the fast method to a file, so that it can be read instead of
   being set up again. The file is about as large as the grid
   (``info["grid_bytes"]``), is in the binary format of the machine, and does not hold
   the coefficients.

.. py:function:: read(filename, radial=None, nthreads=0)

   Reads a model and returns a :py:class:`LoadedModel`. ``radial`` must be the radial
   functions of a model that was set up with them.

.. py:class:: LoadedModel

   A model read from a file. Calling it with ``(r, theta, phi)`` gives an array of shape
   ``(nout,) + shape`` with the outputs in the order of the model that was written;
   :py:meth:`Model.rows` and ``info`` are available.

.. _python-options:

Options
=======

All model classes accept these keyword arguments:

``eps``
   Requested accuracy, relative to the largest magnitude of each field component (a
   scalar, or a horizontal vector or tensor; see :ref:`accuracy`). Default
   :math:`10^{-10}`.
``method``
   ``"auto"`` (default), ``"fast"`` or ``"direct"``.
``npts_hint``
   Expected number of evaluation points. With ``method="auto"``, the method with the
   smaller total time (set-up plus evaluation) is chosen; 0 means unknown (see
   :ref:`choosing-method`).
``nthreads``
   OpenMP threads; 0 means all.
``sigma``, ``w``, ``rank_tol``
   Expert overrides of the oversampling factor, the window width and the radial SVD
   tolerance; 0 means automatic.
``chunk``
   Points per internal chunk; 0 means the default.
``theta_min``, ``theta_max``
   Colatitude range of the evaluation points in radians (``theta_max = 0`` means
   :math:`\pi`); see :ref:`python-regions`.

An unknown keyword raises ``TypeError``.

Coefficients
============

.. py:function:: coef_vector(C, S, lmin=None, lmax=None, mmax=None)

   Builds a coefficient vector from 2-D arrays ``C[l, m]`` (cosine) and ``S[l, m]``
   (sine). ``S`` may be ``None``. By default ``lmin`` is 0 if ``C[0, 0]`` is nonzero and
   1 otherwise.

.. py:function:: idx(l, m, lmin, mmax)

   Index of coefficient :math:`(l, m)` in a coefficient vector. :math:`m < 0` gives the
   sine coefficient of order :math:`|m|`. This is the same as :c:macro:`potsynth_idx`.

.. py:function:: nlm(lmin, lmax, mmax=None)

   Length of a coefficient vector.

Exceptions
==========

.. py:exception:: PotsynthError

   An error reported by the C library. The attribute ``status`` holds the status code.

.. py:exception:: DomainError

   Subclass of :py:exc:`PotsynthError` and :py:exc:`ValueError`. Raised when an
   evaluation point lies outside the valid domain: :math:`\theta \notin [0,\pi]`,
   :math:`r \le 0`, or, for the fast method, :math:`r \notin [r_{\min}, r_{\max}]`.

Threads
=======

A model object holds one evaluation workspace, and calls on the same object are
serialized by a lock. Each call already uses all the threads given by ``nthreads``.
To evaluate from several Python threads at once, create one model object per thread.
