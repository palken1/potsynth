# Changelog

All notable changes to potsynth are listed here. Versions follow
[semantic versioning](https://semver.org); before 1.0 the interface may still change
between minor versions.

## 0.1.0 (unreleased)

Terminology as in the paper (Alken and Burchill): a *functional* is one linear functional
of the potential requested from the core, and an *irreducible component* is the set of
functionals whose accuracy is controlled together. The front-ends keep `quantities` for
the physical quantities they compute.

First public release.

### Core

- Fast synthesis of internal potential-field models at scattered points: after a one-off
  set-up, the cost per point does not depend on the maximum degree L. Accuracy can be set
  from about 1e-4 down to the double-precision rounding floor.
- Optimized direct sum (fused single-pass kernels) for low degrees and small point sets;
  `POTSYNTH_METHOD_AUTO` chooses between the two, including the set-up cost when
  `npts_hint` gives the expected number of points.
- Generic functionals (linear functionals of the potential), including the theta operators
  for gradient tensors. Results are exact at both poles: no division by sin(theta).
- Radial basis weighted by exact mean squares over the sphere (vector spherical
  harmonics): functionals form irreducible components (scalar, horizontal vector, spin-2 tensor part),
  and `eps` is relative to the largest magnitude of each irreducible component, e.g. B_theta and
  B_phi relative to max |B_h|. One radial basis bounds the error of every irreducible component: its
  rank is raised until each irreducible component's own error is within tolerance, so an irreducible component much
  smaller than the others at high degree (e.g. the horizontal gravity gradients of
  EGM2008) is resolved as well as they are.
- The radial step bounds the maximum of each irreducible component's error over the sphere, not its
  root mean square: by the addition theorem a degree-l part is at most sqrt(2l+1) times
  its rms anywhere. The rank is the smallest whose bound, measured for every rank at once
  on a fine set of radii (between the radial nodes too), is within eps/2 (default
  `rank_tol`, was 0.1 eps); the Chebyshev degree is raised by about a quarter while that
  lowers the rank. Below the rounding floor of the bound the basis is taken at the knee
  of the bound, and `err_radial` reports it. This fixes errors of up to 50 eps for EGM2008
  on Earth's topographic surface (shells reaching inside r = a); satellite shells need
  one or two more radial functions.
- Window width chosen per irreducible component and shell: the table measured for the magnetic field
  is a lower bound, and the width is raised where an irreducible component has its power near the band
  edge (where the window's aliasing error is largest; e.g. the horizontal gravity
  gradients) or the shell reaches inside r = a (u > 1), so that every irreducible component meets
  `eps`. Requests like the magnetic field at satellite altitude keep the table's width.
  After the set-up, an automatically chosen width is checked against the direct sum at
  the points where its error is likely to be largest, and widened if any irreducible component
  exceeds `eps`.
- Read-only model workspace shared between threads; evaluation does not allocate memory.
  OpenMP parallelism (optional).
- Models with several coefficient sets (`potsynth_alloc_sets`): an output is a sum of
  terms, each one operator applied to one set, and all outputs share one grid. Fields
  that are not gradients of one potential, and many fields at once, are evaluated with
  one window sum per point. `potsynth_alloc` is the case of one set.
- Radial functions of the caller in place of the powers of a/r (`potsynth_radial`): a
  function that returns the radial basis functions that are nonzero at a radius, for
  fields tabulated at radial nodes with any linear interpolation rule, splines or
  analytic profiles. Only the functions that are nonzero at a point are summed. The
  library has no interpolation of its own; the tests have examples.
- Field components for the accuracy can be given by the caller (`component` and `weight`
  of an output); `ncomp` in `potsynth_info`.
- `potsynth_eval_rows`: the outputs of a point next to each other (rows of a design
  matrix), for both methods and all front-ends.
- Model files (`potsynth_write`, `potsynth_read`): a model of the fast method is written
  after its set-up and read instead of being set up again.
- The window check after the set-up tests each component against a direct sum of its
  own functionals, so its cost does not grow with the number of fields; the copy of the
  coefficients is freed after the set-up of a fast model.
- `err_radial` and `err_check` in `potsynth_info`: the error bound of the radial basis
  and the largest error of the window check after set-up, so that a tolerance below
  the rounding floor (or a non-default `sigma`) no longer passes unnoticed.
- In a model with several coefficient sets, the parts of a functional with the same
  operator are combined before their power is taken, so that sets that nearly cancel
  (e.g. the difference of two epochs) get the radial basis and window of their sum, not
  of each set.
- Window check of a model with more than 16 functionals: 16 candidate nodes per
  irreducible component (was as few as 2 per functional), each component checked against
  the direct sum of its own functionals at its own points.
- Radial functions of the caller are called once per point: validation checks the radius
  against [rmin, rmax] of the model, where the functions must be valid (they were called
  in validation and again in the evaluation).
- The FFTs of the set-up use the header-only C++ pocketfft (BSD, vendored in
  `extern/pocketfft_cpp`, called from one C++ file) instead of the C pocketfft: the many
  transforms of an order or a grid row are done in one call, vectorized across them, and
  the phi pass transforms each row in place. Building needs a C++11 compiler too.
- Robustness after a code review: the radial step returns the best basis of its
  refinement passes and never returns without one (after 8 passes it could); out-of-memory in the radial
  check and the window check is an error instead of a pass; the range check of the
  caller's radial functions cannot overflow; the window is no longer widened again when
  a widening did not halve the error; `potsynth_vsh_nout` returns 0 for invalid flags;
  `potsynth_read` accepts NULL for a model without radial dependence (a VSH model on one
  sphere; the model file records this); pkg-config names libc++ only for Apple's clang.
- The radial step computes each power u^(l+p) once per degree instead of once per radial
  function, and copies the coefficients of functionals with the same radial power.
- Internal cleanup: one table of operator terms built at allocation and used by both
  methods (it was built three times); one routine for the weight rows of a term (it was
  written in both methods); shared helpers for repeated checks and formulas; "layer"
  for one functional and radial function of the grid throughout (was channel/pair);
  debugging output removed.
- Set-up skips the theta-derivative column of the Legendre recursion when no functional
  uses it (potentials and radial components alone).
- Robustness: functions that need the coefficients of a model refuse a fast workspace,
  whose coefficients are freed after set-up (they read freed memory);
  the free functions accept NULL, as `free()` does; colatitudes up to 1e-6 rad
  outside [0, pi] are accepted (pi in single precision), with the smooth continuation of
  the field across the pole; the direct method refuses a point where (a/r)^lmax would overflow
  (it gave NaN); CMake looks for Homebrew's libomp on macOS; `sigma` below 1 is rejected (it overflowed the grid); `potsynth_read`
  checks the sizes in the file header (window inside the stored rows and columns, no
  overflow); a radial function of the caller that fails after the points were validated
  gives NaN for those points and an error, in C and Python (the exception of a Python
  function is raised by the call in which it occurred); a window check that cannot be
  made fails the set-up instead of passing it.
- Colatitude range of the evaluation points (`theta_min`, `theta_max` of the
  parameters): the fast method stores only the grid rows of the range, with identical
  results; for data in a polar cap or a region.
- Interface for the release: `potsynth_version()` and `potsynth_abi()` report the
  library's version and the layout of its public structures (`POTSYNTH_ABI`), and the
  Python package refuses a library of another layout; `potsynth_parameters` and
  `potsynth_info` end in reserved space (zeroed by `potsynth_default_parameters`), so
  that fields can be added without changing their size;
  `potsynth_eval_set_threads` sets the threads of one evaluation workspace (0: those of
  the model).

### Front-ends

- Magnetic (`potsynth_mag.h`): B and its gradient tensor grad B from Schmidt
  semi-normalized Gauss coefficients; Swarm SHC reader (e.g. LCS-1), with truncation at
  a given degree as the `.gfc` reader.
- Gravity (`potsynth_grav.h`): potential V, acceleration g and gravity gradient tensor
  T from fully normalized coefficients; ICGEM `.gfc` reader (EGM2008, EIGEN, XGM, ...).
- Vector spherical harmonics (`potsynth_vsh.h`): any number of vector fields
  sum f_j(r) [a Y r + b Psi + c Phi] with radial functions of the caller; the field,
  its surface divergence and the radial component of its curl; Schmidt or fully
  normalized coefficients. The model is an ordinary workspace. `potsynth_vsh_alloc`
  takes its arguments in the order of the other constructors (coefficients, radial
  range, radial functions, quantities, parameters).
- Geodetic helpers: WGS84/GRS80 ellipsoids, geodetic to geocentric coordinates, normal
  field zonal coefficients, normal gravity. The "Gravity recipes" chapter of the
  documentation shows how to compute disturbances, anomalies, height anomalies,
  deflections and geoid heights.

### Python

- `potsynth` package (ctypes and NumPy): `MagModel`, `GravModel`, the generic `Model`,
  readers, geodetic helpers and tensor rotations.
- `VshModel`, with radial functions of the caller (`Radial`) given as a Python function
  or as a compiled one (ctypes, numba, Cython); `SetsModel` with `Output` and `Term`
  for general models with several coefficient sets; the options of all models
  (`eps`, `method`, ..., `theta_min`, `theta_max`) are keyword arguments handled in
  one place, and an unknown one is an error;
  `rows()` (the outputs of a point next to each other), the colatitude range
  (`theta_min`, `theta_max`) and model files (`write()`, `potsynth.read()`) for all models.
- The same attributes and defaults in all models: `outputs` names the outputs (also of
  `Model` and `SetsModel`, with `names=`), `quantities` holds the requested quantities
  of the front-ends, `rmin` and `rmax` default to None (direct method only, or the
  range of the radial functions), and a coefficient vector of all degrees up to a
  higher degree is truncated at `lmax`. `write()` keeps the output names next to the
  model file (`filename.json`), and `read()` restores them. `SetsModel` needs `a`
  unless radial functions are given (it defaulted to 1).
- One model can be evaluated from several Python threads at once: each thread gets its
  own evaluation workspace; `eval_threads` sets their OpenMP threads.
- `potsynth.__version__`, and a check at import that the library has the layout of the
  bindings.
- Wheels for Linux (x86-64, aarch64), macOS (arm64, x86-64) and Windows (AMD64). The
  x86-64 wheels contain an AVX2/FMA build and a portable build of the library and select
  one at import time.

### Platforms

- Tested on Linux (GCC, x86-64 and aarch64), macOS (Apple clang, arm64) and Windows
  (MSVC), as static and shared libraries, and with AddressSanitizer, UBSan and
  ThreadSanitizer.
