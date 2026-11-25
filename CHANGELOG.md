# Changelog

All notable changes to potsynth are listed here. Versions follow
[semantic versioning](https://semver.org); before 1.0 the interface may still change
between minor versions.

## 0.1.0 (unreleased)

First public release.

### Core

- Fast synthesis of internal potential-field models at scattered points: after a one-off
  set-up, the cost per point does not depend on the maximum degree L. Accuracy can be set
  from about 1e-4 down to the double-precision rounding floor.
- Optimized direct sum (fused single-pass kernels) for low degrees and small point sets;
  `POTSYNTH_METHOD_AUTO` chooses between the two, including the set-up cost when
  `npts_hint` gives the expected number of points.
- Generic channels (linear functionals of the potential), including the theta operators
  for gradient tensors. Results are exact at both poles: no division by sin(theta).
- Radial basis weighted by exact mean squares over the sphere (vector spherical
  harmonics): channels form quantities (scalar, horizontal vector, spin-2 tensor part),
  and `eps` is relative to the largest magnitude of each quantity, e.g. B_theta and
  B_phi relative to max |B_h|. One radial basis bounds the error of every quantity: its
  rank is raised until each quantity's own error is within tolerance, so a quantity much
  smaller than the others at high degree (e.g. the horizontal gravity gradients of
  EGM2008) is resolved as well as they are.
- The radial step bounds the maximum of each quantity's error over the sphere, not its
  root mean square: by the addition theorem a degree-l part is at most sqrt(2l+1) times
  its rms anywhere. The rank and the Chebyshev degree are chosen with this bound within
  eps/2 (default `rank_tol`, was 0.1 eps), checked between the radial nodes too, with
  rounding-floor tests at high accuracy (applied only near the rounding level, so that
  thick shells at high degree, e.g. EGM2008 on 0-2000 km, get their full rank). This fixes errors of up to 50 eps for EGM2008
  on Earth's topographic surface (shells reaching inside r = a); satellite shells need
  one or two more radial functions.
- Window width chosen per quantity and shell: the table measured for the magnetic field
  is a lower bound, and the width is raised where a quantity has its power near the band
  edge (where the window's aliasing error is largest; e.g. the horizontal gravity
  gradients) or the shell reaches inside r = a (u > 1), so that every quantity meets
  `eps`. Requests like the magnetic field at satellite altitude keep the table's width.
  After the set-up, an automatically chosen width is checked against the direct sum at
  the points where its error is likely to be largest, and widened if any quantity
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
  own channels, so its cost does not grow with the number of fields; the copy of the
  coefficients is freed after the set-up of a fast model.
- `err_radial` and `err_check` in `potsynth_info`: the error bound of the radial basis
  and the largest error of the window check after set-up, so that a tolerance below
  the rounding floor (or a non-default `sigma`) no longer passes unnoticed.
- In a model with several coefficient sets, the parts of a channel with the same
  operator are combined before their power is taken, so that sets that nearly cancel
  (e.g. the difference of two epochs) get the radial basis and window of their sum, not
  of each set.
- Robustness: the free functions accept NULL, as `free()` does; colatitudes up to 1e-6 rad
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

### Front-ends

- Magnetic (`potsynth_mag.h`): B and its gradient tensor grad B from Schmidt
  semi-normalized Gauss coefficients; Swarm SHC reader (e.g. LCS-1).
- Gravity (`potsynth_grav.h`): potential V, acceleration g and gravity gradient tensor
  T from fully normalized coefficients; ICGEM `.gfc` reader (EGM2008, EIGEN, XGM, ...).
- Vector spherical harmonics (`potsynth_vsh.h`): any number of vector fields
  sum f_j(r) [a Y r + b Psi + c Phi] with radial functions of the caller; the field,
  its surface divergence and the radial component of its curl; Schmidt or fully
  normalized coefficients. The model is an ordinary workspace.
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
- Wheels for Linux (x86-64, aarch64), macOS (arm64, x86-64) and Windows (AMD64). The
  x86-64 wheels contain an AVX2/FMA build and a portable build of the library and select
  one at import time.

### Platforms

- Tested on Linux (GCC, x86-64 and aarch64), macOS (Apple clang, arm64) and Windows
  (MSVC), as static and shared libraries, and with AddressSanitizer, UBSan and
  ThreadSanitizer.
