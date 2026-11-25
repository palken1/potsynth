"""
Generic potential-field synthesis: channels of an internal potential.
"""
import ctypes
import os
import threading
from dataclasses import dataclass

import numpy as np

from . import _lib
from ._lib import lib

_METHODS = {"auto": _lib.METHOD_AUTO, "fast": _lib.METHOD_FAST, "direct": _lib.METHOD_DIRECT}
_METHOD_NAMES = {v: k for k, v in _METHODS.items()}
_THETA_OPS = {"P": _lib.THETA_P, "dP": _lib.THETA_DP, "Psin": _lib.THETA_PSIN,
              "S2A": _lib.THETA_S2A, "S2B": _lib.THETA_S2B}


class PotsynthError(RuntimeError):
    """Error reported by the potsynth library."""

    def __init__(self, status, what="", reason=""):
        self.status = status
        detail = f"{reason} ({_lib.strerror(status)})" if reason else _lib.strerror(status)
        super().__init__(f"{what}{': ' if what else ''}{detail}")


class DomainError(PotsynthError, ValueError):
    """An evaluation point lies outside the valid domain (theta outside [0, pi],
    r <= 0, or, for the fast method, r outside [rmin, rmax])."""


def _raise(status, what):
    """Raise the first error recorded by the library error handler (or status)."""
    status, reason = _lib.last_error(status)
    if status == _lib.POTSYNTH_EDOM:
        raise DomainError(status, what, reason)
    raise PotsynthError(status, what, reason)


def nlm(lmin, lmax, mmax=None):
    """Length of a coefficient vector for degrees lmin..lmax and orders up to mmax."""
    mmax = lmax if mmax is None else mmax
    return int(lib.potsynth_nlm(lmin, lmax, mmax))


def idx(l, m, lmin, mmax):
    """Index of coefficient (l, m) in a coefficient vector, stored degree by degree
    with orders -M..M within each degree (M = min(l, mmax)): m >= 0 indexes the
    cosine coefficient, m < 0 the sine coefficient of order |m|."""
    if lmin <= mmax:
        if l <= mmax:
            return l * (l + 1) + m - lmin * lmin
        return l * (2 * mmax + 1) + m - mmax * mmax - lmin * lmin
    return (l - lmin) * (2 * mmax + 1) + m + mmax


def coef_vector(C, S, lmin=None, lmax=None, mmax=None):
    """Build a coefficient vector from arrays C[l, m] (cosine) and S[l, m] (sine),
    for example Gauss coefficients g[l, m], h[l, m].

    lmin defaults to 0 if C[0, 0] != 0 and to 1 otherwise; lmax to the largest
    degree in the arrays; mmax to lmax.
    """
    C = np.asarray(C, dtype=float)
    S = np.zeros_like(C) if S is None else np.asarray(S, dtype=float)
    L = C.shape[0] - 1
    lmax = L if lmax is None else lmax
    lmin = (0 if C[0, 0] != 0.0 else 1) if lmin is None else lmin
    mmax = lmax if mmax is None else mmax
    out = np.zeros(nlm(lmin, lmax, mmax))
    for l in range(lmin, lmax + 1):
        M = min(l, mmax)
        out[idx(l, 0, lmin, mmax)] = C[l, 0]
        for m in range(1, M + 1):
            out[idx(l, m, lmin, mmax)] = C[l, m]
            out[idx(l, -m, lmin, mmax)] = S[l, m]
    return out


@dataclass
class Channel:
    """One output quantity of the core model,

        f(r, theta, phi) = sum_{l,m} rfac[l] (a/r)^(l + rpow) T[Pbar_l^m(cos theta)]
                                     (d/dphi)^dphi [C_lm cos(m phi) + S_lm sin(m phi)]

    theta: "P" (Pbar), "dP" (d Pbar / d theta), "Psin" (Pbar / sin theta, m >= 1),
           "S2A" or "S2B" (the second-derivative operators of the gradient tensor).
    rfac: per-degree factors, length lmax + 1, or None for 1.
    """
    rpow: int
    theta: str = "P"
    dphi: int = 0
    rfac: object = None


def _params(eps, method, nthreads, sigma, w, rank_tol, npts_hint=0, theta_min=0.0, theta_max=0.0):
    p = lib.potsynth_default_parameters()
    p.npts_hint = npts_hint
    p.theta_min = theta_min
    p.theta_max = theta_max
    p.eps = eps
    p.method = _METHODS[method]
    p.nthreads = nthreads
    p.sigma = sigma
    p.w = w
    p.rank_tol = rank_tol
    return p


class Radial:
    """Radial functions of a model, in place of the powers of a/r.

    A model with n radial functions f_j(r) has n coefficient sets for each of its
    fields, and its fields are sum_j f_j(r) x (the field of set j). The radial
    functions are those of the caller: interpolation between radial nodes of any
    kind, splines, analytic profiles, or anything else that is linear in the
    coefficients.

    Parameters
    ----------
    n : number of radial functions
    func : the functions, as one of

        - a Python function ``func(r)`` that returns ``(first, values)``: the values
          of the functions first, first + 1, ... that are nonzero at the radius r (all
          others are zero), or None if r is outside of the range of the functions.
          It is called for every evaluation point, one call at a time, so it is slow:
          suited to trying things out and to thousands of points, not millions;
        - a compiled function with the C signature
          ``int f(double r, size_t *first, size_t *n, double *f, double *df, void *params)``
          (see potsynth_radial_fn in the C documentation), given as a ctypes function
          pointer, as an object with an ``address`` attribute (a ``numba.cfunc``), or
          as an integer address (e.g. of a Cython ``cdef`` function). It runs at the
          speed of C and may be called from several threads at once.

    params : for a compiled function, what it receives as ``params``: None, an
        integer address, a ctypes pointer or structure, or a NumPy array (its data)
    rmin, rmax : range of radii in which the functions are defined; the default radius
        range of a model that uses them
    """

    def __init__(self, n, func, params=None, rmin=None, rmax=None):
        if n < 1:
            raise ValueError("n must be at least 1")
        self.n = int(n)
        self.rmin, self.rmax = rmin, rmax
        self._keep = [func, params]        # referenced by the C structure below
        self._error = None                 # exception raised by a Python function

        if isinstance(func, ctypes._CFuncPtr):
            address = ctypes.cast(func, ctypes.c_void_p).value
        elif hasattr(func, "address"):
            address = int(func.address)
        elif isinstance(func, int):
            address = func
        elif callable(func):
            self._callback = _lib.RADIAL_FN_T(self._wrap(func))
            address = ctypes.cast(self._callback, ctypes.c_void_p).value
            params = None
        else:
            raise TypeError("func must be a function, a ctypes function pointer or an address")

        if params is None:
            paddr = None
        elif isinstance(params, int):
            paddr = params
        elif isinstance(params, np.ndarray):
            paddr = params.ctypes.data
        elif isinstance(params, ctypes._Pointer):
            paddr = ctypes.cast(params, ctypes.c_void_p).value
        else:
            paddr = ctypes.addressof(params)

        self._c = _lib.Radial_t(self.n, address, paddr)

    def _wrap(self, func):
        """C callback for a Python function func(r) -> (first, values) or None"""
        n = self.n

        def callback(r, first_p, n_p, f_p, df_p, params):
            try:
                res = func(r)
                if res is None:
                    return _lib.POTSYNTH_EDOM
                first, values = res
                values = np.atleast_1d(np.asarray(values, dtype=float))
                if first < 0 or values.size < 1 or first + values.size > n:
                    raise ValueError(f"radial function returned functions {first}..{first + values.size - 1} of {n}")
                first_p[0] = int(first)
                n_p[0] = values.size
                for j in range(values.size):
                    f_p[j] = values[j]
                return _lib.POTSYNTH_SUCCESS
            except Exception as e:          # an exception cannot cross the C library
                if self._error is None:
                    self._error = e
                return _lib.POTSYNTH_EDOM

        return callback


class _Evaluator:
    """Owns a C workspace and an evaluation workspace; evaluates its channels."""

    def __init__(self, ptr, nchan, chunk, needs_r=True):
        self._ptr = ptr
        self._nchan = nchan
        self._needs_r = needs_r
        _lib.clear_error()
        self._ew = lib.potsynth_eval_alloc(chunk, ptr)
        self._lock = threading.Lock()      # one eval workspace per object
        if not self._ew:
            lib.potsynth_free(ptr)
            self._ptr = None
            raise MemoryError("potsynth_eval_alloc failed")

    def __del__(self):
        if getattr(self, "_ew", None):
            lib.potsynth_eval_free(self._ew)
            self._ew = None
        if getattr(self, "_ptr", None):
            lib.potsynth_free(self._ptr)
            self._ptr = None

    def _points(self, r, theta, phi):
        """Flattened contiguous point arrays (r None for a model without radial
        dependence) and their common shape."""
        if r is None:
            if self._needs_r:
                raise ValueError("this model needs the radii of the points")
            theta, phi = np.broadcast_arrays(np.asarray(theta, float), np.asarray(phi, float))
        else:
            r, theta, phi = np.broadcast_arrays(np.asarray(r, float), np.asarray(theta, float),
                                                np.asarray(phi, float))
            r = np.ascontiguousarray(r.ravel())
        shape = theta.shape
        return r, np.ascontiguousarray(theta.ravel()), np.ascontiguousarray(phi.ravel()), shape

    def _evaluate(self, r, theta, phi, rows=False):
        r, theta, phi, shape = self._points(r, theta, phi)
        n = theta.size
        out = np.empty((n, self._nchan) if rows else (self._nchan, n))
        ptr = lambda a: None if a is None else a.ctypes.data_as(_lib.c_double_p)
        fn = lib.potsynth_eval_rows if rows else lib.potsynth_eval
        with self._lock:
            _lib.clear_error()
            status = fn(n, ptr(r), ptr(theta), ptr(phi), ptr(out), self._ew)
        if status:
            radial = getattr(self, "_radial", None)
            if radial is not None and radial._error is not None:
                error, radial._error = radial._error, None
                raise error
            _raise(status, "potsynth_eval_rows" if rows else "potsynth_eval")
        return out.reshape(shape + (self._nchan,) if rows else (self._nchan,) + shape)

    def rows(self, r, theta, phi):
        """Evaluate with the outputs of a point next to each other: an array of shape
        (shape of the points) + (number of outputs,). With the outputs as the basis
        functions of a fit, these are the rows of its design matrix. The values are
        those of calling the model."""
        return self._evaluate(r, theta, phi, rows=True)

    def write(self, filename):
        """Write the model (fast method only) to a file, to be read with
        potsynth.read() instead of being set up again. Radial functions are not
        stored; give them again when reading."""
        _lib.clear_error()
        status = lib.potsynth_write(os.fsencode(filename), self._ptr)
        if status:
            _raise(status, "potsynth_write")

    @property
    def info(self):
        """Parameters chosen at set-up (method, window width, radial rank, grid size,
        number of field components...)."""
        i = _lib.Info_t()
        lib.potsynth_get_info(ctypes.byref(i), self._ptr)
        d = {name: getattr(i, name) for name, _ in _lib.Info_t._fields_}
        d["method"] = _METHOD_NAMES.get(d["method"], d["method"])
        return d


class LoadedModel(_Evaluator):
    """A model read from a file (potsynth.read()). Calling it with (r, theta, phi)
    returns an array of shape (number of outputs,) + shape of the points, in the
    order of the outputs of the model that was written."""

    def __init__(self, filename, radial=None, nthreads=0, chunk=0):
        self._radial = radial          # keeps the C structures of the radial functions alive
        _lib.clear_error()
        ptr = lib.potsynth_read(os.fsencode(filename), None if radial is None else ctypes.byref(radial._c),
                                nthreads)
        if not ptr:
            _raise(_lib.POTSYNTH_EINVAL, "potsynth_read")
        i = _lib.Info_t()
        lib.potsynth_get_info(ctypes.byref(i), ptr)
        # the radii may be omitted for a model without radial dependence; the library checks
        super().__init__(ptr, int(i.nout), chunk, needs_r=False)

    def __call__(self, r, theta, phi):
        return self._evaluate(r, theta, phi)


def read(filename, radial=None, nthreads=0):
    """Read a model written with the write() method of a model. radial: the radial
    functions of a model that was set up with them (the same as at set-up)."""
    return LoadedModel(filename, radial=radial, nthreads=nthreads)


class Model(_Evaluator):
    """Generic internal potential-field model with user-defined channels.

    Parameters
    ----------
    coef : 1-D array, 4-pi fully normalised coefficients in the layout of idx()
           (see coef_vector()), length nlm(lmin, lmax, mmax)
    lmin, lmax, mmax : degree range and maximum order (mmax defaults to lmax)
    channels : list of Channel
    a : reference radius (same unit as r)
    rmin, rmax : radial range of the evaluation points (required by the fast method)
    eps : requested accuracy relative to the largest magnitude of each quantity
          (a scalar channel, or a horizontal vector or tensor part; see the documentation)
    method : "auto", "fast" or "direct"
    npts_hint : expected number of evaluation points; with "auto", the method with the
                smaller total time (set-up + evaluation) is chosen (0 = unknown)
    nthreads : OpenMP threads (0 = all)
    sigma, w, rank_tol : expert overrides (0 = automatic)
    chunk : points per internal chunk (0 = default)
    theta_min, theta_max : colatitude range of the evaluation points in radians
          (theta_max = 0: pi). The fast method stores only the part of its grid that
          the range needs; points outside of it raise DomainError

    Calling the model with (r, theta, phi) returns an array of shape
    (len(channels),) + broadcast shape of the inputs; rows(r, theta, phi) returns the
    same values with shape (broadcast shape) + (len(channels),).
    """

    def __init__(self, coef, lmin, lmax, channels, a, rmin=0.0, rmax=0.0, mmax=None, eps=1e-10,
                 method="auto", nthreads=0, sigma=0.0, w=0, rank_tol=0.0, chunk=0, npts_hint=0,
                 theta_min=0.0, theta_max=0.0):
        mmax = lmax if mmax is None else mmax
        coef = np.ascontiguousarray(coef, dtype=float)
        if coef.size != nlm(lmin, lmax, mmax):
            raise ValueError(f"coef has length {coef.size}, expected {nlm(lmin, lmax, mmax)}")
        if lmax > _lib.LMAX:
            raise ValueError(f"lmax > {_lib.LMAX} is not supported")

        chans = (_lib.Channel_t * len(channels))()
        self._rfac = []
        for i, ch in enumerate(channels):
            chans[i].rpow = ch.rpow
            chans[i].theta = _THETA_OPS[ch.theta]
            chans[i].dphi = ch.dphi
            if ch.rfac is not None:
                rf = np.ascontiguousarray(ch.rfac, dtype=float)
                if rf.size < lmax + 1:
                    raise ValueError("Channel.rfac must have length lmax + 1")
                self._rfac.append(rf)
                chans[i].rfac = rf.ctypes.data_as(_lib.c_double_p)

        p = _params(eps, method, nthreads, sigma, w, rank_tol, npts_hint, theta_min, theta_max)
        _lib.clear_error()
        ptr = lib.potsynth_alloc(lmin, lmax, mmax, coef.ctypes.data_as(_lib.c_double_p), a,
                                 rmin, rmax, len(channels), chans, ctypes.byref(p))
        if not ptr:
            _raise(_lib.POTSYNTH_EINVAL, "potsynth_alloc")
        super().__init__(ptr, len(channels), chunk)

    def __call__(self, r, theta, phi):
        return self._evaluate(r, theta, phi)
