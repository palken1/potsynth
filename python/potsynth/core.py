"""
Generic potential-field synthesis: functionals of an internal potential.
"""
import ctypes
import json
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


def _reset_radial(radial):
    """Forget an exception of the radial functions from an earlier call."""
    if radial is not None:
        radial._error = None


def _raise_radial(radial, status, what):
    """Raise the exception that the caller's radial functions raised in this call, if
    any (it cannot cross the C library), else the library error."""
    if radial is not None and radial._error is not None:
        error, radial._error = radial._error, None
        raise error
    _raise(status, what)


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
class Functional:
    """One functional of the core model,

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


# options of every model: name -> default. They are described once, in the documentation of
# Model, and accepted as keyword arguments by all model classes
_OPTIONS = {"eps": 1e-10, "method": "auto", "nthreads": 0, "sigma": 0.0, "w": 0, "rank_tol": 0.0,
            "chunk": 0, "npts_hint": 0, "theta_min": 0.0, "theta_max": 0.0}

OPTIONS_DOC = """
    Options (keyword arguments of every model class):

    eps : requested accuracy relative to the largest magnitude of each irreducible component
          (a scalar, or a horizontal vector or tensor; see the documentation), default 1e-10
    method : "auto" (default), "fast" or "direct"
    npts_hint : expected number of evaluation points; with "auto", the method with the
                smaller total time (set-up + evaluation) is chosen (0 = unknown)
    nthreads : OpenMP threads (0 = all)
    sigma, w, rank_tol : expert overrides of the oversampling factor, window width and
                radial SVD tolerance (0 = automatic)
    chunk : points per internal chunk (0 = default)
    theta_min, theta_max : colatitude range of the evaluation points in radians
          (theta_max = 0: pi). The fast method stores only the part of its grid that
          the range needs; points outside of it raise DomainError
"""



def _options(options):
    """C parameters and chunk size from the keyword options of a model"""
    unknown = set(options) - set(_OPTIONS)
    if unknown:
        raise TypeError(f"unknown option(s) {sorted(unknown)}; options are {sorted(_OPTIONS)}")
    o = dict(_OPTIONS, **options)
    if o["method"] not in _METHODS:
        raise ValueError(f"method must be one of {sorted(_METHODS)}")
    p = lib.potsynth_default_parameters()
    p.eps = o["eps"]
    p.method = _METHODS[o["method"]]
    p.nthreads = o["nthreads"]
    p.sigma = o["sigma"]
    p.w = o["w"]
    p.rank_tol = o["rank_tol"]
    p.npts_hint = o["npts_hint"]
    p.theta_min = o["theta_min"]
    p.theta_max = o["theta_max"]
    return p, o["chunk"]


def _lmax_of(n, lmin):
    """Largest degree of a coefficient vector of length n (orders up to the degree)"""
    lmax = lmin
    while nlm(lmin, lmax, lmax) < n:
        lmax += 1
    if nlm(lmin, lmax, lmax) != n:
        raise ValueError(f"{n} coefficients do not form a complete set of degrees {lmin}..")
    return lmax


def _truncate(coef, lmin, lmax, mmax):
    """The coefficients of degrees lmin..lmax and orders up to mmax from vectors (the last
    axis) of a complete set of degrees lmin..L, L >= lmax, or None if coef is not one"""
    n = coef.shape[-1]
    try:
        L = _lmax_of(n, lmin)
    except ValueError:
        return None
    if L < lmax:
        return None
    sel = [np.arange(idx(l, -min(l, mmax), lmin, L), idx(l, min(l, mmax), lmin, L) + 1)
           for l in range(lmin, lmax + 1)]
    return coef[..., np.concatenate(sel)] if sel else coef[..., :0]


def _coefficients(coef, lmin, lmax, mmax, name="coef"):
    """Coefficient vector, lmax and mmax from a vector in the layout of idx() or a tuple
    (C, S) of arrays indexed [l, m]. A complete vector of a higher degree is truncated
    to lmax and mmax."""
    if isinstance(coef, tuple):
        C, S = coef
        lmax = np.asarray(C).shape[0] - 1 if lmax is None else lmax
        mmax = lmax if mmax is None else mmax
        coef = coef_vector(C, S, lmin=lmin, lmax=lmax, mmax=mmax)
    else:
        coef = np.asarray(coef, dtype=float)
        lmax = _lmax_of(coef.shape[-1], lmin) if lmax is None else lmax
        mmax = lmax if mmax is None else mmax
    if coef.shape[-1] > nlm(lmin, lmax, mmax):
        t = _truncate(coef, lmin, lmax, mmax)
        coef = coef if t is None else t
    coef = np.ascontiguousarray(coef, dtype=float)
    if coef.shape[-1] != nlm(lmin, lmax, mmax):
        raise ValueError(f"{name} has length {coef.shape[-1]}, expected {nlm(lmin, lmax, mmax)} "
                         f"(or a complete set of degrees {lmin}..L with L >= {lmax})")
    if lmax > _lib.LMAX:
        raise ValueError(f"lmax > {_lib.LMAX} is not supported")
    return coef, lmax, mmax


def _range(rmin, rmax):
    """rmin and rmax for the C library (0: not given)"""
    return (0.0 if rmin is None else float(rmin)), (0.0 if rmax is None else float(rmax))


def _quantity_list(quantities):
    """A list or comma-separated string of quantities as a tuple"""
    if isinstance(quantities, str):
        quantities = [q.strip() for q in quantities.split(",") if q.strip()]
    return tuple(quantities)


def _names(names, n, what):
    """Output names: those given, or f0, f1, ..."""
    if names is None:
        return tuple(f"f{i}" for i in range(n))
    names = tuple(str(x) for x in names)
    if len(names) != n:
        raise ValueError(f"names has {len(names)} entries for {n} {what}")
    return names


def _sidecar(filename):
    return os.fsdecode(filename) + ".json"


def _rfac(rfac, lmax, keep, what):
    """C pointer to per-degree factors (None: 1); the array is kept alive in keep"""
    if rfac is None:
        return None
    rf = np.ascontiguousarray(rfac, dtype=float)
    if rf.ndim != 1 or rf.size < lmax + 1:
        raise ValueError(f"{what}.rfac must have length lmax + 1")
    keep.append(rf)
    return rf.ctypes.data_as(_lib.c_double_p)


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

        - None, with n = 1: no radial dependence (fields on the sphere; the radii of
          the points are then not used and may be None);

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
        self._sphere = func is None        # no radial dependence: fields on the sphere

        if func is None:
            if self.n != 1:
                raise ValueError("func=None (no radial dependence) needs n = 1")
            self.rmin = 1.0 if rmin is None else rmin
            self.rmax = 1.0 if rmax is None else rmax
            address = None
        elif isinstance(func, ctypes._CFuncPtr):
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


class _ModelHandle:
    """Owns a C model workspace; freed when the last user (the model object and the
    evaluation workspaces of its threads) is gone."""

    def __init__(self, ptr):
        self.ptr = ptr

    def __del__(self):
        if getattr(self, "ptr", None):
            lib.potsynth_free(self.ptr)
            self.ptr = None


class _EvalWorkspace:
    """One C evaluation workspace on a model (one per Python thread)."""

    def __init__(self, model, chunk):
        self._model = model                # keeps the model alive while this exists
        _lib.clear_error()
        self.ptr = lib.potsynth_eval_alloc(chunk, model.ptr)
        self.nthreads = 0                  # 0: those of the model
        if not self.ptr:
            raise MemoryError("potsynth_eval_alloc failed")

    def __del__(self):
        if getattr(self, "ptr", None):
            lib.potsynth_eval_free(self.ptr)
            self.ptr = None


class _Evaluator:
    """Owns a C model workspace and evaluates its functionals. Each Python thread that
    calls the model gets its own evaluation workspace, so one model can be evaluated
    from several threads at once."""

    def __init__(self, ptr, nfunc, chunk, needs_r=True):
        self._model = _ModelHandle(ptr)
        self._ptr = ptr
        self._nfunc = nfunc
        self._needs_r = needs_r
        self._chunk = chunk
        self._eval_threads = 0
        self._local = threading.local()
        self._workspace()                  # fails here, not at the first call, if out of memory

    @property
    def eval_threads(self):
        """OpenMP threads of each evaluation call; 0 (default) means those given at
        set-up (``nthreads``). Set 1 when calling the model from several Python
        threads at once, so that the threads are not multiplied."""
        return self._eval_threads

    @eval_threads.setter
    def eval_threads(self, n):
        n = int(n)
        if n < 0:
            raise ValueError("eval_threads must be >= 0")
        self._eval_threads = n

    def _workspace(self):
        ws = getattr(self._local, "ws", None)
        if ws is None:
            ws = _EvalWorkspace(self._model, self._chunk)
            self._local.ws = ws
        if ws.nthreads != self._eval_threads:
            _lib.clear_error()
            status = lib.potsynth_eval_set_threads(ws.ptr, self._eval_threads)
            if status:
                _raise(status, "potsynth_eval_set_threads")
            ws.nthreads = self._eval_threads
        return ws

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
        out = np.empty((n, self._nfunc) if rows else (self._nfunc, n))
        ptr = lambda a: None if a is None else a.ctypes.data_as(_lib.c_double_p)
        fn = lib.potsynth_eval_rows if rows else lib.potsynth_eval
        radial = getattr(self, "_radial", None)
        ws = self._workspace()
        _lib.clear_error()
        _reset_radial(radial)
        status = fn(n, ptr(r), ptr(theta), ptr(phi), ptr(out), ws.ptr)
        if status:
            _raise_radial(radial, status, "potsynth_eval_rows" if rows else "potsynth_eval")
        return out.reshape(shape + (self._nfunc,) if rows else (self._nfunc,) + shape)

    def rows(self, r, theta, phi):
        """Evaluate with the outputs of a point next to each other: an array of shape
        (shape of the points) + (number of outputs,). With the outputs as the basis
        functions of a fit, these are the rows of its design matrix. The values are
        those of calling the model."""
        return self._evaluate(r, theta, phi, rows=True)

    def write(self, filename):
        """Write the model (fast method only) to a file, to be read with
        potsynth.read() instead of being set up again, and the names of its outputs
        to filename + ".json". Radial functions are not stored; give them again when
        reading."""
        _lib.clear_error()
        status = lib.potsynth_write(os.fsencode(filename), self._ptr)
        if status:
            _raise(status, "potsynth_write")
        meta = {"potsynth": _lib.VERSION, "model": type(self).__name__, "outputs": list(self.outputs)}
        if getattr(self, "quantities", None) is not None:
            meta["quantities"] = list(self.quantities)
        if getattr(self, "nexp", None) is not None:
            meta["nexp"] = self.nexp
        with open(_sidecar(filename), "w") as f:
            json.dump(meta, f, indent=1)

    @property
    def info(self):
        """Parameters chosen at set-up (method, window width, radial rank, grid size,
        number of irreducible components...)."""
        i = _lib.Info_t()
        lib.potsynth_get_info(ctypes.byref(i), self._ptr)
        d = {name: getattr(i, name) for name, _ in _lib.Info_t._fields_ if name != "reserved"}
        d["method"] = _METHOD_NAMES.get(d["method"], d["method"])
        return d


class LoadedModel(_Evaluator):
    """A model read from a file (potsynth.read()). Calling it with (r, theta, phi)
    returns an array of shape (number of outputs,) + shape of the points, in the
    order of the outputs of the model that was written. outputs names them (from the
    file filename + ".json" written with the model; f0, f1, ... without it), model is
    the class of the model that was written and quantities its requested quantities
    (None if unknown)."""

    def __init__(self, filename, radial=None, nthreads=0, chunk=0):
        self._radial = radial          # keeps the C structures of the radial functions alive
        _lib.clear_error()
        _reset_radial(radial)
        ptr = lib.potsynth_read(os.fsencode(filename), None if radial is None else ctypes.byref(radial._c),
                                nthreads)
        if not ptr:
            _raise_radial(radial, _lib.POTSYNTH_EINVAL, "potsynth_read")
        i = _lib.Info_t()
        lib.potsynth_get_info(ctypes.byref(i), ptr)
        # the radii may be omitted for a model without radial dependence; the library checks
        super().__init__(ptr, int(i.nout), chunk, needs_r=False)
        meta = {}
        if os.path.exists(_sidecar(filename)):
            with open(_sidecar(filename)) as f:
                meta = json.load(f)
        nexp = meta.get("nexp", 1)
        outputs = meta.get("outputs")
        if outputs is None or len(outputs) * nexp != int(i.nout):
            outputs, nexp = None, 1
        self.outputs = _names(outputs, int(i.nout) // nexp, "outputs")
        self.model = meta.get("model") if outputs is not None else None
        q = meta.get("quantities") if outputs is not None else None
        self.quantities = None if q is None else tuple(q)
        if nexp > 1:
            self.nexp = nexp

    def __call__(self, r, theta, phi):
        return self._evaluate(r, theta, phi)


def read(filename, radial=None, nthreads=0):
    """Read a model written with the write() method of a model. radial: the radial
    functions of a model that was set up with them (the same as at set-up)."""
    return LoadedModel(filename, radial=radial, nthreads=nthreads)


class Model(_Evaluator):
    """Generic internal potential-field model with user-defined functionals.

    Parameters
    ----------
    coef : 1-D array, 4-pi fully normalised coefficients in the layout of idx()
           (see coef_vector()), length nlm(lmin, lmax, mmax)
    lmin, lmax, mmax : degree range and maximum order (mmax defaults to lmax)
    functionals : list of Functional
    a : reference radius (same unit as r)
    rmin, rmax : radial range of the evaluation points (required by the fast method;
           default: none, i.e. direct method only)
    names : names of the functionals (default f0, f1, ...), as model.outputs
    """ + OPTIONS_DOC + """
    Calling the model with (r, theta, phi) returns an array of shape
    (len(functionals),) + broadcast shape of the inputs; rows(r, theta, phi) returns the
    same values with shape (broadcast shape) + (len(functionals),).
    """

    def __init__(self, coef, lmin, lmax, functionals, a, rmin=None, rmax=None, mmax=None, names=None,
                 **options):
        names = _names(names, len(functionals), "functionals")
        coef, lmax, mmax = _coefficients(coef, lmin, lmax, mmax)
        rmin, rmax = _range(rmin, rmax)
        p, chunk = _options(options)

        funcs = (_lib.Functional_t * len(functionals))()
        self._rfac = []
        for i, ch in enumerate(functionals):
            if ch.theta not in _THETA_OPS:
                raise ValueError(f"unknown operator {ch.theta!r}; use {sorted(_THETA_OPS)}")
            funcs[i].rpow = ch.rpow
            funcs[i].theta = _THETA_OPS[ch.theta]
            funcs[i].dphi = ch.dphi
            funcs[i].rfac = _rfac(ch.rfac, lmax, self._rfac, "Functional")

        _lib.clear_error()
        ptr = lib.potsynth_alloc(lmin, lmax, mmax, coef.ctypes.data_as(_lib.c_double_p), a,
                                 rmin, rmax, len(functionals), funcs, ctypes.byref(p))
        if not ptr:
            _raise(_lib.POTSYNTH_EINVAL, "potsynth_alloc")
        super().__init__(ptr, len(functionals), chunk)
        self.outputs = names

    def __call__(self, r, theta, phi):
        return self._evaluate(r, theta, phi)


@dataclass
class Term:
    """One term of an Output of a SetsModel,

        fac * rfac[l] * T[Pbar_l^m(cos theta)] (d/dphi)^dphi [C_lm cos(m phi) + S_lm sin(m phi)]

    with the coefficients C, S of the coefficient set `set` (with radial functions: of
    the sets set, set + 1, ..., one for each radial function). theta and rfac as for
    Functional.
    """
    set: int
    theta: str = "P"
    dphi: int = 0
    fac: float = 1.0
    rfac: object = None


@dataclass
class Output:
    """One output of a SetsModel: (a/r)^(l + rpow) times the sum of its terms (with
    radial functions of the caller: their sum, and rpow is not used).

    The terms must form one irreducible component: all with the same operator, or the two
    operators ("dP", dphi=0) and ("Psin", dphi=1) of a horizontal vector, or, with a
    component and a weight given, ("S2A", dphi=0) and ("S2B", dphi=1) of a horizontal
    tensor. component: outputs with the same value > 0 form one component, and eps is
    relative to the largest magnitude of each component (0: chosen by the library);
    weight: the weight of the output in the magnitude of its component (0: chosen).
    """
    terms: list
    rpow: int = 0
    component: int = 0
    weight: float = 0.0


class SetsModel(_Evaluator):
    """Model with several coefficient sets: each output is a sum of terms, and each
    term applies one operator to one set. All outputs share one grid of the fast
    method, so many fields, or fields that are not gradients of one potential, are
    evaluated with one window sum per point.

    Parameters
    ----------
    coef : array of shape (nset, nlm), 4-pi fully normalised coefficients of each set
           in the layout of idx(), nlm = nlm(lmin, lmax, mmax)
    lmin, lmax, mmax : degree range and maximum order (lmax from nlm if None; mmax
           defaults to lmax)
    outputs : list of Output
    a : reference radius (same unit as r); required unless radial functions are given
           (then not used)
    rmin, rmax : radial range of the evaluation points (required by the fast method;
           default with radial functions: their range)
    radial : radial functions of the caller (Radial), in place of the powers of a/r,
           or None
    names : names of the outputs (default f0, f1, ...), as model.outputs
    """ + OPTIONS_DOC + """
    Calling the model with (r, theta, phi) returns an array of shape
    (len(outputs),) + broadcast shape of the inputs; rows() as for Model.
    """

    def __init__(self, coef, lmin, lmax, outputs, a=None, rmin=None, rmax=None, mmax=None, radial=None,
                 names=None, **options):
        names = _names(names, len(outputs), "outputs")
        if a is None:
            if radial is None:
                raise ValueError("the reference radius a is required without radial functions")
            a = 1.0
        coef = np.asarray(coef, dtype=float)
        if coef.ndim != 2:
            raise ValueError("coef must have shape (nset, nlm)")
        coef, lmax, mmax = _coefficients(coef, lmin, lmax, mmax)
        nset = coef.shape[0]
        p, chunk = _options(options)

        if radial is not None:
            rmin = radial.rmin if rmin is None else rmin
            rmax = radial.rmax if rmax is None else rmax
            if rmin is None or rmax is None:
                raise ValueError("rmin and rmax must be given (here or in the Radial object)")
        rmin, rmax = _range(rmin, rmax)

        self._radial = radial              # keeps the C structures of the radial functions alive
        self._rfac = []
        outs = (_lib.Output_t * len(outputs))()
        self._terms = []
        for i, o in enumerate(outputs):
            if not o.terms:
                raise ValueError("an Output needs at least one Term")
            terms = (_lib.Term_t * len(o.terms))()
            for k, t in enumerate(o.terms):
                if t.theta not in _THETA_OPS:
                    raise ValueError(f"unknown operator {t.theta!r}; use {sorted(_THETA_OPS)}")
                terms[k].set = t.set
                terms[k].fac = t.fac
                terms[k].rfac = _rfac(t.rfac, lmax, self._rfac, "Term")
                terms[k].theta = _THETA_OPS[t.theta]
                terms[k].dphi = t.dphi
            self._terms.append(terms)
            outs[i].rpow = o.rpow
            outs[i].nterm = len(o.terms)
            outs[i].term = terms
            outs[i].component = o.component
            outs[i].weight = o.weight

        _lib.clear_error()
        _reset_radial(radial)
        ptr = lib.potsynth_alloc_sets(lmin, lmax, mmax, nset, coef.ctypes.data_as(_lib.c_double_p), a, rmin, rmax,
                                      None if radial is None else ctypes.byref(radial._c),
                                      len(outputs), outs, ctypes.byref(p))
        if not ptr:
            _raise_radial(radial, _lib.POTSYNTH_EINVAL, "potsynth_alloc_sets")
        super().__init__(ptr, len(outputs), chunk, needs_r=not (radial is not None and radial.n == 1 and
                                                                radial._sphere))
        self.outputs = names

    def __call__(self, r, theta, phi):
        return self._evaluate(r, theta, phi)
