"""
Gravity front-end: potential V, gravitational acceleration g = grad V and gravity
gradient tensor T = grad grad V for 4-pi fully normalized coefficients (ICGEM .gfc
models such as EGM2008).
"""
import collections
import ctypes
import os
from dataclasses import dataclass

import numpy as np

from . import _lib
from ._lib import lib
from .core import OPTIONS_DOC, _coefficients, _Evaluator, _options, _raise, nlm

GRAV_V = 0x1
GRAV_GR = 0x2
GRAV_G = 0x4
GRAV_T = 0x8

_QUANTITIES = {"V": GRAV_V, "g_r": GRAV_GR, "g": GRAV_G, "T": GRAV_T}


@dataclass
class GfcModel:
    """Static gravity field model read from an ICGEM .gfc file.

    coef holds the 4-pi fully normalized coefficients, indexed by
    idx(l, m, lmin, mmax) (C at m >= 0, S at m < 0); coefficients absent from the
    file are 0. GM is in m^3/s^2 and a in m."""
    lmin: int
    lmax: int
    mmax: int
    GM: float
    a: float
    modelname: str
    tide_system: str
    coef: np.ndarray


def read_gfc(filename, lmax=0):
    """Read a static gravity field model from an ICGEM .gfc file.

    lmax truncates the model (0 = all degrees of the file). Unnormalized files are
    converted to fully normalized coefficients; the tide system is reported but not
    converted. Time-variable models are not supported. Returns a GfcModel."""
    _lib.clear_error()
    ptr = lib.potsynth_gfc_read(os.fsencode(filename), lmax)
    if not ptr:
        _raise(_lib.POTSYNTH_EINVAL, "potsynth_gfc_read")
    try:
        s = ptr.contents
        n = nlm(s.lmin, s.lmax, s.mmax)
        coef = np.ctypeslib.as_array(s.coef, shape=(n,)).copy()
        return GfcModel(int(s.lmin), int(s.lmax), int(s.mmax), float(s.GM), float(s.a),
                        s.modelname.decode(errors="replace"), s.tide_system.decode(errors="replace"),
                        coef)
    finally:
        lib.potsynth_gfc_free(ptr)


def _flags(quantities):
    if isinstance(quantities, str):
        quantities = [q.strip() for q in quantities.split(",")]
    flags = 0
    for q in quantities:
        if q not in _QUANTITIES:
            raise ValueError(f"unknown quantity {q!r}; use 'V', 'g_r', 'g' or 'T'")
        flags |= _QUANTITIES[q]
    if flags == 0:
        raise ValueError("no quantities requested")
    names = (["V"] if flags & GRAV_V else []) + \
            (["g_r", "g_theta", "g_phi"] if flags & GRAV_G else ["g_r"] if flags & GRAV_GR else []) + \
            (["T_rr", "T_rt", "T_rp", "T_tt", "T_tp", "T_pp"] if flags & GRAV_T else [])
    return flags, names


class GravModel(_Evaluator):
    """Gravity field model (EGM2008, EIGEN, XGM, GOCO...): gravitational potential,
    acceleration and gradient tensor of the model coefficients (no centrifugal term,
    no normal field).

    Parameters
    ----------
    coef : 4-pi fully normalized coefficients, either a 1-D vector in the layout of idx()
           (C at m >= 0, S at m < 0), or a tuple (C, S) of 2-D arrays C[l, m], S[l, m]
    GM : gravitational constant times mass in m^3/s^2
    a : reference radius in m
    lmax : maximum degree (default: from coef)
    lmin : minimum degree (default 0)
    mmax : maximum order (default lmax)
    quantities : any of "V" (potential, m^2/s^2), "g_r" (radial acceleration, m/s^2),
                 "g" (acceleration g_r, g_theta, g_phi, m/s^2) and "T" (gradient tensor
                 T_rr, T_rt, T_rp, T_tt, T_tp, T_pp, 1/s^2; 1 E = 1e-9/s^2), as a list or a
                 comma-separated string (default "V,g")
    rmin, rmax : radial range in m of the evaluation points (required by the fast
                 method; default: none, i.e. direct method only)
    """ + OPTIONS_DOC + """
    Calling the model with (r [m], theta [rad, colatitude], phi [rad]) returns a named
    tuple with the requested quantities in the order V, g_r, g_theta, g_phi, T_rr, ...,
    T_pp, each with the broadcast shape of the inputs. g = grad V: g_r is negative (towards the center)
    and g_theta is positive southwards.
    """

    def __init__(self, coef, GM, a, lmax=None, lmin=0, mmax=None, quantities="V,g", rmin=0.0, rmax=0.0,
                 **options):
        flags, names = _flags(quantities)
        coef, lmax, mmax = _coefficients(coef, lmin, lmax, mmax)
        p, chunk = _options(options)

        _lib.clear_error()
        ptr = lib.potsynth_grav_alloc(lmin, lmax, mmax, coef.ctypes.data_as(_lib.c_double_p),
                                      GM, a, rmin, rmax, flags, ctypes.byref(p))
        if not ptr:
            _raise(_lib.POTSYNTH_EINVAL, "potsynth_grav_alloc")
        super().__init__(ptr, len(names), chunk)
        self.quantities = tuple(names)
        self._result = collections.namedtuple("GravField", names)

    @classmethod
    def from_gfc(cls, filename, lmax=0, quantities="V,g", rmin=0.0, rmax=0.0, **kwargs):
        """Build a model from an ICGEM .gfc file (see read_gfc()), with GM and the
        reference radius of the file; lmax truncates the model (0 = all degrees). The
        other arguments are as for GravModel."""
        g = read_gfc(filename, lmax)
        return cls(g.coef, g.GM, g.a, lmax=g.lmax, lmin=g.lmin, mmax=g.mmax,
                   quantities=quantities, rmin=rmin, rmax=rmax, **kwargs)

    def __call__(self, r, theta, phi):
        return self._result(*self._evaluate(r, theta, phi))
