"""
Vector spherical harmonic front-end: vector fields that need not be gradients of
a potential.
"""
import ctypes

import numpy as np

from . import _lib
from ._lib import lib
from .core import OPTIONS_DOC, _coefficients, _Evaluator, _options, _raise

_NORMS = {"full": _lib.VSH_NORM_FULL, "schmidt": _lib.VSH_NORM_SCHMIDT}
_QUANTITIES = {"V": (_lib.VSH_V, ("V_r", "V_theta", "V_phi")),
               "V_r": (_lib.VSH_VR, ("V_r",)),
               "V_h": (_lib.VSH_VH, ("V_theta", "V_phi")),
               "div": (_lib.VSH_DIV, ("div",)),
               "curl": (_lib.VSH_CURL, ("curl",))}
_ORDER = ("V_r", "V_theta", "V_phi", "div", "curl")


def _flags(quantities):
    """Bit flags and output names (in the order of the library) of a list of quantities."""
    if isinstance(quantities, str):
        quantities = [q.strip() for q in quantities.split(",") if q.strip()]
    flags, names = 0, set()
    for q in quantities:
        if q not in _QUANTITIES:
            raise ValueError(f"unknown quantity {q!r}; choose from {sorted(_QUANTITIES)}")
        flags |= _QUANTITIES[q][0]
        names.update(_QUANTITIES[q][1])
    if not flags:
        raise ValueError("no quantities requested")
    return flags, tuple(n for n in _ORDER if n in names)


def _first_given(a, b, c):
    """The first of a, b, c that is given (only its last dimension is used)"""
    return next(x for x in (a, b, c) if x is not None)


class VshModel(_Evaluator):
    """Vector fields given by vector spherical harmonic coefficients,

        V(r, theta, phi) = sum_j f_j(r) sum_lm [ a^j_lm Y_lm rhat + b^j_lm Psi_lm + c^j_lm Phi_lm ]

    with the real harmonics Y_lm = P_l^m(cos theta) cos(m phi) for m >= 0 and
    P_l^|m|(cos theta) sin(|m| phi) for m < 0 (the layout of idx()), and

        Psi_lm = dY/dtheta thetahat + (1/sin theta) dY/dphi phihat,   Phi_lm = rhat x Psi_lm.

    A model holds any number of fields, which share one grid of the fast method.

    Parameters
    ----------
    a, b, c : coefficients of Y rhat, Psi and Phi, each None (zero) or an array of shape
        (nexp, nrad, nlm), with nexp the number of fields, nrad the number of radial
        functions and nlm = nlm(lmin, lmax, mmax). (nexp, nlm) is accepted for
        nrad = 1, and (nlm,) for one field
    lmax, lmin, mmax : degree range and maximum order (lmax is found from nlm if None;
        mmax defaults to lmax)
    norm : "schmidt" (Schmidt semi-normalized P_l^m) or "full" (4-pi fully normalized)
    quantities : any of "V" (V_r, V_theta, V_phi), "V_r", "V_h" (V_theta, V_phi), "div"
        (surface divergence of the tangential part on the unit sphere; divide by r)
        and "curl" (radial component of its surface curl; divide by r), as a list or
        a comma-separated string. The outputs of each field are ordered V_r, V_theta,
        V_phi, div, curl
    radial : radial functions f_j (a Radial object), or None for fields on the
        sphere; the radii of the points are then not used and may be None
    rmin, rmax : radius range of the evaluation points, within the range of the radial
        functions (default: rmin and rmax of the Radial object)
    """ + OPTIONS_DOC + """
    eps is relative to the largest magnitude of each component of each field (V_r; the
    horizontal vector; the divergence; the curl).

    Calling the model with (r, theta, phi) returns an array of shape
    (nexp, nout) + shape of the points, with nout the number of outputs per field
    (model.outputs names them). rows(r, theta, phi) returns the same values with
    shape (shape of the points) + (nexp * nout,), field after field: with the fields as
    the basis functions of a fit, the rows of its design matrix.
    """

    def __init__(self, a=None, b=None, c=None, lmax=None, lmin=1, mmax=None, norm="schmidt",
                 quantities="V", radial=None, rmin=None, rmax=None, **options):
        flags, names = _flags(quantities)
        if norm not in _NORMS:
            raise ValueError("norm must be 'schmidt' or 'full'")
        nrad = 1 if radial is None else radial.n
        given = [x for x in (a, b, c) if x is not None]
        if not given:
            raise ValueError("at least one of a, b, c must be given")

        def shaped(x):
            """(nexp, nrad, nlm), C-contiguous"""
            if x is None:
                return None
            x = np.asarray(x, dtype=float)
            if x.ndim == 1:
                x = x[None, None, :]
            elif x.ndim == 2 and nrad == 1:
                x = x[:, None, :]
            if x.ndim != 3 or x.shape[1] != nrad:
                raise ValueError(f"coefficients must have shape (nexp, {nrad}, nlm)")
            return np.ascontiguousarray(x)

        a, b, c = shaped(a), shaped(b), shaped(c)
        shapes = {x.shape for x in (a, b, c) if x is not None}
        if len(shapes) != 1:
            raise ValueError("a, b and c must have the same shape")
        nexp = shapes.pop()[0]
        _, lmax, mmax = _coefficients(_first_given(a, b, c), lmin, lmax, mmax, "coefficients")
        p, chunk = _options(options)

        if radial is not None:
            rmin = radial.rmin if rmin is None else rmin
            rmax = radial.rmax if rmax is None else rmax
            if rmin is None or rmax is None:
                raise ValueError("rmin and rmax must be given (here or in the Radial object)")
        else:
            rmin = rmax = 0.0

        self._radial = radial           # keeps the C structures of the radial functions alive
        ptr = lambda x: None if x is None else x.ctypes.data_as(_lib.c_double_p)
        _lib.clear_error()
        w_ptr = lib.potsynth_vsh_alloc(lmin, lmax, mmax, nexp, ptr(a), ptr(b), ptr(c), _NORMS[norm],
                                       flags, rmin, rmax,
                                       None if radial is None else ctypes.byref(radial._c),
                                       ctypes.byref(p))
        if not w_ptr:
            _raise(_lib.POTSYNTH_EINVAL, "potsynth_vsh_alloc")
        super().__init__(w_ptr, nexp * len(names), chunk, needs_r=radial is not None)
        self.nexp = nexp
        self.outputs = names

    def __call__(self, r, theta, phi):
        out = self._evaluate(r, theta, phi)
        return out.reshape((self.nexp, len(self.outputs)) + out.shape[1:])
