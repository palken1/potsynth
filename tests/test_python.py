"""
Tests of the Python interface (python/potsynth) against the reference data.

Usage: POTSYNTH_LIBRARY=/path/to/libpotsynth.so python test_python.py

With POTSYNTH_TEST_INSTALLED=1 the installed package (e.g. from a wheel) is
tested instead of the source tree.
"""
import ctypes
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
if os.environ.get("POTSYNTH_TEST_INSTALLED", "0") in ("", "0"):
    sys.path.insert(0, os.path.join(HERE, "..", "python"))
import potsynth  # noqa: E402

nfail = 0


def check(cond, msg):
    global nfail
    if not cond:
        print("FAIL:", msg)
        nfail += 1


def read_ref(path):
    with open(path, "rb") as f:
        L, n = np.fromfile(f, "<i8", 2)
        a, rmin, rmax = np.fromfile(f, "<f8", 3)
        gh = np.fromfile(f, "<f8", potsynth.nlm(1, int(L), int(L)))
        r, th, ph = (np.fromfile(f, "<f8", int(n)) for _ in range(3))
        B = np.fromfile(f, "<f8", 3 * int(n)).reshape(3, int(n))
    return int(L), a, rmin, rmax, gh, r, th, ph, B


def rel_err(B, Bref):
    return np.max(np.abs(np.asarray(B) - Bref)) / np.max(np.abs(Bref))


def test_reference():
    for name in ("ref_L13_s1.bin", "ref_L185_s3.bin"):
        L, a, rmin, rmax, gh, r, th, ph, Bref = read_ref(os.path.join(HERE, "data", name))
        for method, eps, tol in (("direct", 1e-10, 1e-13), ("fast", 1e-10, 1e-10),
                                 ("fast", 1e-14, 1e-13)):
            m = potsynth.MagModel(gh, lmax=L, a=a, rmin=rmin, rmax=rmax, eps=eps, method=method)
            err = rel_err(np.array(m(r, th, ph)), Bref)
            check(err <= tol, f"{name} {method} eps={eps}: {err:.2e}")
            print(f"  {name} {method:6s} eps={eps:.0e} err {err:.2e}")


def test_api():
    L, a, rmin, rmax, gh, r, th, ph, Bref = read_ref(os.path.join(HERE, "data", "ref_L60_s2.bin"))
    m = potsynth.MagModel(gh, a=a, rmin=rmin, rmax=rmax, method="fast")
    check(m.info["lmax"] == L and m.info["method"] == "fast", f"info {m.info}")

    # (g, h) arrays give the same model as the vector
    g = np.zeros((L + 1, L + 1))
    h = np.zeros((L + 1, L + 1))
    for l in range(1, L + 1):
        for k in range(l + 1):
            g[l, k] = gh[potsynth.idx(l, k, 1, L)]
            if k:
                h[l, k] = gh[potsynth.idx(l, -k, 1, L)]
    m2 = potsynth.MagModel((g, h), a=a, rmin=rmin, rmax=rmax, method="fast")
    check(np.array_equal(np.array(m(r, th, ph)), np.array(m2(r, th, ph))), "(g, h) input differs")

    # broadcasting and shapes
    rr = np.full((4, 5), 0.5 * (rmin + rmax))
    Br, Bt, Bp = m(rr, np.linspace(0, np.pi, 5), 1.0)
    check(Br.shape == (4, 5) and Bt.shape == (4, 5), "broadcast shape")
    Br0, _, _ = m(rmin, 0.0, 0.0)
    check(np.ndim(Br0) == 0 and np.isfinite(Br0), "scalar input")

    # NaN in, NaN out; domain errors raise
    Br, _, _ = m([rmin, np.nan], [1.0, 1.0], [0.0, 0.0])
    check(np.isfinite(Br[0]) and np.isnan(Br[1]), "NaN handling")
    try:
        m(rmax * 1.01, 1.0, 0.0)
        check(False, "r > rmax accepted")
    except potsynth.DomainError:
        pass
    try:
        potsynth.MagModel(gh[:-1], lmax=L)
        check(False, "short coefficient vector accepted")
    except ValueError:
        pass

    # NEC conversion
    N, E, C = potsynth.rtp_to_nec(1.0, 2.0, 3.0)
    check((N, E, C) == (-2.0, 3.0, -1.0), "rtp_to_nec")

    # generic core: potential V and dV/dr agree between fast and direct
    nn = np.arange(L + 1, dtype=float)
    funcs = [potsynth.Functional(rpow=1, theta="P"), potsynth.Functional(rpow=2, theta="P", rfac=nn + 1)]
    coef = potsynth.coef_vector(g / np.sqrt(2 * nn + 1)[:, None], h / np.sqrt(2 * nn + 1)[:, None])
    vf = potsynth.Model(coef, 1, L, funcs, a, rmin, rmax, method="fast", eps=1e-12)(r, th, ph)
    vd = potsynth.Model(coef, 1, L, funcs, a, method="direct")(r, th, ph)
    check(rel_err(vf, vd) < 1e-11, f"core fast vs direct {rel_err(vf, vd):.2e}")
    # dV/dr functional equals B_r of the magnetic model
    check(rel_err(vd[1], Bref[0]) < 1e-13, f"core dV/dr vs B_r {rel_err(vd[1], Bref[0]):.2e}")
    # SHC reader: file -> model equals vector -> model; time-dependent files rejected
    import tempfile
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "t.shc")
        with open(f, "w") as fp:
            fp.write("# test\n 1 2 1 1 1\n 2015.0\n")
            for ll in range(1, 3):
                for mm in range(-ll, ll + 1):
                    fp.write(f"{ll} {mm} {0.5 * ll + 0.1 * mm}\n")
        shc = potsynth.read_shc(f)
        check(shc.lmin == 1 and shc.lmax == 2 and shc.epoch == 2015.0, "read_shc header")
        check(shc.gh[potsynth.idx(2, -1, 1, 2)] == 0.9, "read_shc coefficient h21")
        m1 = potsynth.MagModel.from_shc(f, method="direct")
        m2 = potsynth.MagModel(shc.gh, lmax=2, method="direct")
        check(np.array_equal(np.array(m1(7000.0, 1.0, 2.0)), np.array(m2(7000.0, 1.0, 2.0))),
              "MagModel.from_shc")
        with open(f, "w") as fp:
            fp.write("1 1 2 2 1\n2015.0 2020.0\n1 0 1 1\n1 1 1 1\n1 -1 1 1\n")
        try:
            potsynth.read_shc(f)
            check(False, "time-dependent SHC accepted")
        except potsynth.PotsynthError as e:
            check(e.status == 4 and "time-dependent" in str(e), "time-dependent SHC error")

    # npts_hint: few points -> direct, many -> fast (L = 60)
    gh60 = np.ones(potsynth.nlm(1, 60, 60))
    m_few = potsynth.MagModel(gh60, rmin=rmin, rmax=rmax, npts_hint=1000)
    m_many = potsynth.MagModel(gh60, rmin=rmin, rmax=rmax, npts_hint=10_000_000)
    check(m_few.info["method"] == "direct" and m_many.info["method"] == "fast",
          f"npts_hint selection {m_few.info['method']} {m_many.info['method']}")

    print("  API checks done")


def test_grav():
    # gfc reader and gravity model: point mass + J2 against closed forms, from_gfc,
    # quantity selection, time-variable files rejected
    import tempfile
    GM, a, C20 = 3.986004415e14, 6378136.3, -0.484165143790e-03
    head = ("free text\nproduct_type gravity_field\nmodelname TEST\n"
            f"earth_gravity_constant {GM:.10e}\nradius {a:.10e}\nmax_degree 2\n"
            "tide_system zero_tide\nend_of_head\n")
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "t.gfc")
        with open(f, "w") as fp:
            fp.write(head + "gfc 0 0 1.0D+00 0.0\n" + f"gfc 2 0 {C20:.15e} 0.0\n")
        g = potsynth.read_gfc(f)
        check(g.lmax == 2 and g.GM == GM and g.a == a and g.modelname == "TEST"
              and g.tide_system == "zero_tide", "read_gfc header")
        check(g.coef[potsynth.idx(2, 0, 0, 2)] == C20 and g.coef[potsynth.idx(1, 0, 0, 2)] == 0.0,
              "read_gfc coefficients")

        r = np.array([6.3e6, 6.5e6, 6.9e6, 6.4e6])
        th = np.array([0.0, 0.7, 2.0, np.pi])
        ph = np.array([0.0, 1.0, 4.0, 2.0])
        x, q = np.cos(th), a / r
        P2, dP2 = np.sqrt(5) * 0.5 * (3 * x * x - 1), -np.sqrt(5) * 3 * x * np.sin(th)
        V = GM / r * (1 + q * q * C20 * P2)
        gr = -GM / r ** 2 * (1 + 3 * q * q * C20 * P2)
        gt = GM / r ** 2 * q * q * C20 * dP2
        for method in ("direct", "fast"):
            m = potsynth.GravModel.from_gfc(f, rmin=6.3e6, rmax=6.9e6, method=method)
            out = m(r, th, ph)
            tol = 1e-13 if method == "direct" else 1e-9
            check(m.quantities == ("V", "g_r", "g_theta", "g_phi"), "GravModel quantities")
            check(rel_err(out.V, V) < tol, f"GravModel V ({method})")
            check(rel_err(np.array([out.g_r, out.g_theta, out.g_phi]),
                          np.array([gr, gt, 0 * gr])) < tol, f"GravModel g ({method})")

        m = potsynth.GravModel(g.coef, GM, a, quantities=["g_r"], method="direct")
        check(m.quantities == ("g_r",) and rel_err(m(r, th, ph).g_r, gr) < 1e-13,
              "GravModel g_r only")
        try:
            potsynth.GravModel(g.coef, GM, a, quantities="W")
            check(False, "unknown quantity accepted")
        except ValueError:
            pass

        with open(f, "w") as fp:
            fp.write(head + "gfc 0 0 1.0 0.0\ngfct 2 0 -4.8e-4 0.0 0.0 0.0 20050101.0\n")
        try:
            potsynth.read_gfc(f)
            check(False, "time-variable gfc accepted")
        except potsynth.PotsynthError as e:
            check(e.status == 4, "time-variable gfc error")

    print("  gravity checks done")


def test_hessian():
    # gradient tensor functionals (S2A, S2B) against the 40-digit reference values
    with open(os.path.join(HERE, "data", "hess_L24_s12.bin"), "rb") as f:
        L, n = (int(v) for v in np.fromfile(f, "<i8", 2))
        a = np.fromfile(f, "<f8", 1)[0]
        coef = np.fromfile(f, "<f8", potsynth.nlm(0, L, L))
        r, th, ph = (np.fromfile(f, "<f8", n) for _ in range(3))
        T = np.fromfile(f, "<f8", 6 * n).reshape(6, n)
    l = np.arange(L + 1.0)
    funcs = [potsynth.Functional(3, "P", 0, (l + 1) * (l + 2) / a**2),
             potsynth.Functional(3, "dP", 0, -(l + 2) / a**2),
             potsynth.Functional(3, "Psin", 1, -(l + 2) / a**2),
             potsynth.Functional(3, "S2A", 0, np.full(L + 1, 1 / a**2)),
             potsynth.Functional(3, "S2B", 1, np.full(L + 1, 1 / a**2))]
    m = potsynth.Model(coef, 0, L, funcs, a)
    check(m.info["method"] == "direct", "no radial range: direct method")
    Trr, Trt, Trp, D, Ttp = m(r, th, ph)
    Tm = np.array([Trr, Trt, Trp, 0.5 * (-Trr + D), Ttp, 0.5 * (-Trr - D)])
    check(rel_err(Tm, T) < 1e-13, f"Hessian functionals vs reference ({rel_err(Tm, T):.2e})")
    for method, tol in (("direct", 1e-13), ("fast", 1e-11)):
        m = potsynth.Model(coef, 0, L, funcs, a, rmin=1.0, rmax=1.3, method=method, eps=1e-12)
        Trr, Trt, Trp, D, Ttp = m(r, th, ph)
        Tm = np.array([Trr, Trt, Trp, 0.5 * (-Trr + D), Ttp, 0.5 * (-Trr - D)])
        check(rel_err(Tm, T) < tol, f"Hessian functionals vs reference, {method} ({rel_err(Tm, T):.2e})")

    # front-ends: GravModel "T" is the reference tensor for GM = a = 1
    for method in ("direct", "fast"):
        g = potsynth.GravModel(coef, 1.0, 1.0, lmax=L, quantities="T", rmin=1.0, rmax=1.3,
                               method=method, eps=1e-12)
        out = g(r, th, ph)
        check(g.quantities == ("T_rr", "T_rt", "T_rp", "T_tt", "T_tp", "T_pp"), "GravModel T names")
        check(rel_err(np.array(out), T) < (1e-13 if method == "direct" else 1e-11),
              f"GravModel T vs reference, {method} ({rel_err(np.array(out), T):.2e})")

    # MagModel gradB = -T of the same potential (gh = C sqrt(2l+1), GM = a^2)
    am = 1.3
    cm = coef[potsynth.nlm(0, 0, 0):]
    gh = cm * np.concatenate([np.full(2 * l + 1, np.sqrt(2 * l + 1.0)) for l in range(1, L + 1)])
    mm = potsynth.MagModel(gh, lmax=L, a=am, quantities="B,gradB", method="direct")
    gm = potsynth.GravModel(cm, am * am, am, lmax=L, lmin=1, quantities="T", method="direct")
    fm, fg = mm(r * am, th, ph), gm(r * am, th, ph)
    check(mm.quantities[:3] == ("B_r", "B_theta", "B_phi") and len(mm.quantities) == 9,
          "MagModel B,gradB names")
    check(rel_err(np.array(fm[3:]), -np.array(fg)) < 1e-13, "MagModel gradB = -GravModel T")
    check(mm.info["nout"] == 9 and mm.info["nfunc"] == 8, "MagModel B,gradB info nout/nfunc")
    Br, Bt, Bp = potsynth.MagModel(gh, lmax=L, a=am, method="direct")(r * am, th, ph)
    check(np.array_equal(Br, fm.B_r) and np.array_equal(Bp, fm.B_phi), "MagModel B unchanged with gradB")

    # frame helpers: NEC mapping; Cartesian components at the poles independent of phi
    nec = potsynth.tensor_rtp_to_nec(1.0, 2.0, 3.0, 4.0, 5.0, 6.0)
    check([float(v) for v in nec] == [4.0, -5.0, 2.0, 6.0, -3.0, 1.0], "tensor_rtp_to_nec mapping")
    phis = np.array([0.0, 1.0, 2.5, -2.0])
    for th0 in (0.0, np.pi):
        tp = potsynth.GravModel(coef, 1.0, 1.0, lmax=L, quantities="T", method="direct")(1.2, th0, phis)
        X = np.array(potsynth.tensor_rtp_to_xyz(th0, phis, *tp))
        check(np.abs(X - X[:, :1]).max() < 1e-14 * np.abs(X).max(), f"Cartesian tensor at theta={th0:g}")
    M = potsynth.tensor_matrix(*T[:, :5])
    check(M.shape == (5, 3, 3) and np.allclose(M, np.swapaxes(M, -1, -2)), "tensor_matrix shape, symmetry")
    check(np.abs(np.trace(M, axis1=-2, axis2=-1)).max() < 1e-13 * np.abs(M).max(), "reference tensor traceless")

    print("  gradient tensor checks done")


def test_geodesy():
    # helpers: published WGS84 values, copy semantics, shapes
    w = potsynth.WGS84
    check(w.a == 6378137.0 and abs(1 / w.f - 298.257223563) < 1e-9, "WGS84 constants")
    C = potsynth.normal_zonals(10)
    check(abs(C[2] / -0.484166774985e-3 - 1) < 1e-12 and abs(C[4] / 0.790303733511e-6 - 1) < 1e-11,
          "WGS84 normal zonals")
    g = potsynth.normal_gravity(np.radians([0.0, 90.0]))
    check(abs(g[0] - 9.7803253359) < 1e-10 and abs(g[1] - 9.8321849378) < 1e-10, "WGS84 normal gravity")
    r, th = potsynth.geodetic_to_geocentric(np.zeros((2, 3)), 100.0)
    check(r.shape == (2, 3) and np.allclose(r, w.a + 100.0) and np.allclose(th, np.pi / 2),
          "geodetic_to_geocentric shape and equator")
    c = np.zeros(potsynth.nlm(0, 4, 4))
    c[0] = 1.0
    cT = potsynth.subtract_normal_field(c, w.GM, w.a)
    check(c[0] == 1.0 and abs(cT[0]) < 1e-16 and abs(cT[potsynth.idx(2, 0, 0, 4)] - (-C[2])) < 1e-18,
          "subtract_normal_field returns a modified copy")

    # recipe signs against the definitions (finite differences of T in latitude/longitude)
    rng = np.random.default_rng(5)
    L = 30
    coef = rng.normal(0, 1e-6, potsynth.nlm(0, L, L)) / 10
    m = potsynth.GravModel(coef, w.GM, w.a, quantities="V,g", method="direct")
    lat, lon, rr = 0.6, 1.1, 6.4e6
    th, h = np.pi / 2 - lat, 1e-6
    T, gr, gt, gp = (float(v) for v in m(rr, th, lon))
    dTdlat = (m(rr, th - h, lon).V - m(rr, th + h, lon).V) / (2 * h)
    dTdlon = (m(rr, th, lon + h).V - m(rr, th, lon - h).V) / (2 * h)
    gam = 9.8
    xi_def, eta_def = -dTdlat / (rr * gam), -dTdlon / (rr * gam * np.cos(lat))
    check(abs(gt / gam - xi_def) < 1e-8 * abs(xi_def), "recipe: xi = g_theta / gamma")
    check(abs(-gp / gam - eta_def) < 1e-8 * abs(eta_def), "recipe: eta = -g_phi / gamma")
    dTdr = (m(rr + 1.0, th, lon).V - m(rr - 1.0, th, lon).V) / 2.0
    check(abs(-gr - (-dTdr)) < 1e-7 * abs(dTdr), "recipe: delta g = -g_r")

    print("  geodesy checks done")


_RADIAL_C = r"""
#include <stddef.h>
/* barycentric interpolation; params: n, then n nodes, then n weights (doubles) */
int radial_bary(double r, size_t * first, size_t * n, double * f, double * df, void * params)
{
  const double * p = (const double *) params;
  const size_t nn = (size_t) p[0];
  const double * x = p + 1, * w = p + 1 + nn;
  double sum = 0.0;
  size_t j, hit = nn;

  (void) df;
  if (!(r >= x[0] && r <= x[nn - 1]))
    return 3; /* POTSYNTH_EDOM */

  for (j = 0; j < nn; ++j)
    if (r == x[j])
      hit = j;

  for (j = 0; j < nn; ++j)
    {
      f[j] = (hit < nn) ? (j == hit ? 1.0 : 0.0) : w[j] / (r - x[j]);
      sum += f[j];
    }

  for (j = 0; j < nn; ++j)
    f[j] /= sum;

  *first = 0;
  *n = nn;
  return 0;
}
"""


def _compiled_radial(tmpdir):
    """A compiled radial function as a ctypes function pointer, or None without a C compiler"""
    import shutil
    import subprocess
    cc = shutil.which("cc") or shutil.which("gcc") or shutil.which("clang")
    if cc is None:
        return None
    src = os.path.join(tmpdir, "potsynth_test_radial.c")
    lib = os.path.join(tmpdir, "potsynth_test_radial.so")
    with open(src, "w") as f:
        f.write(_RADIAL_C)
    if subprocess.call([cc, "-O2", "-shared", "-fPIC", "-o", lib, src]) != 0:
        return None
    dll = ctypes.CDLL(lib)
    _compiled_radial.keep = dll
    return dll.radial_bary


def test_vsh():
    """Vector spherical harmonic models, radial functions, rows, colatitude range, model files"""
    import tempfile
    rng = np.random.default_rng(1)
    L, N = 20, 400
    n = potsynth.nlm(1, L)
    th = np.arccos(rng.uniform(-1.0, 1.0, N))
    ph = rng.uniform(0.0, 2.0 * np.pi, N)
    a, b, c = (rng.normal(size=(3, n)) * 0.95 ** np.arange(n)[None, :] for _ in range(3))

    # fields on the sphere: fast against direct, all quantities; no radii
    m = potsynth.VshModel(a, b, c, quantities="V,div,curl", method="fast")
    d = potsynth.VshModel(a, b, c, quantities="V,div,curl", method="direct")
    V, Vd = m(None, th, ph), d(None, th, ph)
    check(V.shape == (3, 5, N) and m.outputs == ("V_r", "V_theta", "V_phi", "div", "curl"), f"vsh shape {V.shape}")
    check(m.info["ncomp"] == 12 and m.info["method"] == "fast", f"vsh info {m.info}")
    err = max(np.max(np.abs(V[:, q] - Vd[:, q])) / np.max(np.abs(Vd[:, q])) for q in range(5))
    check(err < 3e-10, f"vsh fast vs direct: {err:.2e}")

    # against the generic model: V_r = sum a Y; V_theta = b dP - c Psin d/dphi (Schmidt -> full)
    full = np.concatenate([np.full(2 * l + 1, 1.0 / np.sqrt(2 * l + 1)) for l in range(1, L + 1)])
    ch = [potsynth.Functional(0, "P"), potsynth.Functional(0, "dP"), potsynth.Functional(0, "Psin", dphi=1)]
    one = np.ones(N)
    ga = potsynth.Model(a[0] * full, 1, L, ch[:1], 1.0, method="direct")(one, th, ph)
    gb = potsynth.Model(b[0] * full, 1, L, ch[1:], 1.0, method="direct")(one, th, ph)
    gc = potsynth.Model(c[0] * full, 1, L, ch[1:], 1.0, method="direct")(one, th, ph)
    ref = np.array([ga[0], gb[0] - gc[1], gb[1] + gc[0]])
    err = np.max(np.abs(Vd[0, :3] - ref)) / np.max(np.abs(ref))
    check(err < 1e-13, f"vsh against generic functionals: {err:.2e}")

    # rows: the same values, point after point
    R = m.rows(None, th, ph)
    check(R.shape == (N, 15) and np.array_equal(R.reshape(N, 3, 5), np.moveaxis(V, -1, 0)), "vsh rows")
    mm = potsynth.MagModel(a[0], lmax=L, rmin=6700.0, rmax=6900.0, method="fast")
    rr = rng.uniform(6700.0, 6900.0, N)
    check(np.array_equal(mm.rows(rr, th, ph).T, np.array(mm(rr, th, ph))), "magnetic rows")

    # parts of the field
    h = potsynth.VshModel(None, b, c, quantities="V_h", method="direct")
    check(h.outputs == ("V_theta", "V_phi") and np.array_equal(h(None, th, ph), Vd[:, 1:3]), "vsh horizontal part")

    # radial functions of the caller: coefficient sets p(r_j) a give p(r) V_0
    nodes = np.array([6700.0, 6780.0, 6850.0, 6900.0])
    wts = np.array([1.0 / np.prod(np.delete(nodes[j] - nodes, j)) for j in range(4)])
    p = lambda r: 2.0 + 0.3 * (r - 6800.0) - 1e-3 * (r - 6800.0) ** 2
    dp = lambda r: 0.3 - 2e-3 * (r - 6800.0)
    v0 = potsynth.VshModel(a[:1], None, None, quantities="V_r", method="direct")(None, th, ph)[0, 0]

    def bary(r):
        """barycentric interpolation between the nodes: all four functions"""
        if not nodes[0] <= r <= nodes[-1]:
            return None
        d = r - nodes
        if np.any(d == 0.0):
            return 0, (d == 0.0).astype(float)
        f = wts / d
        return 0, f / f.sum()

    def hermite(r):
        """cubic Hermite interpolation: value and derivative functions of two nodes"""
        if not nodes[0] <= r <= nodes[-1]:
            return None
        k = min(int(np.searchsorted(nodes, r, side="right")) - 1, nodes.size - 2)
        h = nodes[k + 1] - nodes[k]
        t = (r - nodes[k]) / h
        return 2 * k, (2 * t**3 - 3 * t**2 + 1, h * (t**3 - 2 * t**2 + t), -2 * t**3 + 3 * t**2, h * (t**3 - t**2))

    A = {"bary": a[:1, None, :] * p(nodes)[None, :, None],
         "hermite": a[:1, None, :] * np.stack([p(nodes), dp(nodes)], axis=1).ravel()[None, :, None]}
    models = {}
    for kind, func, nthreads in (("bary", bary, 1), ("bary", bary, 2), ("hermite", hermite, 2)):
        rad = potsynth.Radial(A[kind].shape[1], func, rmin=nodes[0], rmax=nodes[-1])
        mr = potsynth.VshModel(A[kind], None, None, quantities="V_r", radial=rad, method="fast", eps=1e-12,
                               nthreads=nthreads)
        v = mr(rr, th, ph)[0, 0]
        err = np.max(np.abs(v - p(rr) * v0)) / np.max(np.abs(p(rr) * v0))
        check(err < 1e-11, f"Python radial function {kind}, {nthreads} threads: {err:.2e}")
        models[kind] = (mr, rad, v)
        try:
            mr(rr + 500.0, th, ph)
            check(False, f"radial {kind}: radius outside accepted")
        except potsynth.DomainError:
            pass
        try:
            mr(None, th, ph)
            check(False, f"radial {kind}: missing radii accepted")
        except ValueError:
            pass

    # an exception in a Python radial function is raised by the evaluation
    def broken(r):
        if r > 6850.0:
            raise ZeroDivisionError("radial function failed")
        return bary(r)

    rad = potsynth.Radial(4, broken, rmin=nodes[0], rmax=6840.0)
    mb = potsynth.VshModel(A["bary"], None, None, quantities="V_r", radial=rad, method="fast")
    try:
        mb(np.full(3, 6860.0), th[:3], ph[:3])
        check(False, "exception of a radial function not raised")
    except ZeroDivisionError:
        pass
    check(np.allclose(mb(rr[:5] * 0 + 6800.0, th[:5], ph[:5])[0, 0], p(6800.0) * v0[:5], rtol=1e-8),
          "model usable after an exception")

    # an exception on the second pass over the points (after the library validated them)
    # is raised by that call, not by a later one
    calls = [0]

    def late(r):
        calls[0] += 1
        if calls[0] > 4:
            raise ArithmeticError("radial function failed late")
        return bary(r)

    rad = potsynth.Radial(4, late, rmin=nodes[0], rmax=nodes[-1])
    calls[0] = -10**9
    ml = potsynth.VshModel(A["bary"], None, None, quantities="V_r", radial=rad, method="fast", nthreads=1)
    calls[0] = 0
    try:
        ml(np.full(3, 6800.0), th[:3], ph[:3])
        check(False, "late exception of a radial function not raised")
    except ArithmeticError:
        pass
    calls[0] = -10**9
    check(np.allclose(ml(rr[:5] * 0 + 6800.0, th[:5], ph[:5])[0, 0], p(6800.0) * v0[:5], rtol=1e-8),
          "no stale exception after a late failure")

    # a compiled radial function (needs a C compiler; skipped without one)
    compiled = _compiled_radial(tempfile.gettempdir())
    if compiled is None:
        print("  no C compiler: compiled radial function not tested")
    else:
        par = np.concatenate([[4.0], nodes, wts])
        for give in ("ctypes", "address"):
            func = compiled if give == "ctypes" else ctypes.cast(compiled, ctypes.c_void_p).value
            rad = potsynth.Radial(4, func, params=par, rmin=nodes[0], rmax=nodes[-1])
            mc = potsynth.VshModel(A["bary"], None, None, quantities="V_r", radial=rad, method="fast", eps=1e-12,
                                   nthreads=2)
            v = mc(rr, th, ph)[0, 0]
            err = np.max(np.abs(v - p(rr) * v0)) / np.max(np.abs(p(rr) * v0))
            check(err < 1e-11, f"compiled radial function ({give}): {err:.2e}")
            try:
                mc(rr + 500.0, th, ph)
                check(False, "compiled radial function: radius outside accepted")
            except potsynth.DomainError:
                pass

    try:
        import numba
    except ImportError:
        numba = None
    if numba is not None:
        from numba import carray, cfunc, types
        sig = types.intc(types.double, types.CPointer(types.uintp), types.CPointer(types.uintp),
                         types.CPointer(types.double), types.CPointer(types.double), types.voidptr)
        xn, wn = tuple(nodes), tuple(wts)

        @cfunc(sig)
        def nb_bary(r, first, n, f, df, params):
            if not (r >= xn[0] and r <= xn[3]):
                return 3
            out = carray(f, 4)
            s = 0.0
            hit = -1
            for j in range(4):
                if r == xn[j]:
                    hit = j
            for j in range(4):
                out[j] = (1.0 if j == hit else 0.0) if hit >= 0 else wn[j] / (r - xn[j])
                s += out[j]
            for j in range(4):
                out[j] /= s
            first[0] = 0
            n[0] = 4
            return 0

        rad = potsynth.Radial(4, nb_bary, rmin=nodes[0], rmax=nodes[-1])
        mn = potsynth.VshModel(A["bary"], None, None, quantities="V_r", radial=rad, method="fast", eps=1e-12)
        v = mn(rr, th, ph)[0, 0]
        err = np.max(np.abs(v - p(rr) * v0)) / np.max(np.abs(p(rr) * v0))
        check(err < 1e-11, f"numba radial function: {err:.2e}")

    # model files
    mr, rad, v = models["bary"]
    name = os.path.join(tempfile.gettempdir(), "potsynth_test_python_model.bin")
    mr.write(name)
    m2 = potsynth.read(name, radial=rad)
    check(np.array_equal(m2(rr, th, ph)[0], v) and m2.info["rho"] == rad.n, "model file")
    try:
        potsynth.read(name)
        check(False, "model file read without its radial functions")
    except potsynth.PotsynthError:
        pass
    os.remove(name)
    try:
        d.write(name)
        check(False, "direct model written")
    except potsynth.PotsynthError:
        pass

    # colatitude range
    cap = potsynth.MagModel(a[0], lmax=L, rmin=6700.0, rmax=6900.0, method="fast", theta_max=0.8)
    inside = th < 0.8
    check(cap.info["nrows"] < mm.info["nrows"], "colatitude range: fewer rows")
    check(np.array_equal(np.array(cap(rr[inside], th[inside], ph[inside])),
                         np.array(mm(rr[inside], th[inside], ph[inside]))), "colatitude range: same values")
    try:
        cap(6800.0, 1.0, 0.0)
        check(False, "colatitude range: point outside accepted")
    except potsynth.DomainError:
        pass

    print("  vsh, radial functions, rows, model files and colatitude range checks done")


def test_sets():
    """General models with several coefficient sets, options handling"""
    rng = np.random.default_rng(3)
    L, N = 16, 300
    n = potsynth.nlm(1, L)
    th = np.arccos(rng.uniform(-1.0, 1.0, N))
    ph = rng.uniform(0.0, 2.0 * np.pi, N)
    rr = rng.uniform(6700.0, 6900.0, N)
    coef = rng.normal(size=(2, n)) * 0.9 ** np.arange(n)[None, :]

    # 2 x set 0 - 3 x set 1 under P, and a horizontal vector of the two sets, against functionals
    outs = [potsynth.Output([potsynth.Term(0, "P", fac=2.0), potsynth.Term(1, "P", fac=-3.0)], rpow=1),
            potsynth.Output([potsynth.Term(0, "dP"), potsynth.Term(1, "Psin", dphi=1, fac=-1.0)], rpow=1),
            potsynth.Output([potsynth.Term(0, "Psin", dphi=1), potsynth.Term(1, "dP")], rpow=1)]
    ch = [potsynth.Functional(1, "P"), potsynth.Functional(1, "dP"), potsynth.Functional(1, "Psin", dphi=1)]
    a = 6371.2
    c0 = potsynth.Model(coef[0], 1, L, ch, a, method="direct")(rr, th, ph)
    c1 = potsynth.Model(coef[1], 1, L, ch, a, method="direct")(rr, th, ph)
    ref = np.array([2.0 * c0[0] - 3.0 * c1[0], c0[1] - c1[2], c0[2] + c1[1]])
    for method, tol in (("direct", 1e-13), ("fast", 3e-11)):
        m = potsynth.SetsModel(coef, 1, L, outs, a=a, rmin=6700.0, rmax=6900.0, method=method, eps=1e-11)
        v = m(rr, th, ph)
        err = max(np.max(np.abs(v[k] - ref[k])) / np.max(np.abs(ref[k])) for k in range(3))
        check(err < tol and m.info["ncomp"] == 2, f"SetsModel {method}: {err:.2e}, {m.info['ncomp']} components")
        check(np.array_equal(m.rows(rr, th, ph).T, v), f"SetsModel {method} rows")

    # a field on the sphere (no radial dependence): radii not needed
    sph = potsynth.Radial(1, None)
    ms = potsynth.SetsModel(coef, 1, L, outs[:1], radial=sph, method="direct")
    v = ms(None, th, ph)[0]
    one = np.ones(N)
    r1 = 2.0 * potsynth.Model(coef[0], 1, L, ch[:1], 1.0, method="direct")(one, th, ph)[0] \
        - 3.0 * potsynth.Model(coef[1], 1, L, ch[:1], 1.0, method="direct")(one, th, ph)[0]
    check(np.max(np.abs(v - r1)) < 1e-12 * np.max(np.abs(r1)), "SetsModel on the sphere")

    # errors: mixed component, unknown option, unknown operator
    for bad, kind in ((lambda: potsynth.SetsModel(coef, 1, L, [potsynth.Output([potsynth.Term(0, "P"),
                                                                                  potsynth.Term(1, "dP")])],
                                                   method="direct"), potsynth.PotsynthError),
                      (lambda: potsynth.MagModel(coef[0], lmax=L, rmin=6700.0, rmax=6900.0, epsilon=1e-6), TypeError),
                      (lambda: potsynth.SetsModel(coef, 1, L, [potsynth.Output([potsynth.Term(0, "Q")])]), ValueError)):
        try:
            bad()
            check(False, f"invalid model accepted ({kind.__name__})")
        except kind:
            pass

    print("  sets and options checks done")


def test_version():
    # the Python package and the C library carry the same version
    import re
    with open(os.path.join(HERE, "..", "include", "potsynth.h")) as f:
        h = f.read()
    v = ".".join(re.search(rf"POTSYNTH_VERSION_{k}\s+(\d+)", h).group(1) for k in ("MAJOR", "MINOR", "PATCH"))
    check(v == potsynth.__version__, f"version: header {v}, package {potsynth.__version__}")
    print(f"  library: {potsynth.library_path}")


def test_former_names():
    # Channel and info["nchan"] are the former names of Functional and info["nfunc"]
    check(potsynth.Channel is potsynth.Functional, "former name Channel")
    m = potsynth.MagModel(np.array([0.0, 1.0, 2.0]), lmax=1, quantities="B", method="direct")
    check(m.info["nchan"] == m.info["nfunc"] == 3, "former name info['nchan']")


if __name__ == "__main__":
    test_version()
    test_reference()
    test_api()
    test_grav()
    test_hessian()
    test_geodesy()
    test_vsh()
    test_sets()
    test_former_names()
    print("FAILED" if nfail else "all passed")
    sys.exit(1 if nfail else 0)
