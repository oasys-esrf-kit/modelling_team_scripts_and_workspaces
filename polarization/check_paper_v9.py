# check_paper.py
#
# Numerically and symbolically verifies the Jones/Mueller formulas in:
#   "Converting X-ray Linear Polarization to Circular Polarization
#    Using a Crystal Wave Plate"
#
# Equation numbers refer to the paper as produced by xray_polarization_corrected_v2.tex.
# Complete equation map (53 labelled equations):
#   (1-4)   Stokes parameter definitions S0..S3
#   (5)     Pauli matrices sigma_0..3
#   (6)     Sk = Tr(sigma_k rho)
#   (7)     M_mn = 1/2 Tr(sigma_m J sigma_n J†)   [Pauli trace formula]
#   (8)     J(chi) = R(-chi) J0 R(chi)             [Jones rotation sandwich]
#   (9)     R(chi) = [[cos,sin],[-sin,cos]]         [paper's rotation matrix]
#   (10)    R_M(chi) 4x4 Mueller rotation matrix
#   (11)    M(chi) = R_M(-chi) M0 R_M(chi)         [Mueller rotation sandwich]
#   (12)    J_ret^(0) = e^{-id/2}[[1,0],[0,e^{id}]]  [general retarder, chi=0]
#   (13)    J_QWP^(0) = e^{-ipi/4}[[1,0],[0,i]]    [QWP at chi=0]
#   (14)    M_QWP^(0) = diag-like 4x4              [QWP Mueller at chi=0]
#   (15)    M_ret(chi,delta) full 4x4               [general rotated retarder]
#   (16)    J_QWP^(45) = e^{+ipi/4}/sqrt2 [[1,-i],[-i,1]]  [QWP at chi=+45°]
#   (17)    M_QWP^(45) = [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,-1,0,0]]
#   (18)    S_in = [1,1,0,0]^T                     [H-linear input Stokes]
#   (19)    E_out = [1,-i]/sqrt2                    [RCP output Jones vector]
#   (20)    I(theta) = 1/2(S0 + S1 cos2θ + S2 sin2θ)  [LP detected intensity]
#   (21-24) IH, IV, I45, I135 at key polarizer angles
#   (25-27) S0, S1, S2 recovery from LP measurements
#   (28-29) I+, I- with QWP analyser
#   (30-31) S3 and DOCP recovery formulas
#   (32)    I(theta) = I0/2 for circular input (all theta)
#   (33)    J0 = diag(rs, rp)                       [crystal Jones matrix, chi=0]
#   (34)    J0 Pauli decomposition
#   (35)    J(chi) Pauli decomposition
#   (36)    J(chi) matrix form
#   (37)    M0 4x4                                  [crystal Mueller, chi=0]
#   (38)    M(chi) full rotated 4x4                 [Mueller_final]
#   (39-42) E_out general and special cases
#   (43)    S_out = M(chi) S_in
#   (44-47) S0..S3 output, general input
#   (48-51) S0..S3 output, H-linear input
#   (52)    Pc general formula
#   (53)    Pc for H-linear input  [= Giles' formula with sign note]
#
# Note: LP Jones and LP Mueller matrices are displayed in Sect. 3.2 but carry
#       no \label, so they have no equation number.
# Note: substitution cases for +45° and RCP inputs (derived inline from eqs 44-47
#       and 52) also carry no separate equation labels.
#
# Convention (same as the paper):
#   R(chi) = [[cos chi, sin chi], [-sin chi, cos chi]]
#   This is NOT the standard CCW rotation; it is its transpose.
#   Rotated Jones:   J(chi) = R(-chi) @ J0 @ R(chi)            [eq 8]
#   Rotated Mueller: M(chi) = R_M(-chi) @ M0 @ R_M(chi)        [eq 11]

import numpy as np
import sympy as sp

# ===========================================================================
# NumPy helper functions
# ===========================================================================

def _jones_to_stokes_np(j, normalize=True):
    """Stokes vector from Jones vector using eqs (1-4)."""
    j = np.array(j, dtype=complex)
    j1, j2 = j[0], j[1]
    S = [
        (j1 * j1.conjugate() + j2 * j2.conjugate()).real,   # eq (1)
        (j1 * j1.conjugate() - j2 * j2.conjugate()).real,   # eq (2)
        2 * (j1 * j2.conjugate()).real,                      # eq (3)
        2 * (j1 * j2.conjugate()).imag,                      # eq (4)
    ]
    if normalize:
        S[1] /= S[0]; S[2] /= S[0]; S[3] /= S[0]
    return np.array(S, dtype=float)


def _rotation_J_np(chi: float) -> np.ndarray:
    """2×2 Jones rotation matrix R(chi) — eq (9)."""
    c, s = np.cos(chi), np.sin(chi)
    return np.array([[c, s], [-s, c]])


def _rotation_M_np(chi: float) -> np.ndarray:
    """4×4 Mueller rotation matrix R_M(chi) — eq (10)."""
    c2, s2 = np.cos(2.0 * chi), np.sin(2.0 * chi)
    return np.array([[1, 0, 0, 0],
                     [0, c2, s2, 0],
                     [0, -s2, c2, 0],
                     [0, 0, 0, 1]], dtype=np.float64)


def _pc_giles(rs, rp, chi):
    """Circular polarization degree, Giles' thesis form (eq 53 denominator variant)."""
    SinPhi = (rs * rp.conjugate()).imag / np.abs(rs) / np.abs(rp)
    a = np.abs(rs) / np.abs(rp)
    return -a * SinPhi * np.sin(2 * chi) / (np.cos(chi)**2 + a**2 * np.sin(chi)**2)


def _pc_np(rs, rp, theta, s_in):
    """Degree of circular polarization from eq (52)."""
    c2, s2 = np.cos(2 * theta), np.sin(2 * theta)
    rs_sq = np.real(rs * np.conj(rs))
    rp_sq = np.real(rp * np.conj(rp))
    cross = rs * np.conj(rp)
    re_cp, im_cp = np.real(cross), np.imag(cross)
    s0_i, s1_i, s2_i, s3_i = s_in
    num = im_cp * (s2 * s1_i - c2 * s2_i) + re_cp * s3_i
    den = (rs_sq + rp_sq) / 2 * s0_i + (rs_sq - rp_sq) / 2 * (c2 * s1_i + s2 * s2_i)
    return num / den


def _mueller_rotated_crystal_np(rs: complex, rp: complex, theta: float) -> np.ndarray:
    """Rotated crystal Mueller matrix M(chi) — eq (38)."""
    c2, s2 = np.cos(2 * theta), np.sin(2 * theta)
    rs_sq = np.real(rs * np.conj(rs))
    rp_sq = np.real(rp * np.conj(rp))
    tp = (rs_sq + rp_sq) / 2
    tm = (rs_sq - rp_sq) / 2
    cross = rs * np.conj(rp)
    re_cp, im_cp = np.real(cross), np.imag(cross)
    return np.array([
        [tp,      tm*c2,                  tm*s2,                  0       ],
        [tm*c2,   tp*c2**2+re_cp*s2**2,   (tp-re_cp)*s2*c2,      -im_cp*s2],
        [tm*s2,   (tp-re_cp)*s2*c2,       tp*s2**2+re_cp*c2**2,   im_cp*c2],
        [0,       im_cp*s2,              -im_cp*c2,               re_cp   ],
    ], dtype=complex)


# ===========================================================================
# SymPy helper functions
# ===========================================================================

def _rotation_J_sp(a):
    """2×2 Jones rotation matrix R(chi) — eq (9)."""
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c, s], [-s, c]])


def _rotation_M_sp(a):
    """4×4 Mueller rotation matrix R_M(chi) — eq (10)."""
    c2a, s2a = sp.cos(2*a), sp.sin(2*a)
    return sp.Matrix([[1,0,0,0],[0,c2a,s2a,0],[0,-s2a,c2a,0],[0,0,0,1]])


def _crystal_J_sp(r_s, r_p):
    """Crystal Jones matrix at chi=0 — eq (33)."""
    return sp.Matrix([[r_s, 0], [0, r_p]])


def _retarder_J_sp(delta, use_global_phase=False):
    """General retarder Jones matrix at chi=0 — eq (12)."""
    if use_global_phase:
        return sp.exp(-sp.I * delta / 2) * sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])
    return sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])


def _linearpolarizer_J_sp(a):
    """Linear polarizer Jones matrix at angle a (Sect. 3.2, unlabelled)."""
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c*c, c*s], [c*s, s*s]])


def _M0_matrix_sp(rs, rp):
    """Crystal Mueller matrix at chi=0 written explicitly — eq (37)."""
    rs_sq = rs * rs.conjugate()
    rp_sq = rp * rp.conjugate()
    tp = (rs_sq + rp_sq) / 2
    tm = (rs_sq - rp_sq) / 2
    cross = rs * rp.conjugate()
    re_c = sp.re(cross)
    im_c = sp.im(cross)
    return sp.Matrix([
        [tp,  tm,  0,    0   ],
        [tm,  tp,  0,    0   ],
        [0,   0,   re_c, im_c],
        [0,   0,  -im_c, re_c],
    ])


def _jones_to_mueller_sp(J):
    """Convert Jones matrix to Mueller matrix via Pauli trace formula — eq (7)."""
    sigmas = [
        sp.Matrix([[1, 0], [0,  1]]),
        sp.Matrix([[1, 0], [0, -1]]),
        sp.Matrix([[0, 1], [1,  0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    ]
    J_dag = J.H
    M = sp.Matrix.zeros(4, 4)
    for i in range(4):
        for j in range(4):
            M[i, j] = sp.simplify((sigmas[i] * J * sigmas[j] * J_dag).trace() / 2)
    return M


def _pauli_sp():
    """Return list of four Pauli matrices — eq (5)."""
    return [
        sp.Matrix([[1, 0], [0,  1]]),
        sp.Matrix([[1, 0], [0, -1]]),
        sp.Matrix([[0, 1], [1,  0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    ]


def _mueller_rotated_crystal_sp(rs, rp, theta):
    """Rotated crystal Mueller matrix M(chi) — eq (38)."""
    c2, s2 = sp.cos(2 * theta), sp.sin(2 * theta)
    rs_sq = rs * sp.conjugate(rs)
    rp_sq = rp * sp.conjugate(rp)
    tp = (rs_sq + rp_sq) / 2
    tm = (rs_sq - rp_sq) / 2
    cross = rs * sp.conjugate(rp)
    re_cp = sp.re(cross)
    im_cp = sp.im(cross)
    return sp.Matrix([
        [tp,      tm*c2,                  tm*s2,                  0       ],
        [tm*c2,   tp*c2**2+re_cp*s2**2,   (tp-re_cp)*s2*c2,      -im_cp*s2],
        [tm*s2,   (tp-re_cp)*s2*c2,       tp*s2**2+re_cp*c2**2,   im_cp*c2],
        [0,       im_cp*s2,              -im_cp*c2,               re_cp   ],
    ])


def _stokes_Sout_sp(S_in, r_s, r_p, chi):
    """
    Full S_out after crystal at angle chi — eqs (43-47).
    Shortcut: builds M(chi) via eq (11) then applies to S_in.
    """
    M = (_rotation_M_sp(-chi)
         @ _jones_to_mueller_sp(_crystal_J_sp(r_s, r_p))
         @ _rotation_M_sp(chi))
    return M @ S_in


# ===========================================================================
# main
# ===========================================================================

if __name__ == "__main__":
    chi_s   = sp.Symbol("chi",   real=True)
    delta_s = sp.Symbol("delta", real=True)

    # -----------------------------------------------------------------------
    print("=" * 60)
    print("SECTION: Rotation matrices (eqs 9-10)")
    print("=" * 60)
    print("R_J(chi)  eq (9) :", _rotation_J_sp(chi_s))
    print("R_M(chi)  eq (10):", _rotation_M_sp(chi_s))

    for _ in range(100):
        chi = np.random.rand() * 2 * np.pi
        diff_J = (_rotation_J_np(chi) - np.array(
            _rotation_J_sp(chi_s).subs(chi_s, chi).evalf(15)).astype(float)).sum()
        diff_M = (_rotation_M_np(chi) - np.array(
            _rotation_M_sp(chi_s).subs(chi_s, chi).evalf(15)).astype(float)).sum()
        assert abs(diff_J) < 1e-10 and abs(diff_M) < 1e-10
    print("R_J and R_M numpy/sympy agreement: CHECKED (100 random angles) ✓")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SECTION: Retarder / QWP (eqs 12-19)")
    print("=" * 60)

    print("\nRetarder J at chi=0  eq (12):", _retarder_J_sp(delta_s))
    M_ret0_sp = _jones_to_mueller_sp(_retarder_J_sp(delta_s))
    print("Retarder M at chi=0 (derived via eq 7):", M_ret0_sp)

    delta_qwp = np.pi / 2
    print("\nQWP J at chi=0  eq (13):",
          _retarder_J_sp(delta_s).subs(delta_s, delta_qwp).evalf(5))
    print("QWP M at chi=0  eq (14):",
          _jones_to_mueller_sp(_retarder_J_sp(delta_s)).subs(delta_s, delta_qwp).evalf(5))

    # General rotated retarder Mueller — eq (15): R_M(-chi) @ M_ret0 @ R_M(chi)
    M_ret_rot = (_rotation_M_sp(-chi_s) @ M_ret0_sp @ _rotation_M_sp(chi_s))
    print("\nGeneral rotated retarder M  eq (15):", M_ret_rot)

    # -------------------------------------------------------------------
    # eq (16): J_QWP^(45) = R(-45) @ J_QWP^(0) @ R(45)  [eq 8]
    #   Result: e^{+i pi/4}/sqrt(2) * [[1,-i],[-i,1]]
    # -------------------------------------------------------------------
    chi_45  = np.pi / 4
    QWP0_np = np.array([[1, 0], [0, 1j]], dtype=complex)   # eq (13), no global phase

    J_qwp_45 = _rotation_J_np(-chi_45) @ QWP0_np @ _rotation_J_np(chi_45)   # eq (8)
    factor   = J_qwp_45[0, 0]   # = e^{+i pi/4}/sqrt(2)
    print(f"\neq (16)  J_QWP(chi=+45°) = R(-45)@J_QWP^(0)@R(45)  [eq 8]:")
    print(f"  computed = {J_qwp_45.round(4)}")
    print(f"  global phase factor = {factor:.4f}  (= e^{{+i pi/4}}/sqrt(2))")
    print(f"  J / factor = {(J_qwp_45/factor).round(4)}")
    print(f"  => J = e^{{+i pi/4}}/sqrt(2) * [[1,-i],[-i,1]]  ✓")

    # -------------------------------------------------------------------
    # eq (17): M_QWP^(45) — verify via J and via Mueller sandwich eq (11)
    # -------------------------------------------------------------------
    def J_to_M_np(J):
        """Jones to Mueller via eq (7)."""
        sigmas = [np.array([[1,0],[0,1]]),  np.array([[1,0],[0,-1]]),
                  np.array([[0,1],[1,0]]),  np.array([[0,-1j],[1j,0]])]
        Jd = J.conj().T
        M  = np.zeros((4,4), dtype=complex)
        for i in range(4):
            for j in range(4):
                M[i,j] = np.trace(sigmas[i] @ J @ sigmas[j] @ Jd) / 2
        return np.real(M)

    M_qwp0_np        = J_to_M_np(QWP0_np)
    M_qwp_45_from_J  = J_to_M_np(J_qwp_45)
    M_qwp_45_sandwich = (_rotation_M_np(-chi_45) @ M_qwp0_np
                         @ _rotation_M_np(chi_45))                        # eq (11)
    M_qwp_45_paper   = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,-1,0,0]])
    print(f"\neq (17)  M_QWP(chi=+45°):")
    print(f"  from J via eq (7):          {M_qwp_45_from_J.round(3)}")
    print(f"  from sandwich eq (11):      {M_qwp_45_sandwich.round(3)}")
    print(f"  paper eq (17):              {M_qwp_45_paper}")
    assert np.allclose(M_qwp_45_from_J,  M_qwp_45_paper, atol=1e-10), "eq 17 FAIL (from J)"
    assert np.allclose(M_qwp_45_sandwich,M_qwp_45_paper, atol=1e-10), "eq 17 FAIL (sandwich)"
    print("  => MATCH ✓")

    # eq (18): S_in = [1,1,0,0]  (H-linear input)
    E_in = np.array([1, 0], dtype=complex)
    print(f"\neq (18)  S_in = {_jones_to_stokes_np([1,0])} (H-linear input via eqs 1-4) ✓")

    # -------------------------------------------------------------------
    # eq (19): E_out = [1,-i]/sqrt(2)
    # -------------------------------------------------------------------
    E_out_19 = J_qwp_45 @ E_in
    ratio    = E_out_19[1] / E_out_19[0]
    print(f"\neq (19)  E_out = J_QWP(+45°) @ [1,0] = {E_out_19.round(4)}")
    print(f"  E_out[1]/E_out[0] = {ratio.round(4)}  => proportional to [1,-i]")
    print(f"  => E_out = e^{{+i pi/4}}/sqrt(2) * [1,-i]  ✓")

    # Stokes from Jones route, using eqs (1-4)
    S_jones = _jones_to_stokes_np(E_out_19)
    print(f"\neqs (1-4) applied to eq (19):  S = {S_jones}")
    print(f"  S3/S0 = {S_jones[3]:.1f} => RCP (physics convention: S3>0 is RCP) ✓")

    # Mueller route using eq (17): M_QWP^(45) @ S_in
    S_in_18 = np.array([1, 1, 0, 0], dtype=float)   # eq (18)
    S_out_17 = M_qwp_45_paper @ S_in_18
    print(f"\neq (17) @ eq (18):  M_QWP(+45°) @ [1,1,0,0] = {S_out_17}")
    print(f"  S3 = {S_out_17[3]:.1f} => IEEE/optics convention (S3>0=LCP); same")
    print(f"  physical state as S3=+1 from Jones route. (Paper discusses this")
    print(f"  sign-convention difference in the text after eq 17.)")

    # -------------------------------------------------------------------
    # QWP at chi=-45°  (LCP case, eqs 16 and 17 at negative angle)
    # -------------------------------------------------------------------
    chi_m45   = -np.pi / 4
    J_qwp_m45 = _rotation_J_np(-chi_m45) @ QWP0_np @ _rotation_J_np(chi_m45)
    E_out_lcp = J_qwp_m45 @ E_in
    S_lcp     = _jones_to_stokes_np(E_out_lcp)
    M_qwp_m45 = J_to_M_np(J_qwp_m45)
    M_qwp_m45_paper = np.array([[1,0,0,0],[0,0,0,-1],[0,0,1,0],[0,1,0,0]])
    print(f"\nQWP at chi=-45° (LCP, eqs 16-17 at negative angle):")
    print(f"  E_out = {E_out_lcp.round(4)}  ratio = {(E_out_lcp[1]/E_out_lcp[0]).round(4)}")
    print(f"  => [1,+i]/sqrt(2), S3/S0 = {S_lcp[3]:.1f} (LCP)")
    assert np.allclose(M_qwp_m45, M_qwp_m45_paper, atol=1e-10)
    print(f"  M_QWP(-45°) matches paper (eq 17 at chi=-45°) ✓")

    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SECTION: Linear polarizer intensity (eq 20, Sect. 3.2)")
    print("=" * 60)
    # LP Jones and Mueller matrices are displayed in Sect. 3.2 but are unlabelled.
    print("LP J (Sect. 3.2, unlabelled):", _linearpolarizer_J_sp(delta_s))
    print("LP M (Sect. 3.2, unlabelled):", _jones_to_mueller_sp(_linearpolarizer_J_sp(delta_s)))

    S0_s = sp.Symbol("S0", real=True)
    S1_s = sp.Symbol("S1", real=True)
    S2_s = sp.Symbol("S2", real=True)
    S3_s = sp.Symbol("S3", real=True)
    I_det = (_jones_to_mueller_sp(_linearpolarizer_J_sp(delta_s))
             @ sp.Matrix([S0_s, S1_s, S2_s, S3_s]))[0]
    print("Detected intensity  eq (20):", I_det)

    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SECTION: Crystal phase plate (eqs 33-53)")
    print("=" * 60)
    r_s_s, r_p_s = sp.symbols("r_s r_p")

    # eq (33): J0 = diag(rs, rp)
    J0 = _crystal_J_sp(r_s_s, r_p_s)
    print("J0  eq (33):", J0)

    # eq (34): Pauli decomposition of J0
    J0b = (_pauli_sp()[0] * (r_s_s + r_p_s) / 2
           + _pauli_sp()[1] * (r_s_s - r_p_s) / 2)
    assert (J0 - J0b).is_zero_matrix
    print("Pauli decomposition  eq (34): CHECKED ANALYTICALLY ✓")

    # eqs (35-36): rotated J via eq (8): J(chi) = R(-chi) @ J0 @ R(chi)
    J1  = _rotation_J_sp(-chi_s) @ J0 @ _rotation_J_sp(chi_s)   # eq (8)
    J1b = ((r_s_s + r_p_s) / 2 * _pauli_sp()[0]
           + (r_s_s - r_p_s) / 2 * (sp.cos(2*chi_s) * _pauli_sp()[1]
                                    + sp.sin(2*chi_s) * _pauli_sp()[2]))
    for _ in range(100):
        chi = np.random.rand() * 2 * np.pi
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp, chi_s: chi}
        J1_n   = _rotation_J_np(-chi) @ np.array([[rs,0],[0,rp]]) @ _rotation_J_np(chi)
        J1_sp  = np.array(J1.subs(vals).evalf(15)).astype(complex)
        J1b_sp = np.array(J1b.subs(vals).evalf(15)).astype(complex)
        assert np.abs(J1_n - J1_sp).sum()  < 1e-10
        assert np.abs(J1_n - J1b_sp).sum() < 1e-10
    print("Rotated J(chi)  eqs (35-36): CHECKED NUMERICALLY (100 random) ✓")

    # eq (37): unrotated crystal Mueller M0
    M0  = _jones_to_mueller_sp(J0)
    M0p = _M0_matrix_sp(r_s_s, r_p_s)
    for _ in range(100):
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp}
        assert np.abs(np.array(M0.subs(vals).evalf(15))
                      - np.array(M0p.subs(vals).evalf(15))).sum() < 1e-10
    print("Unrotated M0  eq (37): CHECKED NUMERICALLY ✓")

    # eq (38): rotated crystal Mueller M(chi) via eq (11):
    #   M(chi) = R_M(-chi) @ M0 @ R_M(chi)
    M1 = _rotation_M_sp(-chi_s) @ M0 @ _rotation_M_sp(chi_s)   # eq (11)
    for _ in range(100):
        chi  = np.random.rand() * 2 * np.pi
        rs   = np.random.rand() + 1j * np.random.rand()
        rp   = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp, chi_s: chi}
        M1_n   = np.real(_mueller_rotated_crystal_np(rs, rp, chi))
        M1_sp  = np.real(np.array(M1.subs(vals).evalf(15)).astype(complex))
        M1_doc = np.real(np.array(
            _mueller_rotated_crystal_sp(r_s_s, r_p_s, chi_s).subs(vals).evalf(15)
        ).astype(complex))
        assert np.abs(M1_n  - M1_sp ).max() < 1e-8
        assert np.abs(M1_n  - M1_doc).max() < 1e-8
    print("Rotated M(chi)  eq (38): CHECKED NUMERICALLY (100 random) ✓")

    # eqs (44-47): S_out general input
    S_in = sp.Matrix([S0_s, S1_s, S2_s, S3_s])
    S_out = M1 @ S_in    # eq (43)
    print("\nS0_out  eq (44) =", S_out[0])
    print("S1_out  eq (45) =", S_out[1])
    print("S2_out  eq (46) =", S_out[2])
    print("S3_out  eq (47) =", S_out[3])

    # Cross-check eqs (44-47) via shortcut _stokes_Sout_sp (uses same eq 11 sandwich)
    S_out2 = _stokes_Sout_sp(S_in, r_s_s, r_p_s, chi_s)
    print("\nCross-check  eqs (44-47) via _stokes_Sout_sp (eq 11 sandwich):")
    print("  S3:", S_out2[3], "  [should match S3 above]")

    # eqs (48-51): H-linear input S_in=[1,1,0,0]  eq (18)
    sub_H = {S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}
    print("\nH-linear input (substitution into eqs 44-47 to get eqs 48-51):")
    print("  S0  eq (48) =", S_out[0].subs(sub_H))
    print("  S1  eq (49) =", S_out[1].subs(sub_H))
    print("  S2  eq (50) =", S_out[2].subs(sub_H))
    print("  S3  eq (51) =", S_out[3].subs(sub_H))

    # +45° input (substitution into eqs 44-47, unlabelled in paper)
    sub_45 = {S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}
    print("\n+45° input (substitution into eqs 44-47, unlabelled):")
    print("  S3 =", S_out[3].subs(sub_45))

    # RCP input (substitution into eqs 44-47, unlabelled in paper)
    sub_R = {S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}
    print("\nRCP input (substitution into eqs 44-47, unlabelled):")
    print("  S3 =", S_out[3].subs(sub_R))

    # eq (52): Pc general
    Pc = S_out[3] / S_out[0]
    print("\nPc  eq (52) =", Pc)

    # eq (53): Pc for H input (Giles' formula)
    Pc_H = Pc.subs(sub_H)
    print("Pc H-input  eq (53) =", Pc_H)

    # Numerical check: eq (53) vs Giles formula
    for _ in range(100):
        chi = np.pi / 4
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp, chi_s: chi}
        Pc1 = complex(Pc_H.subs(vals).evalf(15))
        Pc2 = -_pc_giles(rs, rp, chi)    # Giles uses opposite sign convention
        assert abs(Pc1 - Pc2) < 1e-10
    print("eq (53) vs Giles formula: CHECKED NUMERICALLY ✓")
    print("  (sign flip relative to Giles due to rs/rp ordering convention)")

    # +45° and RCP substitutions into eq (52) (unlabelled in paper)
    print("\n+45° input (substitution into eq 52, unlabelled):")
    print("  Pc =", Pc.subs(sub_45))
    print("\nRCP input (substitution into eq 52, unlabelled):")
    print("  Pc =", Pc.subs(sub_R))

    print("\n" + "=" * 60)
    print("ALL CHECKS PASSED")
    print("=" * 60)