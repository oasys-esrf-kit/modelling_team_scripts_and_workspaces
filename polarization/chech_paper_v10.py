# check_paper_v10.py
#
# Numerically and symbolically verifies the Jones/Mueller formulas in:
#   "Converting X-ray Linear Polarization to Circular Polarization
#    Using a Crystal Wave Plate"
#
# Equation numbers refer to main10.tex (v10). Map of the equations checked:
#   (1)     physical wave  E(t,r) = Re[(Ex x + Ey y) e^{-i(wt-kz)}]
#   (2)     real field components  |E|cos(tau - phi)
#   (3)     Jones vector (complex envelope of eq 1)
#   (4-7)   Stokes definitions; S3 = 2 Im(Ex* Ey)
#   (8)     Pauli matrices sigma_0..3
#   (9)     rho = E E† = 1/2 sum Sk sigma_k        [density matrix]
#   (10)    Sk = Tr(sigma_k rho)
#   (11)    rho = 1/2 [[S0+S1, S2-iS3],[S2+iS3, S0-S1]]
#   (14)    M_mn = 1/2 Tr(sigma_m J sigma_n J†)   [Pauli trace formula]
#   (16)    J(chi) = R(-chi) J0 R(chi)             [Jones rotation sandwich]
#   (17)    R(chi) = [[cos,sin],[-sin,cos]]        [passive rotation matrix]
#   (22)    R_M(chi) 4x4 Mueller rotation matrix
#   (23)    M(chi) = R_M(-chi) M0 R_M(chi)         [Mueller rotation sandwich]
#   (24)    J_ret^(0) = e^{-id/2}[[1,0],[0,e^{id}]]  [general retarder, chi=0]
#   (25)    M_ret^(0): delta-rotation in the (S2,S3) plane   [NEW in v10]
#   (26)    J_QWP^(0) = e^{-ipi/4}[[1,0],[0,i]]    [QWP at chi=0]
#   (27)    M_QWP^(0)                              [QWP Mueller at chi=0]
#   (28)    M_ret(chi,delta) full 4x4              [general rotated retarder]
#   (29)    J_QWP^(45) = (1/sqrt2)[[1,-i],[-i,1]]  [QWP at chi=+45, incl. e^{-ipi/4}]
#   (30)    M_QWP^(45) = [[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,-1,0,0]]
#   (32)    S_in = [1,1,0,0]^T                     [H-linear input Stokes]
#   (34)    E_out = [1,-i]/sqrt2, S3 = -1          [RCP output Jones vector]
#   (39)    rho_RCP = 1/2 [[1, i],[-i, 1]]         [NEW in v10]
#   (40-44) chi = -45: E_out = [1,+i]/sqrt2, S3 = +1 (LCP)
#   (45)    rho_LCP = 1/2 [[1,-i],[ i, 1]]         [NEW in v10]
#   Table 2 basic states: H, +45, V, RCP, LCP      [NEW in v10]
#   (48)    I(theta) = 1/2(S0 + S1 cos2t + S2 sin2t)  [LP detected intensity]
#   (49-52) IH, IV, I45, I135 at key polarizer angles
#   (53-55) S0, S1, S2 recovery from LP measurements
#   (56-57) I+, I- with QWP analyser
#   (58-59) S3 and DOCP recovery formulas
#   (61)    I(theta) = I0/2 for circular input (all theta)
#   Table 5: RCP beam (1,0,0,-1): I+ = 0, I- = I0  [NEW in v10]
#   (68)    Re(rs rp*) =  |rs||rp| cos(dphi),
#           Im(rs rp*) = -|rs||rp| sin(dphi), dphi = phi_p - phi_s  [NEW in v10]
#   (69)    J0 = diag(rs, rp)                      [crystal Jones matrix, chi=0]
#   (70)    J0 Pauli decomposition
#   (71)    J(chi) Pauli decomposition
#   (72)    J(chi) matrix form (corrected +sin2chi off-diagonals)
#   (73)    M0 4x4                                 [crystal Mueller, chi=0]
#   (75)    M(chi) full rotated 4x4                [Mueller_final]
#   (77-80) E_out general / Pauli / H input / V input (corrected signs)
#   (81)    S_out = M(chi) S_in
#   (82-85) S0..S3 output, general input
#   (86-89) S0..S3 output, H-linear input
#   (99)    Pc general formula
#   (100)   Pc for H-linear input  [= -(Giles' formula): opposite handedness conv.]
#
# Convention (same as the paper, Sec. 2.1):
#   time factor e^{-i(wt-kz)};  S3 = 2 Im(Ex* Ey);  sigma_3 = [[0,-i],[i,0]]
#   anchor: (1,+i)/sqrt2 -> S3=+1 (LCP);  (1,-i)/sqrt2 -> S3=-1 (RCP)
#   R(chi) = [[cos chi, sin chi], [-sin chi, cos chi]]   (passive; transpose of
#   the active CCW matrix).  Positive chi rotates an element's axis from +x
#   toward +y.
#   Rotated Jones:   J(chi) = R(-chi) @ J0 @ R(chi)            [eq 16]
#   Rotated Mueller: M(chi) = R_M(-chi) @ M0 @ R_M(chi)        [eq 23]

import numpy as np
import sympy as sp

# ===========================================================================
# NumPy helper functions
# ===========================================================================

def _jones_to_stokes_np(j, normalize=True):
    """Stokes vector from Jones vector using eqs (4-7).

    NOTE (v10): S3 = 2 Im(Ex* Ey), NOT 2 Im(Ex Ey*).  This is the unique
    choice consistent with Sk = Tr(sigma_k rho) and sigma_3 = [[0,-i],[i,0]]
    (eqs 5-8 of the paper), hence with every Mueller matrix derived via the
    Pauli trace formula.  With this convention (1,i)/sqrt2 has S3=+1 and is
    LCP (counterclockwise looking toward the source, Detlefs et al. 2012);
    (1,-i)/sqrt2 has S3=-1 and is RCP."""
    j = np.array(j, dtype=complex)
    j1, j2 = j[0], j[1]
    S = [
        (j1 * j1.conjugate() + j2 * j2.conjugate()).real,   # eq (4)
        (j1 * j1.conjugate() - j2 * j2.conjugate()).real,   # eq (5)
        2 * (j1.conjugate() * j2).real,                      # eq (6)
        2 * (j1.conjugate() * j2).imag,                      # eq (7)
    ]
    if normalize:
        S[1] /= S[0]; S[2] /= S[0]; S[3] /= S[0]
    return np.array(S, dtype=float)


def _rotation_J_np(chi: float) -> np.ndarray:
    """2×2 Jones rotation matrix R(chi) — eq (17)."""
    c, s = np.cos(chi), np.sin(chi)
    return np.array([[c, s], [-s, c]])


def _rotation_M_np(chi: float) -> np.ndarray:
    """4×4 Mueller rotation matrix R_M(chi) — eq (22)."""
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
    """Degree of circular polarization from eq (99)."""
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
    """Rotated crystal Mueller matrix M(chi) — eq (75)."""
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
    """2×2 Jones rotation matrix R(chi) — eq (17)."""
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c, s], [-s, c]])


def _rotation_M_sp(a):
    """4×4 Mueller rotation matrix R_M(chi) — eq (22)."""
    c2a, s2a = sp.cos(2*a), sp.sin(2*a)
    return sp.Matrix([[1,0,0,0],[0,c2a,s2a,0],[0,-s2a,c2a,0],[0,0,0,1]])


def _crystal_J_sp(r_s, r_p):
    """Crystal Jones matrix at chi=0 — eq (69)."""
    return sp.Matrix([[r_s, 0], [0, r_p]])


def _retarder_J_sp(delta, use_global_phase=False):
    """General retarder Jones matrix at chi=0 — eq (24)."""
    if use_global_phase:
        return sp.exp(-sp.I * delta / 2) * sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])
    return sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])


def _linearpolarizer_J_sp(a):
    """Linear polarizer Jones matrix at angle a (Sect. 3.2, unlabelled)."""
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c*c, c*s], [c*s, s*s]])


def _M0_matrix_sp(rs, rp):
    """Crystal Mueller matrix at chi=0 written explicitly — eq (73)."""
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
    """Convert Jones matrix to Mueller matrix via Pauli trace formula — eq (14)."""
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
    """Return list of four Pauli matrices — eq (8)."""
    return [
        sp.Matrix([[1, 0], [0,  1]]),
        sp.Matrix([[1, 0], [0, -1]]),
        sp.Matrix([[0, 1], [1,  0]]),
        sp.Matrix([[0, -sp.I], [sp.I, 0]]),
    ]


def _mueller_rotated_crystal_sp(rs, rp, theta):
    """Rotated crystal Mueller matrix M(chi) — eq (75)."""
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
    Full S_out after crystal at angle chi — eqs (81-85).
    Shortcut: builds M(chi) via eq (23) then applies to S_in.
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
    print("R_J(chi)  eq (17) :", _rotation_J_sp(chi_s))
    print("R_M(chi)  eq (22):", _rotation_M_sp(chi_s))

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

    print("\nRetarder J at chi=0  eq (24):", _retarder_J_sp(delta_s))
    M_ret0_sp = _jones_to_mueller_sp(_retarder_J_sp(delta_s))
    print("Retarder M at chi=0 (derived via eq 7):", M_ret0_sp)

    delta_qwp = np.pi / 2
    print("\nQWP J at chi=0  eq (26):",
          _retarder_J_sp(delta_s).subs(delta_s, delta_qwp).evalf(5))
    print("QWP M at chi=0  eq (27):",
          _jones_to_mueller_sp(_retarder_J_sp(delta_s)).subs(delta_s, delta_qwp).evalf(5))

    # General rotated retarder Mueller — eq (28): R_M(-chi) @ M_ret0 @ R_M(chi)
    M_ret_rot = (_rotation_M_sp(-chi_s) @ M_ret0_sp @ _rotation_M_sp(chi_s))
    print("\nGeneral rotated retarder M  eq (28):", M_ret_rot)

    # -------------------------------------------------------------------
    # eq (29): J_QWP^(45) = R(-45) @ J_QWP^(0) @ R(45)  [eq 8]
    #   Result: e^{+i pi/4}/sqrt(2) * [[1,-i],[-i,1]]
    # -------------------------------------------------------------------
    chi_45  = np.pi / 4
    QWP0_np = np.array([[1, 0], [0, 1j]], dtype=complex)   # eq (26), no global phase

    J_qwp_45 = _rotation_J_np(-chi_45) @ QWP0_np @ _rotation_J_np(chi_45)   # eq (16)
    factor   = J_qwp_45[0, 0]   # = e^{+i pi/4}/sqrt(2)
    print(f"\neq (16)  J_QWP(chi=+45°) = R(-45)@J_QWP^(0)@R(45)  [eq 8]:")
    print(f"  computed = {J_qwp_45.round(4)}")
    print(f"  global phase factor = {factor:.4f}  (= e^{{+i pi/4}}/sqrt(2))")
    print(f"  J / factor = {(J_qwp_45/factor).round(4)}")
    print(f"  => J = e^{{+i pi/4}}/sqrt(2) * [[1,-i],[-i,1]]  ✓")

    # -------------------------------------------------------------------
    # eq (30): M_QWP^(45) — verify via J and via Mueller sandwich eq (23)
    # -------------------------------------------------------------------
    def J_to_M_np(J):
        """Jones to Mueller via eq (14)."""
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
                         @ _rotation_M_np(chi_45))                        # eq (23)
    M_qwp_45_paper   = np.array([[1,0,0,0],[0,0,0,1],[0,0,1,0],[0,-1,0,0]])
    print(f"\neq (17)  M_QWP(chi=+45°):")
    print(f"  from J via eq (14):          {M_qwp_45_from_J.round(3)}")
    print(f"  from sandwich eq (23):      {M_qwp_45_sandwich.round(3)}")
    print(f"  paper eq (30):              {M_qwp_45_paper}")
    assert np.allclose(M_qwp_45_from_J,  M_qwp_45_paper, atol=1e-10), "eq 17 FAIL (from J)"
    assert np.allclose(M_qwp_45_sandwich,M_qwp_45_paper, atol=1e-10), "eq 17 FAIL (sandwich)"
    print("  => MATCH ✓")

    # eq (32): S_in = [1,1,0,0]  (H-linear input)
    E_in = np.array([1, 0], dtype=complex)
    print(f"\neq (18)  S_in = {_jones_to_stokes_np([1,0])} (H-linear input via eqs 1-4) ✓")

    # -------------------------------------------------------------------
    # eq (34): E_out = [1,-i]/sqrt(2)
    # -------------------------------------------------------------------
    E_out_19 = J_qwp_45 @ E_in
    ratio    = E_out_19[1] / E_out_19[0]
    print(f"\neq (19)  E_out = J_QWP(+45°) @ [1,0] = {E_out_19.round(4)}")
    print(f"  E_out[1]/E_out[0] = {ratio.round(4)}  => proportional to [1,-i]")
    assert abs(ratio + 1j) < 1e-12

    # Stokes from Jones route, using eqs (4-7)
    S_jones = _jones_to_stokes_np(E_out_19)
    print(f"\neqs (1-4) applied to eq (34):  S = {S_jones}")
    print(f"  S3/S0 = {S_jones[3]:.1f} => RCP (S3=-1 is RCP; S3=+1 is LCP) ✓")
    assert abs(S_jones[3] + 1) < 1e-12

    # Mueller route using eq (30): M_QWP^(45) @ S_in
    S_in_18 = np.array([1, 1, 0, 0], dtype=float)   # eq (32)
    S_out_17 = M_qwp_45_paper @ S_in_18
    print(f"\neq (17) @ eq (32):  M_QWP(+45°) @ [1,1,0,0] = {S_out_17}")
    assert np.allclose(S_out_17, S_jones, atol=1e-12), "Jones and Mueller routes disagree!"
    print(f"  Jones route and Mueller route now agree exactly: S3 = -1 (RCP) ✓")

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
    print(f"  => [1,+i]/sqrt(2), S3/S0 = {S_lcp[3]:.1f} (LCP: S3=+1)")
    assert abs(E_out_lcp[1]/E_out_lcp[0] - 1j) < 1e-12 and abs(S_lcp[3] - 1) < 1e-12
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
    print("Detected intensity  eq (48):", I_det)

    # -----------------------------------------------------------------------
    print("\n" + "=" * 60)
    print("SECTION: Crystal phase plate (eqs 33-53)")
    print("=" * 60)
    r_s_s, r_p_s = sp.symbols("r_s r_p")

    # eq (69): J0 = diag(rs, rp)
    J0 = _crystal_J_sp(r_s_s, r_p_s)
    print("J0  eq (69):", J0)

    # eq (70): Pauli decomposition of J0
    J0b = (_pauli_sp()[0] * (r_s_s + r_p_s) / 2
           + _pauli_sp()[1] * (r_s_s - r_p_s) / 2)
    assert (J0 - J0b).is_zero_matrix
    print("Pauli decomposition  eq (70): CHECKED ANALYTICALLY ✓")

    # eqs (71-72): rotated J via eq (16): J(chi) = R(-chi) @ J0 @ R(chi)
    J1  = _rotation_J_sp(-chi_s) @ J0 @ _rotation_J_sp(chi_s)   # eq (16)
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
    print("Rotated J(chi)  eqs (71-72): CHECKED NUMERICALLY (100 random) ✓")

    # eq (73): unrotated crystal Mueller M0
    M0  = _jones_to_mueller_sp(J0)
    M0p = _M0_matrix_sp(r_s_s, r_p_s)
    for _ in range(100):
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp}
        assert np.abs(np.array(M0.subs(vals).evalf(15))
                      - np.array(M0p.subs(vals).evalf(15))).sum() < 1e-10
    print("Unrotated M0  eq (73): CHECKED NUMERICALLY ✓")

    # eq (75): rotated crystal Mueller M(chi) via eq (23):
    #   M(chi) = R_M(-chi) @ M0 @ R_M(chi)
    M1 = _rotation_M_sp(-chi_s) @ M0 @ _rotation_M_sp(chi_s)   # eq (23)
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
    print("Rotated M(chi)  eq (75): CHECKED NUMERICALLY (100 random) ✓")

    # eqs (82-85): S_out general input
    S_in = sp.Matrix([S0_s, S1_s, S2_s, S3_s])
    S_out = M1 @ S_in    # eq (81)
    print("\nS0_out  eq (82) =", S_out[0])
    print("S1_out  eq (83) =", S_out[1])
    print("S2_out  eq (84) =", S_out[2])
    print("S3_out  eq (85) =", S_out[3])

    # Cross-check eqs (82-85) via shortcut _stokes_Sout_sp (uses same eq 11 sandwich)
    S_out2 = _stokes_Sout_sp(S_in, r_s_s, r_p_s, chi_s)
    print("\nCross-check  eqs (82-85) via _stokes_Sout_sp (eq 11 sandwich):")
    print("  S3:", S_out2[3], "  [should match S3 above]")

    # eqs (86-89): H-linear input S_in=[1,1,0,0]  eq (32)
    sub_H = {S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}
    print("\nH-linear input (substitution into eqs 44-47 to get eqs 48-51):")
    print("  S0  eq (86) =", S_out[0].subs(sub_H))
    print("  S1  eq (87) =", S_out[1].subs(sub_H))
    print("  S2  eq (88) =", S_out[2].subs(sub_H))
    print("  S3  eq (89) =", S_out[3].subs(sub_H))

    # +45° input (substitution into eqs 44-47, unlabelled in paper)
    sub_45 = {S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}
    print("\n+45° input (substitution into eqs 44-47, unlabelled):")
    print("  S3 =", S_out[3].subs(sub_45))

    # RCP input (substitution into eqs 44-47, unlabelled in paper)
    sub_R = {S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}
    print("\nRCP input (substitution into eqs 44-47, unlabelled):")
    print("  S3 =", S_out[3].subs(sub_R))

    # eq (99): Pc general
    Pc = S_out[3] / S_out[0]
    print("\nPc  eq (99) =", Pc)

    # eq (100): Pc for H input (Giles' formula)
    Pc_H = Pc.subs(sub_H)
    print("Pc H-input  eq (100) =", Pc_H)

    # Numerical check: eq (100) vs Giles formula
    for _ in range(100):
        chi = np.pi / 4
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        vals = {r_s_s: rs, r_p_s: rp, chi_s: chi}
        Pc1 = complex(Pc_H.subs(vals).evalf(15))
        Pc2 = -_pc_giles(rs, rp, chi)    # Giles uses opposite handedness convention for Pc
        assert abs(Pc1 - Pc2) < 1e-10
    print("eq (100) vs Giles formula: CHECKED NUMERICALLY ✓")
    print("  Giles' expression equals -Pc of eq (100): opposite handedness")
    print("  convention for the circular polarization rate, no physical difference.")

    # +45° and RCP substitutions into eq (99) (unlabelled in paper)
    print("\n+45° input (substitution into eq 52, unlabelled):")
    print("  Pc =", Pc.subs(sub_45))
    print("\nRCP input (substitution into eq 52, unlabelled):")
    print("  Pc =", Pc.subs(sub_R))

    print("\n" + "=" * 60)
    print("ALL CHECKS PASSED")
    print("=" * 60)
# ===========================================================================
# NEW in v10: convention cross-checks (added after coauthor discussion)
# ===========================================================================

def _convention_checks():
    print("\n" + "=" * 60)
    print("SECTION: Convention cross-checks (v10)")
    print("=" * 60)

    # (a) eq (7) vs eq (10): S3 = 2 Im(Ex* Ey) must equal Tr(sigma_3 rho)
    sigma3 = np.array([[0, -1j], [1j, 0]])
    for _ in range(50):
        E = np.random.randn(2) + 1j * np.random.randn(2)
        rho = np.outer(E, E.conj())
        s3_def   = 2 * (E[0].conjugate() * E[1]).imag
        s3_trace = np.real(np.trace(sigma3 @ rho))
        assert abs(s3_def - s3_trace) < 1e-12
    print("(a) S3 = 2 Im(Ex* Ey) == Tr(sigma_3 rho): eqs (7) and (10) consistent ✓")

    # (b) handedness anchors (Detlefs et al. 2012, footnote 1)
    S_lcp = _jones_to_stokes_np([1/np.sqrt(2),  1j/np.sqrt(2)])
    S_rcp = _jones_to_stokes_np([1/np.sqrt(2), -1j/np.sqrt(2)])
    assert abs(S_lcp[3] - 1) < 1e-12 and abs(S_rcp[3] + 1) < 1e-12
    print("(b) (1,+i)/sqrt2 -> S3=+1 (LCP);  (1,-i)/sqrt2 -> S3=-1 (RCP) ✓")

    # (c) element-rotation convention: eq (16) sandwich puts a polarizer
    #     'rotated by +theta' with transmission axis (cos t, sin t),
    #     so the detected intensity is (S0 + S1 cos2t + S2 sin2t)/2  [eq 20].
    #     The opposite sandwich R(chi) J R(-chi) gives axis (cos t, -sin t)
    #     and hence  -S2 sin2t  (the sign observed in CD's notebook).
    LP0 = np.array([[1, 0], [0, 0]], dtype=complex)
    for _ in range(50):
        th = np.random.rand() * np.pi
        t_axis = np.array([np.cos(th), np.sin(th)])
        LP_paper   = _rotation_J_np(-th) @ LP0 @ _rotation_J_np(th)
        LP_flipped = _rotation_J_np(th) @ LP0 @ _rotation_J_np(-th)
        assert np.allclose(LP_paper, np.outer(t_axis, t_axis), atol=1e-12)
        assert np.allclose(LP_flipped,
                           np.outer([np.cos(th), -np.sin(th)],
                                    [np.cos(th), -np.sin(th)]), atol=1e-12)
        # intensity through paper-convention polarizer for random input
        E = np.random.randn(2) + 1j * np.random.randn(2)
        S = _jones_to_stokes_np(E, normalize=False)
        I_direct = abs(t_axis @ E) ** 2
        I_eq20   = 0.5 * (S[0] + S[1]*np.cos(2*th) + S[2]*np.sin(2*th))
        assert abs(I_direct - I_eq20) < 1e-10
    print("(c) eq (16) sandwich == polarizer with axis (cos t, sin t);")
    print("    eq (48) intensity formula verified; opposite sandwich = axis at -t ✓")

    # (d) chi -> -chi parity of the H-input outputs (Sec. 5 discussion):
    #     P1 even, P2 and P3 odd.
    for _ in range(50):
        chi = np.random.rand() * np.pi
        rs  = np.random.rand() + 1j * np.random.rand()
        rp  = np.random.rand() + 1j * np.random.rand()
        Sp = np.real(_mueller_rotated_crystal_np(rs, rp,  chi) @ [1, 1, 0, 0])
        Sm = np.real(_mueller_rotated_crystal_np(rs, rp, -chi) @ [1, 1, 0, 0])
        assert abs(Sp[1] - Sm[1]) < 1e-12          # P1 even
        assert abs(Sp[2] + Sm[2]) < 1e-12          # P2 odd
        assert abs(Sp[3] + Sm[3]) < 1e-12          # P3 odd
    print("(d) chi -> -chi: P1 even, P2/P3 odd for H input (basis of Sec. 5) ✓")


if __name__ == "__main__":
    _convention_checks()

# ===========================================================================
# NEW in v10: full-paper re-check against main10 (new formulas, new numbering)
# ===========================================================================

def _full_paper_recheck_v10():
    print("\n" + "=" * 60)
    print("SECTION: Full-paper re-check vs main10")
    print("=" * 60)
    rng = np.random.default_rng(0)
    sig = [np.eye(2), np.diag([1, -1]),
           np.array([[0, 1], [1, 0]]), np.array([[0, -1j], [1j, 0]])]

    def J2M(J):
        return np.array([[np.real(np.trace(sig[m] @ J @ sig[n] @ J.conj().T)) / 2
                          for n in range(4)] for m in range(4)])

    # --- eqs (1)-(2): real field components from the complex envelope
    for _ in range(20):
        Ex, Ey = rng.normal(size=2) + 1j * rng.normal(size=2)
        tau = rng.uniform(0, 2 * np.pi)     # tau = wt - kz
        assert abs(np.real(Ex * np.exp(-1j * tau))
                   - abs(Ex) * np.cos(tau - np.angle(Ex))) < 1e-12
        assert abs(np.real(Ey * np.exp(-1j * tau))
                   - abs(Ey) * np.cos(tau - np.angle(Ey))) < 1e-12
    print("eqs (1)-(2): Re[E e^{-i tau}] == |E| cos(tau - phi) ✓")

    # --- eqs (9)-(11): rho = E E† = 1/2 sum Sk sigma_k
    #                       = 1/2 [[S0+S1, S2-iS3],[S2+iS3, S0-S1]]
    for _ in range(20):
        E = rng.normal(size=2) + 1j * rng.normal(size=2)
        rho = np.outer(E, E.conj())
        S = _jones_to_stokes_np(E, normalize=False)
        rho_sum = 0.5 * sum(S[k] * sig[k] for k in range(4))
        rho_11  = 0.5 * np.array([[S[0] + S[1], S[2] - 1j * S[3]],
                                  [S[2] + 1j * S[3], S[0] - S[1]]])
        assert np.allclose(rho, rho_sum, atol=1e-12)
        assert np.allclose(rho, rho_11,  atol=1e-12)
    print("eqs (9)-(11): rho = E E† = ½ΣSkσk = ½[[S0+S1,S2-iS3],[S2+iS3,S0-S1]] ✓")

    # --- eq (25): aligned Mueller matrix of the general retarder
    for _ in range(20):
        d = rng.uniform(0, 2 * np.pi)
        Jr = np.exp(-1j * d / 2) * np.diag([1, np.exp(1j * d)])   # eq (24)
        M25 = np.array([[1, 0, 0, 0],
                        [0, 1, 0, 0],
                        [0, 0, np.cos(d), -np.sin(d)],
                        [0, 0, np.sin(d),  np.cos(d)]])
        assert np.allclose(J2M(Jr), M25, atol=1e-12)
    # delta = pi/2 recovers M_QWP^(0), eq (27)
    M27 = np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    assert np.allclose(J2M(np.exp(-1j*np.pi/4)*np.diag([1, 1j])), M27, atol=1e-12)
    print("eq (25): M_ret^(0) = delta-rotation in (S2,S3); delta=pi/2 -> eq (27) ✓")

    # --- eqs (34)-(45) and Table 2: the five basic states,
    #     Jones vs Stokes vs density matrix
    states = {
        "H  (0 deg)":  (np.array([1, 0]),               [1,  1, 0,  0]),
        "+45 linear":  (np.array([1, 1]) / np.sqrt(2),  [1,  0, 1,  0]),
        "V  (90 deg)": (np.array([0, 1]),               [1, -1, 0,  0]),
        "RCP":         (np.array([1, -1j]) / np.sqrt(2), [1, 0, 0, -1]),
        "LCP":         (np.array([1,  1j]) / np.sqrt(2), [1, 0, 0, +1]),
    }
    for name, (E, S_expected) in states.items():
        S = _jones_to_stokes_np(E)
        assert np.allclose(S, S_expected, atol=1e-12), name
        rho = np.outer(E, E.conj())
        rho_11 = 0.5 * np.array([[S[0] + S[1], S[2] - 1j * S[3]],
                                 [S[2] + 1j * S[3], S[0] - S[1]]])
        assert np.allclose(rho, rho_11, atol=1e-12), name
    rho_rcp = 0.5 * np.array([[1,  1j], [-1j, 1]])                 # eq (39)
    rho_lcp = 0.5 * np.array([[1, -1j], [ 1j, 1]])                 # eq (45)
    assert np.allclose(np.outer(states["RCP"][0], states["RCP"][0].conj()),
                       rho_rcp, atol=1e-12)
    assert np.allclose(np.outer(states["LCP"][0], states["LCP"][0].conj()),
                       rho_lcp, atol=1e-12)
    print("Table 2 + eqs (39),(45): H/+45/V/RCP/LCP consistent in all three")
    print("  representations; rho_RCP = ½[[1,i],[-i,1]], rho_LCP = ½[[1,-i],[i,1]] ✓")

    # --- Table 5 / eqs (56)-(61): polarimetry of the RCP beam (1,0,0,-1)
    M_qwp_45 = np.array([[1, 0, 0, 0], [0, 0, 0, 1],
                         [0, 0, 1, 0], [0, -1, 0, 0]], dtype=float)  # eq (30)
    M_qwp_m45 = np.array([[1, 0, 0, 0], [0, 0, 0, -1],
                          [0, 0, 1, 0], [0, 1, 0, 0]], dtype=float)
    def M_LP(th):
        c2, s2 = np.cos(2 * th), np.sin(2 * th)
        return 0.5 * np.array([[1, c2, s2, 0],
                               [c2, c2**2, s2 * c2, 0],
                               [s2, s2 * c2, s2**2, 0],
                               [0, 0, 0, 0]])
    S_rcp_beam = np.array([1, 0, 0, -1], dtype=float)   # I0 = 1
    for th in np.linspace(0, np.pi, 13):
        assert abs((M_LP(th) @ S_rcp_beam)[0] - 0.5) < 1e-12       # eq (61)
    I_plus  = (M_LP(0) @ M_qwp_45  @ S_rcp_beam)[0]
    I_minus = (M_LP(0) @ M_qwp_m45 @ S_rcp_beam)[0]
    assert abs(I_plus - 0.0) < 1e-12 and abs(I_minus - 1.0) < 1e-12
    assert abs((I_plus - I_minus) - S_rcp_beam[3]) < 1e-12          # eq (58)
    print("eqs (56)-(61)/Table 5: RCP beam -> I(theta)=I0/2 for all theta,")
    print("  I+ = 0, I- = I0, S3 = I+ - I- = -I0, DOCP = -1 ✓")

    # --- eq (68): cross-product vs Delta-phi identity
    for _ in range(20):
        rs = rng.uniform(0.1, 1) * np.exp(1j * rng.uniform(-np.pi, np.pi))
        rp = rng.uniform(0.1, 1) * np.exp(1j * rng.uniform(-np.pi, np.pi))
        dphi = np.angle(rp) - np.angle(rs)                          # phi_p - phi_s
        cross = rs * np.conj(rp)
        assert abs(cross.real - abs(rs) * abs(rp) * np.cos(dphi)) < 1e-9
        assert abs(cross.imag + abs(rs) * abs(rp) * np.sin(dphi)) < 1e-9
    print("eq (68): Re(rs rp*) = |rs||rp|cos(dphi), Im(rs rp*) = -|rs||rp|sin(dphi) ✓")

    # --- eqs (77), (79), (80): corrected output Jones vectors
    for _ in range(20):
        rs = rng.normal() + 1j * rng.normal()
        rp = rng.normal() + 1j * rng.normal()
        chi = rng.uniform(0, 2 * np.pi)
        c2, s2 = np.cos(2 * chi), np.sin(2 * chi)
        Jchi = _rotation_J_np(-chi) @ np.diag([rs, rp]) @ _rotation_J_np(chi)
        E0 = rng.normal(size=2) + 1j * rng.normal(size=2)
        # eq (77) general
        E77 = (rs + rp) / 2 * E0 + (rs - rp) / 2 * np.array(
            [c2 * E0[0] + s2 * E0[1], s2 * E0[0] - c2 * E0[1]])
        assert np.allclose(Jchi @ E0, E77, atol=1e-10)
        # eq (79) H input, eq (80) V input
        E79 = np.array([rs * np.cos(chi)**2 + rp * np.sin(chi)**2,
                        (rs - rp) * np.sin(chi) * np.cos(chi)])
        E80 = np.array([(rs - rp) * np.sin(chi) * np.cos(chi),
                        rs * np.sin(chi)**2 + rp * np.cos(chi)**2])
        assert np.allclose(Jchi @ [1, 0], E79, atol=1e-10)
        assert np.allclose(Jchi @ [0, 1], E80, atol=1e-10)
        # cross-check eq (79) Stokes against Mueller route, eqs (86)-(89)
        S_j = _jones_to_stokes_np(Jchi @ [1, 0], normalize=False)
        S_m = np.real(_mueller_rotated_crystal_np(rs, rp, chi) @ [1, 1, 0, 0])
        assert np.allclose(S_j, S_m, atol=1e-9)
    print("eqs (77),(79),(80): corrected Jones outputs match sandwich; Jones and")
    print("  Mueller routes agree for the crystal at arbitrary chi ✓")

    # --- Table 7 (unified summary): ideal QWP crystal at chi=45
    r = 0.7
    Jq = _rotation_J_np(-np.pi/4) @ np.diag([r, 1j * r]) @ _rotation_J_np(np.pi/4)
    # Pauli decomposition r(1+i)/2 sigma_0 + r(1-i)/2 sigma_2
    assert np.allclose(Jq, r*(1+1j)/2 * sig[0] + r*(1-1j)/2 * sig[2], atol=1e-12)
    S_out = _jones_to_stokes_np(Jq @ [1, 0], normalize=False)
    assert np.allclose(S_out, [r**2, 0, 0, -r**2], atol=1e-12)
    print("Table 7: crystal QWP (rp = i rs) at chi=45: J = r(1+i)/2 σ0 + r(1-i)/2 σ2,")
    print("  S_out = I0 r^2 (1,0,0,-1) (RCP) ✓")

    print("\nFULL-PAPER RE-CHECK vs main10: ALL PASSED")


if __name__ == "__main__":
    _full_paper_recheck_v10()

