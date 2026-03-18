# _np means numoy (numeric)
# _sp means sympy (symbolic)

import numpy as np
import sympy as sp


# ===========================================================================
# numpy auxiliar functions
# ============================================================================
def _jones_to_stokes_np(j, normalize=True):
    j = np.array(j, dtype=complex)
    j1 = j[0]
    j2 = j[1]

    S = [
        (j1 * j1.conjugate() + j2 * j2.conjugate()).real,
        (j1 * j1.conjugate() - j2 * j2.conjugate()).real,
        2 * (j1 * j2.conjugate()).real,
        2 * (j1 * j2.conjugate()).imag
    ]

    if normalize:
        S[1] /= S[0]
        S[2] /= S[0]
        S[3] /= S[0]
    return np.array(S, dtype=float)

def _rotation_J_np(chi: float) -> np.ndarray:
    """2x2 Jones rotation matrix R_J(chi)."""
    c, s = np.cos(chi), np.sin(chi)
    return np.array([[c, s], [-s, c]]) # as in the paper!

def _rotation_M_np(chi: float) -> np.ndarray:
    """4x4 Mueller rotation matrix R_M(chi) — rotates (S1,S2) by 2*chi."""
    c2, s2 = np.cos(2.0 * chi), np.sin(2.0 * chi)
    return np.array([[1,   0,   0, 0],
                     [0,  c2,  s2, 0],
                     [0, -s2,  c2, 0],
                     [0,   0,   0, 1]], dtype=np.float64)

def _pc_giles(rs, rp, chi): # circular polarization factor as in Giles thesis
    SinPhi = (rs * rp.conjugate()).imag / np.abs(rs) / np.abs(rp)
    a = np.abs(rs) / np.abs(rp)
    return - a * SinPhi * np.sin(2 * chi) / (np.cos(chi)**2 + a**2 * np.sin(chi)**2)

def _pc_np(rs, rp, theta, s_in):
    """
    Calculates Circular Polarizance (Pc) using NumPy with equation (89) in our paper.

    Parameters:
    rs, rp (complex): Complex reflection/transmission coefficients.
    theta (float or ndarray): Rotation angle(s) in radians.
    s_in (list or ndarray): Input Stokes vector [S0, S1, S2, S3].
    """
    # 1. Trig components (supports arrays of theta)
    c2 = np.cos(2 * theta)
    s2 = np.sin(2 * theta)

    # 2. Complex term calculations
    # Magnitude squares: |rs|^2 and |rp|^2
    rs_sq = np.real(rs * np.conj(rs))
    rp_sq = np.real(rp * np.conj(rp))

    # Cross terms: rs * rp*
    cross_term = rs * np.conj(rp)
    re_cp = np.real(cross_term)
    im_cp = np.imag(cross_term)

    # 3. Input Stokes components
    s0_i, s1_i, s2_i, s3_i = s_in

    # 4. Numerator calculation
    # Im(rs*rp*)(s2*S1 - c2*S2) + Re(rs*rp*)*S3
    num = im_cp * (s2 * s1_i - c2 * s2_i) + re_cp * s3_i

    # 5. Denominator calculation (Total Output Intensity S0_out)
    term_plus = (rs_sq + rp_sq) / 2
    term_minus = (rs_sq - rp_sq) / 2
    den = term_plus * s0_i + term_minus * (c2 * s1_i + s2 * s2_i)

    return num / den

def mueller_rotated_crystal_document_np(
        rs: complex, rp: complex, theta: float) -> np.ndarray:
    """
    Returns the general rotated Mueller matrix M_chi numerically (numpy).
    theta: the rotation angle
    rs, rp: complex reflection/transmission coefficients
    """
    # 1. Define trig shortcuts
    c2 = np.cos(2 * theta)
    s2 = np.sin(2 * theta)

    # 2. Define fundamental complex terms
    rs_sq = rs * np.conjugate(rs)
    rp_sq = rp * np.conjugate(rp)

    # Using your derived identities
    term_plus = (rs_sq + rp_sq) / 2
    term_minus = (rs_sq - rp_sq) / 2

    cross_term = rs * sp.conjugate(rp)
    re_cp = np.real(cross_term)
    im_cp = np.imag(cross_term)

    # 3. Construct the matrix row by row
    row0 = [
        term_plus,
        term_minus * c2,
        term_minus * s2,
        0
    ]

    row1 = [
        term_minus * c2,
        term_plus * c2 ** 2 + re_cp * s2 ** 2,
        (term_plus - re_cp) * s2 * c2,
        -im_cp * s2
    ]

    row2 = [
        term_minus * s2,
        (term_plus - re_cp) * s2 * c2,
        term_plus * s2 ** 2 + re_cp * c2 ** 2,
        im_cp * c2
    ]

    row3 = [
        0,
        im_cp * s2,
        -im_cp * c2,
        re_cp
    ]

    M_chi = np.array([row0, row1, row2, row3], dtype=complex)

    return M_chi

# ===========================================================================
# sympy auxiliar functions
# ============================================================================
def _rotation_J_sp(a): # Rotation matrix 2x2 (for Jones)
    c, s = sp.cos(a), sp.sin(a)
    # return sp.Matrix([[c, -s], [s, c]])
    return sp.Matrix([[c, s], [-s, c]]) # as in the paper

def _rotation_M_sp(a): # Rotation matrix 4x4 (for Mueller)
    c2a, s2a = sp.cos(2*a), sp.sin(2*a)
    return sp.Matrix([[1,0,0,0],[0,c2a,s2a,0],[0,-s2a,c2a,0],[0,0,0,1]])

def _crystal_J_sp(r_s, r_p): # crystal plate Jones matrix (unrotated)
    return sp.Matrix([[r_s, 0], [0, r_p]])

def _retarder_J_sp(delta, use_constant=False): # phase plate Jones matrix (unrotated, any retardation)
    if use_constant:
        return sp.exp(-sp.I * delta / 2) * sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])
    else:
        return sp.Matrix([[1, 0], [0, sp.exp(sp.I * delta)]])

def _linearpolarizer_J_sp(a): # Jones matrix for the linear polarizer at angle θ to the horizontal (eq 36)
    c, s = sp.cos(a), sp.sin(a)
    return sp.Matrix([[c*c,c*s],[c*s,s*s]])

def _M0_matrix_sp(rs, rp): # Crystal (unrotated) Mueller matrix as written in the paper (eq 62)
    """
    Returns the symbolic unrotated 4x4 Mueller matrix M0 based on complex
    reflection/transmission coefficients rs and rp.
    """
    # Define terms for clarity
    # rs * rs* and rp * rp* are essentially the squared magnitudes
    rs_sq = rs * rs.conjugate()
    rp_sq = rp * rp.conjugate()

    term_plus = (rs_sq + rp_sq) / 2
    term_minus = (rs_sq - rp_sq) / 2

    # Off-diagonal complex interactions
    cross_term = rs * rp.conjugate()
    re_cross = sp.re(cross_term)
    im_cross = sp.im(cross_term)

    # Construct the matrix
    M0 = sp.Matrix([
        [term_plus, term_minus, 0, 0],
        [term_minus, term_plus, 0, 0],
        [0, 0, re_cross, im_cross],
        [0, 0, -im_cross, re_cross]
    ])

    return M0

def _stokes_Sout_sp(S_in_s, r_s_s, r_p_s, chi_s):
    # a shortcut for full Sout after crystal plate, useful for experimental check
    # eqs 72-75
    M1 = _rotation_M_sp(-chi_s) @ _jones_to_mueller_sp(_crystal_J_sp(r_s_s, r_p_s)) @ _rotation_M_sp(chi_s)
    return M1 @ S_in_s

def _jones_to_mueller_sp(J):
    """
    Converts a Jones matrix J to a Mueller matrix M using a Pauli basis:
    Sigma[0]: Identity
    Sigma[1]: Diagonal (Standard Z)
    Sigma[2]: Off-diagonal real (Standard X)
    Sigma[3]: Off-diagonal imaginary (Standard Y)
    """
    # 1. Define your specific basis
    I_unit = sp.I
    sigmas = [
        sp.Matrix([[1, 0], [0, 1]]),  # Sigma 0
        sp.Matrix([[1, 0], [0, -1]]),  # Sigma 1 (Diagonal)
        sp.Matrix([[0, 1], [1, 0]]),  # Sigma 2 (Real off-diagonal)
        sp.Matrix([[0, -I_unit], [I_unit, 0]])  # Sigma 3 (Imaginary off-diagonal)
    ]

    # 2. Pre-calculate the Hermitian conjugate (Adjoint) of J
    J_dag = J.H

    # 3. Build the 4x4 Mueller Matrix
    M = sp.Matrix.zeros(4, 4)

    for i in range(4):
        for j in range(4):
            # Formula: 1/2 * Tr( Sigma_i * J * Sigma_j * J_dagger )
            # Note: order is Sigma_i * J * Sigma_j * J_dag
            # This ensures M maps S_in (j) to S_out (i)
            stat = sigmas[i] * J * sigmas[j] * J_dag
            M[i, j] = sp.simplify(stat.trace() / 2)

    return M

def _pauli_sp():
    return [sp.Matrix([[1, 0], [0,  1]]),
         sp.Matrix([[1, 0], [0, -1]]),
         sp.Matrix([[0, 1], [1,  0]]),
         sp.Matrix([[0, -sp.I], [sp.I, 0]])]

def _mueller_rotated_crystal_document_sp(rs, rp, theta):
    """
    Returns the general rotated Mueller matrix M_chi symbolically.
    theta: the rotation angle
    rs, rp: complex reflection/transmission coefficients
    """
    # 1. Define trig shortcuts
    c2 = sp.cos(2 * theta)
    s2 = sp.sin(2 * theta)

    # 2. Define fundamental complex terms
    rs_sq = rs * sp.conjugate(rs)
    rp_sq = rp * sp.conjugate(rp)

    # Using your derived identities
    term_plus = (rs_sq + rp_sq) / 2
    term_minus = (rs_sq - rp_sq) / 2

    cross_term = rs * sp.conjugate(rp)
    re_cp = sp.re(cross_term)
    im_cp = sp.im(cross_term)

    # 3. Construct the matrix row by row
    row0 = [
        term_plus,
        term_minus * c2,
        term_minus * s2,
        0
    ]

    row1 = [
        term_minus * c2,
        term_plus * c2 ** 2 + re_cp * s2 ** 2,
        (term_plus - re_cp) * s2 * c2,
        -im_cp * s2
    ]

    row2 = [
        term_minus * s2,
        (term_plus - re_cp) * s2 * c2,
        term_plus * s2 ** 2 + re_cp * c2 ** 2,
        im_cp * c2
    ]

    row3 = [
        0,
        im_cp * s2,
        -im_cp * c2,
        re_cp
    ]

    M_chi = sp.Matrix([row0, row1, row2, row3])

    return M_chi

if __name__ == "__main__":

    chi_s = sp.Symbol("chi", real=True)
    delta_s = sp.symbols("delta", real=True)

    if True: # rotations, retarders and polarizers

        print("Rotation matrix J (eq 14)", _rotation_J_sp(chi_s))
        print("Rotation matrix M (eq 15)", _rotation_M_sp(chi_s))

        # check numeric vs symbolic values
        for i in range(100):
            chi = np.random.rand() * 2 * np.pi
            # print( _rotation_J_np(chi=chi) ,  _rotation_J_sp(chi_s).subs({chi_s: chi}) )
            # print( _rotation_M_np(chi=chi) ,  _rotation_M_sp(chi_s).subs({chi_s: chi}))
            assert ( ( _rotation_J_np(chi=chi) -  np.array(_rotation_J_sp(chi_s).subs({chi_s: chi}).evalf(15)).astype(np.float64) ).sum() < 1e-10 )
            assert ( ( _rotation_M_np(chi=chi) -  np.array(_rotation_M_sp(chi_s).subs({chi_s: chi}).evalf(15)).astype(np.float64) ).sum() < 1e-10 )



        J0 = _rotation_J_sp(chi_s) # sp.Matrix([[j00, j01], [j10, j11]])

        M0 = _jones_to_mueller_sp(J0)
        print("Rotation M from J (eq 15): ", _rotation_M_sp(chi_s))



        #
        # retarders
        #

        print("Retarder/plate J (eq 17): ", _retarder_J_sp(delta_s))
        print("Retarder/plate M (eq **): ", _jones_to_mueller_sp(_retarder_J_sp(delta_s)))

        delta = np.pi/2
        print("retarder/plate QWP J (eq 18): ", _retarder_J_sp(delta_s).subs({delta_s:delta}).evalf(5))
        print("retarder/plate QWP M (eq 19): ", _jones_to_mueller_sp(_retarder_J_sp(delta_s)).subs({delta_s: delta}).evalf(5))

        RMrot = _rotation_M_sp(-chi_s) @ _jones_to_mueller_sp(_retarder_J_sp(delta_s)) @ _rotation_M_sp(chi_s)
        print("retarder/plate rotated M (eq 20): ", RMrot)

        # eq 21 (Jones)
        delta = np.pi/2
        chi = np.pi/4
        PJrot = _rotation_J_np(-chi) @ \
              _retarder_J_sp(delta_s, use_constant=False).subs({delta_s:delta}).evalf(15) @ \
              _rotation_J_np(chi)
        print("retarder/plate QWP J (eq 21) : ", PJrot)
        print("retarder/plate QWP M from sandwiched J (eq 22) : ", _jones_to_mueller_sp(PJrot))
        print("retarder/plate QWP M from sandwich (eq 22b): ", RMrot.subs({chi_s:chi, delta_s:delta}).evalf(15))

        print("J to S [1,0]): (eq 24)", _jones_to_stokes_np([1,0]))
        print("J circular pol (eq 25-26)", PJrot @ np.array([1,0], dtype=complex))
        print("S circular pol (eq 27-30)", _jones_to_stokes_np(PJrot @ np.array([1, 0], dtype=complex)))
        print("S circular pol (eq 31)", RMrot.subs({chi_s:chi, delta_s:delta}).evalf(15) @ np.array([1, 1, 0, 0], dtype=float))

        # reverse sign
        delta = np.pi / 2
        chi = -np.pi / 4
        PJrot = _rotation_J_np(-chi) @ \
                _retarder_J_sp(delta_s, use_constant=False).subs({delta_s: delta}).evalf(15) @ \
                _rotation_J_np(chi)
        print("retarder/plate QWP J (eq 34) : ", PJrot)
        print("J circular pol (eq 34)", PJrot @ np.array([1, 0], dtype=complex))
        print("J to S: (eq 34)", _jones_to_stokes_np(PJrot @ np.array([1, 0], dtype=complex)))
        print("retarder/plate QWP M (eq 35) : ", _jones_to_mueller_sp(PJrot))


        # summary
        print("Summary retarders: ")
        print("    retarder/plate J: ", _retarder_J_sp(delta_s))
        print("    retarder/plate M: ", _jones_to_mueller_sp(_retarder_J_sp(delta_s)))

        #
        # section: 3 Experimental Polarization Analysis
        #
        print("\n")
        print("Linear Polarizer J (eq 36): ", _linearpolarizer_J_sp(delta_s))
        print("Linear Polarizer M (eq 36): ", _jones_to_mueller_sp(_linearpolarizer_J_sp(delta_s)))
        S0_s = sp.Symbol("S0", real=True)
        S1_s = sp.Symbol("S1", real=True)
        S2_s = sp.Symbol("S2", real=True)
        S3_s = sp.Symbol("S3", real=True)

        Sout = _jones_to_mueller_sp(_linearpolarizer_J_sp(delta_s)) @ sp.Matrix([S0_s, S1_s, S2_s, S3_s])
        print("Intensity (eq 38): ", Sout[0])


    if True: # crystals
        print("\n---crystals---")
        #
        # section 4 Crystal Phase Plate with General Amplitudes rs, rp
        #
        r_s_s, r_p_s   = sp.symbols("r_s r_p")

        J0 = _crystal_J_sp(r_s_s, r_p_s)
        print("J crystal plate unrotated (eq 58): ", J0)

        # check eq 59
        J0b = _pauli_sp()[0] * (r_s_s + r_p_s) / 2 + _pauli_sp()[1] * (r_s_s - r_p_s) / 2
        assert( (J0 - J0b).is_zero_matrix)

        J1 = _rotation_J_sp(-chi_s) @ _crystal_J_sp(r_s_s, r_p_s) @ _rotation_J_sp(chi_s)
        print("J crystal plate rotated (eq 60-61) : ", J1)

        # check eq 60
        J1b = (r_s_s + r_p_s) / 2 *  _pauli_sp()[0] + \
              (r_s_s - r_p_s) / 2 * (sp.cos(2 * chi_s) * _pauli_sp()[1] + sp.sin(2 * chi_s) * _pauli_sp()[2])

        print( J1 - J1b)
        if (J1 - J1b).is_zero_matrix:
            print("J crystal plate rotated (eq 61) CHECKED ANALYTICALLY")
        else: # check numeric
            for i in range(100):
                chi = np.random.rand() * 2 * np.pi
                r_s = np.random.rand() + np.random.rand() * 1j
                r_p = np.random.rand() + np.random.rand() * 1j
                values = {r_s_s: r_s, r_p_s: r_p, chi_s: chi}

                J1_subs = J1.subs(values).evalf(15)
                J1b_subs = J1b.subs(values).evalf(15)
                J1_numeric = _rotation_J_np(-chi) @ np.array([[r_s, 0], [0, r_p]], dtype=complex) @ _rotation_J_sp(chi)
                diff1 = np.abs(J1_subs - J1_numeric).sum()
                diff2 = np.abs(J1b_subs - J1_numeric).sum()
                diff3 = np.abs(J1b_subs - J1_subs).sum()
                assert (diff1 < 1e-10)
                assert (diff2 < 1e-10)
                assert (diff3 < 1e-10)
            print("J crystal plate rotated (eq 61) CHECKED NUMERICALLY!")



        #
        # final formula Mueller matrix
        #

        # unrotated
        print("\nUnrotated:")
        J0 = _crystal_J_sp(r_s_s, r_p_s)
        M0 = _jones_to_mueller_sp(J0)
        # sp.print_latex(M0)
        print("M0 from unrotated J0 (replace eq 62)", M0)

        #
        M0p = _M0_matrix_sp(r_s_s, r_p_s)
        print("M0 from unrotated J0 written in the paper (eq 62)", M0p)
        # check eq 62
        print( (M0 - M0p))
        if (M0 - M0p).is_zero_matrix:
            print("M0 from unrotated J0 written in the paper (eq 62) CHECKED ANALYTICALLY")
        else:  # check numeric
            for i in range(100):
                r_s = np.random.rand() + np.random.rand() * 1j
                r_p = np.random.rand() + np.random.rand() * 1j
                values = {r_s_s: r_s, r_p_s: r_p}

                M0_subs = M0.subs(values).evalf(15)
                M0p_subs = M0p.subs(values).evalf(15)
                diff1 = np.abs(M0_subs - M0p_subs).sum()
                assert (diff1 < 1e-10)
            print("M0 from unrotated J0 written in the paper (eq 62) CHECKED NUMERICALLY!")

        #
        # rotated
        #

        print("\nRotated:")
        M1 = _rotation_M_sp(-chi_s) @ M0 @ _rotation_M_sp(chi_s)
        print("M1 from sandwich M0 (eq 64)", M1)
        # sp.print_latex(M1)

        M1b = _jones_to_mueller_sp(J1)
        print("M1 from rotated J1 (eq 64)", M1b)

        print( M1 - M1b)
        if (M1 - M1b).is_zero_matrix:
            print("M1 from sandwich M0 (eq 64) CHECKED ANALYTICALLY")
        else: # check numeric
            for i in range(100):
                chi = np.random.rand() * 2 * np.pi
                r_s = np.random.rand() + np.random.rand() * 1j
                r_p = np.random.rand() + np.random.rand() * 1j
                values = {r_s_s: r_s, r_p_s: r_p, chi_s: chi}

                M1_subs    = M1.subs(values).evalf(15)
                M1b_subs   = M1b.subs(values).evalf(15)
                # M1_numeric = mueller_rotated_crystal_document_np(r_s, r_p, chi)
                M1_numeric = _mueller_rotated_crystal_document_sp(r_s_s, r_p_s, chi_s).subs(values).evalf(15)
                diff1 = np.abs(J1_subs - J1_numeric).sum()
                diff2 = np.abs(J1b_subs - J1_numeric).sum()
                diff3 = np.abs(J1b_subs - J1_subs).sum()
                # print(diff1, diff2, diff3)
                assert (diff1 < 1e-10)
                assert (diff2 < 1e-10)
                assert (diff3 < 1e-10)
            print("M1 from sandwich M0 (eq 64) CHECKED NUMERICALLY!")

            #
            # degree of pol, P1, P2, P3
            #

            S0_s = sp.Symbol("S0", real=True)
            S1_s = sp.Symbol("S1", real=True)
            S2_s = sp.Symbol("S2", real=True)
            S3_s = sp.Symbol("S3", real=True)
            # Create the input Stokes vector (S_in)
            S_in = sp.Matrix([S0_s, S1_s, S2_s, S3_s])

            S_out = M1 @ S_in # or
            print("\n")
            print("S0 (eq 72) = ", S_out[0])
            print("S1 (eq 73) = ", S_out[1])
            print("S2 (eq 74) = ", S_out[2])
            print("S3 (eq 75) = ", S_out[3])

            S_out2 = _stokes_Sout_sp(S_in, r_s_s, r_p_s, chi_s) # a shortcut, useful for experimental check
            print("\n")
            print("S0 ** new ** (eq 72) = ", S_out2[0])
            print("S1 ** new ** (eq 73) = ", S_out2[1])
            print("S2 ** new ** (eq 74) = ", S_out2[2])
            print("S3 ** new ** (eq 75) = ", S_out2[3])


            print("\n")
            print("S0 for linear polarizarion 1100 (eq 76) = ", S_out[0].subs({S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}))
            print("S1 for linear polarizarion 1100 (eq 77) = ", S_out[1].subs({S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}))
            print("S2 for linear polarizarion 1100 (eq 78) = ", S_out[2].subs({S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}))
            print("S3 for linear polarizarion 1100 (eq 79) = ", S_out[3].subs({S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}))

            print("\n")
            print("S0 for 45 de linear polarizarion 1100 (eq 81) = ", S_out[0].subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))
            print("S1 for 45 de linear polarizarion 1100 (eq 82) = ", S_out[1].subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))
            print("S2 for 45 de linear polarizarion 1100 (eq 83) = ", S_out[2].subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))
            print("S3 for 45 de linear polarizarion 1100 (eq 84) = ", S_out[3].subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))

            print("\n")
            print("S0 for right-circular polarizarion 1100 (eq 85) = ", S_out[0].subs({S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}))
            print("S1 for right-circular polarizarion 1100 (eq 86) = ", S_out[1].subs({S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}))
            print("S2 for right-circular polarizarion 1100 (eq 87) = ", S_out[2].subs({S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}))
            print("S3 for right-circular polarizarion 1100 (eq 88) = ", S_out[3].subs({S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}))


            Pc = S_out[3]/ S_out[0]
            print("\n")
            print("Pc (eq 89) = ", Pc)
            print("Pc for linear polarizarion 1100 (eq 90<--Giles) = ",         Pc.subs({S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0}))
            print("Pc for 45 de linear polarizarion 1100 (eq 91) = ",   Pc.subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))
            print("Giles: ")
            sp.print_latex(Pc.subs({S0_s: 1, S1_s: 0, S2_s: 1, S3_s: 0}))


            for i in range(100):
                chi = np.pi / 4 # np.random.rand() * 2 * np.pi
                r_s = np.random.rand() + np.random.rand() * 1j
                r_p = np.random.rand() + np.random.rand() * 1j
                values = {S0_s: 1, S1_s: 1, S2_s: 0, S3_s: 0, r_s_s: r_s, r_p_s: r_p, chi_s: chi}

                Pc1 = Pc.subs(values).evalf(15)
                Pc2 = - _pc_giles(r_s, r_p, chi)
                assert( np.abs(Pc1 - Pc2) < 1e-10)
            print("Giles formula (eq 90) CHECKED NUMERICALLY!  *** MINUS SIGN CHANGED !! ***")

            print("Pc for right-circular polarizarion 1100 (eq 92) = ", Pc.subs({S0_s: 1, S1_s: 0, S2_s: 0, S3_s: 1}))

    # print(" 15: ", _rotation_J_np(np.radians(15)) @ np.array([1, 1], dtype=float))
    # print("-15: ", _rotation_J_np(np.radians(-15)) @ np.array([1,1], dtype=float))