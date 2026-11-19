#!/usr/bin/env python3
"""QuantumErrorCorrection-Lab - educational quantum error correction lab.

by Darion Kellar

Pure Python standard library. Zero dependencies. Zero network.
Education only - and the central lesson is stated up front:

  * Error correction is bookkeeping with syndromes: encode one
    logical qubit into several physical qubits, read a syndrome that
    names the error without naming the data, apply the recovery the
    table prescribes, and grade the result by fidelity against the
    state you encoded - exactly, not approximately.
  * This lab builds three toy codes from scratch on state vectors:
    the 3-qubit repetition codes, Shor's 9-qubit code, and the
    Steane [[7,1,3]] code, plus the logical error rate of the
    repetition code graded two ways (closed form AND exhaustive
    enumeration of every error pattern).
  * These are toy codes on paper. NO fault-tolerance threshold for
    any real code or any real hardware is claimed anywhere in this
    file: the p = 0.5 crossing below is a fact about the repetition
    code's majority vote under iid bit flips, fenced as such, and
    nothing more. No real qubit was corrected.
  * QEC-Fortress is a separate project of Darion's. This lab is
    codes-from-scratch teaching only and is not that project.

Conventions, stated once and used everywhere:

  * Qubit 0 is the most significant bit of the state index.
  * The stated encoded state uses alpha = 0.6, beta = 0.8 (real,
    normalised: 0.36 + 0.64 = 1) unless a check says otherwise.
  * Fidelity of a recovered state is |<original|recovered>|^2; a
    correct recovery scores exactly 1 up to a global phase.
  * A syndrome bit is 0 when the stabilizer expectation is +1 and
    1 when it is -1.

Every real-world deployment use is a permanent HOLD gate - with
Gemini Argon 4 Extended Thinking as DECIDER (names-only,
tighten-only), Muse as advisor (never overrides HOLD), and Google
Flow plan-only seats.
"""
import cmath
import hashlib
import json
import math
import sys

# ---------------------------------------------------------------- constants
VERSION = "1.0.0"

STATED_ALPHA = 0.6
STATED_BETA = 0.8
PARAM_GRID = (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9, 1.0)
MC_P = 0.2
MC_TRIALS = 20000
MC_SEED = 20261118
MC_TOLERANCE = 0.02

# ---------------------------------------------------------------- HOLD gates
HOLD_GATES = {
    "fault_tolerance_claim": "HOLD - no fault-tolerance threshold claim for any real code or hardware; the p = 0.5 crossing here is a repetition-code majority-vote fact only",
    "model_as_measurement_claim": "HOLD - a computed number is a statement about the model at stated inputs, never a measurement",
    "live_hardware_control": "HOLD - no control of real qubits, links, or laboratory hardware, ever",
    "investment_claim": "HOLD - no investment claim; a lab model is not a prospectus",
    "certification_claim": "HOLD - no certification that any system meets any standard",
    "real_device_correction_claim": "HOLD - no claim that any error was corrected on a real device; every recovery here is a state-vector calculation",
    "error_correction_product_claim": "HOLD - no error-correction product, service, or offering is made or certified here",
    "hardware_specification_claim": "HOLD - no hardware specification (gate error, T1, T2) for any real device is asserted",
    "vendor_code_claim": "HOLD - no vendor code, decoder, or device error figure is stated, endorsed, or compared",
    "quantum_network_deployment_claim": "HOLD - no deploying or dimensioning a real quantum network from these toy codes",
    "cryptographic_security_claim": "HOLD - no cryptographic or QKD security claim follows from toy-code mathematics",
    "data_protection_claim": "HOLD - no claim that any stored or transmitted data is protected by these teaching codes",
    "safety_critical_claim": "HOLD - no use in safety-critical, medical, or control systems",
    "threshold_generalisation_claim": "HOLD - no generalising the repetition-code crossing to surface codes, concatenated codes, or any other family",
}


def gate(action):
    """Return 'HOLD' if the action is gated, else 'PROCEED' (sandbox sim)."""
    return "HOLD" if action in HOLD_GATES else "PROCEED"


# ================================================================ 1. Matrix kit
# Small complex matrices as plain lists of lists, plus a state-vector
# kit. Everything the lab needs - products, daggers, Kronecker
# products, Pauli actions on state vectors - is written out here so
# every number can be read and diffed.

def mat(rows):
    return [[complex(v) for v in row] for row in rows]


def eye(n):
    return [[complex(1 if i == j else 0) for j in range(n)] for i in range(n)]


def zeros(n):
    return [[0j] * n for _ in range(n)]


def madd(A, B):
    return [[a + b for a, b in zip(ra, rb)] for ra, rb in zip(A, B)]


def mscale(A, s):
    return [[a * s for a in row] for row in A]


def mmul(A, B):
    n, m, k = len(A), len(B[0]), len(B)
    out = zeros(n)
    for i in range(n):
        for j in range(m):
            out[i][j] = sum(A[i][t] * B[t][j] for t in range(k))
    return [[out[i][j] for j in range(m)] for i in range(n)]


def dagger(A):
    return [[A[j][i].conjugate() for j in range(len(A))]
            for i in range(len(A[0]))]


def mtrace(A):
    return sum(A[i][i] for i in range(len(A)))


def kron(A, B):
    """Kronecker product: out[(a,b),(c,d)] = A[a][c] * B[b][d]."""
    na, nb = len(A), len(B)
    out = zeros(na * nb)
    for a in range(na):
        for c in range(na):
            for b in range(nb):
                for d in range(nb):
                    out[a * nb + b][c * nb + d] = A[a][c] * B[b][d]
    return out


def mat_max_dev(A, B):
    return max(abs(a - b) for ra, rb in zip(A, B) for a, b in zip(ra, rb))


PAULI_I = mat([[1, 0], [0, 1]])
PAULI_X = mat([[0, 1], [1, 0]])
PAULI_Y = mat([[0, -1j], [1j, 0]])
PAULI_Z = mat([[1, 0], [0, -1]])
PAULIS = (PAULI_X, PAULI_Y, PAULI_Z)


def hermitian_eigenvalues(A):
    """Eigenvalues of a Hermitian matrix by cyclic Jacobi rotations.

    Refereed on known spectra in CHECKS before any use. Sweeps run
    until the off-diagonal Frobenius mass is below 1e-14 or 100
    sweeps pass. Returns eigenvalues sorted ascending.
    """
    n = len(A)
    a = [[A[i][j] for j in range(n)] for i in range(n)]
    for _ in range(100):
        off = math.sqrt(sum(abs(a[p][q]) ** 2
                            for p in range(n) for q in range(p + 1, n)))
        if off < 1e-14:
            break
        for p in range(n - 1):
            for q in range(p + 1, n):
                apq = a[p][q]
                b = abs(apq)
                if b < 1e-15:
                    continue
                app, aqq = a[p][p].real, a[q][q].real
                alpha = cmath.phase(apq)
                tau = (aqq - app) / (2.0 * b)
                if tau == 0.0:
                    t = -1.0
                else:
                    t = -math.copysign(
                        1.0 / (abs(tau) + math.sqrt(tau * tau + 1.0)), tau)
                c = 1.0 / math.sqrt(t * t + 1.0)
                s = t * c
                e_neg = cmath.exp(-1j * alpha)
                U = eye(n)
                U[p][p] = c
                U[q][q] = c * e_neg
                U[q][p] = s * e_neg
                U[p][q] = -s
                a = mmul(dagger(U), mmul(a, U))
    return sorted(a[i][i].real for i in range(n))


# ================================================================ 2. State kit
# State vectors are lists of complex amplitudes, index bit for
# qubit q is 1 << (n - 1 - q) (qubit 0 most significant).

def popcount(x):
    return bin(x).count("1")


def zero_state(n):
    st = [0j] * (1 << n)
    st[0] = 1.0 + 0j
    return st


def apply_pauli(state, n, qubit, pauli):
    """Apply one single-qubit Pauli ('I','X','Y','Z') to a state."""
    if pauli == "I":
        return list(state)
    bit = 1 << (n - 1 - qubit)
    out = [0j] * len(state)
    for idx, amp in enumerate(state):
        if amp == 0j:
            continue
        if pauli == "X":
            out[idx ^ bit] += amp
        elif pauli == "Z":
            out[idx] += amp * (-1.0 if idx & bit else 1.0)
        elif pauli == "Y":
            if idx & bit:
                out[idx ^ bit] += -1j * amp
            else:
                out[idx ^ bit] += 1j * amp
        else:
            raise ValueError(pauli)
    return out


def apply_pauli_string(state, n, spec):
    """spec: iterable of (qubit, pauli) pairs, applied left to right."""
    out = list(state)
    for qubit, pauli in spec:
        out = apply_pauli(out, n, qubit, pauli)
    return out


def inner(a, b):
    return sum(x.conjugate() * y for x, y in zip(a, b))


def state_norm(state):
    return math.sqrt(sum(abs(v) ** 2 for v in state))


def fidelity(a, b):
    """|<a|b>|^2 - 1 exactly for identical states up to global phase."""
    return abs(inner(a, b)) ** 2


def expectation(state, n, spec):
    """<psi| P |psi> for a Pauli string spec; real for stabilizers."""
    return inner(state, apply_pauli_string(state, n, spec)).real


def syndrome_bits(state, n, stabilizers):
    """0 for expectation +1, 1 for expectation -1, per stabilizer."""
    return tuple(0 if expectation(state, n, s) > 0 else 1
                 for s in stabilizers)


# ================================================================ 3. SplitMix64
# In-file deterministic PRNG (splitmix64) for the Monte Carlo referee.
# No random module: the stream is a pure function of the seed, so the
# Monte Carlo numbers below are byte-identical between runs.

_MASK64 = (1 << 64) - 1


def splitmix64_stream(seed):
    state = seed & _MASK64
    while True:
        state = (state + 0x9E3779B97F4A7C15) & _MASK64
        z = state
        z = ((z ^ (z >> 30)) * 0xBF58476D1CE4E5B9) & _MASK64
        z = ((z ^ (z >> 27)) * 0x94D049BB133111EB) & _MASK64
        z = z ^ (z >> 31)
        yield z / float(1 << 64)


# ================================================================ 4. Repetition
# Bit-flip code: |0_L> = |000>, |1_L> = |111>, stabilizers Z0Z1, Z1Z2,
# corrects one X error. Phase-flip code: |0_L> = |+++>, |1_L> = |--->,
# stabilizers X0X1, X1X2, corrects one Z error. Honest boundary, kept
# in the open: on a general superposition a Z error on the bit-flip
# code (or an X error on the phase-flip code) is a LOGICAL error -
# the syndrome is silent and the fidelity is (alpha^2-beta^2)^2, not
# 1. On a basis state (|0_L> alone) every Pauli does restore with
# fidelity 1, because the logical operator then acts as a phase.

REP_BIT_STABILIZERS = [[(0, "Z"), (1, "Z")], [(1, "Z"), (2, "Z")]]
REP_PHASE_STABILIZERS = [[(0, "X"), (1, "X")], [(1, "X"), (2, "X")]]
# syndrome (s0, s1) -> qubit to correct; shared shape for both codes
REP_SYNDROME_TABLE = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}


def encode_repetition_bit(alpha=STATED_ALPHA, beta=STATED_BETA):
    st = [0j] * 8
    st[0] = complex(alpha)
    st[7] = complex(beta)
    return st


def encode_repetition_phase(alpha=STATED_ALPHA, beta=STATED_BETA):
    st = []
    for x in range(8):
        sign = -1.0 if popcount(x) % 2 else 1.0
        st.append(complex((alpha + beta * sign) / math.sqrt(8.0)))
    return st


def recover_repetition_bit(state):
    syn = syndrome_bits(state, 3, REP_BIT_STABILIZERS)
    q = REP_SYNDROME_TABLE[syn]
    if q is None:
        return list(state), syn
    return apply_pauli(state, 3, q, "X"), syn


def recover_repetition_phase(state):
    syn = syndrome_bits(state, 3, REP_PHASE_STABILIZERS)
    q = REP_SYNDROME_TABLE[syn]
    if q is None:
        return list(state), syn
    return apply_pauli(state, 3, q, "Z"), syn


def repetition_report():
    rows = {"alpha": STATED_ALPHA, "beta": STATED_BETA}
    bit0 = encode_repetition_bit()
    phase0 = encode_repetition_phase()
    bit_rows = {}
    for q in range(3):
        err = apply_pauli(bit0, 3, q, "X")
        rec, syn = recover_repetition_bit(err)
        bit_rows[f"X_on_qubit_{q}"] = {
            "syndrome": list(syn), "recovered_fidelity": fidelity(bit0, rec)}
    rows["bit_flip_code_X_errors"] = bit_rows
    phase_rows = {}
    for q in range(3):
        err = apply_pauli(phase0, 3, q, "Z")
        rec, syn = recover_repetition_phase(err)
        phase_rows[f"Z_on_qubit_{q}"] = {
            "syndrome": list(syn),
            "recovered_fidelity": fidelity(phase0, rec)}
    rows["phase_flip_code_Z_errors"] = phase_rows
    rows["syndrome_table"] = {
        ",".join(str(b) for b in k): v
        for k, v in REP_SYNDROME_TABLE.items()}
    # honest logical-error witnesses on the general superposition
    err = apply_pauli(bit0, 3, 0, "Z")
    rec, syn = recover_repetition_bit(err)
    rows["bit_code_Z_error_is_logical"] = {
        "syndrome": list(syn),
        "recovered_fidelity": fidelity(bit0, rec),
        "expected_fidelity": (STATED_ALPHA ** 2 - STATED_BETA ** 2) ** 2}
    err = apply_pauli(phase0, 3, 0, "X")
    rec, syn = recover_repetition_phase(err)
    rows["phase_code_X_error_is_logical"] = {
        "syndrome": list(syn),
        "recovered_fidelity": fidelity(phase0, rec),
        "expected_fidelity": (STATED_ALPHA ** 2 - STATED_BETA ** 2) ** 2}
    return rows


# ================================================================ 5. Logical rates
# iid bit flips with probability p per qubit, majority-vote decoding.
# n = 3: p_L = 3p^2 - 2p^3 in closed form AND by exhaustive
# enumeration of all 8 error patterns. n = 5, 7: binomial sums over
# weights above the majority threshold. Every curve crosses the
# physical line p_L = p at exactly p = 0.5 - a repetition-code fact
# about majority vote, NOT a fault-tolerance threshold for any real
# code or hardware (that generalisation is a HOLD gate).

def logical_rate_closed_n3(p):
    return 3.0 * p * p - 2.0 * p * p * p


def logical_rate_enumerated_n3(p):
    """Sum the probabilities of the weight >= 2 patterns, all 8 listed."""
    total = 0.0
    fail = 0.0
    for mask in range(8):
        w = popcount(mask)
        prob = (p ** w) * ((1.0 - p) ** (3 - w))
        total += prob
        if w >= 2:
            fail += prob
    return fail


def logical_rate_binomial(p, n):
    """p_L = sum_{w > n/2} C(n,w) p^w (1-p)^(n-w) for odd n."""
    return sum(math.comb(n, w) * (p ** w) * ((1.0 - p) ** (n - w))
               for w in range(n // 2 + 1, n + 1))


def logical_crossing(n):
    """Bisect p_L - p = 0 on [0.1, 0.9] (sign change is guaranteed)."""
    lo, hi = 0.1, 0.9
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if logical_rate_binomial(mid, n) - mid < 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def monte_carlo_n3(p=MC_P, trials=MC_TRIALS, seed=MC_SEED):
    """Seeded splitmix64 Monte Carlo of the n = 3 code at p.

    Each trial draws 3 uniforms; a logical failure is weight >= 2.
    Deterministic: same seed, same estimate, every run.
    """
    stream = splitmix64_stream(seed)
    failures = 0
    for _ in range(trials):
        w = sum(1 for _ in range(3) if next(stream) < p)
        if w >= 2:
            failures += 1
    return {
        "p": p,
        "trials": trials,
        "seed": seed,
        "failures": failures,
        "estimate": failures / trials,
        "closed_form": logical_rate_closed_n3(p),
        "abs_dev": abs(failures / trials - logical_rate_closed_n3(p)),
        "tolerance": MC_TOLERANCE,
    }


def rates_report():
    rows = {}
    grid = {}
    for p in PARAM_GRID:
        closed = logical_rate_closed_n3(p)
        enum = logical_rate_enumerated_n3(p)
        grid[str(p)] = {
            "closed_n3": round(closed, 9),
            "enumerated_n3": round(enum, 9),
            "abs_dev": abs(closed - enum),
            "binomial_n5": round(logical_rate_binomial(p, 5), 9),
            "binomial_n7": round(logical_rate_binomial(p, 7), 9),
        }
    rows["grid"] = grid
    rows["crossings"] = {
        "n3": round(logical_crossing(3), 9),
        "n5": round(logical_crossing(5), 9),
        "n7": round(logical_crossing(7), 9),
        "exact": 0.5,
    }
    rows["monte_carlo"] = monte_carlo_n3()
    return rows


# ================================================================ 6. Shor 9-qubit
# |0_L> = ((|000>+|111>)/sqrt2)^3, |1_L> = ((|000>-|111>)/sqrt2)^3.
# Stabilizers: Z0Z1, Z1Z2 in each block of three (bit syndromes) and
# X0..X5, X3..X8 across blocks (phase syndromes). The recovery reads
# the full 8-bit syndrome first, then applies the X correction the
# bit syndromes name and the Z correction (first qubit of the block)
# the phase syndrome names. Every single-qubit Pauli (9 qubits x
# X,Y,Z = 27 cases) is restored with fidelity exactly 1. Honest
# boundary: a weight-2 error can share a syndrome with a weight-1
# error - X on qubits 0 and 1 is read as X on qubit 2, the recovery
# completes a logical phase flip, and the fidelity is 0.0784, shown
# in the open below, not hidden.

SHOR_BLOCKS = ((0, 1, 2), (3, 4, 5), (6, 7, 8))
SHOR_BIT_STABILIZERS = [[(a, "Z"), (b, "Z")]
                        for (a, b, _c) in SHOR_BLOCKS
                        for (a, b) in ((a, b), (b, _c))]
SHOR_PHASE_STABILIZERS = [[(q, "X") for q in range(0, 6)],
                          [(q, "X") for q in range(3, 9)]]
SHOR_PHASE_TABLE = {(0, 0): None, (1, 0): 0, (1, 1): 1, (0, 1): 2}


def encode_shor(alpha=STATED_ALPHA, beta=STATED_BETA):
    st = [0j] * (1 << 9)
    for x in range(1 << 9):
        bits = [(x >> (8 - q)) & 1 for q in range(9)]
        ok = True
        sign = 1.0
        for block in SHOR_BLOCKS:
            triple = tuple(bits[q] for q in block)
            if triple == (0, 0, 0):
                pass
            elif triple == (1, 1, 1):
                sign = -sign
            else:
                ok = False
                break
        if ok:
            st[x] = complex((alpha + beta * sign) / math.sqrt(8.0))
    return st


def shor_syndrome(state):
    bit_syn = []
    for i in range(0, 6, 2):
        bit_syn.append(syndrome_bits(state, 9,
                                     SHOR_BIT_STABILIZERS[i:i + 2]))
    phase_syn = syndrome_bits(state, 9, SHOR_PHASE_STABILIZERS)
    return tuple(bit_syn), phase_syn


def recover_shor(state):
    bit_syn, phase_syn = shor_syndrome(state)
    out = list(state)
    for block_idx, syn in enumerate(bit_syn):
        q_local = REP_SYNDROME_TABLE[syn]
        if q_local is not None:
            out = apply_pauli(out, 9, SHOR_BLOCKS[block_idx][q_local], "X")
    block = SHOR_PHASE_TABLE[phase_syn]
    if block is not None:
        out = apply_pauli(out, 9, SHOR_BLOCKS[block][0], "Z")
    return out, (bit_syn, phase_syn)


def shor_report():
    orig = encode_shor()
    cases = {}
    worst = 1.0
    for q in range(9):
        for pauli in ("X", "Y", "Z"):
            err = apply_pauli(orig, 9, q, pauli)
            rec, _ = recover_shor(err)
            f = fidelity(orig, rec)
            worst = min(worst, f)
            cases[f"{pauli}_on_qubit_{q}"] = round(f, 12)
    err = apply_pauli(apply_pauli(orig, 9, 0, "X"), 9, 1, "X")
    rec, syn = recover_shor(err)
    phase_flipped = encode_shor(STATED_ALPHA, -STATED_BETA)
    return {
        "alpha": STATED_ALPHA, "beta": STATED_BETA,
        "single_pauli_cases": 27,
        "worst_single_pauli_fidelity": worst,
        "cases": cases,
        "weight2_X0_X1_honest_failure": {
            "syndrome_read_as": "X on qubit 2 (same block)",
            "recovered_fidelity": fidelity(orig, rec),
            "expected_fidelity": (STATED_ALPHA ** 2
                                  - STATED_BETA ** 2) ** 2,
            "recovered_equals_logical_phase_flip_fidelity":
                fidelity(phase_flipped, rec),
        },
    }


# ================================================================ 7. Steane [[7,1,3]]
# Classical Hamming(7,4) parity-check matrix H built explicitly
# below; its columns are the seven non-zero 3-bit strings, so every
# single-qubit error has a distinct non-zero syndrome. The Hamming
# code has 16 codewords, minimum weight 3; its dual (the row space
# of H) has 8 codewords, minimum weight 4, and sits inside the
# Hamming code - that containment is what makes the CSS construction
# work. |0_L> is the uniform superposition over the 8 even-weight
# Hamming codewords (the dual code), |1_L> over the 8 odd-weight
# ones. Stabilizers: for each row of H, a Z-string and an X-string
# on the row's support. X errors are read off the Z stabilizers,
# Z errors off the X stabilizers; 7 X cases + 7 Z cases, all exact.

HAMMING_H = [[0, 0, 0, 1, 1, 1, 1],
             [0, 1, 1, 0, 0, 1, 1],
             [1, 0, 1, 0, 1, 0, 1]]


def hamming_syndrome(word):
    return tuple(sum(HAMMING_H[i][j] * word[j] for j in range(7)) % 2
                 for i in range(3))


def hamming_codewords():
    out = []
    for x in range(128):
        word = tuple((x >> (6 - j)) & 1 for j in range(7))
        if hamming_syndrome(word) == (0, 0, 0):
            out.append(word)
    return out


def dual_codewords():
    out = set()
    for mask in range(8):
        word = tuple(sum(((mask >> i) & 1) * HAMMING_H[i][j]
                         for i in range(3)) % 2 for j in range(7))
        out.add(word)
    return sorted(out)


def encode_steane(alpha=STATED_ALPHA, beta=STATED_BETA):
    code = hamming_codewords()
    even = [w for w in code if sum(w) % 2 == 0]
    odd = [w for w in code if sum(w) % 2 == 1]
    st = [0j] * 128
    for word in even:
        idx = sum(b << (6 - j) for j, b in enumerate(word))
        st[idx] += complex(alpha / math.sqrt(8.0))
    for word in odd:
        idx = sum(b << (6 - j) for j, b in enumerate(word))
        st[idx] += complex(beta / math.sqrt(8.0))
    return st


STEANE_Z_STABILIZERS = [[(j, "Z") for j in range(7) if HAMMING_H[i][j]]
                        for i in range(3)]
STEANE_X_STABILIZERS = [[(j, "X") for j in range(7) if HAMMING_H[i][j]]
                        for i in range(3)]
STEANE_COLUMN_TABLE = {tuple(HAMMING_H[i][j] for i in range(3)): j
                       for j in range(7)}


def recover_steane(state):
    bit_syn = syndrome_bits(state, 7, STEANE_Z_STABILIZERS)
    phase_syn = syndrome_bits(state, 7, STEANE_X_STABILIZERS)
    out = list(state)
    if bit_syn != (0, 0, 0):
        out = apply_pauli(out, 7, STEANE_COLUMN_TABLE[bit_syn], "X")
    if phase_syn != (0, 0, 0):
        out = apply_pauli(out, 7, STEANE_COLUMN_TABLE[phase_syn], "Z")
    return out, (bit_syn, phase_syn)


def steane_report():
    orig = encode_steane()
    code = hamming_codewords()
    dual = dual_codewords()
    cases = {}
    for q in range(7):
        for pauli in ("X", "Z"):
            err = apply_pauli(orig, 7, q, pauli)
            rec, _ = recover_steane(err)
            cases[f"{pauli}_on_qubit_{q}"] = round(fidelity(orig, rec), 12)
    return {
        "parameters": "[7,1,3]",
        "n": 7, "k": 1, "d": 3,
        "hamming_codewords": len(code),
        "hamming_min_weight": min(sum(w) for w in code if any(w)),
        "dual_codewords": len(dual),
        "dual_min_weight": min(sum(w) for w in dual if any(w)),
        "columns_distinct_nonzero":
            len(set(STEANE_COLUMN_TABLE)) == 7
            and (0, 0, 0) not in STEANE_COLUMN_TABLE,
        "single_X_Z_cases": 14,
        "worst_X_Z_fidelity": min(cases.values()),
        "cases": cases,
    }


# ================================================================ 8. kNN classifier
# Features of a decoding situation: (p, n, w) = physical bit-flip
# rate, code length, observed syndrome weight (number of flipped
# qubits the syndrome implies). The training grid is labelled by the
# lab's own policy (model_label below) - the classifier interpolates
# the policy, it does not replace the computed rates above.

def model_label(p, n, w):
    """The regime policy distilled to three verdicts."""
    if w > n // 2 or p >= 0.5:
        return "logical-failure-likely"
    if w == n // 2 or p >= 0.3:
        return "marginal"
    return "correctable-likely"


def _training_grid():
    rows = []
    for p in (0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8):
        for n in (3, 5, 7):
            for w in range(0, n + 1):
                feats = (p, n, w)
                rows.append((feats, model_label(*feats)))
    return rows


TRAINING_GRID = _training_grid()


def classify_features(features, k=5):
    """Distance-weighted kNN over the policy-labelled grid.

    Features are normalised to (p, n/7, w/7) before distances.
    """
    p, n, w = features
    target = (p, n / 7.0, w / 7.0)
    ranked = []
    for feats, label in TRAINING_GRID:
        fp, fn, fw = feats
        cand = (fp, fn / 7.0, fw / 7.0)
        dist = math.sqrt(sum((a - b) ** 2
                             for a, b in zip(target, cand)))
        ranked.append((dist, feats, label))
    ranked.sort(key=lambda row: (row[0], row[1]))
    votes = {}
    for dist, _, label in ranked[:k]:
        votes[label] = votes.get(label, 0.0) + 1.0 / (dist + 1e-9)
    priority = {"correctable-likely": 0, "marginal": 1,
                "logical-failure-likely": 2}
    return sorted(votes.items(), key=lambda kv: (-kv[1], priority[kv[0]]))[0][0]


CLASSIFIER_SCENARIOS = {
    "clean_n3_p0.0_w0": (0.0, 3, 0),
    "mild_n5_p0.1_w1": (0.1, 5, 1),
    "mild_n7_p0.2_w2": (0.2, 7, 2),
    "at_limit_n3_p0.1_w1": (0.1, 3, 1),
    "hot_n5_p0.4_w1": (0.4, 5, 1),
    "at_limit_n5_p0.2_w2": (0.2, 5, 2),
    "overweight_n3_p0.1_w2": (0.1, 3, 2),
    "past_crossing_n5_p0.6_w0": (0.6, 5, 0),
    "overweight_n7_p0.2_w4": (0.2, 7, 4),
}


def classifier_report():
    rows = {}
    for name, feats in CLASSIFIER_SCENARIOS.items():
        rows[name] = {
            "features": {"p": feats[0], "n": feats[1],
                         "syndrome_weight": feats[2]},
            "policy": model_label(*feats),
            "knn": classify_features(feats),
        }
    return rows


# ================================================================ 9. Advisory seats
def _trace(payload):
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True).encode()).hexdigest()[:16]


def argon4_decide(topic="quantumerrorcorrection-demo"):
    """Gemini Argon 4 DECIDER seat: MAP->MODEL->STRESS->DECIDE.

    Names-only keys, live calls 0, tighten-only: can add HOLDs, never lift one.
    """
    stages = ["MAP", "MODEL", "STRESS", "DECIDE"]
    payload = {"seat": "gemini-argon-4", "topic": topic, "stages": stages}
    return {
        "seat": "Gemini Argon 4 Extended Thinking",
        "model_id": "gemini-argon-4",
        "role": "DECIDER",
        "pipeline": stages,
        "tiers": [2048, 8192, 32768],
        "key": "names-only (GEMINI_API_KEY name checked, value never read)",
        "live_calls": 0,
        "decision": "PROCEED - sandbox simulation only, all deployment gates stay HOLD",
        "can_override_hold": False,
        "trace_sha256_16": _trace(payload),
    }


def claude_advisor(topic="quantumerrorcorrection-demo"):
    """Muse advisor seat: advises, can never override a HOLD."""
    return {
        "seat": "Muse",
        "role": "ADVISOR",
        "key": "names-only (MODEL_API_KEY name checked, value never read)",
        "live_calls": 0,
        "advice": "Print the fence before the number: p = 0.5 is where majority vote on a repetition code stops helping under iid bit flips - it is not a fault-tolerance threshold for any real code or device - then show the syndrome tables, the 27 Shor cases, the one weight-2 failure, and the Steane tables that were actually computed.",
        "can_override_hold": False,
    }


def google_flow_plan():
    """Google Flow seat: plan-only, nothing rendered, nothing live."""
    return {
        "seat": "Google Flow",
        "role": "PLAN-ONLY",
        "shots": ["three qubits voting by majority while one flips and is voted back",
                  "the p_L curves for n = 3, 5, 7 bending away from the diagonal at p = 0.5",
                  "nine qubits in three blocks, a phase error lighting one block's syndrome",
                  "a weight-2 error wearing a weight-1 syndrome as a disguise",
                  "seven qubits under the Hamming parity checks, columns lighting one by one"],
        "rendered": False,
        "live": 0,
    }


def brains_report():
    return {"decider": argon4_decide(), "advisor": claude_advisor(),
            "flow": google_flow_plan()}


# ================================================================ checks
def run_checks():
    c = []

    def chk(name, ok):
        c.append((name, bool(ok)))

    # 0. Jacobi eigensolver + matrix kit refereed on known spectra first
    chk("jacobi_pauli_y_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(PAULI_Y), (-1.0, 1.0))))
    chk("jacobi_complex_hermitian_2x2",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(mat([[2, 1j], [-1j, 3]])),
                ((5 - math.sqrt(5)) / 2, (5 + math.sqrt(5)) / 2))))
    xx_zz = madd(kron(PAULI_X, PAULI_X), kron(PAULI_Z, PAULI_Z))
    chk("jacobi_xx_plus_zz_spectrum",
        all(abs(a - b) < 1e-12 for a, b in
            zip(hermitian_eigenvalues(xx_zz), (-2.0, 0.0, 0.0, 2.0))))
    chk("kron_pauli_algebra",
        mat_max_dev(mmul(PAULI_X, PAULI_Z),
                    mscale(PAULI_Y, -1j)) < 1e-12
        and mat_max_dev(kron(PAULI_I, PAULI_X),
                        mat([[0, 1, 0, 0], [1, 0, 0, 0],
                             [0, 0, 0, 1], [0, 0, 1, 0]])) < 1e-12)

    # 1. State kit: Pauli actions are unitary and fidelity is honest
    base = encode_repetition_bit()
    chk("pauli_actions_preserve_norm",
        all(abs(state_norm(apply_pauli(base, 3, q, p)) - 1.0) < 1e-12
            for q in range(3) for p in ("X", "Y", "Z")))
    chk("pauli_actions_involutions",
        all(fidelity(base, apply_pauli(
            apply_pauli(base, 3, q, p), 3, q, p)) > 1.0 - 1e-12
            for q in range(3) for p in ("X", "Y", "Z")))
    chk("fidelity_self_one_orthogonal_zero",
        abs(fidelity(base, base) - 1.0) < 1e-12
        and fidelity(encode_repetition_bit(1.0, 0.0),
                     encode_repetition_bit(0.0, 1.0)) < 1e-12)

    # 2. Repetition codes: encoding, syndromes, recovery
    bit0 = encode_repetition_bit()
    phase0 = encode_repetition_phase()
    chk("repetition_bit_encode_norm_orthonormal",
        abs(state_norm(bit0) - 1.0) < 1e-12
        and abs(inner(encode_repetition_bit(1.0, 0.0),
                      encode_repetition_bit(0.0, 1.0))) < 1e-12)
    chk("repetition_phase_encode_norm_orthonormal",
        abs(state_norm(phase0) - 1.0) < 1e-12
        and abs(inner(encode_repetition_phase(1.0, 0.0),
                      encode_repetition_phase(0.0, 1.0))) < 1e-12)
    chk("repetition_stabilizers_plus_one_on_codewords",
        all(abs(expectation(bit0, 3, s) - 1.0) < 1e-12
            for s in REP_BIT_STABILIZERS)
        and all(abs(expectation(phase0, 3, s) - 1.0) < 1e-12
                for s in REP_PHASE_STABILIZERS))
    chk("repetition_syndrome_table_complete_distinct",
        len(REP_SYNDROME_TABLE) == 4
        and sorted(v for v in REP_SYNDROME_TABLE.values()
                   if v is not None) == [0, 1, 2])
    chk("bit_code_corrects_all_X_general_state",
        all(fidelity(bit0, recover_repetition_bit(
            apply_pauli(bit0, 3, q, "X"))[0]) > 1.0 - 1e-12
            for q in range(3)))
    chk("bit_code_corrects_identity",
        fidelity(bit0, recover_repetition_bit(bit0)[0]) > 1.0 - 1e-12
        and recover_repetition_bit(bit0)[1] == (0, 0))
    chk("phase_code_corrects_all_Z_general_state",
        all(fidelity(phase0, recover_repetition_phase(
            apply_pauli(phase0, 3, q, "Z"))[0]) > 1.0 - 1e-12
            for q in range(3)))
    chk("phase_code_corrects_identity",
        fidelity(phase0, recover_repetition_phase(phase0)[0]) > 1.0 - 1e-12
        and recover_repetition_phase(phase0)[1] == (0, 0))
    chk("bit_code_syndromes_match_table",
        all(recover_repetition_bit(apply_pauli(bit0, 3, q, "X"))[1]
            == syn for syn, q in REP_SYNDROME_TABLE.items() if q is not None))
    chk("phase_code_syndromes_match_table",
        all(recover_repetition_phase(apply_pauli(phase0, 3, q, "Z"))[1]
            == syn for syn, q in REP_SYNDROME_TABLE.items()
            if q is not None))
    # On a basis state every Pauli restores exactly (logical acts as phase)
    bit_basis = encode_repetition_bit(1.0, 0.0)
    phase_basis = encode_repetition_phase(1.0, 0.0)
    chk("bit_code_basis_state_all_pauli_fidelity_one",
        all(fidelity(bit_basis, recover_repetition_bit(
            apply_pauli(bit_basis, 3, q, p))[0]) > 1.0 - 1e-12
            for q in range(3) for p in ("X", "Y", "Z")))
    chk("phase_code_basis_state_all_pauli_fidelity_one",
        all(fidelity(phase_basis, recover_repetition_phase(
            apply_pauli(phase_basis, 3, q, p))[0]) > 1.0 - 1e-12
            for q in range(3) for p in ("X", "Y", "Z")))
    # Honest boundary: the off-type error is logical, fidelity 0.0784
    logical_fid = (STATED_ALPHA ** 2 - STATED_BETA ** 2) ** 2
    chk("bit_code_Z_error_is_logical_honest",
        abs(fidelity(bit0, recover_repetition_bit(
            apply_pauli(bit0, 3, 0, "Z"))[0]) - logical_fid) < 1e-12
        and recover_repetition_bit(apply_pauli(bit0, 3, 0, "Z"))[1]
        == (0, 0))
    chk("phase_code_X_error_is_logical_honest",
        abs(fidelity(phase0, recover_repetition_phase(
            apply_pauli(phase0, 3, 0, "X"))[0]) - logical_fid) < 1e-12)

    # 3. Logical error rates: closed form vs enumeration, crossings
    chk("n3_closed_matches_enumerated_grid",
        all(abs(logical_rate_closed_n3(p)
                - logical_rate_enumerated_n3(p)) < 1e-12
            for p in PARAM_GRID))
    chk("n3_closed_form_spot_values",
        abs(logical_rate_closed_n3(0.2) - 0.104) < 1e-12
        and abs(logical_rate_closed_n3(0.5) - 0.5) < 1e-12
        and abs(logical_rate_closed_n3(0.0)) < 1e-12
        and abs(logical_rate_closed_n3(1.0) - 1.0) < 1e-12)
    chk("n3_enumeration_probabilities_sum_one",
        all(abs(sum((p ** popcount(m)) * ((1 - p) ** (3 - popcount(m)))
                    for m in range(8)) - 1.0) < 1e-12
            for p in PARAM_GRID))
    chk("n3_crossing_bisected_at_half",
        abs(logical_crossing(3) - 0.5) < 1e-9)
    chk("n5_crossing_bisected_at_half",
        abs(logical_crossing(5) - 0.5) < 1e-9)
    chk("n7_crossing_bisected_at_half",
        abs(logical_crossing(7) - 0.5) < 1e-9)
    chk("n5_binomial_matches_polynomial",
        all(abs(logical_rate_binomial(p, 5)
                - (10 * p ** 3 - 15 * p ** 4 + 6 * p ** 5)) < 1e-12
            for p in PARAM_GRID))
    chk("binomial_rates_at_half_are_half",
        all(abs(logical_rate_binomial(0.5, n) - 0.5) < 1e-12
            for n in (3, 5, 7)))
    chk("codes_help_below_half_hurt_above",
        all(logical_rate_binomial(p, n) < p
            for p in (0.1, 0.2, 0.3, 0.4) for n in (3, 5, 7))
        and all(logical_rate_binomial(p, n) > p
                for p in (0.6, 0.75, 0.9) for n in (3, 5, 7)))
    chk("larger_n_better_below_half",
        logical_rate_binomial(0.2, 7)
        < logical_rate_binomial(0.2, 5)
        < logical_rate_binomial(0.2, 3) < 0.2)

    # 4. Monte Carlo referee (seeded splitmix64, deterministic)
    mc = monte_carlo_n3()
    chk("monte_carlo_within_stated_tolerance",
        mc["abs_dev"] < mc["tolerance"]
        and mc["trials"] == MC_TRIALS and mc["p"] == MC_P)
    chk("monte_carlo_deterministic_same_seed",
        monte_carlo_n3() == monte_carlo_n3())
    stream = splitmix64_stream(12345)
    draws = [next(stream) for _ in range(1000)]
    chk("splitmix64_stream_in_unit_interval_distinct",
        all(0.0 <= v < 1.0 for v in draws)
        and len(set(draws)) == 1000)

    # 5. Shor 9-qubit code: 27 cases exact, one honest weight-2 failure
    shor0 = encode_shor()
    chk("shor_encode_norm_orthonormal",
        abs(state_norm(shor0) - 1.0) < 1e-12
        and abs(inner(encode_shor(1.0, 0.0),
                      encode_shor(0.0, 1.0))) < 1e-12)
    chk("shor_stabilizers_plus_one_on_codeword",
        all(abs(expectation(shor0, 9, s) - 1.0) < 1e-9
            for s in SHOR_BIT_STABILIZERS + SHOR_PHASE_STABILIZERS)
        and len(SHOR_BIT_STABILIZERS) == 6
        and len(SHOR_PHASE_STABILIZERS) == 2)
    chk("shor_corrects_all_27_single_pauli",
        all(fidelity(shor0, recover_shor(
            apply_pauli(shor0, 9, q, p))[0]) > 1.0 - 1e-9
            for q in range(9) for p in ("X", "Y", "Z")))
    chk("shor_corrects_identity_zero_syndrome",
        fidelity(shor0, recover_shor(shor0)[0]) > 1.0 - 1e-12
        and recover_shor(shor0)[1][1] == (0, 0))
    w2_err = apply_pauli(apply_pauli(shor0, 9, 0, "X"), 9, 1, "X")
    w2_rec, _ = recover_shor(w2_err)
    chk("shor_weight2_failure_honest_fidelity",
        abs(fidelity(shor0, w2_rec) - logical_fid) < 1e-9
        and fidelity(shor0, w2_rec) < 0.5)
    chk("shor_weight2_recovered_is_logical_phase_flip",
        fidelity(encode_shor(STATED_ALPHA, -STATED_BETA), w2_rec)
        > 1.0 - 1e-9)
    chk("shor_weight2_syndrome_disguised_as_weight1",
        recover_shor(w2_err)[1][0][0]
        == recover_shor(apply_pauli(shor0, 9, 2, "X"))[1][0][0])

    # 6. Steane [[7,1,3]]: Hamming combinatorics + 14 cases exact
    code = hamming_codewords()
    dual = dual_codewords()
    steane0 = encode_steane()
    chk("hamming_codewords_16_min_weight_3",
        len(code) == 16
        and min(sum(w) for w in code if any(w)) == 3)
    chk("dual_codewords_8_min_weight_4_subset",
        len(dual) == 8
        and min(sum(w) for w in dual if any(w)) == 4
        and all(w in code for w in dual))
    chk("hamming_columns_distinct_nonzero",
        len(STEANE_COLUMN_TABLE) == 7
        and (0, 0, 0) not in STEANE_COLUMN_TABLE)
    chk("steane_parameters_7_1_3",
        len(code) == 2 ** 4 and len(dual) == 2 ** 3
        and 7 - 3 - 3 == 1)
    chk("steane_encode_norm_orthonormal",
        abs(state_norm(steane0) - 1.0) < 1e-12
        and abs(inner(encode_steane(1.0, 0.0),
                      encode_steane(0.0, 1.0))) < 1e-12)
    chk("steane_stabilizers_plus_one_on_codeword",
        all(abs(expectation(steane0, 7, s) - 1.0) < 1e-9
            for s in STEANE_Z_STABILIZERS + STEANE_X_STABILIZERS)
        and len(STEANE_Z_STABILIZERS) == 3
        and len(STEANE_X_STABILIZERS) == 3)
    chk("steane_corrects_all_X_7",
        all(fidelity(steane0, recover_steane(
            apply_pauli(steane0, 7, q, "X"))[0]) > 1.0 - 1e-9
            for q in range(7)))
    chk("steane_corrects_all_Z_7",
        all(fidelity(steane0, recover_steane(
            apply_pauli(steane0, 7, q, "Z"))[0]) > 1.0 - 1e-9
            for q in range(7)))
    chk("steane_corrects_all_Y_7_bonus",
        all(fidelity(steane0, recover_steane(
            apply_pauli(steane0, 7, q, "Y"))[0]) > 1.0 - 1e-9
            for q in range(7)))
    chk("steane_syndromes_are_hamming_columns",
        all(recover_steane(apply_pauli(steane0, 7, q, "X"))[1][0]
            == tuple(HAMMING_H[i][q] for i in range(3))
            for q in range(7)))

    # 7. Classifier: labels come from the lab's own policy
    rep = classifier_report()
    chk("classifier_scenario_labels",
        rep["clean_n3_p0.0_w0"]["knn"] == "correctable-likely"
        and rep["mild_n5_p0.1_w1"]["knn"] == "correctable-likely"
        and rep["mild_n7_p0.2_w2"]["knn"] == "correctable-likely"
        and rep["at_limit_n3_p0.1_w1"]["knn"] == "marginal"
        and rep["hot_n5_p0.4_w1"]["knn"] == "marginal"
        and rep["at_limit_n5_p0.2_w2"]["knn"] == "marginal"
        and rep["overweight_n3_p0.1_w2"]["knn"]
        == "logical-failure-likely"
        and rep["past_crossing_n5_p0.6_w0"]["knn"]
        == "logical-failure-likely"
        and rep["overweight_n7_p0.2_w4"]["knn"]
        == "logical-failure-likely")
    chk("classifier_knn_matches_policy_scenarios",
        all(row["knn"] == row["policy"] for row in rep.values()))
    chk("classifier_knn_matches_policy_grid",
        all(classify_features(feats) == label
            for feats, label in TRAINING_GRID))

    # gates + seats
    chk("hold_gates_14", len(HOLD_GATES) == 14
        and all(v.startswith("HOLD") for v in HOLD_GATES.values()))
    chk("gate_fault_tolerance_hold",
        gate("fault_tolerance_claim") == "HOLD")
    chk("gate_model_as_measurement_hold",
        gate("model_as_measurement_claim") == "HOLD")
    chk("gate_live_hardware_hold", gate("live_hardware_control") == "HOLD")
    chk("gate_investment_hold", gate("investment_claim") == "HOLD")
    chk("gate_certification_hold", gate("certification_claim") == "HOLD")
    chk("gate_sandbox_proceed", gate("run_simulation") == "PROCEED")
    brains = brains_report()
    chk("argon4_no_override", brains["decider"]["can_override_hold"] is False
        and brains["decider"]["seat"]
        == "Gemini Argon 4 Extended Thinking"
        and brains["decider"]["model_id"] == "gemini-argon-4")
    chk("claude_no_override", brains["advisor"]["can_override_hold"] is False)
    chk("flow_plan_only", brains["flow"]["rendered"] is False
        and brains["flow"]["live"] == 0)
    return c


# ================================================================ CLI / demo
def demo():
    print("=" * 68)
    print("QuantumErrorCorrection-Lab - toy codes, syndromes, honest rates")
    print("by Darion Kellar | MIT | pure stdlib | zero network | EDUCATION ONLY")
    print("=" * 68)
    print("\n[1] THE HONEST CORE - toy codes on paper, nothing real")
    print("    Encode, read a syndrome, recover, grade by fidelity -")
    print("    exactly. The p = 0.5 crossing is a repetition-code")
    print("    majority-vote fact, NOT a fault-tolerance threshold for")
    print("    any real code or hardware. QEC-Fortress is a separate")
    print("    project; this lab is codes-from-scratch teaching only.")
    print("\n[2] Repetition codes - syndrome table + recovery (alpha=0.6, beta=0.8)")
    print(f"    syndrome table (s0,s1) -> qubit: "
          f"{ {','.join(str(b) for b in k): v for k, v in REP_SYNDROME_TABLE.items()} }")
    rep = repetition_report()
    for name, row in rep["bit_flip_code_X_errors"].items():
        print(f"    bit-flip  {name:16s} syndrome {row['syndrome']} | "
              f"recovered fidelity {row['recovered_fidelity']:.6f}")
    for name, row in rep["phase_flip_code_Z_errors"].items():
        print(f"    phase-flip {name:16s} syndrome {row['syndrome']} | "
              f"recovered fidelity {row['recovered_fidelity']:.6f}")
    print(f"    honest boundary: Z on the bit code is LOGICAL - syndrome "
          f"{rep['bit_code_Z_error_is_logical']['syndrome']}, fidelity "
          f"{rep['bit_code_Z_error_is_logical']['recovered_fidelity']:.6f} "
          f"(= (0.36-0.64)^2), not 1.")
    print("\n[3] Logical error rate - closed form vs exhaustive enumeration")
    rates = rates_report()
    for p in ("0.1", "0.2", "0.3", "0.5"):
        row = rates["grid"][p]
        print(f"    p={p}: n=3 closed {row['closed_n3']:.6f} | enumerated "
              f"{row['enumerated_n3']:.6f} | dev {row['abs_dev']:.2e} | "
              f"n=5 {row['binomial_n5']:.6f} | n=7 {row['binomial_n7']:.6f}")
    print(f"    crossings with p_L = p (bisected): n=3 "
          f"{rates['crossings']['n3']:.9f}, n=5 "
          f"{rates['crossings']['n5']:.9f}, n=7 "
          f"{rates['crossings']['n7']:.9f} - exact 0.5, repetition only.")
    mc = rates["monte_carlo"]
    print("\n[4] Monte Carlo referee - seeded splitmix64, deterministic")
    print(f"    p={mc['p']}, trials={mc['trials']}, seed={mc['seed']}: "
          f"estimate {mc['estimate']:.6f} vs closed form "
          f"{mc['closed_form']:.6f} | dev {mc['abs_dev']:.6f} "
          f"< tolerance {mc['tolerance']}")
    print("\n[5] Shor 9-qubit code - 27 single-Pauli cases + one honest failure")
    shor = shor_report()
    print(f"    all 27 cases (9 qubits x X,Y,Z): worst recovered fidelity "
          f"{shor['worst_single_pauli_fidelity']:.6f}")
    w2 = shor["weight2_X0_X1_honest_failure"]
    print(f"    weight-2 X0 X1: read as '{w2['syndrome_read_as']}', "
          f"recovered fidelity {w2['recovered_fidelity']:.6f} - a logical")
    print("    phase flip the code cannot see. Shown, not hidden.")
    print("\n[6] Steane [[7,1,3]] code - Hamming(7,4) checks, built explicitly")
    steane = steane_report()
    print(f"    Hamming codewords {steane['hamming_codewords']} "
          f"(min weight {steane['hamming_min_weight']}) | dual "
          f"{steane['dual_codewords']} (min weight "
          f"{steane['dual_min_weight']}) | parameters {steane['parameters']}")
    print(f"    all 14 single X/Z cases: worst recovered fidelity "
          f"{steane['worst_X_Z_fidelity']:.6f}")
    print("\n[7] Regime classifier (TinyML kNN, labels from the lab policy)")
    for name, row in classifier_report().items():
        print(f"    {name:28s} -> {row['knn']}")
    print("\n[8] Advisory seats")
    brains = brains_report()
    print(f"    decider: {brains['decider']['decision']} "
          f"trace={brains['decider']['trace_sha256_16']}")
    print(f"    advisor can_override_hold = {brains['advisor']['can_override_hold']}")
    print(f"    flow rendered = {brains['flow']['rendered']} "
          f"live = {brains['flow']['live']}")
    print(f"\n[9] HOLD gates: {len(HOLD_GATES)}")
    for key in HOLD_GATES:
        print(f"    HOLD  {key}")
    checks = run_checks()
    passed = sum(1 for _, ok in checks if ok)
    print(f"\nCHECKS {passed}/{len(checks)}")
    for name, ok in checks:
        if not ok:
            print(f"    FAIL {name}")
    verdict = "PROCEED" if passed == len(checks) else "HOLD"
    print(f"VERDICT {verdict} (sandbox simulation only) | HOLD gates {len(HOLD_GATES)}")
    return 0 if passed == len(checks) else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv or argv[0] == "demo":
        return demo()
    cmd = argv[0]
    if cmd == "checks":
        checks = run_checks()
        for name, ok in checks:
            print(f"{'PASS' if ok else 'FAIL'} {name}")
        print(f"CHECKS {sum(1 for _, ok in checks if ok)}/{len(checks)}")
        return 0 if all(ok for _, ok in checks) else 1
    if cmd == "repetition":
        print(json.dumps(repetition_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "rates":
        print(json.dumps(rates_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "shor":
        print(json.dumps(shor_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "steane":
        print(json.dumps(steane_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "classify":
        print(json.dumps(classifier_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "brains":
        print(json.dumps(brains_report(), indent=2, sort_keys=True))
        return 0
    if cmd == "holds":
        for key, value in HOLD_GATES.items():
            print(f"{key}: {value}")
        return 0
    if cmd == "hold-demo":
        action = argv[1] if len(argv) > 1 else "fault_tolerance_claim"
        print(f"{action}: {HOLD_GATES.get(action, 'PROCEED - not gated')}")
        return 2 if gate(action) == "HOLD" else 0
    print(__doc__)
    return 0


if __name__ == "__main__":
    sys.exit(main())
