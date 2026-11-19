"""Tests for QuantumErrorCorrection-Lab. Every test grades a computed
number against a closed form or an independent implementation in
main.py. Tests import ONLY main.py - no fixtures, no binary files."""
import math
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import main as qec

TOL = 1e-9


def test_jacobi_known_spectra():
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(qec.hermitian_eigenvalues(qec.PAULI_Y), (-1.0, 1.0)))
    xx_zz = qec.madd(qec.kron(qec.PAULI_X, qec.PAULI_X),
                     qec.kron(qec.PAULI_Z, qec.PAULI_Z))
    assert all(abs(a - b) < 1e-12 for a, b in
               zip(qec.hermitian_eigenvalues(xx_zz), (-2.0, 0.0, 0.0, 2.0)))


def test_pauli_actions_unitary():
    base = qec.encode_repetition_bit()
    for q in range(3):
        for p in ("X", "Y", "Z"):
            out = qec.apply_pauli(base, 3, q, p)
            assert abs(qec.state_norm(out) - 1.0) < 1e-12
            back = qec.apply_pauli(out, 3, q, p)
            assert qec.fidelity(base, back) > 1 - 1e-12


def test_repetition_bit_corrects_X():
    orig = qec.encode_repetition_bit()
    for q in range(3):
        rec, syn = qec.recover_repetition_bit(
            qec.apply_pauli(orig, 3, q, "X"))
        assert qec.fidelity(orig, rec) > 1 - 1e-12
        assert syn != (0, 0)


def test_repetition_phase_corrects_Z():
    orig = qec.encode_repetition_phase()
    for q in range(3):
        rec, syn = qec.recover_repetition_phase(
            qec.apply_pauli(orig, 3, q, "Z"))
        assert qec.fidelity(orig, rec) > 1 - 1e-12
        assert syn != (0, 0)


def test_repetition_syndrome_table():
    assert len(qec.REP_SYNDROME_TABLE) == 4
    orig = qec.encode_repetition_bit()
    for syn, q in qec.REP_SYNDROME_TABLE.items():
        if q is not None:
            _, got = qec.recover_repetition_bit(
                qec.apply_pauli(orig, 3, q, "X"))
            assert got == syn


def test_repetition_basis_all_pauli():
    for encode, recover in (
            (qec.encode_repetition_bit, qec.recover_repetition_bit),
            (qec.encode_repetition_phase, qec.recover_repetition_phase)):
        basis = encode(1.0, 0.0)
        for q in range(3):
            for p in ("X", "Y", "Z"):
                rec, _ = recover(qec.apply_pauli(basis, 3, q, p))
                assert qec.fidelity(basis, rec) > 1 - 1e-12


def test_repetition_off_type_is_logical_honest():
    orig = qec.encode_repetition_bit()
    rec, syn = qec.recover_repetition_bit(qec.apply_pauli(orig, 3, 0, "Z"))
    assert syn == (0, 0)
    assert abs(qec.fidelity(orig, rec) - 0.0784) < 1e-12


def test_n3_closed_vs_enumerated():
    for p in qec.PARAM_GRID:
        assert abs(qec.logical_rate_closed_n3(p)
                   - qec.logical_rate_enumerated_n3(p)) < 1e-12


def test_n3_spot_values():
    assert abs(qec.logical_rate_closed_n3(0.2) - 0.104) < 1e-12
    assert abs(qec.logical_rate_closed_n3(0.5) - 0.5) < 1e-12


def test_crossings_at_half():
    for n in (3, 5, 7):
        assert abs(qec.logical_crossing(n) - 0.5) < 1e-9
        assert abs(qec.logical_rate_binomial(0.5, n) - 0.5) < 1e-12


def test_n5_polynomial():
    for p in qec.PARAM_GRID:
        expect = 10 * p ** 3 - 15 * p ** 4 + 6 * p ** 5
        assert abs(qec.logical_rate_binomial(p, 5) - expect) < 1e-12


def test_codes_help_below_half():
    for n in (3, 5, 7):
        assert qec.logical_rate_binomial(0.2, n) < 0.2
        assert qec.logical_rate_binomial(0.8, n) > 0.8


def test_monte_carlo_deterministic_and_close():
    mc = qec.monte_carlo_n3()
    assert mc == qec.monte_carlo_n3()
    assert mc["trials"] == 20000
    assert mc["abs_dev"] < mc["tolerance"]


def test_splitmix64_stream():
    stream = qec.splitmix64_stream(999)
    draws = [next(stream) for _ in range(500)]
    assert all(0.0 <= v < 1.0 for v in draws)
    assert len(set(draws)) == 500


def test_shor_encode_and_stabilizers():
    orig = qec.encode_shor()
    assert abs(qec.state_norm(orig) - 1.0) < 1e-12
    for s in qec.SHOR_BIT_STABILIZERS + qec.SHOR_PHASE_STABILIZERS:
        assert abs(qec.expectation(orig, 9, s) - 1.0) < 1e-9


def test_shor_corrects_all_27():
    orig = qec.encode_shor()
    for q in range(9):
        for p in ("X", "Y", "Z"):
            rec, _ = qec.recover_shor(qec.apply_pauli(orig, 9, q, p))
            assert qec.fidelity(orig, rec) > 1 - 1e-9


def test_shor_weight2_honest_failure():
    orig = qec.encode_shor()
    err = qec.apply_pauli(qec.apply_pauli(orig, 9, 0, "X"), 9, 1, "X")
    rec, _ = qec.recover_shor(err)
    assert abs(qec.fidelity(orig, rec) - 0.0784) < 1e-9
    flipped = qec.encode_shor(qec.STATED_ALPHA, -qec.STATED_BETA)
    assert qec.fidelity(flipped, rec) > 1 - 1e-9


def test_hamming_combinatorics():
    code = qec.hamming_codewords()
    dual = qec.dual_codewords()
    assert len(code) == 16
    assert min(sum(w) for w in code if any(w)) == 3
    assert len(dual) == 8
    assert min(sum(w) for w in dual if any(w)) == 4
    assert all(w in code for w in dual)
    assert len(qec.STEANE_COLUMN_TABLE) == 7


def test_steane_encode_and_stabilizers():
    orig = qec.encode_steane()
    assert abs(qec.state_norm(orig) - 1.0) < 1e-12
    for s in qec.STEANE_Z_STABILIZERS + qec.STEANE_X_STABILIZERS:
        assert abs(qec.expectation(orig, 7, s) - 1.0) < 1e-9


def test_steane_corrects_X_and_Z_14():
    orig = qec.encode_steane()
    for q in range(7):
        for p in ("X", "Z"):
            rec, _ = qec.recover_steane(qec.apply_pauli(orig, 7, q, p))
            assert qec.fidelity(orig, rec) > 1 - 1e-9


def test_steane_syndromes_are_columns():
    orig = qec.encode_steane()
    for q in range(7):
        _, (bit_syn, _) = qec.recover_steane(
            qec.apply_pauli(orig, 7, q, "X"))
        assert bit_syn == tuple(qec.HAMMING_H[i][q] for i in range(3))


def test_classifier_scenarios():
    rep = qec.classifier_report()
    assert rep["clean_n3_p0.0_w0"]["knn"] == "correctable-likely"
    assert rep["at_limit_n3_p0.1_w1"]["knn"] == "marginal"
    assert rep["overweight_n3_p0.1_w2"]["knn"] == "logical-failure-likely"
    assert all(row["knn"] == row["policy"] for row in rep.values())


def test_classifier_grid_consistency():
    for feats, label in qec.TRAINING_GRID:
        assert qec.classify_features(feats) == label


def test_hold_gates():
    assert len(qec.HOLD_GATES) == 14
    assert qec.gate("fault_tolerance_claim") == "HOLD"
    assert qec.gate("model_as_measurement_claim") == "HOLD"
    assert qec.gate("live_hardware_control") == "HOLD"
    assert qec.gate("investment_claim") == "HOLD"
    assert qec.gate("certification_claim") == "HOLD"
    assert qec.gate("run_simulation") == "PROCEED"


def test_seats_never_override():
    brains = qec.brains_report()
    assert brains["decider"]["can_override_hold"] is False
    assert brains["decider"]["model_id"] == "gemini-argon-4"
    assert brains["advisor"]["can_override_hold"] is False
    assert brains["flow"]["rendered"] is False
    assert brains["flow"]["live"] == 0


def test_checks_all_pass():
    checks = qec.run_checks()
    assert all(ok for _, ok in checks), [n for n, ok in checks if not ok]
    assert len(checks) == 64


def test_cli_hold_demo_exit_2():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "hold-demo"],
        capture_output=True, text=True)
    assert proc.returncode == 2
    assert "HOLD" in proc.stdout


def test_cli_checks_command():
    proc = subprocess.run(
        [sys.executable, str(Path(__file__).resolve().parent.parent
                             / "main.py"), "checks"],
        capture_output=True, text=True)
    assert proc.returncode == 0
    assert "CHECKS 64/64" in proc.stdout


def test_demo_deterministic():
    main_py = str(Path(__file__).resolve().parent.parent / "main.py")
    a = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    b = subprocess.run([sys.executable, main_py, "demo"],
                       capture_output=True, text=True)
    assert a.stdout == b.stdout and a.returncode == 0
