"""Small helpers for code-capacity simulations of CSS codes.

Codes are given by binary parity-check matrices H_X (X-type stabilisers) and
H_Z (Z-type stabilisers).  We study X errors only: they are detected by H_Z and
cause a logical failure if the residual error anticommutes with a logical Z.
"""
import numpy as np


# ---------------------------------------------------------------- GF(2) algebra
def gf2_rref(M):
    """Row-reduce a binary matrix. Returns (R, pivot_columns)."""
    R = (np.array(M, dtype=np.uint8) % 2).copy()
    rows, cols = R.shape
    pivots, r = [], 0
    for c in range(cols):
        if r == rows:
            break
        hits = np.nonzero(R[r:, c])[0]
        if hits.size == 0:
            continue
        p = r + hits[0]
        R[[r, p]] = R[[p, r]]
        others = np.nonzero(R[:, c])[0]
        others = others[others != r]
        R[others] ^= R[r]
        pivots.append(c)
        r += 1
    return R[:r], pivots


def gf2_rank(M):
    return len(gf2_rref(M)[1])


def gf2_nullspace(M):
    """Basis (as rows) of {x : M x = 0 mod 2}."""
    M = np.array(M, dtype=np.uint8) % 2
    n = M.shape[1]
    R, piv = gf2_rref(M)
    free = [c for c in range(n) if c not in piv]
    basis = []
    for f in free:
        v = np.zeros(n, dtype=np.uint8)
        v[f] = 1
        for i, p in enumerate(piv):
            v[p] = R[i, f]
        basis.append(v)
    return np.array(basis, dtype=np.uint8).reshape(-1, n)


def logical_z_operators(H_X, H_Z):
    """Z-type logicals: kernel of H_X modulo the row space of H_Z."""
    ker = gf2_nullspace(H_X)
    # Greedy: keep kernel vectors that increase the rank beyond rowspace(H_Z).
    base_rank = gf2_rank(H_Z)
    chosen = []
    current = H_Z.copy() % 2
    for v in ker:
        trial = np.vstack([current, v])
        if gf2_rank(trial) > gf2_rank(current):
            chosen.append(v)
            current = trial
    L = np.array(chosen, dtype=np.uint8)
    assert gf2_rank(current) - base_rank == L.shape[0]
    return L


# ---------------------------------------------------------------- code families
def cyclic_shift(size):
    return np.roll(np.eye(size, dtype=np.uint8), 1, axis=1)


def bivariate_bicycle(l, m, a_terms, b_terms):
    """Bivariate bicycle code (Bravyi et al., Nature 2024).

    a_terms / b_terms are lists of ('x', power) or ('y', power), e.g. the gross
    code uses A = x^3 + y + y^2 and B = y^3 + x + x^2 with l=12, m=6.
    """
    x = np.kron(cyclic_shift(l), np.eye(m, dtype=np.uint8))
    y = np.kron(np.eye(l, dtype=np.uint8), cyclic_shift(m))
    def poly(terms):
        out = np.zeros((l * m, l * m), dtype=np.uint8)
        for var, power in terms:
            base = x if var == "x" else y
            out ^= np.linalg.matrix_power(base, power).astype(np.uint8) % 2
        return out
    A, B = poly(a_terms), poly(b_terms)
    H_X = np.hstack([A, B]) % 2
    H_Z = np.hstack([B.T, A.T]) % 2
    return H_X.astype(np.uint8), H_Z.astype(np.uint8)


def gross_code():
    """The [[144, 12, 12]] 'gross' code from IBM's qLDPC roadmap."""
    return bivariate_bicycle(12, 6, [("x", 3), ("y", 1), ("y", 2)],
                             [("y", 3), ("x", 1), ("x", 2)])


def repetition_checks(d):
    H = np.zeros((d - 1, d), dtype=np.uint8)
    for i in range(d - 1):
        H[i, i] = H[i, i + 1] = 1
    return H


def planar_surface_code(d):
    """[[d^2 + (d-1)^2, 1, d]] surface code as a hypergraph product of repetition codes."""
    H1 = repetition_checks(d)
    m, n = H1.shape
    H_X = np.hstack([np.kron(H1, np.eye(n, dtype=np.uint8)),
                     np.kron(np.eye(m, dtype=np.uint8), H1.T)]) % 2
    H_Z = np.hstack([np.kron(np.eye(n, dtype=np.uint8), H1),
                     np.kron(H1.T, np.eye(m, dtype=np.uint8))]) % 2
    return H_X.astype(np.uint8), H_Z.astype(np.uint8)


def code_parameters(H_X, H_Z):
    n = H_X.shape[1]
    assert not ((H_X.astype(int) @ H_Z.T.astype(int)) % 2).any(), "checks do not commute"
    k = n - gf2_rank(H_X) - gf2_rank(H_Z)
    return n, k


# ---------------------------------------------------------------- Monte Carlo
def bposd_batch(decoder):
    """Wrap an ldpc BpOsdDecoder so it decodes a batch of syndromes."""
    return lambda S: np.array([decoder.decode(s) for s in S], dtype=np.uint8)


def code_capacity_failures(H_Z, L_Z, p, decode_batch, rng,
                           target_failures=200, max_shots=200_000, chunk=5_000):
    """Monte Carlo estimate of the logical failure probability under i.i.d. X noise.

    Samples until `target_failures` logical failures or `max_shots` shots.
    Returns (failure_rate, failures, shots).
    """
    n = H_Z.shape[1]
    HZt, LZt = H_Z.T.astype(np.uint8), L_Z.T.astype(np.uint8)
    shots = fails = 0
    while fails < target_failures and shots < max_shots:
        k = min(chunk, max_shots - shots)
        E = (rng.random((k, n)) < p).astype(np.uint8)
        S = (E @ HZt) % 2
        R = E ^ decode_batch(S.astype(np.uint8))
        assert not ((R @ HZt) % 2).any(), "decoder output inconsistent with syndrome"
        fails += int(((R @ LZt) % 2).any(axis=1).sum())
        shots += k
    return fails / shots, fails, shots
