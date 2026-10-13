#!/usr/bin/env python3
"""Exact, independent verification accompanying the all-degree proof in proof.tex.

Python >=3.10 and SymPy >=1.12. No network access and no floating-point roots.
The finite calculations audit, but do not substitute for, the mathematical proof.
All coefficient arrays have ASCENDING powers; exact rationals serialize as strings.
Run: python verify.py --output verification_report.json
"""
from __future__ import annotations

import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import platform
import runpy
import sys
from typing import Sequence

import sympy as sp

HERE = Path(__file__).resolve().parent
X = sp.Symbol("x")


def require(test: bool, message: str) -> None:
    """Use explicit checks so Python -O cannot disable verification."""
    if not test:
        raise ArithmeticError(message)


def trim(p: Sequence) -> list:
    q = list(p)
    while len(q) > 1 and q[-1] == 0:
        q.pop()
    return q


def mul(p: Sequence, q: Sequence) -> list:
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i+j] += a*b
    return trim(out)


def scaled_add(p: list, q: Sequence, scale=1) -> None:
    if len(p) < len(q):
        p.extend([0] * (len(q)-len(p)))
    for i, value in enumerate(q):
        p[i] += scale*value


@lru_cache(None)
def rising(length: int, shift=0) -> tuple:
    require(length >= 0, "negative rising-factorial length")
    p = [1]
    for j in range(length):
        p = mul(p, [shift+j, 1])
    return tuple(p)


@lru_cache(None)
def kernel(n: int, k: int) -> tuple:
    p = list(mul(rising(k), rising(n-k)))
    scaled_add(p, mul(rising(k, 1), rising(n-k, -1)), -1)
    return tuple(trim([comb(n, k)*v for v in p]))


def direct_transform_weights(w: Sequence) -> list:
    n = len(w)-1
    p = [0]
    for k in range(n+1):
        scaled_add(p, kernel(n, k), w[k])
    return trim(p)


def input_coefficients(parameters: Sequence, c=1) -> list:
    p = [Fraction(c)]
    for a in parameters:
        require(a > 0, "input parameters must be strictly positive")
        p = mul(p, [1, a])
    require(c > 0, "c must be strictly positive")
    return p


def symmetric_bernstein_to_q(w: Sequence) -> list:
    """Find Q exactly from sum C(n,k)w_k t^k(1-t)^(n-k)=Q(4t(1-t))."""
    n = len(w)-1
    require(list(w) == list(w)[::-1], "weights must be palindromic")
    b = [0] * (n+1)
    for k, wk in enumerate(w):
        for h in range(n-k+1):
            b[k+h] += comb(n,k)*wk*comb(n-k,h)*(-1)**h
    q = [0]*(n//2+1)
    for j in range(n//2, -1, -1):
        q[j] = Fraction(b[2*j], (-4)**j)
        for h in range(j+1):
            b[j+h] -= q[j]*4**j*comb(j,h)*(-1)**h
    require(all(v == 0 for v in b), "Bernstein-to-Q remainder is nonzero")
    return trim(q)


def beta_transform(g: Sequence, beta, r: int) -> list:
    require(len(trim(g))-1 <= r, "input degree exceeds r")
    out = [0]
    for j, gj in enumerate(g):
        term = mul(rising(j), rising(r-j, beta+j))
        scaled_add(out, term, gj)
    return trim(out)


def derivative(p: Sequence) -> list:
    return trim([j*p[j] for j in range(1, len(p))] or [0])


def sympy_poly(p: Sequence) -> sp.Poly:
    coeffs = [sp.Rational(v.numerator, v.denominator)
              if isinstance(v, Fraction) else sp.Rational(v) for v in p]
    return sp.Poly.from_list(coeffs[::-1], X, domain=sp.QQ)


def negative_count(p: Sequence) -> int:
    """Exact Sturm counts, including multiplicities via square-free factorization."""
    polynomial = sympy_poly(p)
    require(not polynomial.is_zero, "zero polynomial in root count")
    require(polynomial.eval(0) != 0, "zero endpoint in strict-negative root count")
    return sum(int(mult)*int(factor.count_roots(-sp.oo, 0))
               for factor, mult in polynomial.sqf_list()[1])


def serialize(p: Sequence) -> list[str]:
    return [str(v) for v in p]


def audit_one(parameters: Sequence, c=1, label="", archived=None) -> dict:
    a = [Fraction(v) for v in parameters]
    c = Fraction(c)
    n = len(a)
    f = input_coefficients(a, c)
    w = [f[k]*f[n-k] for k in range(n+1)]
    v = [comb(n,k)*w[k] for k in range(n+1)]
    q = symmetric_bernstein_to_q(w)
    g = derivative(q)
    m = n//2
    r = m-1
    D = list(rising((n-1)//2, 1))
    R = beta_transform(g, Fraction(3,2), r)
    P = direct_transform_weights(w)
    rhs = [2**n * value for value in mul(D,R)]
    require(P == rhs, f"factorization identity failed at {label}")
    require(len(P)-1 == n-2, f"degree failed at {label}")
    require(all(value > 0 for value in P), f"coefficient positivity failed at {label}")
    require(q[0] == w[0], f"q0 identity failed at {label}")
    require(q[1] == Fraction(n,4)*(w[1]-w[0]) > 0, f"q1 failed at {label}")
    require(all(value > 0 for value in q), f"Q coefficients failed at {label}")
    require(len(R)-1 == r, f"residual degree failed at {label}")
    counts = {}
    for name, p in [("W",w), ("V",v), ("Q",q), ("Q_prime",g), ("R",R), ("P",P)]:
        counts[name] = negative_count(p)
        require(counts[name] == len(p)-1, f"negative root count failed: {label}/{name}")
    require(sp.gcd(sympy_poly(R), sympy_poly(derivative(R))).degree() == 0,
            f"residual simplicity failed at {label}")
    if archived is not None:
        require(f == archived['f'], f"archived f mismatch: {label}")
        original = [archived['content']*v for v in archived['P_primitive_ascending']]
        require(P == original, f"archived P mismatch: {label}")
        require(counts['P'] == archived['negative_roots_with_multiplicity'],
                f"archived count mismatch: {label}")
    return {
        "label": label, "n": n, "c": str(c), "root_parameters": serialize(a),
        "f_ascending": serialize(f), "Q_ascending": serialize(q),
        "R_ascending": serialize(R), "P_ascending": serialize(P),
        "negative_roots_with_multiplicity": counts,
        "residual_square_free": True,
    }


def audit_linear_identity() -> dict:
    """Verify the identity for every basis vector of symmetric w, n=3,...,40."""
    records = []
    for n in range(3,41):
        m = n//2
        for k in range(m+1):
            w = [int(j == k or j == n-k) for j in range(n+1)]
            q = symmetric_bernstein_to_q(w)
            R = beta_transform(derivative(q), Fraction(3,2), m-1)
            rhs = [2**n*v for v in mul(rising((n-1)//2,1),R)]
            require(direct_transform_weights(w) == trim(rhs),
                    f"linear identity basis failure at n={n},k={k}")
        records.append({"n": n, "independent_symmetric_basis_vectors": m+1,
                        "all_zero_remainders": True})
    return {"records": records, "total_basis_vectors": sum(r['independent_symmetric_basis_vectors'] for r in records)}


def audit_recurrence() -> dict:
    beta, a = sp.symbols('beta a')
    for d in range(9):
        b = sp.symbols('b0:'+str(d+1))
        R = sum(b[j]*sp.rf(X,j)*sp.rf(X+beta+j,d-j) for j in range(d+1))
        new_b = [b[0]] + [b[j]+a*b[j-1] for j in range(1,d+1)] + [a*b[d]]
        S = sum(new_b[j]*sp.rf(X,j)*sp.rf(X+beta+j,d+1-j) for j in range(d+2))
        require(sp.expand(S-(X+beta+d)*R-a*X*R.subs(X,X+1)) == 0,
                f"symbolic recurrence failure at d={d}")
    return {"d_values": list(range(9)), "generic_beta_a_and_input_coefficients": True,
            "all_zero_remainders": True}


def additional_inputs():
    for n in (3,4,5,6,8,9,12,16,20):
        yield "equal_parameters", [Fraction(1)]*n, Fraction(3,7)
        yield "heterogeneous_rational", [Fraction(i+1,n-i+1) for i in range(n)], Fraction(1)
        scale = [Fraction(1,10000), Fraction(1,100), Fraction(1), Fraction(100), Fraction(10000)]
        yield "multiscale", [scale[i%len(scale)] for i in range(n)], Fraction(7,5)
        yield "near_coincident", [1+Fraction(i-n//2,10**6) for i in range(n)], Fraction(11,13)
    for n in (24,32):
        yield "higher_degree_equal", [Fraction(1)]*n, Fraction(1)


def audit_degree_padding() -> list:
    """Exercise the degree-deficient cases of the transform lemma independently."""
    results = []
    for beta in [Fraction(1,100),Fraction(1,2),Fraction(3,2),Fraction(5)]:
        for r in range(8):
            for d in range(r+1):
                g = input_coefficients([Fraction(j+1,j+2) for j in range(d)])
                out = beta_transform(g,beta,r)
                small = beta_transform(g,beta,d)
                require(out == mul(rising(r-d,beta+d),small), "padding identity failed")
                require(negative_count(out) == r, "padding negative root count failed")
                require(sp.gcd(sympy_poly(out),sympy_poly(derivative(out))).degree() == 0,
                        "padding residual not square-free")
                results.append({"beta":str(beta),"r":r,"input_degree":d,
                                "negative_root_count":r,"square_free":True})
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, default=HERE/'reproduced_verification_report.json')
    args = parser.parse_args()
    source = HERE/'data/archived_168_records.json'
    source_bytes = source.read_bytes()
    expected_hash = 'ea75209aa8ebad900b74368c64f4b32d03cff364b1751e14c1dd65779f37fbe3'
    require(hashlib.sha256(source_bytes).hexdigest() == expected_hash, "supplied data checksum mismatch")
    archived = json.loads(source_bytes)
    require(len(archived) == 168, "unexpected archive record count")
    counts = Counter(row['n'] for row in archived)
    require(counts == {n:12 for n in range(3,17)}, "unexpected archive n coverage")
    # The supplied checker uses asserts; reject optimized mode before executing it.
    require(sys.flags.optimize == 0, "run without -O; the supplied checker uses assertions")
    runpy.run_path(str(HERE/'verify_low_degree.py'),run_name='__main__')
    print('Checking full palindromic linear identity, n=3,...,40.', flush=True)
    identities = audit_linear_identity()
    print('Checking generic recurrence, d=0,...,8.', flush=True)
    recurrence = audit_recurrence()
    replay = []
    for i,row in enumerate(archived):
        replay.append(audit_one(row['root_parameters'],1,f'archive_{i:03d}',row))
        if (i+1)%12 == 0:
            print(f'Archived exact cases verified: {i+1}/168 (through n={row["n"]}).',flush=True)
    new = []
    for family,parameters,c in additional_inputs():
        new.append(audit_one(parameters,c,f'{family}_n{len(parameters)}'))
        print(f'New exact case verified: {new[-1]["label"]}.',flush=True)
    padding = audit_degree_padding()
    report = {
        "status":"PASS", "proof_dependency":"No finite calculation is a premise of the all-degree proof.",
        "arithmetic":"Exact rational coefficients; SymPy square-free factorization and Sturm real-root counts.",
        "environment":{"python":platform.python_version(),"sympy":sp.__version__},
        "supplied_data_sha256":expected_hash,
        "supplied_symbolic_checker":"PASS: n=3,...,7 identities and discriminant identities",
        "linear_identity":identities, "generic_recurrence":recurrence,
        "supplied_cases_replayed":len(replay), "new_exact_cases":len(new),
        "degree_padding_cases":len(padding),
        "archive_replay":replay, "additional_cases":new, "padding_cases":padding,
    }
    args.output.parent.mkdir(parents=True,exist_ok=True)
    args.output.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(f'PASS. Full exact data and checks written to {args.output}.',flush=True)


if __name__ == '__main__':
    main()
