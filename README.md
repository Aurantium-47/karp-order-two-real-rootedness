# Karp's order-two rising-factorial transform

**Author: Aurantium.** The [complete English proof](proof.pdf) and [editable TeX](proof.tex) establish the following positive, full-support statement.

For every integer `n>=3`, real `c>0` and real `a_1,...,a_n>0`, with repetitions allowed, write

```text
F(t) = c product_(i=1)^n (1+a_i t) = sum_(k=0)^n f_k t^k.
P_n(x) = sum_(k=0)^n binom(n,k) f_k f_(n-k)
         [(x)_k (x)_(n-k) - (x+1)_k (x-1)_(n-k)],
```

where `(x)_k` is the rising factorial and `(x)_0=1`. Then `P_n` is nonzero, has degree `n-2`, and all its zeros are strictly negative and real.

More precisely, let `m=floor(n/2)`, `epsilon=n-2m`, and `r=m-1`. The proof gives

```text
P_n(x) = 2^n (x+1)_(m-1+epsilon) R_n(x).
```

The residual has degree `r`, positive leading coefficient, and simple negative roots in `[1/2-m,0)`, with consecutive gaps at least one. For `n=3` the residual is a positive constant. **Simplicity is asserted for the residual only:** the two factors can share a root, so `P_n` itself may have repeated roots. The symbol `r` in this factorization is a residual degree, not the Toeplitz determinant order, which is two throughout the target.

The argument is uniform in all positive real parameters. It preserves the exact binomial weights and the shifts `x+1,x-1`. Zero boundary parameters, zero padding below the declared degree, arbitrary shifts, and higher-order Toeplitz determinants are outside the theorem. The paper also proves a sufficient condition for positive symmetric weights whose ordinary generating polynomial is negative-rooted; it does not assert the same result for arbitrary positive or merely log-concave weights.

The proof includes coefficient convolution, reciprocal root pairing, a complete beta-moment factorization, a root-spacing induction, degree padding and the boundary cases. Karp's original degree/coefficient-positivity result and the exact source conjectures are identified in [REFERENCES.md](REFERENCES.md). No worldwide-priority claim is made.

## Files and reproduction

- `proof.tex`, `proof.pdf`: the complete natural-language mathematical proof.
- `verify.py`, `verify_low_degree.py`: exact rational/SymPy calculations and the retained low-degree helper.
- `data/archived_168_records.json`: archived inputs and exact outputs for twelve parameter vectors at each `n=3,...,16`.
- `data/verification_report.json`: the delivered historical report, including complete rational inputs and output coefficients.
- `requirements.txt`: the recorded SymPy version; use Python 3.10 or newer.
- `check_package.py`, `PUBLIC_FILES.txt`, `SHA256SUMS.json`: file-integrity checks, not mathematical verification.

With the dependencies already available, run without `-O`:

```sh
python verify.py --output reproduced_verification_report.json
```

The complete suite was not rerun for this edition. Its stored counts and the role of the auxiliary checker are detailed in [REPRODUCIBILITY.md](REPRODUCIBILITY.md). The code checks negative-root counts and residual square-freeness; the universal mesh and root-location assertions come from the proof, not from those checks.

```sh
python check_package.py
```

This second command checks only sizes and hashes. The mathematical results do not follow from a hash match or from finite tests. This work was prepared with AI assistance; no independent human peer-review or proof-assistant certification is claimed. See [RIGHTS.md](RIGHTS.md) for attribution and reuse status.
