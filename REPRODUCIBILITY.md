# Exact records and execution boundaries

Use Python 3.10 or newer with the SymPy version in `requirements.txt` (1.14.0). The delivered historical report records Python 3.13.5 and SymPy 1.14.0; this is not a claim that the complete suite ran in the preparation environment.

From the package root, with assertions enabled:

```sh
python verify.py --output reproduced_verification_report.json
```

The main verifier reads `data/archived_168_records.json` and checks its fixed SHA256 before using it. It invokes the included `verify_low_degree.py`, which uses assertions; therefore optimized Python mode (`-O`) is rejected. No original input archive, source PDF, network access or external data directory is needed. Use a separate output filename to preserve the delivered report in `data/verification_report.json`.

The retained historical coverage is:

| Component | Recorded range |
|---|---|
| Archived exact examples | 168 cases, twelve parameter vectors at each `n=3,...,16` |
| Additional exact examples | 38 cases: four families at `n=3,4,5,6,8,9,12,16,20`, plus equal parameters at `n=24,32` |
| Linear factorization identity | 437 symmetric-weight basis vectors for `n=3,...,40` |
| Generic transform recurrence | `d=0,...,8`, symbolic beta, factor parameter and coefficients |
| Degree padding | 144 cases; four positive beta values, `0<=r<=7`, `0<=d<=r` |
| Low-degree helper | `n=3,...,7` identities and the associated symbolic discriminant identities |

All polynomial arrays use ascending powers. Rational values are exact strings such as `3/2`. The report includes full input parameters and coefficients rather than floating-point root approximations. Square-free factorization and Sturm counts restore multiplicities explicitly. A positive constant has zero roots and is treated separately.

The verifier checks `P=2^n D R`, degree, negative-root counts and residual square-freeness. It does **not** certify universal root spacing or the residual interval bounds; these are established in the manuscript. Neither 206 positive examples nor a finite list of basis-vector identities proves the theorem for every n and every positive real parameter.

The complete mathematical suite was not rerun for this edition. The saved `status` fields and hashes identify delivered records; they do not independently authenticate past execution. Legacy report keys beginning `supplied_` identify the archived-data and low-degree-helper checks; their mathematical schema has been retained so that the report data remain unchanged.

For a newly reproduced JSON report, compare all mathematical fields to `data/verification_report.json`; the `environment` field can differ with the actual interpreter and dependency versions. Do not discard other fields or replace exact comparisons with rounded numerical root lists.

```sh
python check_package.py
```

This performs only file existence, byte-size and SHA256 checks against `SHA256SUMS.json`, which excludes itself. It neither imports nor runs either mathematical verifier. The manuscript source is standalone; rebuilding its PDF requires the listed standard LaTeX packages. PDF timestamps can change the PDF hash after a new build, while the source hash identifies the mathematical text.
