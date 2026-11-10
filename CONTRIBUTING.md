# Contributing, Template Reuse & Editorial Co-Stewardship (`mizopl/hadamard-maxdet`)

Thank you for contributing to or building upon the **Hadamard Maximal Determinant Catalog ($1 \le n \le 120$)**. Every mathematical claim, witness matrix, and structural tag (`*`) in this repository is machine-checked from first principles with zero external runtime dependencies.

---

## 1. Submitting a New Lower-Bound Record or Optimality Proof

1. **Verify your candidate matrix locally** (accepts `+/-` sign grids, JSON 2D arrays `[[1, -1, ...]]`, CSV/whitespace `±1` matrices, or order-$(n-1)$ `{0, 1}` binary matrices):
   ```bash
   python3 verify.py --check path/to/your_matrix.txt --sympy
   ```
   Confirm that the output prints `VERDICT : NEW RECORD (+...% strictly above catalog best!)`.
2. **Open a Pull Request or GitHub Issue** containing:
   - the explicit $n \times n$ `{-1, +1}` matrix file (`+/-` format preferred),
   - a concise **algebraic construction formula / generator description** (e.g., circulant first rows `(a, b)` and block `Type`, Goethals–Seidel / Williamson 4-block generators, or Hadamard switching coordinates), and
   - the bibliographic reference / author attribution to appear in `catalog/hadamard_maxdet_n1_to_120.json` and `README.md`.

---

## 2. Submitting an Alternative Non-Isomorphic Record-Tying Construction

Structured algebraic constructions that **tie** an existing open record—especially **block-circulant (`*`)** or closed-form generator families replacing an unstructured local-search grid—are warmly welcomed:
1. Run `python3 verify.py --check path/to/your_matrix.txt` and inspect both:
   - `Structure (*) : ...` (verified over $\mathbb{F}_2$ up to row/column sign flips by [`verifiers/circulant_structure.py`](verifiers/circulant_structure.py)), and
   - `Gram invariant: ...` (confirming whether your matrix is provably Hadamard-inequivalent to the archived witness).
2. Include both the `+/-` matrix and its compact generator formula in your PR/Issue.

---

## 3. Architecture of the `verifiers/` Suite (Reusing This Repo as a Template)

All verification and maintenance code lives openly in [`verify.py`](verify.py) and [`verifiers/`](verifiers/) so other open-problem catalogs can reuse this exact-verification pattern:

| File | Role | Standalone CLI |
| :--- | :--- | :--- |
| [`verify.py`](verify.py) | Master 5-stage repository verifier (`--sympy`, `--corrupt`, `--check FILE`) | `python3 verify.py` |
| [`verifiers/determinant_and_bounds.py`](verifiers/determinant_and_bounds.py) | Exact fraction-free Bareiss integer determinant + rational Hadamard / Barba / Ehlich–Wojtas / Ehlich($s$) bounds | `python3 verifiers/determinant_and_bounds.py matrices/literature_witnesses/n051.txt` |
| [`verifiers/circulant_structure.py`](verifiers/circulant_structure.py) | Exact $\mathbb{F}_2$ bitmask Gaussian-elimination solver certifying 1-, 2-, and 4-block circulant/back-circulant cores (`*`) under arbitrary row/column sign gauges | `python3 verifiers/circulant_structure.py --inspect matrices/literature_witnesses/n051.txt` |
| [`verifiers/candidate_audit.py`](verifiers/candidate_audit.py) | Flexible external matrix parser (`+/-`, JSON, CSV, `{0,1}` binary) + Hadamard-equivalence Gram invariant fingerprinter | `python3 verifiers/candidate_audit.py candidate.json` |
| [`verifiers/admin_update_catalog.py`](verifiers/admin_update_catalog.py) | Maintainer CLI for ingesting verified new records into JSON/CSV/`matrices/`/`README.md` and refreshing `SHA256SUMS` | `python3 verifiers/admin_update_catalog.py --refresh-hashes` |

---

## 4. Maintainer Workflow & Open Invitation for Editorial Co-Stewardship

Co-maintainers can ingest a verified record, update `catalog/hadamard_maxdet_n1_to_120.json`, `.csv`, and both `README.md` leaderboard tables (including automatic `block_circulant_witness` `*` detection), recompute `SHA256SUMS`, and re-run the full verification suite in a single command:

```bash
python3 verifiers/admin_update_catalog.py \
  --ingest path/to/verified_candidate.txt \
  --holder "Author Name(s)" \
  --year 2026 \
  --ref-key Auth2026 \
  --citation "Author Name(s), Title, Journal / arXiv (2026)."
```

**Editorial Co-Stewardship:** This repository was created to keep the classical `1 <= n <= 120` maximal-determinant tables permanently live, machine-readable, and cryptographically verified following the retirement of Indiana University's legacy `~maxdet` web hosting. Researchers active in Hadamard matrices, $D$-optimal designs, or integer sequence curation who would like **co-maintainer / editorial steward access** (or wish to mirror the catalog under an academic institution) are warmly invited to open an Issue or contact `@mizopl`.
