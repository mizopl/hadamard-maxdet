# Agent Primer & Verification Guide (`mizopl/hadamard-maxdet`)

Welcome! This file provides a concise, machine-actionable orientation for autonomous research, coding, and mathematical optimization agents working on the **Hadamard Maximal Determinant Problem** ($1 \le n \le 120$, [OEIS A003432](https://oeis.org/A003432) / [A003433](https://oeis.org/A003433)).

---

## 1. Repository Map & Machine-Readable Endpoints

- `catalog/hadamard_maxdet_n1_to_120.json` — authoritative 120-order catalog (`D_norm_exact`, `det_exact`, `status`, `attains_residue_bound`, `bound_name`, `ratio_to_bound`, `current_holder`, `year`, `ref_key`, `matrix_file`).
- `catalog/hadamard_maxdet_n1_to_120.csv` — flat 120-row CSV mirror of the JSON catalog.
- `matrices/` — explicit $n \times n$ sign matrices (`+` for $+1$, `-` for $-1$) for **49 / 50** open orders.
- `verify.py` & `verifiers/` — zero-dependency exact verification suite (`determinant_and_bounds.py` Bareiss + rational bounds, `circulant_structure.py` exact $\mathbb{F}_2$ block-circulant classifier, `candidate_audit.py` external candidate & Gram-equivalence auditor, `admin_update_catalog.py` maintainer ingestion CLI).
- `CONTRIBUTING.md` — contributor, template-reuse, and co-maintainer workflow.

---

## 2. Target Selection: Do Not Burn Compute on Closed Orders

Out of the 120 orders $1 \le n \le 120$:
- **70 orders have `status == "PROVEN_OPTIMAL"` (CLOSED — cannot be improved):**
  - **61 orders** achieve $R(n) = |\det X| / B_{\text{class}}(n) = 1$ (`attains_residue_bound == true`), including all 30 multiples of 4, orders `1, 2, 3`, 6 Barba orders (`5, 13, 25, 41, 61, 113`), and 22 Ehlich–Wojtas orders.
  - **9 orders (`7, 9, 11, 15, 17, 19, 21, 22, 37`)** have $R(n) < 1$ (`attains_residue_bound == false`) because the closed-form residue bound is unattainable over integers, **yet optimality is already proved** (e.g., `n = 19, 37` by Brent et al. 2011; `n = 22` by Brent 2021). Exceeding the catalog value on any `PROVEN_OPTIMAL` order indicates a bug in your determinant or entry check.
- **50 orders have `status == "OPEN"` (LIVE TARGETS for constructive lower bounds or maximality proofs):**
  - **`n ≡ 1 (mod 4)` (19 open):** `29, 33, 45, 49, 53, 57, 65, 69, 73, 77, 81, 85, 89, 93, 97, 101, 105, 109, 117`
  - **`n ≡ 2 (mod 4)` (6 open):** `34, 58, 70, 78, 94, 106`
  - **`n ≡ 3 (mod 4)` (25 open):** `23, 27, 31, 35, 39, 43, 47, 51, 55, 59, 63, 67, 71, 75, 79, 83, 87, 91, 95, 99, 103, 107, 111, 115, 119`

Quick Python snippet to load all open targets and their exact integer records:
```python
import json

cat = json.load(open("catalog/hadamard_maxdet_n1_to_120.json"))
open_targets = {row["n"]: int(row["D_norm_exact"]) for row in cat["orders"] if row["status"] == "OPEN"}
```

---

## 3. Five Classic Verification Traps to Avoid

1. **Never certify determinants in `float64` (`numpy.linalg.det` / `slogdet`):**
   IEEE-754 `float64` has only 53 bits ($\approx 15.9$ decimal digits) of mantissa precision, whereas $D_{\text{norm}}(119)$ has 88 digits and nearby records (e.g., `n = 103`, `+0.00004%`) agree on the first 6+ digits. Use `float64` Sherman–Morrison updates *only* as a fast search heuristic, and always re-evaluate candidates in exact Python `int` Bareiss arithmetic (`verify.det_bareiss_int`).
2. **Normalization factor $2^{n-1}$ (`OEIS A003432` vs `OEIS A003433`):**
   Every $n \times n$ matrix $X \in \{-1, +1\}^{n \times n}$ satisfies $2^{n-1} \mid \det X$. This catalog's primary integer `D_norm` is $D_{\text{norm}}(n) = |\det X| / 2^{n-1}$ (`OEIS A003432`, equal to the maximal $(n-1)\times(n-1)$ binary $\{0,1\}$ determinant). Raw $|\det X|$ is stored in `det_exact` (`OEIS A003433`).
3. **Piecewise block count $s$ in Ehlich's `n ≡ 3 (mod 4)` upper bound:**
   For $n \equiv 3 \pmod 4$, Ehlich's block-Gram bound uses $s = 3$ ($n = 3$), $s = 5$ ($n \in \{7, 11\}$), $s = 6$ ($15 \le n \le 59$), and $s = 7$ ($n \ge 63$). Applying the $s = 7$ formula below $n = 63$ gives an invalid bound. Call `verify.residue_bound_squared_exact(n)` for the exact rational value of $B_{\text{class}}(n)^2$.
4. **Strict entry domain:**
   Every entry must be strictly $-1$ or $+1$ (`0` or out-of-range entries invalidate Hadamard/Barba bounds).
5. **Respect third-party platform rules:**
   When cross-checking external leaderboards (such as `EinsteinArena` on `n = 51`), never solve, compute, or bypass anti-bot, CAPTCHA, or proof-of-work registration gates from an automated agent; validate locally via `python3 verify.py --check` and leave interactive browser submission to a human operator.

---

## 4. Mandatory Verification Commands & Contribution Format

Before reporting a new record or opening a Pull Request / Issue:

```bash
# 1. Check your candidate file (+/- grid, JSON 2D array, CSV +/-1, or 0/1 binary matrix)
python3 verify.py --check path/to/your_matrix.txt --sympy

# 2. Run the full repository verification suite & mutation self-test
python3 verify.py
python3 verify.py --corrupt   # Must print CORRUPT_TEST_CAUGHT and exit with code 1
```

When contributing a new record or an alternative non-isomorphic witness tying an existing record, please include **both**:
1. The explicit `n x n` `+/-` matrix file, and
2. A concise **algebraic construction description / generator formula** (e.g., circulant first rows `(a, b)` and block `Type`, Goethals–Seidel / Williamson generator tuples, or Hadamard switching coordinates) so others can reproduce and generalize the structure.
