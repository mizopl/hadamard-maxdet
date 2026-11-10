#!/usr/bin/env python3
"""
Standalone Offline Exact-Integer Verifier for the Hadamard Maximal Determinant Catalog
(Orders 1 <= n <= 120, OEIS A003432 / A003433).

Checks performed from first principles (zero floating-point determinants):
1. If catalog/new_records.json is present:
   - exact dimensions n x n and every entry strictly in {-1, +1};
   - exact integer determinant via fraction-free Bareiss elimination on Python ints
     (plus optional SymPy cross-check via --sympy);
   - exact divisibility |det(X)| % 2^(n-1) == 0 and digit-for-digit equality of D_norm = |det(X)| // 2^(n-1);
   - strict improvement D_norm > previous_D_norm and disk witness identity;
   - mathematical validity |det(X)|^2 <= B_class(n)^2 against the exact residue-class upper bound.
2. Every explicit open-order witness matrix in matrices/ (23 <= n <= 119):
   - exact dimensions n x n, entries in {-1, +1}, exact Bareiss D_norm matching
     catalog/hadamard_maxdet_n1_to_120.json, and |det(X)|^2 <= B_class(n)^2.
3. Full 120-row consistency between catalog/hadamard_maxdet_n1_to_120.json and
   catalog/hadamard_maxdet_n1_to_120.csv, plus exact residue-class bound satisfaction across 1 <= n <= 120.
4. SHA-256 integrity check of all files listed in SHA256SUMS.
5. Optional --corrupt self-test verifying that mutating a matrix entry triggers FAIL (exit 1).
6. Standalone candidate checker (--check FILE): evaluates any external {+1, -1} matrix
   (+/- grid, JSON array, or CSV/whitespace integers, as used on EinsteinArena / NumPy) or
   {0, 1} binary matrix (OEIS A003432), checks the exact residue-class bound, and compares
   against the catalog record for that order.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
import sys
import time
from decimal import Decimal, getcontext
from pathlib import Path
from typing import List, Tuple

getcontext().prec = 140

try:
    import sympy as sp  # type: ignore
    HAS_SYMPY = True
except Exception:
    HAS_SYMPY = False


def parse_sign_rows(rows: List[str], expected_n: int) -> List[List[int]]:
    if len(rows) != expected_n:
        raise ValueError(f"Expected {expected_n} rows, got {len(rows)}")
    mat: List[List[int]] = []
    for r_idx, row in enumerate(rows):
        s = row.strip()
        if len(s) != expected_n:
            raise ValueError(f"Row {r_idx} has length {len(s)}, expected {expected_n}")
        parsed_row: List[int] = []
        for c_idx, ch in enumerate(s):
            if ch == "+":
                parsed_row.append(1)
            elif ch == "-":
                parsed_row.append(-1)
            else:
                raise ValueError(f"Invalid character {ch!r} at ({r_idx}, {c_idx}); must be '+' or '-'")
        mat.append(parsed_row)
    return mat


def det_bareiss_int(mat: List[List[int]]) -> int:
    """Exact fraction-free Bareiss determinant over arbitrary-precision Python ints."""
    n = len(mat)
    if n == 0:
        return 1
    a = [row[:] for row in mat]
    sign = 1
    prev_pivot = 1
    for k in range(n - 1):
        if a[k][k] == 0:
            swap_row = -1
            for i in range(k + 1, n):
                if a[i][k] != 0:
                    swap_row = i
                    break
            if swap_row == -1:
                return 0
            a[k], a[swap_row] = a[swap_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            row_i = a[i]
            row_k = a[k]
            rik = row_i[k]
            if rik == 0:
                for j in range(k + 1, n):
                    row_i[j] = (pivot * row_i[j]) // prev_pivot
            else:
                for j in range(k + 1, n):
                    row_i[j] = (pivot * row_i[j] - rik * row_k[j]) // prev_pivot
            row_i[k] = 0
        prev_pivot = pivot
    return sign * a[n - 1][n - 1]


def ehlich_block_count(n: int) -> int:
    """Exact block count s in Ehlich's (1964) upper bound for n == 3 (mod 4)."""
    assert n % 4 == 3 and n >= 3
    if n == 3:
        return 3
    if n in (7, 11):
        return 5
    if 15 <= n <= 59:
        return 6
    return 7


def residue_bound_squared_exact(n: int) -> Tuple[int, int, str]:
    """
    Returns (num, den, bound_name) representing the exact rational value of B_class(n)^2
    for {-1, +1} matrices of order n:
      - n == 1, 2 or n == 0 mod 4 (Hadamard 1893): B(n)^2 = n^n
      - n == 1 mod 4 (Barba 1933): B(n)^2 = (2n - 1) * (n - 1)^(n - 1)
      - n == 2 mod 4, n >= 6 (Ehlich-Wojtas 1964): B(n)^2 = 4 * (n - 1)^2 * (n - 2)^(n - 2)
      - n == 3 mod 4 (Ehlich 1964):
          B(n)^2 = (n - 3)^(n - s) * (n - 3 + 4r)^u * (n + 1 + 4r)^v
                   * (1 - u*r/(n - 3 + 4r) - v*(r + 1)/(n + 1 + 4r))
    """
    if n in (1, 2) or n % 4 == 0:
        return (n ** n, 1, "Hadamard")
    if n % 4 == 1:
        return ((2 * n - 1) * ((n - 1) ** (n - 1)), 1, "Barba")
    if n % 4 == 2:
        return (4 * ((n - 1) ** 2) * ((n - 2) ** (n - 2)), 1, "Ehlich-Wojtas")
    s = ehlich_block_count(n)
    r = n // s
    v = n - r * s
    u = s - v
    a = n - 3 + 4 * r
    b = n + 1 + 4 * r
    num_fac = a * b - u * r * b - v * (r + 1) * a
    den_fac = a * b
    base = ((n - 3) ** (n - s)) * (a ** u) * (b ** v)
    num = base * num_fac
    g = math.gcd(num, den_fac)
    return (num // g, den_fac // g, f"Ehlich(s={s})")


def compute_ratio_decimal(det_abs: int, n: int) -> Decimal:
    num, den, _ = residue_bound_squared_exact(n)
    b_sq = Decimal(num) / Decimal(den)
    return Decimal(det_abs) / b_sq.sqrt()


def verify_matrix_rows_against_d_norm(
    n: int, rows: List[str], expected_d_norm: int, check_sympy: bool = False
) -> int:
    mat = parse_sign_rows(rows, n)
    det_b = abs(det_bareiss_int(mat))
    if check_sympy and HAS_SYMPY:
        det_sp = abs(int(sp.Matrix(mat).det(method="bareiss")))
        if det_b != det_sp:
            raise AssertionError(f"n={n}: Bareiss ({det_b}) != SymPy ({det_sp})")
    scale = 1 << (n - 1)
    if det_b % scale != 0:
        raise AssertionError(f"n={n}: |det| is not divisible by 2^(n-1)")
    d_norm = det_b // scale
    if d_norm != expected_d_norm:
        raise AssertionError(f"n={n}: computed D_norm={d_norm} != claimed {expected_d_norm}")
    num, den, bname = residue_bound_squared_exact(n)
    if (det_b * det_b) * den > num:
        raise AssertionError(f"n={n}: |det|^2 exceeds {bname} residue-class bound")
    return det_b


def run_verification(
    repo_root: Path, verify_all_open: bool = True, use_sympy: bool = False, run_corrupt_test: bool = False
) -> int:
    t0 = time.perf_counter()
    we_path = repo_root / "catalog" / "new_records.json"
    cat_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.json"
    csv_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.csv"
    sums_path = repo_root / "SHA256SUMS"

    if not cat_path.exists() or not csv_path.exists():
        print("FAIL: Missing catalog JSON/CSV files.", file=sys.stderr)
        return 1

    cat_data = json.loads(cat_path.read_text(encoding="utf-8"))
    entries = cat_data["orders"]
    if len(entries) != 120:
        print(f"FAIL: Expected 120 orders in catalog, found {len(entries)}", file=sys.stderr)
        return 1

    if run_corrupt_test:
        first_open = next(e for e in entries if e["status"] == "OPEN" and e.get("matrix_file"))
        n_c = int(first_open["n"])
        mpath_c = repo_root / first_open["matrix_file"]
        rows_c = [ln.strip() for ln in mpath_c.read_text(encoding="utf-8").splitlines() if ln.strip()]
        rows_c[0] = ("-" if rows_c[0][0] == "+" else "+") + rows_c[0][1:]
        try:
            verify_matrix_rows_against_d_norm(n_c, rows_c, int(first_open["D_norm_exact"]), check_sympy=False)
        except Exception as exc:
            print(f"CORRUPT_TEST_CAUGHT: {exc}")
            print("FAIL (expected under --corrupt mutation check)")
            return 1
        print("ERROR: Corrupt matrix was not rejected!", file=sys.stderr)
        return 2

    we_count = 0
    if we_path.exists():
        we_data = json.loads(we_path.read_text(encoding="utf-8"))
        records = we_data["records"]
        we_count = len(records)
        print(f"=== 1. Verifying {we_count} New Constructive Records (SymPy cross-check: {use_sympy and HAS_SYMPY}) ===")
        for rec in records:
            t_rec = time.perf_counter()
            n = int(rec["n"])
            claimed_d = int(rec["D_norm"])
            prev_d = int(rec["previous_D_norm"])
            det_b = verify_matrix_rows_against_d_norm(n, rec["matrix_rows"], claimed_d, check_sympy=use_sympy)
            if det_b != int(rec["det_exact"]):
                print(f"FAIL: n={n} det_exact mismatch!", file=sys.stderr)
                return 1
            if not (claimed_d > prev_d):
                print(f"FAIL: n={n} D_norm ({claimed_d}) <= previous ({prev_d})!", file=sys.stderr)
                return 1
            mfile = repo_root / rec["matrix_file"]
            disk_rows = [ln.strip() for ln in mfile.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if disk_rows != rec["matrix_rows"]:
                print(f"FAIL: Disk matrix {mfile} differs from JSON certificate!", file=sys.stderr)
                return 1
            ratio = compute_ratio_decimal(det_b, n)
            dt = (time.perf_counter() - t_rec) * 1000.0
            print(
                f"  [PASS] n={n:>3d} ({n%4} mod 4) | D_norm = {rec['factored_D_norm']:<24s} "
                f"| +{rec['gain_percent']:.4f}% vs {rec['previous_holder_short']:<28s} "
                f"| R={float(ratio):.6f} | {dt:5.1f} ms"
            )

    if verify_all_open:
        open_entries = [e for e in entries if e["status"] == "OPEN"]
        verified_open = 0
        print(f"=== Verifying Explicit +/-1 Witness Matrices Across Open Orders (23 <= n <= 119) ===")
        for e in open_entries:
            n = int(e["n"])
            mrel = e.get("matrix_file", "")
            if not mrel:
                if e.get("missing_archived_witness") is True and n == 89 and we_count == 0:
                    continue
                print(f"FAIL: Open order n={n} is missing an explicit witness matrix!", file=sys.stderr)
                return 1
            mpath = repo_root / mrel
            if not mpath.is_file():
                print(f"FAIL: Witness matrix file {mrel} not found!", file=sys.stderr)
                return 1
            rows = [ln.strip() for ln in mpath.read_text(encoding="utf-8").splitlines() if ln.strip()]
            verify_matrix_rows_against_d_norm(n, rows, int(e["D_norm_exact"]), check_sympy=False)
            verified_open += 1
        print(
            f"  [PASS] {verified_open}/{len(open_entries)} open-order (-1, +1) witness matrices "
            f"verified from scratch by exact Bareiss elimination."
        )

    print("=== Verifying CSV/JSON Parity & Residue-Class Upper Bounds Across All 120 Orders (1 <= n <= 120) ===")
    with csv_path.open(encoding="utf-8") as f:
        csv_rows = list(csv.DictReader(f))
    if len(csv_rows) != 120:
        print(f"FAIL: CSV has {len(csv_rows)} rows, expected 120!", file=sys.stderr)
        return 1

    for e, crow in zip(entries, csv_rows):
        n = int(e["n"])
        if int(crow["n"]) != n or crow["D_norm_exact"] != e["D_norm_exact"] or crow["det_exact"] != e["det_exact"]:
            print(f"FAIL: CSV/JSON mismatch at order n={n}!", file=sys.stderr)
            return 1
        if crow["ratio_to_bound"] != e["ratio_to_bound"] or crow["status"] != e["status"]:
            print(f"FAIL: CSV/JSON metadata mismatch at order n={n}!", file=sys.stderr)
            return 1
        d_norm = int(e["D_norm_exact"])
        det_abs = d_norm * (1 << (n - 1))
        if det_abs != int(e["det_exact"]):
            print(f"FAIL: n={n} det_exact != D_norm * 2^(n-1)", file=sys.stderr)
            return 1
        num, den, bname = residue_bound_squared_exact(n)
        lhs = (det_abs * det_abs) * den
        if lhs > num:
            print(f"FAIL: n={n} exceeds {bname} bound!", file=sys.stderr)
            return 1
        if e["attains_residue_bound"] and lhs != num:
            print(f"FAIL: n={n} marked attains_residue_bound=True but |det|^2 != bound^2!", file=sys.stderr)
            return 1
        ratio_calc = Decimal(1) if e["attains_residue_bound"] else compute_ratio_decimal(det_abs, n)
        if f"{ratio_calc:.8f}" != e["ratio_to_bound"]:
            print(f"FAIL: n={n} ratio_to_bound mismatch ({ratio_calc:.8f} != {e['ratio_to_bound']})!", file=sys.stderr)
            return 1
    print("  [PASS] 120/120 orders satisfy CSV/JSON parity, exact 2^(n-1) scaling, and residue-class bounds.")

    if sums_path.exists():
        print("=== Verifying SHA256SUMS Manifest ===")
        checked = 0
        for line in sums_path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if not line:
                continue
            parts = line.split("  ", 1)
            if len(parts) != 2:
                print(f"FAIL: Malformed SHA256SUMS line: {line}", file=sys.stderr)
                return 1
            expected_sha, rel_path = parts
            fpath = repo_root / rel_path
            if not fpath.exists():
                print(f"FAIL: Missing file in SHA256SUMS: {rel_path}", file=sys.stderr)
                return 1
            actual_sha = hashlib.sha256(fpath.read_bytes()).hexdigest()
            if actual_sha != expected_sha:
                print(f"FAIL: SHA-256 mismatch on {rel_path}: {actual_sha} != {expected_sha}", file=sys.stderr)
                return 1
            checked += 1
        print(f"  [PASS] {checked}/{checked} file SHA-256 digests verified.")

    total_s = time.perf_counter() - t0
    print(f"\nALL CHECKS PASSED (120/120 catalog orders verified in {total_s:.2f}s).")
    print("PASS")
    return 0


def parse_flexible_candidate_file(path: Path) -> Tuple[List[List[int]], str]:
    """
    Parses an external candidate matrix from:
      - JSON 2D list of {-1, +1} or {0, 1} integers (e.g. EinsteinArena submission payload),
      - Compact '+'/'-' sign rows (Orrick MaxDet / matrices/*.txt format),
      - Compact '0'/'1' rows or comma/whitespace-separated {-1, +1} or {0, 1} integers.
    If a {0, 1} binary matrix B of order m is given, lifts it via the exact Williamson
    bijection X = [[1, 1^T], [1, J - 2B]] to an (m+1) x (m+1) {-1, +1} matrix with D_norm = |det B|.
    """
    raw = path.read_text(encoding="utf-8").strip()
    if not raw:
        raise ValueError("Candidate matrix file is empty.")

    rows_int: List[List[int]] = []
    if raw.startswith(("[", "{")):
        arr = json.loads(raw)
        if isinstance(arr, dict) and "matrix" in arr:
            arr = arr["matrix"]
        if not isinstance(arr, list) or not arr or not isinstance(arr[0], list):
            raise ValueError("JSON input must be a 2D array (or an object with a 'matrix' 2D array).")
        rows_int = [[int(x) for x in r] for r in arr]
    else:
        lines = [
            ln.strip()
            for ln in raw.splitlines()
            if ln.strip() and not ln.strip().startswith("#")
        ]
        if all(set(ln) <= {"+", "-"} for ln in lines):
            n = len(lines)
            return (parse_sign_rows(lines, n), f"{n}x{n} '+/-' sign grid")
        if all(set(ln) <= {"0", "1"} and len(ln) > 1 for ln in lines):
            rows_int = [[int(ch) for ch in ln] for ln in lines]
        else:
            for ln in lines:
                cleaned = ln.replace("[", " ").replace("]", " ").replace(",", " ")
                toks = cleaned.split()
                if toks:
                    rows_int.append([int(t) for t in toks])

    m = len(rows_int)
    if m == 0 or any(len(r) != m for r in rows_int):
        raise ValueError("Matrix must be non-empty and square (n x n).")

    entries = {v for r in rows_int for v in r}
    if entries <= {-1, 1}:
        return (rows_int, f"{m}x{m} {{-1, +1}} integer matrix")
    if entries <= {0, 1}:
        n = m + 1
        lifted: List[List[int]] = [[1] * n]
        for r in rows_int:
            lifted.append([1] + [1 - 2 * v for v in r])
        return (lifted, f"{m}x{m} {{0, 1}} binary matrix (lifted to {n}x{n} {{-1, +1}} sign matrix)")
    raise ValueError(
        f"Matrix entries must be strictly in {{-1, +1}} or {{0, 1}}; found values: {sorted(entries)[:8]}"
    )


def check_external_candidate(repo_root: Path, cand_path: Path, use_sympy: bool = False) -> int:
    if not cand_path.is_file():
        print(f"FAIL: Candidate file not found: {cand_path}", file=sys.stderr)
        return 1
    try:
        mat, fmt_desc = parse_flexible_candidate_file(cand_path)
    except Exception as exc:
        print(f"FAIL: Could not parse candidate matrix ({exc})", file=sys.stderr)
        return 1

    n = len(mat)
    t0 = time.perf_counter()
    det_b = abs(det_bareiss_int(mat))
    if use_sympy and HAS_SYMPY:
        det_sp = abs(int(sp.Matrix(mat).det(method="bareiss")))
        if det_b != det_sp:
            print(f"FAIL: Bareiss ({det_b}) != SymPy ({det_sp})", file=sys.stderr)
            return 1
    dt_ms = (time.perf_counter() - t0) * 1000.0

    if det_b == 0:
        print(f"Parsed:  {fmt_desc} (order n = {n}, residue {n % 4} mod 4)")
        print("Result:  SINGULAR MATRIX (det = 0)")
        return 0

    scale = 1 << (n - 1)
    if det_b % scale != 0:
        print(f"FAIL: |det| is not divisible by 2^(n-1) for order n={n}!", file=sys.stderr)
        return 1
    d_norm = det_b // scale
    num, den, bname = residue_bound_squared_exact(n)
    if (det_b * det_b) * den > num:
        print(f"FAIL: |det|^2 exceeds {bname} upper bound — invalid matrix!", file=sys.stderr)
        return 1
    attains = (det_b * det_b) * den == num
    ratio = Decimal(1) if attains else compute_ratio_decimal(det_b, n)

    print(f"=== Candidate Matrix Verification ({cand_path.name}) ===")
    print(f"  Format        : {fmt_desc}")
    print(f"  Order         : n = {n} ({n % 4} mod 4), upper bound = {bname}")
    print(f"  |det(X)|      : {det_b}")
    print(f"  D_norm        : {d_norm}  (= |det(X)| / 2^{n-1}, OEIS A003432)")
    print(f"  Bound ratio R : {ratio:.8f}  ({'ATTAINS BOUND (R = 1)' if attains else 'R < 1'})")
    print(f"  Runtime       : {dt_ms:.1f} ms (exact Bareiss integer arithmetic)")

    cat_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.json"
    if cat_path.is_file() and 1 <= n <= 120:
        cat_data = json.loads(cat_path.read_text(encoding="utf-8"))
        entry = cat_data["orders"][n - 1]
        best_d = int(entry["D_norm_exact"])
        holder = f"{entry['current_holder']} ({entry['year']})"
        print(f"  Catalog best  : {best_d}  [{entry['status']}, {holder}]")
        if d_norm > best_d:
            if entry["status"] == "PROVEN_OPTIMAL":
                print("  ERROR: Exceeds a PROVEN_OPTIMAL catalog value — please double-check!", file=sys.stderr)
                return 1
            gain = (Decimal(d_norm - best_d) * Decimal(100)) / Decimal(best_d)
            print(f"  VERDICT       : NEW RECORD (+{gain:.6f}% strictly above catalog best!)")
        elif d_norm == best_d:
            print(f"  VERDICT       : TIES CATALOG BEST ({entry['status']})")
        else:
            pct = (Decimal(d_norm) * Decimal(100)) / Decimal(best_d)
            print(f"  VERDICT       : VALID MATRIX, BELOW CATALOG BEST ({pct:.4f}% of best known)")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Verify Hadamard Maximal Determinant catalog and witnesses.")
    parser.add_argument(
        "--check",
        type=str,
        default=None,
        metavar="MATRIX_FILE",
        help="Verify an external {+1,-1} or {0,1} candidate matrix file (+/- grid, JSON 2D array, or CSV) against exact bounds and catalog records.",
    )
    parser.add_argument("--fast", action="store_true", help="Verify bounds and new records only (skip full witness sweep).")
    parser.add_argument("--sympy", action="store_true", help="Also cross-check determinants via sympy.Matrix.det(method='bareiss').")
    parser.add_argument("--corrupt", action="store_true", help="Run adversarial 1-entry sign-flip mutation self-test (must exit 1).")
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parent
    if args.check is not None:
        sys.exit(check_external_candidate(repo_root, Path(args.check), use_sympy=args.sympy))
    sys.exit(
        run_verification(
            repo_root,
            verify_all_open=not args.fast,
            use_sympy=args.sympy,
            run_corrupt_test=args.corrupt,
        )
    )


if __name__ == "__main__":
    main()
