#!/usr/bin/env python3
"""
External Candidate Matrix Auditor & Hadamard-Equivalence Invariant Fingerprinter.

Supports four community formats:
  1. Compact '+/-' sign grids (`matrices/*.txt` & Orrick `MaxDet` format),
  2. JSON 2D arrays `[[1, -1, ...], ...]` or `{"matrix": [[...]]}` (EinsteinArena payloads),
  3. CSV or whitespace-separated `{-1, +1}` integer matrices (`numpy.savetxt`),
  4. Binary `{0, 1}` matrices of order m = n - 1 (OEIS A003432 convention, lifted via
     Williamson's bijection X = [[1, 1^T], [1, J - 2B]] with D_norm(m+1) = |det B|).

Also computes an exact Hadamard-equivalence invariant (`hadamard_gram_fingerprint`,
invariant under all row/column permutations and +/-1 row/column sign multiplications)
and runs `classify_block_circulant` so contributors and co-maintainers see immediately
whether a candidate has a block-circulant core (`*`) or is provably non-isomorphic to
an existing record witness.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import sys
import time
from decimal import Decimal
from pathlib import Path
from typing import List, Tuple

sys.dont_write_bytecode = True

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from verifiers.circulant_structure import classify_block_circulant
from verifiers.determinant_and_bounds import (
    HAS_SYMPY,
    compute_ratio_decimal,
    det_bareiss_int,
    parse_sign_rows,
    residue_bound_squared_exact,
)

try:
    import sympy as sp  # type: ignore
except Exception:
    sp = None  # type: ignore


def matrix_to_sign_rows(mat: List[List[int]]) -> List[str]:
    return ["".join("+" if v == 1 else "-" for v in row) for row in mat]


def hadamard_gram_fingerprint(mat: List[List[int]]) -> str:
    """
    Computes a fast, deterministic 16-hex-char SHA-256 fingerprint that is strictly
    invariant under Hadamard equivalence (arbitrary row/column permutations and
    arbitrary +/-1 row/column sign flips):
      - sorted per-row sorted tuples of absolute row-Gram inner products |(X X^T)_{i,j}| (j != i),
      - sorted per-col sorted tuples of absolute col-Gram inner products |(X^T X)_{i,j}| (j != i).
    Distinct fingerprints rigorously prove two matrices are Hadamard-inequivalent (up to
    transposition if row/col multisets are swapped; we sort {row_sig, col_sig} to be
    transpose-invariant as well).
    """
    n = len(mat)
    row_profiles = []
    for i in range(n):
        ri = mat[i]
        degs = sorted(
            abs(sum(ri[k] * mat[j][k] for k in range(n)))
            for j in range(n)
            if j != i
        )
        row_profiles.append(tuple(degs))
    row_profiles.sort()

    col_profiles = []
    for c1 in range(n):
        degs = sorted(
            abs(sum(mat[r][c1] * mat[r][c2] for r in range(n)))
            for c2 in range(n)
            if c2 != c1
        )
        col_profiles.append(tuple(degs))
    col_profiles.sort()

    pair = sorted([repr(row_profiles), repr(col_profiles)])
    payload = ("||".join(pair)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()[:16]


def parse_flexible_candidate_file(path: Path) -> Tuple[List[List[int]], str]:
    """Parses '+/-', JSON 2D array/dict, CSV/whitespace +/-1, or {0,1} binary matrix."""
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
    if use_sympy and HAS_SYMPY and sp is not None:
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

    sign_rows = matrix_to_sign_rows(mat)
    circ_info = classify_block_circulant(sign_rows)
    gram_fp = hadamard_gram_fingerprint(mat)

    print(f"=== Candidate Matrix Verification ({cand_path.name}) ===")
    print(f"  Format        : {fmt_desc}")
    print(f"  Order         : n = {n} ({n % 4} mod 4), upper bound = {bname}")
    print(f"  |det(X)|      : {det_b}")
    print(f"  D_norm        : {d_norm}  (= |det(X)| / 2^{n-1}, OEIS A003432)")
    print(f"  Bound ratio R : {ratio:.8f}  ({'ATTAINS BOUND (R = 1)' if attains else 'R < 1'})")
    print(f"  Structure (*) : {circ_info.summary()}")
    print(f"  Gram invariant: {gram_fp} (Hadamard- & transpose-invariant fingerprint)")
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
            eq_note = ""
            mrel = entry.get("matrix_file", "")
            if mrel and (repo_root / mrel).is_file():
                cat_rows = [
                    ln.strip()
                    for ln in (repo_root / mrel).read_text(encoding="utf-8").splitlines()
                    if ln.strip()
                ]
                cat_fp = hadamard_gram_fingerprint(parse_sign_rows(cat_rows, n))
                if gram_fp != cat_fp:
                    eq_note = f" — PROVABLY NON-ISOMORPHIC to catalog witness (Gram {gram_fp} != {cat_fp})"
                else:
                    eq_note = f" — matches catalog Gram profile ({cat_fp})"
            print(f"  VERDICT       : TIES CATALOG BEST ({entry['status']}){eq_note}")
        else:
            pct = (Decimal(d_norm) * Decimal(100)) / Decimal(best_d)
            print(f"  VERDICT       : VALID MATRIX, BELOW CATALOG BEST ({pct:.4f}% of best known)")
    return 0


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit an external +/-1 or 0/1 candidate matrix against the catalog.")
    parser.add_argument("matrix_file", type=str, help="Path to candidate matrix file.")
    parser.add_argument("--sympy", action="store_true", help="Cross-check Bareiss determinant with SymPy.")
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parent.parent
    sys.exit(check_external_candidate(repo_root, Path(args.matrix_file), use_sympy=args.sympy))


if __name__ == "__main__":
    main()
