#!/usr/bin/env python3
"""
Exact Arbitrary-Precision Integer Determinant & Residue-Class Bound Verifier.

Zero floating-point determinant calls are used anywhere in this module:
- Determinants are computed over Python `int` via fraction-free Bareiss elimination
  (`det_bareiss_int`), with optional independent cross-verification via SymPy (`--sympy`).
- Classical residue-class upper bounds B_class(n)^2 (Hadamard 1893, Barba 1933,
  Ehlich-Wojtas 1964, and Ehlich's piecewise s-block Gram bound 1964) are represented
  as exact reduced rational squares (numerator, denominator) so |det X|^2 <= B_class(n)^2
  is checked by pure integer multiplication.
"""

from __future__ import annotations

import argparse
import math
import sys
from decimal import Decimal, getcontext
from pathlib import Path
from typing import List, Tuple

sys.dont_write_bytecode = True
getcontext().prec = 140

try:
    import sympy as sp  # type: ignore
    HAS_SYMPY = True
except Exception:
    HAS_SYMPY = False


def parse_sign_rows(rows: List[str], expected_n: int) -> List[List[int]]:
    """Parses n compact '+/-' strings into an n x n matrix of {-1, +1} Python ints."""
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
    if n % 4 != 3 or n < 3:
        raise ValueError(f"Ehlich block count requires n >= 3 with n == 3 (mod 4), got {n}")
    if n == 3:
        return 3
    if n in (7, 11):
        return 5
    if 15 <= n <= 59:
        return 6
    return 7


def residue_bound_squared_exact(n: int) -> Tuple[int, int, str]:
    """
    Returns (num, den, bound_name) representing the exact reduced rational value
    num / den = B_class(n)^2 for n x n {-1, +1} matrices:
      - n in {1, 2} or n == 0 mod 4 (Hadamard 1893): B(n)^2 = n^n
      - n == 1 mod 4, n >= 5 (Barba 1933):           B(n)^2 = (2n - 1) * (n - 1)^(n - 1)
      - n == 2 mod 4, n >= 6 (Ehlich-Wojtas 1964):   B(n)^2 = 4 * (n - 1)^2 * (n - 2)^(n - 2)
      - n == 3 mod 4, n >= 3 (Ehlich 1964):
          B(n)^2 = (n - 3)^(n - s) * (n - 3 + 4r)^u * (n + 1 + 4r)^v
                   * (1 - u*r/(n - 3 + 4r) - v*(r + 1)/(n + 1 + 4r))
    """
    if n < 1:
        raise ValueError(f"Order n must be >= 1, got {n}")
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
    """High-precision (140-digit Decimal) ratio R(n) = |det X| / B_class(n)."""
    num, den, _ = residue_bound_squared_exact(n)
    b_sq = Decimal(num) / Decimal(den)
    return Decimal(det_abs) / b_sq.sqrt()


def verify_matrix_rows_against_d_norm(
    n: int, rows: List[str], expected_d_norm: int, check_sympy: bool = False
) -> int:
    """Verifies shape, {-1,+1} entries, exact Bareiss D_norm, and residue-class upper bound."""
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


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute exact Bareiss determinant & bound ratio of a +/- sign matrix.")
    parser.add_argument("matrix_file", type=str, help="Path to n x n '+/-' sign matrix file.")
    parser.add_argument("--sympy", action="store_true", help="Cross-check Bareiss determinant with SymPy.")
    args = parser.parse_args()

    p = Path(args.matrix_file)
    rows = [ln.strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
    n = len(rows)
    mat = parse_sign_rows(rows, n)
    det_b = abs(det_bareiss_int(mat))
    if args.sympy and HAS_SYMPY:
        det_sp = abs(int(sp.Matrix(mat).det(method="bareiss")))
        assert det_b == det_sp, f"Bareiss ({det_b}) != SymPy ({det_sp})"
    d_norm = det_b // (1 << (n - 1))
    num, den, bname = residue_bound_squared_exact(n)
    assert (det_b * det_b) * den <= num, f"Exceeds {bname} bound!"
    ratio = compute_ratio_decimal(det_b, n)
    print(f"n = {n} ({n % 4} mod 4) | bound = {bname} | D_norm = {d_norm} | R = {ratio:.8f}")


if __name__ == "__main__":
    main()
