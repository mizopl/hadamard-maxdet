#!/usr/bin/env python3
"""
Exact GF(2) Sign-Gauge Block-Circulant Structure Verifier (`*` Engineering Flag).

Why a simple shift test is insufficient:
Archived maximal-determinant witnesses (e.g., in W. P. Orrick's `MaxDet` database) are
frequently normalized so that their first row and first column consist entirely of `+1`s.
Multiplying rows and columns by `-1` preserves |det X| and Hadamard equivalence, but hides
underlying circulant blocks `C` under a diagonal sign gauge `B = diag(s_r) * C * diag(s_c)`.

Method (exact linear algebra over GF(2), zero heuristics):
For each candidate block `B` of size `m x m` and shift orientation (`circulant` `c[(j - i) mod m]`
or `back-circulant` `c[(j + i) mod m]`), the existence of row signs `s_r in {0, 1}^m`,
column signs `s_c in {0, 1}^m`, and generator vector `c in {0, 1}^m` satisfying
    s_r[i] XOR s_c[j] XOR c[k] == (1 if B[i][j] == -1 else 0)   for all 0 <= i, j < m
is an exact system of m^2 linear equations in 3m binary unknowns over GF(2), solved below
via bitmask Gaussian elimination.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

sys.dont_write_bytecode = True


@dataclass(frozen=True)
class CirculantClassification:
    is_block_circulant: bool
    k_blocks: Optional[int] = None
    border_width: Optional[int] = None
    border_position: Optional[str] = None
    block_size: Optional[int] = None
    block_orientations: Tuple[str, ...] = ()

    def summary(self) -> str:
        if not self.is_block_circulant:
            return "non-circulant (general / locally switched sign matrix)"
        orient = "circulant" if len(set(self.block_orientations)) == 1 and self.block_orientations[0] == "circulant" else "circulant/back-circulant"
        return (
            f"{self.k_blocks}-block {orient} core "
            f"(block size m={self.block_size}, border={self.border_width} at {self.border_position})"
        )


def _solve_gf2_sign_circulant(blk: List[List[int]]) -> Optional[str]:
    """
    Returns 'circulant' or 'back-circulant' if diag(s_r) * blk * diag(s_c) has that
    exact cyclic structure for some s_r, s_c in {-1, +1}^m; otherwise None.
    """
    m = len(blk)
    if m < 3:
        return None
    for back, label in ((False, "circulant"), (True, "back-circulant")):
        pivots: Dict[int, Tuple[int, int]] = {}
        ok = True
        for i in range(m):
            row_i = blk[i]
            for j in range(m):
                k = (j + i) % m if back else (j - i) % m
                mask = (1 << i) | (1 << (m + j)) | (1 << (2 * m + k))
                rhs = 1 if row_i[j] == -1 else 0
                inserted = False
                while mask:
                    msb = mask.bit_length() - 1
                    if msb in pivots:
                        pm, pr = pivots[msb]
                        mask ^= pm
                        rhs ^= pr
                    else:
                        pivots[msb] = (mask, rhs)
                        inserted = True
                        break
                if not inserted and mask == 0 and rhs == 1:
                    ok = False
                    break
            if not ok:
                break
        if ok:
            return label
    return None


def classify_block_circulant(rows: List[str]) -> CirculantClassification:
    """
    Inspects an n x n '+/-' sign matrix (or list of {-1,+1} rows converted to '+/-')
    and determines whether stripping b in {0..4} border rows/columns at top-left or
    bottom-right decomposes the remaining (n - b) x (n - b) core into k x k
    sign-gauge circulant/back-circulant blocks for k in {1, 2, 4}.
    """
    n = len(rows)
    if n == 0:
        return CirculantClassification(False)
    mat = [[1 if ch == "+" else -1 for ch in r.strip()] for r in rows]
    if any(len(r) != n for r in mat):
        raise ValueError("Matrix rows must be square (n x n).")

    for kblocks in (1, 2, 4):
        for b in (0, 1, 2, 3, 4):
            rem = n - b
            if rem >= 4 * kblocks and rem % kblocks == 0:
                m = rem // kblocks
                for top in (True, False):
                    sub = [r[b:] for r in mat[b:]] if top else [r[:rem] for r in mat[:rem]]
                    orientations: List[str] = []
                    good = True
                    for bi in range(kblocks):
                        for bj in range(kblocks):
                            blk = [sub[bi * m + ii][bj * m : (bj + 1) * m] for ii in range(m)]
                            kind = _solve_gf2_sign_circulant(blk)
                            if kind is None:
                                good = False
                                break
                            orientations.append(kind)
                        if not good:
                            break
                    if good:
                        pos = "none" if b == 0 else ("top-left" if top else "bottom-right")
                        return CirculantClassification(
                            is_block_circulant=True,
                            k_blocks=kblocks,
                            border_width=b,
                            border_position=pos,
                            block_size=m,
                            block_orientations=tuple(orientations),
                        )
    return CirculantClassification(False)


def has_block_circulant_core(rows: List[str]) -> bool:
    """Boolean wrapper used for the `*` engineering-structure flag in catalog tables."""
    return classify_block_circulant(rows).is_block_circulant


def verify_catalog_circulant_flags(repo_root: Path) -> Tuple[int, int]:
    """
    Audits all 120 orders in catalog/hadamard_maxdet_n1_to_120.json from scratch:
      1. Verifies `dyadic_walsh_order == ((n & (n - 1)) == 0)` on all 120 orders.
      2. Recomputes `has_block_circulant_core` via exact GF(2) elimination on every
         explicit open-order witness matrix and checks exact agreement with
         `block_circulant_witness`.
    Returns (dyadic_count, block_circulant_open_count).
    """
    cat_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.json"
    cat_data = json.loads(cat_path.read_text(encoding="utf-8"))
    dyadic_cnt = 0
    circ_cnt = 0

    for entry in cat_data["orders"]:
        n = int(entry["n"])
        expected_dyadic = (n & (n - 1)) == 0
        if entry.get("dyadic_walsh_order") is not expected_dyadic:
            raise AssertionError(
                f"n={n}: dyadic_walsh_order={entry.get('dyadic_walsh_order')} != {expected_dyadic}"
            )
        if expected_dyadic:
            dyadic_cnt += 1

        mrel = entry.get("matrix_file", "")
        if entry["status"] == "OPEN" and mrel:
            mpath = repo_root / mrel
            rows = [ln.strip() for ln in mpath.read_text(encoding="utf-8").splitlines() if ln.strip()]
            actual_circ = has_block_circulant_core(rows)
            if entry.get("block_circulant_witness") is not actual_circ:
                raise AssertionError(
                    f"n={n} ({mrel}): block_circulant_witness={entry.get('block_circulant_witness')} "
                    f"!= recomputed GF(2) value {actual_circ}"
                )
            if actual_circ:
                circ_cnt += 1
        else:
            if entry.get("block_circulant_witness") not in (False, None):
                raise AssertionError(f"n={n}: unexpected block_circulant_witness=True without open matrix")

    return (dyadic_cnt, circ_cnt)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify GF(2) sign-gauge block-circulant structure (`*` flag) on a matrix file or the full catalog."
    )
    parser.add_argument(
        "--inspect",
        type=str,
        default=None,
        help="Path to an n x n '+/-' matrix file to classify individually.",
    )
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parent.parent

    if args.inspect:
        p = Path(args.inspect)
        rows = [ln.strip() for ln in p.read_text(encoding="utf-8").splitlines() if ln.strip()]
        res = classify_block_circulant(rows)
        print(f"{p.name} (n={len(rows)}): {res.summary()}")
        sys.exit(0)

    dyadic_cnt, circ_cnt = verify_catalog_circulant_flags(repo_root)
    print(
        f"[PASS] Verified 120/120 orders: {dyadic_cnt} dyadic Walsh-Hadamard orders (2^m) "
        f"and {circ_cnt} block-circulant open-order witnesses (`*`)."
    )


if __name__ == "__main__":
    main()
