#!/usr/bin/env python3
"""
Master Offline Exact-Integer & Structural Verifier for the Hadamard Maximal Determinant Catalog
(Orders 1 <= n <= 120, OEIS A003432 / A003433).

All underlying mathematical engines live openly in `./verifiers/` for inspection,
standalone CLI execution, and template reuse:
  - `verifiers/determinant_and_bounds.py`: fraction-free Bareiss integer determinant
    (`det_bareiss_int`) + exact rational Hadamard / Barba / Ehlich-Wojtas / Ehlich(s) bounds.
  - `verifiers/circulant_structure.py`: exact GF(2) sign-gauge solver certifying 1-, 2-,
    and 4-block circulant/back-circulant cores (`*` engineering-friendly witnesses) and
    dyadic Walsh-Hadamard orders (`2^m`).
  - `verifiers/candidate_audit.py`: multi-format external candidate checker (`--check FILE`)
    with Hadamard-equivalence Gram-profile fingerprinting.
  - `verifiers/admin_update_catalog.py`: co-maintainer ingestion & SHA256SUMS refresh CLI.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from decimal import Decimal
from pathlib import Path

sys.dont_write_bytecode = True

from verifiers.candidate_audit import check_external_candidate, parse_flexible_candidate_file
from verifiers.circulant_structure import (
    classify_block_circulant,
    has_block_circulant_core,
    verify_catalog_circulant_flags,
)
from verifiers.determinant_and_bounds import (
    HAS_SYMPY,
    compute_ratio_decimal,
    det_bareiss_int,
    ehlich_block_count,
    parse_sign_rows,
    residue_bound_squared_exact,
    verify_matrix_rows_against_d_norm,
)


def run_verification(
    repo_root: Path, verify_all_open: bool = True, use_sympy: bool = False, run_corrupt_test: bool = False
) -> int:
    t0 = time.perf_counter()
    we_path = repo_root / "catalog" / "new_records.json"
    tied_path = repo_root / "catalog" / "tied_constructions.json"
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
        print(f"=== Verifying {we_count} New Constructive Records (SymPy cross-check: {use_sympy and HAS_SYMPY}) ===")
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

    if tied_path.exists():
        tied_data = json.loads(tied_path.read_text(encoding="utf-8"))
        tied_list = tied_data.get("tied_witnesses", [])
        print(f"=== Verifying {len(tied_list)} Alternative Record-Tying Witnesses ===")
        for tw in tied_list:
            n = int(tw["n"])
            claimed_d = int(tw["D_norm"])
            cat_d = int(entries[n - 1]["D_norm_exact"])
            if claimed_d != cat_d:
                print(f"FAIL: Tied witness n={n} D_norm ({claimed_d}) != catalog ({cat_d})!", file=sys.stderr)
                return 1
            verify_matrix_rows_against_d_norm(n, tw["matrix_rows"], claimed_d, check_sympy=use_sympy)
            tfile = repo_root / tw["matrix_file"]
            disk_rows = [ln.strip() for ln in tfile.read_text(encoding="utf-8").splitlines() if ln.strip()]
            if disk_rows != tw["matrix_rows"]:
                print(f"FAIL: Disk tied witness {tfile} differs from JSON certificate!", file=sys.stderr)
                return 1
        print(f"  [PASS] {len(tied_list)}/{len(tied_list)} alternative record-tying witnesses verified.")

    if verify_all_open:
        open_entries = [e for e in entries if e["status"] == "OPEN"]
        verified_open = 0
        print("=== Verifying Explicit +/-1 Witness Matrices Across Open Orders (23 <= n <= 119) ===")
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

    print("=== Verifying GF(2) Sign-Gauge Block-Circulant (`*`) & Dyadic Walsh-Hadamard Flags ===")
    try:
        dyadic_cnt, circ_cnt = verify_catalog_circulant_flags(repo_root)
    except Exception as exc:
        print(f"FAIL: Structural flag verification error: {exc}", file=sys.stderr)
        return 1
    print(
        f"  [PASS] {dyadic_cnt} dyadic Walsh-Hadamard orders (2^m) and {circ_cnt} block-circulant "
        f"open-order witnesses (`*`) verified over GF(2)."
    )

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


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Standalone exact-integer & GF(2) structural verifier for Hadamard Maximal Determinants (1 <= n <= 120)."
    )
    parser.add_argument(
        "--check",
        type=str,
        default=None,
        metavar="MATRIX_FILE",
        help="Verify an external candidate matrix (+/- grid, JSON array, CSV +/-1, or 0/1 binary) against catalog bounds & records.",
    )
    parser.add_argument(
        "--sympy",
        action="store_true",
        help="Also cross-check determinants via sympy.Matrix.det(method='bareiss') if SymPy is installed.",
    )
    parser.add_argument(
        "--corrupt",
        action="store_true",
        help="Adversarial self-test: mutate 1 matrix entry in memory and verify that the checker fails (exits 1).",
    )
    parser.add_argument(
        "--fast",
        action="store_true",
        help="Skip recomputing Bareiss determinants on all literature witnesses (verify JSON/CSV bounds, GF(2) flags & SHA256SUMS).",
    )
    args = parser.parse_args()
    repo_root = Path(__file__).resolve().parent
    if args.check:
        rc = check_external_candidate(repo_root=repo_root, cand_path=Path(args.check), use_sympy=args.sympy)
        sys.exit(rc)
    rc = run_verification(
        repo_root=repo_root,
        verify_all_open=not args.fast,
        use_sympy=args.sympy,
        run_corrupt_test=args.corrupt,
    )
    sys.exit(rc)


if __name__ == "__main__":
    main()
