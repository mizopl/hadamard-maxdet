#!/usr/bin/env python3
"""
Maintainer & Co-Admin Catalog Update Utility (`verifiers/admin_update_catalog.py`).

Designed so that any co-maintainer of `hadamard-maxdet` (or anyone adapting this
repository as a template for an exact combinatorial record table) can ingest a newly
verified witness matrix or refresh cryptographic manifests without external tooling:

  # 1. Refresh SHA256SUMS across all tracked repository files:
  python3 verifiers/admin_update_catalog.py --refresh-hashes

  # 2. Ingest a new strictly-improving open-order record into JSON, CSV & matrices/:
  python3 verifiers/admin_update_catalog.py \
      --ingest path/to/candidate.txt \
      --holder "J. Doe & A. Smith" \
      --year 2026 \
      --ref-key DS2026 \
      --citation "J. Doe and A. Smith, Title, Journal / arXiv (2026)."
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import hashlib
import json
import subprocess
import sys
from decimal import Decimal
from pathlib import Path
from typing import List, Tuple

sys.dont_write_bytecode = True

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from verifiers.candidate_audit import matrix_to_sign_rows, parse_flexible_candidate_file
from verifiers.circulant_structure import has_block_circulant_core
from verifiers.determinant_and_bounds import (
    HAS_SYMPY,
    compute_ratio_decimal,
    det_bareiss_int,
    residue_bound_squared_exact,
)

try:
    import sympy as sp  # type: ignore
except Exception:
    sp = None  # type: ignore


CSV_COLUMNS = [
    "n",
    "residue_mod_4",
    "status",
    "attains_residue_bound",
    "bound_name",
    "display_D_norm",
    "scientific_D_norm",
    "ratio_to_bound",
    "current_holder",
    "year",
    "ref_key",
    "previous_holder",
    "gain_percent",
    "matrix_file",
    "D_norm_exact",
    "det_exact",
]


def format_sci_exact(val: int, decimals: int = 4) -> str:
    s = str(val)
    if len(s) <= 8:
        return s
    exp = len(s) - 1
    mant = Decimal(val) / (Decimal(10) ** exp)
    if round(mant, decimals) >= Decimal(10):
        exp += 1
        mant = Decimal(1)
    return f"{mant:.{decimals}f}e+{exp}"


def format_factored_d_norm(n: int, d_norm: int) -> Tuple[str, str]:
    m = n // 4
    if m <= 1 or d_norm <= 10**9:
        return (str(d_norm), str(d_norm))
    max_exp = 2 * m - 1
    e = 0
    rem = d_norm
    while e < max_exp and rem % m == 0:
        rem //= m
        e += 1
    if e >= 6 and len(str(rem)) <= 14:
        fac = f"{rem} * {m}^{e}"
        return (fac, fac)
    sci = format_sci_exact(d_norm, 4)
    if e >= 1:
        return (f"{rem} * {m}^{e}", sci)
    return (str(d_norm), sci)


def format_sci_html(val: int, decimals: int = 4) -> str:
    if val < 1_000_000:
        return f"{val:,}"
    s = str(val)
    exp = len(s) - 1
    mant = Decimal(val) / (Decimal(10) ** exp)
    if round(mant, decimals) >= Decimal(10):
        exp += 1
        mant = Decimal(1)
    return f"{mant:.{decimals}f} × 10<sup>{exp}</sup>"


def format_factored_html(n: int, d_norm: int) -> str:
    if d_norm < 1_000_000_000:
        return f"{d_norm:,}"
    m = n // 4
    if m >= 2:
        max_exp = 2 * m - 1
        e = 0
        rem = d_norm
        while e < max_exp and rem % m == 0:
            rem //= m
            e += 1
        if e >= 6 and len(str(rem)) <= 8:
            rem_s = f"{rem:,}" if rem >= 10_000 else str(rem)
            return f"{rem_s} · {m}<sup>{e}</sup>"
    return "—"


def update_readme_for_entry(
    repo_root: Path, entry: dict, gain_pct: float, ref_key: str, citation: str
) -> None:
    readme_path = repo_root / "README.md"
    if not readme_path.is_file():
        return
    n = int(entry["n"])
    d_int = int(entry["D_norm_exact"])
    sci_cell = format_sci_html(d_int, 4)
    fac_cell = format_factored_html(n, d_int)
    r_cell = f"{float(entry['ratio_to_bound']):.4f}"
    gp_disp = f"+{gain_pct:.5f}%" if gain_pct < 0.01 else f"+{gain_pct:.2f}%"
    circ_star = "*" if entry.get("block_circulant_witness") else ""
    m_rel = entry["matrix_file"]
    st = "OPT" if entry["status"] == "PROVEN_OPTIMAL" else "OPEN"
    n_master = f"{n}*" if entry.get("dyadic_walsh_order") else f"{n}"

    open_row_new = (
        f"| **{n}** | **{sci_cell}** | {f'**{fac_cell}**' if fac_cell != '—' else '—'} "
        f"| **{r_cell}** | **{entry['current_holder']}** ({gp_disp}) | **{entry['year']}** "
        f"| [n{n:03d}.txt]({m_rel}){circ_star} |"
    )
    master_row_new = (
        f"| {n_master} | {entry['residue_mod_4']} | {st} | {sci_cell} | {fac_cell} "
        f"| {r_cell} | {entry['current_holder']} | {entry['year']} | {ref_key} "
        f"| [txt]({m_rel}){circ_star} |"
    )

    lines = readme_path.read_text(encoding="utf-8").splitlines()
    out_lines: List[str] = []
    open_prefix = f"| **{n}** |"
    master_prefix = f"| {n_master} | {entry['residue_mod_4']} |"
    for ln in lines:
        if ln.startswith(open_prefix):
            if entry["status"] == "OPEN":
                out_lines.append(open_row_new)
            # If newly PROVEN_OPTIMAL, omit from the open-orders subtable
        elif ln.startswith(master_prefix):
            out_lines.append(master_row_new)
        else:
            out_lines.append(ln)

    text = "\n".join(out_lines) + "\n"
    if citation and f"**[{ref_key}]**" not in text:
        bib_marker = "</details>\n\n---\n\n<sub>Curated"
        if bib_marker in text:
            text = text.replace(bib_marker, f"- **[{ref_key}]** {citation}\n{bib_marker}", 1)
    readme_path.write_text(text, encoding="utf-8")


def collect_tracked_files(repo_root: Path) -> List[Path]:
    root_docs = [
        "README.md",
        "AGENTS.md",
        "CONTRIBUTING.md",
        "llms.txt",
        "CITATION.cff",
        "LICENSE",
        "verify.py",
    ]
    tracked: List[Path] = [repo_root / name for name in root_docs if (repo_root / name).is_file()]
    ver_dir = repo_root / "verifiers"
    if ver_dir.is_dir():
        tracked.extend(sorted(ver_dir.glob("*.py")))
    assets_dir = repo_root / "assets"
    if assets_dir.is_dir():
        tracked.extend(sorted(p for p in assets_dir.iterdir() if p.is_file() and not p.name.startswith(".")))
    cat_dir = repo_root / "catalog"
    for cname in (
        "hadamard_maxdet_n1_to_120.json",
        "hadamard_maxdet_n1_to_120.csv",
        "new_records.json",
        "tied_constructions.json",
    ):
        cp = cat_dir / cname
        if cp.is_file():
            tracked.append(cp)
    mat_root = repo_root / "matrices"
    if mat_root.is_dir():
        for sub in sorted(mat_root.iterdir()):
            if sub.is_dir():
                tracked.extend(sorted(sub.glob("*.txt")))
    return tracked


def refresh_sha256sums(repo_root: Path) -> int:
    tracked = collect_tracked_files(repo_root)
    lines: List[str] = []
    for fp in tracked:
        digest = hashlib.sha256(fp.read_bytes()).hexdigest()
        rel = fp.relative_to(repo_root).as_posix()
        lines.append(f"{digest}  {rel}")
    (repo_root / "SHA256SUMS").write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def ingest_new_record(
    repo_root: Path,
    cand_path: Path,
    holder: str,
    year: int,
    ref_key: str,
    citation: str,
    subdir: str = "community_records",
) -> None:
    mat, _ = parse_flexible_candidate_file(cand_path)
    n = len(mat)
    if not (1 <= n <= 120):
        raise ValueError(f"Order n={n} is outside catalog range 1..120.")

    det_b = abs(det_bareiss_int(mat))
    if HAS_SYMPY and sp is not None:
        det_sp = abs(int(sp.Matrix(mat).det(method="bareiss")))
        if det_b != det_sp:
            raise AssertionError(f"Bareiss ({det_b}) != SymPy ({det_sp})")

    scale = 1 << (n - 1)
    if det_b == 0 or det_b % scale != 0:
        raise ValueError(f"Invalid determinant scaling on order n={n}.")
    d_norm = det_b // scale

    num, den, bname = residue_bound_squared_exact(n)
    if (det_b * det_b) * den > num:
        raise ValueError(f"Candidate exceeds {bname} upper bound!")

    cat_json_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.json"
    cat_csv_path = repo_root / "catalog" / "hadamard_maxdet_n1_to_120.csv"
    cat_data = json.loads(cat_json_path.read_text(encoding="utf-8"))
    entry = cat_data["orders"][n - 1]
    prev_d = int(entry["D_norm_exact"])

    if entry["status"] == "PROVEN_OPTIMAL":
        raise ValueError(f"Order n={n} is already PROVEN_OPTIMAL ({prev_d}).")
    if d_norm <= prev_d:
        raise ValueError(
            f"Candidate D_norm ({d_norm}) does not strictly exceed catalog record ({prev_d}) on n={n}."
        )

    sign_rows = matrix_to_sign_rows(mat)
    out_dir = repo_root / "matrices" / subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    m_rel = f"matrices/{subdir}/n{n:03d}.txt"
    (repo_root / m_rel).write_text("\n".join(sign_rows) + "\n", encoding="utf-8")

    attains = (det_b * det_b) * den == num
    ratio = Decimal(1) if attains else compute_ratio_decimal(det_b, n)
    fac_exact, fac_disp = format_factored_d_norm(n, d_norm)
    gain_pct = float((Decimal(d_norm - prev_d) * Decimal(100)) / Decimal(prev_d))
    gain_str = f"+{gain_pct:.6f}%" if gain_pct < 0.01 else f"+{gain_pct:.4f}%"

    prev_holder_str = f"{entry['current_holder']} ({entry['year']})"
    entry.update(
        {
            "status": "PROVEN_OPTIMAL" if attains else "OPEN",
            "attains_residue_bound": attains,
            "D_norm_exact": str(d_norm),
            "det_exact": str(det_b),
            "factored_D_norm": fac_exact,
            "display_D_norm": fac_disp,
            "scientific_D_norm": format_sci_exact(d_norm, 4),
            "ratio_to_bound": f"{ratio:.8f}",
            "current_holder": holder,
            "year": year,
            "ref_key": ref_key,
            "previous_D_norm": str(prev_d),
            "previous_holder": prev_holder_str,
            "gain_percent": gain_str,
            "matrix_file": m_rel,
            "block_circulant_witness": has_block_circulant_core(sign_rows),
        }
    )
    entry.pop("missing_archived_witness", None)

    if citation:
        cat_data["references"][ref_key] = citation
    cat_data["generated_utc"] = dt.datetime.now(dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    cat_data["summary"]["proven_optimal_orders"] = sum(
        1 for e in cat_data["orders"] if e["status"] == "PROVEN_OPTIMAL"
    )
    cat_data["summary"]["open_orders"] = sum(
        1 for e in cat_data["orders"] if e["status"] == "OPEN"
    )
    cat_data["summary"]["open_orders_with_explicit_matrix"] = sum(
        1 for e in cat_data["orders"] if e["status"] == "OPEN" and bool(e.get("matrix_file"))
    )

    cat_json_path.write_text(json.dumps(cat_data, indent=2) + "\n", encoding="utf-8")
    with cat_csv_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=CSV_COLUMNS, extrasaction="ignore")
        writer.writeheader()
        for row in cat_data["orders"]:
            writer.writerow(row)

    update_readme_for_entry(repo_root, entry, gain_pct, ref_key, citation)
    hashed = refresh_sha256sums(repo_root)
    print(
        f"Ingested n={n} ({gain_str} over {prev_holder_str}) -> {m_rel}; "
        f"updated JSON, CSV, README.md & {hashed} SHA-256 digests."
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Co-maintainer catalog ingestion & SHA256SUMS refresh tool.")
    parser.add_argument("--refresh-hashes", action="store_true", help="Recompute SHA256SUMS across all tracked files.")
    parser.add_argument("--ingest", type=str, default=None, help="Path to new record matrix file to ingest.")
    parser.add_argument("--holder", type=str, default="", help="Author / record-holder attribution string.")
    parser.add_argument("--year", type=int, default=2026, help="Year of construction.")
    parser.add_argument("--ref-key", type=str, default="", help="Short reference key (e.g. Author2026).")
    parser.add_argument("--citation", type=str, default="", help="Full bibliographic citation for references dict.")
    args = parser.parse_args()

    repo_root = Path(__file__).resolve().parent.parent
    if args.ingest:
        if not args.holder or not args.ref_key:
            parser.error("--ingest requires --holder and --ref-key.")
        ingest_new_record(
            repo_root=repo_root,
            cand_path=Path(args.ingest),
            holder=args.holder,
            year=args.year,
            ref_key=args.ref_key,
            citation=args.citation,
        )
    elif args.refresh_hashes:
        cnt = refresh_sha256sums(repo_root)
        print(f"Refreshed SHA256SUMS ({cnt} files).")
    else:
        parser.print_help()
        sys.exit(0)

    # Run full verification suite immediately after any admin update
    rc = subprocess.call([sys.executable, str(repo_root / "verify.py")])
    sys.exit(rc)


if __name__ == "__main__":
    main()
