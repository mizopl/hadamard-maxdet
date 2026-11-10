"""
Exact-integer verification and structural classification suite for the
Hadamard Maximal Determinant Catalog (1 <= n <= 120, OEIS A003432 / A003433).

Modules:
- determinant_and_bounds: Fraction-free Bareiss integer determinant, exact 2^(n-1)
  scaling, and exact rational upper bounds across all four residue classes mod 4.
- circulant_structure: Exact GF(2) sign-gauge solver certifying 1-, 2-, and 4-block
  circulant / back-circulant cores (`*` engineering-friendly constructions).
- candidate_audit: Multi-format external matrix parser (+/- grids, JSON arrays, CSV,
  and OEIS A003432 {0, 1} binary matrices) plus Hadamard-equivalence Gram invariants.
- admin_update_catalog: Co-maintainer CLI for ingesting verified new records and
  refreshing catalog manifests and SHA256SUMS.
"""

import sys

sys.dont_write_bytecode = True

from .determinant_and_bounds import (
    HAS_SYMPY,
    compute_ratio_decimal,
    det_bareiss_int,
    ehlich_block_count,
    parse_sign_rows,
    residue_bound_squared_exact,
    verify_matrix_rows_against_d_norm,
)
from .circulant_structure import (
    CirculantClassification,
    classify_block_circulant,
    has_block_circulant_core,
    verify_catalog_circulant_flags,
)
from .candidate_audit import (
    check_external_candidate,
    hadamard_gram_fingerprint,
    parse_flexible_candidate_file,
)

__all__ = [
    "HAS_SYMPY",
    "CirculantClassification",
    "check_external_candidate",
    "classify_block_circulant",
    "compute_ratio_decimal",
    "det_bareiss_int",
    "ehlich_block_count",
    "hadamard_gram_fingerprint",
    "has_block_circulant_core",
    "parse_flexible_candidate_file",
    "parse_sign_rows",
    "residue_bound_squared_exact",
    "verify_catalog_circulant_flags",
    "verify_matrix_rows_against_d_norm",
]
