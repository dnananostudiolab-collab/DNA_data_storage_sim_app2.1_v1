"""Diagnostics for the v4 positional rule-based DNA codec.

Run:
    python test_rulebase_positional_v4.py

This checks:
- clean round-trip for R0/R1/R2/Rinf;
- R0 HomoMax=1 across block boundaries;
- deterministic positional permutation metadata;
- sequence-distribution diagnostics on highly repetitive input.
"""
from __future__ import annotations

from collections import Counter
from math import log2
import os

import dna_codec
from dna_mapping import encode_bytes_to_dna, decode_dna_with_design
from dna_design_rules import DNA_DESIGN_PERMUTATION_MODE

RULES = ("R0_B9", "R1_B12", "R2_B15", "RINF_B16")


def entropy(symbols) -> float:
    symbols = list(symbols)
    if not symbols:
        return 0.0
    counts = Counter(symbols)
    n = len(symbols)
    return -sum((count / n) * log2(count / n) for count in counts.values())


def main() -> None:
    # Clean round-trip across several payload lengths.
    for size in (1, 2, 3, 4, 7, 16, 31, 100):
        source = os.urandom(size)
        for rule in RULES:
            dna, _, meta = encode_bytes_to_dna(source, rule)
            decoded, _, decoded_meta = decode_dna_with_design(dna, rule, meta)
            assert decoded == source, f"Round-trip failed: {rule}, {size} bytes"
            assert meta.get("dna_design_permutation_mode") == DNA_DESIGN_PERMUTATION_MODE
            assert decoded_meta.get("dna_design_permutation_mode") == DNA_DESIGN_PERMUTATION_MODE

    # R0 must remain homopolymer-free across independent block boundaries.
    for source in (b"\x00" * 300, b"\xff" * 300, (b"ABC" * 100)):
        dna, _, meta = encode_bytes_to_dna(source, "R0_B9")
        assert dna_codec.longest_homopolymer(dna) == 1
        assert meta.get("rule_boundary_mode") == "self_contained_anchor"
        assert int(meta.get("rule_anchor_dimers_per_block", 0)) == 1

    # Diagnostics on an intentionally low-entropy source.
    source = b"\x00" * 300
    print("\nV4 positional-distribution diagnostics (source = 300 zero bytes)")
    print("Rule      nt     Base-H   Dimer-H  Unique dimers   GC      HomoMax")
    for rule in RULES:
        dna, _, _ = encode_bytes_to_dna(source, rule)
        dimers = [dna[i:i + 2] for i in range(0, len(dna), 2)]
        print(
            f"{rule:9s} {len(dna):5d}  {entropy(dna):7.3f}  {entropy(dimers):7.3f}"
            f"      {len(set(dimers)):2d}        {dna_codec.gc_content(dna):.3f}      "
            f"{dna_codec.longest_homopolymer(dna)}"
        )

    print("\nPASS: positional V4 codec diagnostics completed.")


if __name__ == "__main__":
    main()
