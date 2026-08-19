"""Regression tests for substitution error containment in rule-based DNA designs.

Core invariant:
- R0/R1/R2: one substitution may damage the containing 3-byte source block,
  but bytes before and after that block must remain unchanged.
- Rinf: continuous base-16/nibble encoding; one substitution may damage only
  the source byte containing the affected dimer, never downstream bytes.
"""
from __future__ import annotations

from dna_design_rules import (
    BLOCK_LOCAL_RULE_DESIGNS,
    DNA_DESIGN_RINF,
    RULE_BLOCK_BYTES,
)
from dna_mapping import decode_dna_with_design, encode_bytes_to_dna

BASES = "ACGT"


def _mutate(seq: str, position: int, replacement: str) -> str:
    return seq[:position] + replacement + seq[position + 1:]


def _test_block_local_design(dna_design_name: str, source: bytes) -> tuple[int, int]:
    dna, _, dna_design_meta = encode_bytes_to_dna(source, dna_design_name)
    clean, _, clean_meta = decode_dna_with_design(dna, dna_design_name, dna_design_meta)
    assert clean == source, f"Clean round-trip failed for {dna_design_name}"
    assert clean_meta.get("dna_design_mode") == "block_local"
    assert int(clean_meta.get("rule_block_bytes")) == RULE_BLOCK_BYTES

    nt_per_block = int(dna_design_meta["rule_nt_per_block"])
    block_count = int(dna_design_meta["rule_block_count"])
    tested = detected = 0

    # Exhaust every single-base substitution at every encoded position.
    for position in range(len(dna)):
        target_block_index = position // nt_per_block
        assert target_block_index < block_count
        byte_start = target_block_index * RULE_BLOCK_BYTES
        byte_end = min(byte_start + RULE_BLOCK_BYTES, len(source))

        for replacement in BASES:
            if replacement == dna[position]:
                continue
            tested += 1
            mutated = _mutate(dna, position, replacement)
            decoded, _, decoded_meta = decode_dna_with_design(
                mutated, dna_design_name, dna_design_meta
            )

            # The containing block may be wrong/erased, but no byte outside it
            # may change. This is the core no-propagation invariant.
            assert decoded[:byte_start] == source[:byte_start], (
                f"Backward propagation in {dna_design_name} at nt {position + 1}"
            )
            assert decoded[byte_end:] == source[byte_end:], (
                f"Forward propagation in {dna_design_name} at nt {position + 1}"
            )

            if decoded_meta.get("corrupted_block_count", 0):
                detected += 1
                assert decoded_meta.get("corrupted_blocks") == [target_block_index + 1]
                mask = decoded_meta.get("erasure_block_mask", [])
                assert len(mask) == block_count
                assert mask[target_block_index] == 1

    return tested, detected


def _test_rinf_continuous(source: bytes) -> int:
    dna, _, dna_design_meta = encode_bytes_to_dna(source, DNA_DESIGN_RINF)
    clean, _, clean_meta = decode_dna_with_design(dna, DNA_DESIGN_RINF, dna_design_meta)
    assert clean == source, "Clean round-trip failed for Rinf"
    assert clean_meta.get("dna_design_mode") == "continuous"
    assert clean_meta.get("block_local") is False
    assert len(dna) == len(source) * 4  # 2 dimers / byte, 2 nt / dimer

    tested = 0
    for position in range(len(dna)):
        affected_dimer = position // 2
        affected_byte = affected_dimer // 2
        for replacement in BASES:
            if replacement == dna[position]:
                continue
            tested += 1
            mutated = _mutate(dna, position, replacement)
            decoded, _, decoded_meta = decode_dna_with_design(
                mutated, DNA_DESIGN_RINF, dna_design_meta
            )
            assert len(decoded) == len(source)
            assert decoded[:affected_byte] == source[:affected_byte]
            assert decoded[affected_byte + 1:] == source[affected_byte + 1:]
            assert decoded_meta.get("error_propagation_scope") == "single_dimer_4_bits"
    return tested


def main() -> None:
    source = bytes(range(1, 40))  # 13 complete 3-byte blocks

    for dna_design_name in sorted(BLOCK_LOCAL_RULE_DESIGNS):
        tested, detected = _test_block_local_design(dna_design_name, source)
        print(
            f"PASS {dna_design_name}: {tested} single substitutions; "
            f"0 propagated beyond the 3-byte block; {detected} detected as erasure"
        )

    tested_rinf = _test_rinf_continuous(source)
    print(
        f"PASS {DNA_DESIGN_RINF}: {tested_rinf} single substitutions; "
        "0 propagated beyond the affected source byte; continuous/no blocks"
    )


if __name__ == "__main__":
    main()
