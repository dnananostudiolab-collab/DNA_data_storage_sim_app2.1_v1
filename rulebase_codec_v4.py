from __future__ import annotations

"""Canonical v4 rule-based DNA codec used by the app.

Design invariants
-----------------
* R0/R1/R2 encode independent 3-byte source blocks. Decoder state never crosses
  a block boundary, so a substitution cannot numerically/FSM-propagate into a
  later source block.
* R0 adds one self-contained anchor dimer per block. The anchor is carried in
  the physical DNA itself and is selected to preserve HomoMax=1 across block
  boundaries without making decoding depend on the previous block.
* Rinf is continuous and position-addressed: each 4-bit nibble maps to one dimer
  using a deterministic permutation of all 16 dimers at that absolute dimer
  position. It has no previous-dimer state.
* Positional permutation is deterministic scrambling only; it is not ECC.
"""

import hashlib
import math
import random
from typing import Any, Dict, Iterable, List, Sequence, Tuple

import dna_codec
from dna_design_rules import (
    BLOCK_LOCAL_RULE_DESIGNS,
    DNA_DESIGN_PERMUTATION_MODE,
    DNA_DESIGN_PERMUTATION_SEED,
    DNA_DESIGN_R0,
    DNA_DESIGN_R1,
    DNA_DESIGN_R2,
    DNA_DESIGN_RINF,
    RULE_BASE,
    RULE_BLOCK_BYTES,
    RULE_CODEC_VERSION,
    RULE_ERROR_POLICY,
    RULE_INITIAL_DIMER,
    add_legacy_rule_meta_aliases,
)


def _stable_rng(key: str) -> random.Random:
    digest = hashlib.sha256(key.encode("utf-8")).digest()
    return random.Random(digest)


def _permuted(items: Sequence[str], key: str) -> List[str]:
    out = list(items)
    _stable_rng(key).shuffle(out)
    return out


def _fixed_digits_for_bits(base: int, n_bits: int) -> int:
    k = 1
    limit = 1 << int(n_bits)
    while base ** k < limit:
        k += 1
    return k


def _int_to_fixed_base(value: int, base: int, count: int) -> List[int]:
    if value < 0:
        raise ValueError("value must be non-negative")
    digits = [0] * int(count)
    n = int(value)
    for i in range(count - 1, -1, -1):
        digits[i] = n % base
        n //= base
    if n:
        raise ValueError("value does not fit fixed radix width")
    return digits


def _fixed_base_to_int(digits: Iterable[int], base: int) -> int:
    value = 0
    for d in digits:
        d = int(d)
        if not 0 <= d < base:
            raise ValueError(f"digit {d} out of range for base {base}")
        value = value * base + d
    return value


def _position_key(rule: str, block_index: int, local_position: int, prev: str = "") -> str:
    return (
        f"{DNA_DESIGN_PERMUTATION_SEED}|{DNA_DESIGN_PERMUTATION_MODE}|"
        f"{rule}|block={block_index}|pos={local_position}|prev={prev}"
    )


def _rinf_position_key(position: int) -> str:
    return (
        f"{DNA_DESIGN_PERMUTATION_SEED}|{DNA_DESIGN_PERMUTATION_MODE}|"
        f"{DNA_DESIGN_RINF}|absolute={position}"
    )


def _r0_anchor(block_index: int, previous_physical_last: str) -> str:
    candidates = sorted(dna_codec.TABLE_R0_B9.keys())
    if previous_physical_last:
        safe = [d for d in candidates if d[0] != previous_physical_last and d[0] != d[1]]
        if safe:
            candidates = safe
    candidates = _permuted(
        candidates,
        f"{DNA_DESIGN_PERMUTATION_SEED}|R0-anchor|block={block_index}|left={previous_physical_last or '-'}",
    )
    return candidates[0]


def _common_meta(rule: str, bytes_len: int) -> Dict[str, Any]:
    return {
        "dna_design_name": rule,
        "mapping": rule,  # legacy UI alias
        "rule_name": rule,
        "scheme_name": rule,
        "rule_base": int(RULE_BASE[rule]),
        "rule_initial_dimer": RULE_INITIAL_DIMER,
        "init_dimer": RULE_INITIAL_DIMER,
        "bytes_len": int(bytes_len),
        "dna_design_codec_version": RULE_CODEC_VERSION,
        "dna_design_permutation_mode": DNA_DESIGN_PERMUTATION_MODE,
        "dna_design_permutation_seed": DNA_DESIGN_PERMUTATION_SEED,
        "error_policy": RULE_ERROR_POLICY,
    }


def encode_rulebase(data: bytes, rule: str) -> Tuple[str, Dict[str, Any]]:
    raw = bytes(data or b"")
    rule = str(rule)
    if rule == DNA_DESIGN_RINF:
        out: List[str] = []
        for byte_index, value in enumerate(raw):
            nibbles = ((value >> 4) & 0xF, value & 0xF)
            for nibble_offset, nibble in enumerate(nibbles):
                pos = byte_index * 2 + nibble_offset
                lut = _permuted(dna_codec.DIMERS, _rinf_position_key(pos))
                out.append(lut[nibble])
        meta = _common_meta(rule, len(raw))
        meta.update({
            "dna_design_mode": "continuous",
            "block_local": False,
            "rule_block_bytes": 0,
            "rule_block_count": 0,
            "rule_digits_per_block": 0,
            "rule_nt_per_block": 0,
            "rule_anchor_dimers_per_block": 0,
            "dna_design_position_scope": "absolute_dimer_position",
            "error_propagation_scope": "single_dimer_4_bits",
            "dna_length_nt": len(out) * 2,
        })
        return "".join(out), add_legacy_rule_meta_aliases(meta)

    if rule not in BLOCK_LOCAL_RULE_DESIGNS:
        raise ValueError(f"Unsupported rule design: {rule}")

    base = int(RULE_BASE[rule])
    digits_per_block = _fixed_digits_for_bits(base, RULE_BLOCK_BYTES * 8)
    anchor_count = 1 if rule == DNA_DESIGN_R0 else 0
    nt_per_block = 2 * (digits_per_block + anchor_count)
    blocks = max(1, math.ceil(len(raw) / RULE_BLOCK_BYTES)) if raw else 0
    table = dna_codec.TABLES[rule]

    out: List[str] = []
    previous_physical_last = ""
    for block_index in range(blocks):
        chunk = raw[block_index * RULE_BLOCK_BYTES:(block_index + 1) * RULE_BLOCK_BYTES]
        padded = chunk + b"\x00" * (RULE_BLOCK_BYTES - len(chunk))
        value = int.from_bytes(padded, "big")
        digits = _int_to_fixed_base(value, base, digits_per_block)

        if rule == DNA_DESIGN_R0:
            anchor = _r0_anchor(block_index, previous_physical_last)
            out.append(anchor)
            prev = anchor
        else:
            prev = RULE_INITIAL_DIMER

        for local_position, digit in enumerate(digits):
            if prev not in table:
                raise ValueError(f"No table row for prev dimer {prev} in {rule}")
            allowed = table[prev]
            lut = _permuted(allowed, _position_key(rule, block_index, local_position, prev))
            next_dimer = lut[digit]
            out.append(next_dimer)
            prev = next_dimer

        if out:
            previous_physical_last = out[-1][-1]

    dna = "".join(out)
    meta = _common_meta(rule, len(raw))
    meta.update({
        "dna_design_mode": "block_local",
        "block_local": True,
        "rule_block_bytes": RULE_BLOCK_BYTES,
        "rule_block_count": blocks,
        "rule_digits_per_block": digits_per_block,
        "rule_nt_per_block": nt_per_block,
        "rule_anchor_dimers_per_block": anchor_count,
        "rule_boundary_mode": "self_contained_anchor" if rule == DNA_DESIGN_R0 else "fixed_initial_state_per_block",
        "dna_design_position_scope": "block_index_and_local_dimer_position",
        "error_propagation_scope": f"one_{RULE_BLOCK_BYTES}_byte_source_block",
        "erasure_placeholder_byte": 0,
        "dna_length_nt": len(dna),
    })
    return dna, add_legacy_rule_meta_aliases(meta)


def decode_rulebase(dna: str, rule: str, meta: Dict[str, Any] | None = None) -> Tuple[bytes, Dict[str, Any]]:
    seq = dna_codec.clean_dna_text(dna)
    rule = str(rule)
    meta_in = dict(meta or {})
    bytes_len = int(meta_in.get("bytes_len", 0))

    if rule == DNA_DESIGN_RINF:
        if len(seq) % 4 != 0:
            raise ValueError("R∞ DNA length must be a multiple of 4 nt (two dimers per source byte).")
        n_dimers = len(seq) // 2
        nibbles: List[int] = []
        for pos in range(n_dimers):
            dimer = seq[2 * pos:2 * pos + 2]
            lut = _permuted(dna_codec.DIMERS, _rinf_position_key(pos))
            nibbles.append(lut.index(dimer))  # all canonical dimers are valid in R∞
        if len(nibbles) % 2:
            raise ValueError("R∞ stream contains an odd number of nibble dimers.")
        out = bytes(((nibbles[i] << 4) | nibbles[i + 1]) for i in range(0, len(nibbles), 2))
        if bytes_len:
            out = out[:bytes_len]
        dmeta = _common_meta(rule, bytes_len or len(out))
        dmeta.update({
            "dna_design_mode": "continuous",
            "block_local": False,
            "dna_design_position_scope": "absolute_dimer_position",
            "error_propagation_scope": "single_dimer_4_bits",
            "corrupted_block_count": 0,
            "corrupted_blocks": [],
            "erasure_block_mask": [],
            "erasure_byte_ranges": [],
        })
        return out, add_legacy_rule_meta_aliases(dmeta)

    if rule not in BLOCK_LOCAL_RULE_DESIGNS:
        raise ValueError(f"Unsupported rule design: {rule}")

    base = int(RULE_BASE[rule])
    digits_per_block = int(meta_in.get("rule_digits_per_block") or _fixed_digits_for_bits(base, RULE_BLOCK_BYTES * 8))
    anchor_count = int(meta_in.get("rule_anchor_dimers_per_block", 1 if rule == DNA_DESIGN_R0 else 0))
    nt_per_block = int(meta_in.get("rule_nt_per_block") or 2 * (digits_per_block + anchor_count))
    block_count = int(meta_in.get("rule_block_count") or (math.ceil(bytes_len / RULE_BLOCK_BYTES) if bytes_len else 0))
    if block_count <= 0:
        if not seq:
            return b"", add_legacy_rule_meta_aliases(_common_meta(rule, 0))
        if len(seq) % nt_per_block:
            raise ValueError("Rule DNA length does not align to complete source blocks.")
        block_count = len(seq) // nt_per_block
        bytes_len = block_count * RULE_BLOCK_BYTES

    table = dna_codec.TABLES[rule]
    output = bytearray()
    corrupted_blocks: List[int] = []
    erasure_rows: List[Dict[str, Any]] = []
    mask: List[int] = []

    for block_index in range(block_count):
        start = block_index * nt_per_block
        block_seq = seq[start:start + nt_per_block]
        corrupted_reason = ""
        digits: List[int] = []

        if len(block_seq) != nt_per_block:
            corrupted_reason = "missing_or_wrong_length_block"
        else:
            cursor = 0
            if rule == DNA_DESIGN_R0:
                anchor = block_seq[:2]
                cursor = 2
                if anchor not in table:
                    corrupted_reason = "invalid_r0_anchor"
                    prev = RULE_INITIAL_DIMER
                else:
                    prev = anchor
            else:
                prev = RULE_INITIAL_DIMER

            if not corrupted_reason:
                for local_position in range(digits_per_block):
                    dimer = block_seq[cursor:cursor + 2]
                    cursor += 2
                    allowed = table.get(prev)
                    if not allowed:
                        corrupted_reason = f"missing_table_row_prev_{prev}"
                        break
                    lut = _permuted(allowed, _position_key(rule, block_index, local_position, prev))
                    if dimer not in lut:
                        corrupted_reason = f"invalid_dimer_at_local_position_{local_position}"
                        break
                    digit = lut.index(dimer)
                    digits.append(digit)
                    prev = dimer

        if not corrupted_reason:
            value = _fixed_base_to_int(digits, base)
            if value >= (1 << (RULE_BLOCK_BYTES * 8)):
                corrupted_reason = "impossible_radix_codeword"

        if corrupted_reason:
            block_bytes = b"\x00" * RULE_BLOCK_BYTES
            corrupted_blocks.append(block_index + 1)
            mask.append(1)
            erasure_rows.append({
                "Block": block_index + 1,
                "Byte start": block_index * RULE_BLOCK_BYTES,
                "Byte end": block_index * RULE_BLOCK_BYTES + RULE_BLOCK_BYTES - 1,
                "Reason": corrupted_reason,
            })
        else:
            block_bytes = value.to_bytes(RULE_BLOCK_BYTES, "big")
            mask.append(0)

        output.extend(block_bytes)

    if bytes_len:
        output = output[:bytes_len]

    dmeta = _common_meta(rule, bytes_len or len(output))
    dmeta.update({
        "dna_design_mode": "block_local",
        "block_local": True,
        "rule_block_bytes": RULE_BLOCK_BYTES,
        "rule_block_count": block_count,
        "rule_digits_per_block": digits_per_block,
        "rule_nt_per_block": nt_per_block,
        "rule_anchor_dimers_per_block": anchor_count,
        "rule_boundary_mode": "self_contained_anchor" if rule == DNA_DESIGN_R0 else "fixed_initial_state_per_block",
        "dna_design_position_scope": "block_index_and_local_dimer_position",
        "error_propagation_scope": f"one_{RULE_BLOCK_BYTES}_byte_source_block",
        "error_containment": True,
        "valid_block_count": block_count - len(corrupted_blocks),
        "corrupted_block_count": len(corrupted_blocks),
        "corrupted_blocks": corrupted_blocks,
        "rule_erasure_count": len(corrupted_blocks),
        "rule_erasure_rows": erasure_rows,
        "erasure_block_mask": mask,
        "erasure_byte_ranges": [
            [r["Byte start"], min(r["Byte end"], max(0, (bytes_len or len(output)) - 1))]
            for r in erasure_rows
        ],
        "erasure_placeholder_byte": 0,
        "rule_override_rows": [],
    })
    return bytes(output), add_legacy_rule_meta_aliases(dmeta)
