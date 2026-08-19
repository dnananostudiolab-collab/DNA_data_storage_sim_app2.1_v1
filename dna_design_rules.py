"""Canonical definitions for DNA design names and rule-based framing.

The project uses one vocabulary everywhere:
    dna_design_name  -> selected DNA design (R0/R1/R2/Rinf/SM)
    dna_design_meta  -> metadata required to decode that design

R0/R1/R2 are intentionally block-local to bound substitution error propagation.
V4 applies deterministic position-dependent dimer permutations to all rule designs.
R0 additionally stores a self-contained boundary anchor per block so clean DNA can
retain HomoMax=1 without using the previous decoded block as state. Rinf remains
continuous and uses the absolute dimer position for its 16-way permutation.
"""
from __future__ import annotations

from typing import Any, Dict

DNA_DESIGN_SIMPLE = "Simple Mapping"
DNA_DESIGN_R0 = "R0_B9"
DNA_DESIGN_R1 = "R1_B12"
DNA_DESIGN_R2 = "R2_B15"
DNA_DESIGN_RINF = "RINF_B16"
DNA_DESIGN_PAPER_RINF_P8 = "RINF_P8_PAPER"
DNA_DESIGN_PAPER_R1_P8 = "R1_P8_PAPER"

RULE_BASED_DNA_DESIGNS = (
    DNA_DESIGN_R0,
    DNA_DESIGN_R1,
    DNA_DESIGN_R2,
    DNA_DESIGN_RINF,
)

# Original paper Rulebase-P controls. These are intentionally separate from
# the app rulebase because the paper uses direct 4-bit->dimer P8 mapping,
# not the app block/radix architecture.
PAPER_P8_DNA_DESIGNS = (
    DNA_DESIGN_PAPER_RINF_P8,
    DNA_DESIGN_PAPER_R1_P8,
)

# R0/R1/R2 use the same 3-byte block size on purpose.  This makes the
# containment comparison fair: a substitution can affect at most one 3-byte
# source block, regardless of which constrained rule is selected.
BLOCK_LOCAL_RULE_DESIGNS = frozenset({DNA_DESIGN_R0, DNA_DESIGN_R1, DNA_DESIGN_R2})
CONTINUOUS_RULE_DESIGNS = frozenset({DNA_DESIGN_RINF})

RULE_INITIAL_DIMER = "TA"
RULE_BLOCK_BYTES = 3
RULE_CODEC_VERSION = "rulebase-4.0"
RULE_ERROR_POLICY = "mark_block_erasure_no_repair"

# Deterministic positional permutation used to decorrelate repetitive payloads.
# A fixed default seed keeps the current app deterministic; it can be exposed in
# the UI later without changing the codec architecture.
DNA_DESIGN_PERMUTATION_MODE = "POSITIONAL_V1"
DNA_DESIGN_PERMUTATION_SEED = "app5-rulebase-positional-v1"

RULE_BASE = {
    DNA_DESIGN_R0: 9,
    DNA_DESIGN_R1: 12,
    DNA_DESIGN_R2: 15,
    DNA_DESIGN_RINF: 16,
}

DNA_DESIGN_DISPLAY = {
    DNA_DESIGN_SIMPLE: "SM",
    DNA_DESIGN_R0: "R0",
    DNA_DESIGN_R1: "R1",
    DNA_DESIGN_R2: "R2",
    DNA_DESIGN_RINF: "R∞",
    DNA_DESIGN_PAPER_RINF_P8: "R∞-P8 (Paper)",
    DNA_DESIGN_PAPER_R1_P8: "R1-P8 (Paper)",
}


def is_rule_based_design(dna_design_name: str) -> bool:
    return str(dna_design_name) in RULE_BASED_DNA_DESIGNS


def is_block_local_rule(dna_design_name: str) -> bool:
    return str(dna_design_name) in BLOCK_LOCAL_RULE_DESIGNS


def is_continuous_rule(dna_design_name: str) -> bool:
    return str(dna_design_name) in CONTINUOUS_RULE_DESIGNS


def display_dna_design(dna_design_name: str) -> str:
    return DNA_DESIGN_DISPLAY.get(str(dna_design_name), str(dna_design_name))


def prepare_rule_decode_meta(dna_design_meta: Dict[str, Any] | None, dna_design_name: str) -> Dict[str, Any]:
    """Return decode metadata with the canonical no-propagation policy.

    Legacy metadata keys are accepted so saved sessions from older app versions
    remain decodable.  New code should use the dna_design_* / rule_* keys.
    """
    meta = dict(dna_design_meta or {})
    name = str(dna_design_name)

    # Promote legacy keys into canonical names when needed.
    meta.setdefault("dna_design_name", meta.get("mapping", name))
    meta.setdefault("dna_design_mode", "block_local" if is_block_local_rule(name) else "continuous")
    meta.setdefault("rule_name", meta.get("scheme_name", name))
    meta.setdefault("rule_initial_dimer", meta.get("init_dimer", RULE_INITIAL_DIMER))

    if is_block_local_rule(name):
        meta.setdefault("rule_block_bytes", meta.get("block_bytes", RULE_BLOCK_BYTES))
        meta.setdefault("error_policy", RULE_ERROR_POLICY)
        meta.setdefault("erasure_placeholder_byte", 0)
        meta["invalid_rule_handling"] = "Mark block corrupted (no propagation)"

    return meta


def add_legacy_rule_meta_aliases(meta: Dict[str, Any]) -> Dict[str, Any]:
    """Add read-only compatibility aliases expected by older UI code."""
    out = dict(meta)
    name = str(out.get("dna_design_name", out.get("mapping", "")))
    if name:
        out.setdefault("mapping", name)
        out.setdefault("scheme_name", out.get("rule_name", name))
    if "rule_initial_dimer" in out:
        out.setdefault("init_dimer", out["rule_initial_dimer"])
    if "rule_block_bytes" in out:
        out.setdefault("block_bytes", out["rule_block_bytes"])
    if "rule_digits_per_block" in out:
        out.setdefault("digits_per_block", out["rule_digits_per_block"])
    if "rule_nt_per_block" in out:
        out.setdefault("dna_nt_per_block", out["rule_nt_per_block"])
    if "dna_design_mode" in out:
        out.setdefault("block_local", out["dna_design_mode"] == "block_local")
    return out
