from __future__ import annotations

import os
import tempfile
from pathlib import Path
import sys, types
if "streamlit" not in sys.modules:
    sys.modules["streamlit"] = types.ModuleType("streamlit")

from compression_pipeline import run_compression_benchmark
from dna_mapping import encode_bytes_to_dna, decode_dna_with_mapping
from fragments import prepare_dna_strands
from sequencing_simulator_v2 import simulate_reads, consensus_by_si

RULES = ["Simple Mapping", "R0_B9", "R1_B12", "R2_B15", "RINF_B16"]


def clean_reassemble(rows, n):
    return "".join(r["Payload"] for r in rows)[:n]


def main():
    source = (b"Whole-file DNA storage regression test\n" * 11) + bytes(range(64))

    # Compression core: selected candidate is a self-describing stored byte stream.
    with tempfile.TemporaryDirectory() as td:
        p = Path(td) / "sample.bin"
        p.write_bytes(source)
        best, candidates = run_compression_benchmark(str(p), source)
        assert best.data
        assert candidates
        stored_cases = [("No compression", source), (best.method, best.data)]

        for storage_name, stored in stored_cases:
            for rule in RULES:
                dna, bits, meta = encode_bytes_to_dna(stored, rule)
                decoded, _, dmeta = decode_dna_with_mapping(dna, rule, meta)
                assert decoded == stored, (storage_name, rule, "direct roundtrip")

                rows = prepare_dna_strands(
                    dna,
                    fbr="ACACGACGCTCTTCCGATCT",
                    rbr="AGATCGGAAGAGCACACGTCT",
                    index_len=8,
                    target_total_len=125,
                    add_filler=True,
                )
                clean_dna = clean_reassemble(rows, len(dna))
                assert clean_dna == dna

                reads, stats = simulate_reads(
                    rows,
                    coverage=10,
                    substitution_rate=0.02,
                    insertion_rate=0.0,
                    deletion_rate=0.0,
                    strand_dropout_rate=0.0,
                    seed=19,
                )
                crows, cdna, cstats, _ = consensus_by_si(rows, reads, len(dna))
                assert cstats["all_strands_complete"]
                assert len(cdna) == len(dna)
                recovered, _, _ = decode_dna_with_mapping(cdna, rule, meta)
                assert recovered == stored, (storage_name, rule, "consensus roundtrip", cstats)

    # No ECC means dropout must remain visible and exact payload reconstruction must stop.
    dna, _, _ = encode_bytes_to_dna(source, "RINF_B16")
    rows = prepare_dna_strands(dna, fbr="ACACGACGCTCTTCCGATCT", rbr="AGATCGGAAGAGCACACGTCT", index_len=8, target_total_len=125)
    reads, _ = simulate_reads(rows, coverage=5, strand_dropout_rate=1.0, seed=1)
    _, cdna, cstats, _ = consensus_by_si(rows, reads, len(dna))
    assert not cdna
    assert cstats["missing_strands"] == len(rows)

    print("PASS: whole-file/compression -> rulebase -> strands -> sequencing consensus -> stored_bytes")


if __name__ == "__main__":
    main()
