from __future__ import annotations

from config import MAPPING_OPTIONS
from dna_mapping import decode_dna_with_mapping, encode_bytes_to_dna
from fragments import prepare_dna_strands
from sequencing_simulator_v2 import generate_fastq, simulate_reads


def candidate_matches(dna: str, reference: bytes):
    matches = []
    for rule in MAPPING_OPTIONS:
        meta = {"bytes_len": len(reference)}
        if rule == "Simple Mapping":
            meta["bits_len"] = len(reference) * 8
        try:
            decoded, _, _ = decode_dna_with_mapping(dna, rule, meta)
            if decoded == reference:
                matches.append(rule)
        except Exception:
            pass
    return matches


def main() -> None:
    reference = b"\x89PNG\r\n\x1a\n" + bytes(range(64)) * 3
    for rule in MAPPING_OPTIONS:
        dna, _, _ = encode_bytes_to_dna(reference, rule)
        assert candidate_matches(dna, reference) == [rule]

    dna, _, _ = encode_bytes_to_dna(bytes(range(128)) * 4, "RINF_B16")
    rows = prepare_dna_strands(
        dna,
        fbr="ACACGACGCTCTTCCGATCT",
        rbr="AGATCGGAAGAGCACACGTCT",
        index_len=8,
        target_total_len=125,
        add_filler=True,
    )[:5]

    profiles = [
        ("Normal Q distribution", {"q_mean": 25, "q_std": 3}),
        ("Linear quality decay", {"q_start": 35, "q_end": 10}),
        ("Fixed Q", {"q_fixed": 20}),
    ]
    for profile, kwargs in profiles:
        reads, stats = simulate_reads(
            rows,
            coverage=3,
            seed=7,
            error_probability_model="Phred Q-score driven",
            qscore_profile=profile,
            **kwargs,
        )
        fastq, fq_stats = generate_fastq(reads)
        lines = fastq.splitlines()
        assert len(lines) == 4 * len(reads)
        assert fq_stats["records"] == len(reads)
        assert stats["qscore_profile"] == profile
        for i in range(0, len(lines), 4):
            assert len(lines[i + 1]) == len(lines[i + 3])

    test_consensus_methods()
    print("PASS: FASTQ Q-scores, majority/Q-weighted consensus, and five-rule auto-detection core checks")

def test_consensus_methods() -> None:
    from sequencing_simulator_v2 import consensus_by_si

    strand_rows = [{"No.": 1, "Strand index": "AAAA", "Payload": "C"}]
    # Three low-quality A calls versus two very high-quality C calls.
    # Majority should pick A; Q-score weighted consensus should pick C.
    read_rows = []
    for i in range(3):
        read_rows.append({
            "Read ID": f"lowA_{i}", "Strand ID": 1, "SI": "AAAA",
            "Read payload": "A", "Payload quality string": chr(2 + 33),
        })
    for i in range(2):
        read_rows.append({
            "Read ID": f"highC_{i}", "Strand ID": 1, "SI": "AAAA",
            "Read payload": "C", "Payload quality string": chr(40 + 33),
        })

    _, dna_majority, stats_majority, _ = consensus_by_si(
        strand_rows, read_rows, 1, method="Majority voting"
    )
    _, dna_q, stats_q, _ = consensus_by_si(
        strand_rows, read_rows, 1, method="Q-score weighted consensus"
    )
    assert dna_majority == "A"
    assert dna_q == "C"
    assert stats_majority["method"] == "Majority voting"
    assert stats_q["method"] == "Q-score weighted consensus"

if __name__ == "__main__":
    main()
