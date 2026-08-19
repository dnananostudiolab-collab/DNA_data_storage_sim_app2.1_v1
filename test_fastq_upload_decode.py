from __future__ import annotations

from config import MAPPING_OPTIONS
from dna_mapping import decode_dna_with_mapping, encode_bytes_to_dna
from fragments import prepare_dna_strands
from sequencing_simulator_v2 import (
    consensus_by_si,
    generate_fastq,
    parse_fastq_reads,
    simulate_reads,
)

FBR = "ACACGACGCTCTTCCGATCT"
RBR = "AGATCGGAAGAGCACACGTCT"


def main() -> None:
    reference = b"FASTQ upload decoding regression\n" * 8 + bytes(range(80))
    for rule in MAPPING_OPTIONS:
        dna, _, meta = encode_bytes_to_dna(reference, rule)
        rows = prepare_dna_strands(
            dna, fbr=FBR, rbr=RBR, index_len=8,
            target_total_len=125, add_filler=True,
        )
        reads, _ = simulate_reads(
            rows,
            coverage=10,
            seed=23,
            error_probability_model="Phred Q-score driven",
            qscore_profile="Fixed Q",
            q_fixed=25,
            insertion_rate=0.0,
            deletion_rate=0.0,
            strand_dropout_rate=0.0,
        )
        fastq, _ = generate_fastq(reads)
        uploaded_reads, pstats = parse_fastq_reads(fastq, rows)
        assert pstats["records"] == len(reads)
        assert pstats["matched_reads"] == len(reads)
        assert pstats["unmatched_reads"] == 0
        _, cdna, cstats, _ = consensus_by_si(
            rows, uploaded_reads, len(dna), method="Q-score weighted consensus"
        )
        assert cstats["all_strands_complete"], (rule, cstats)
        recovered, _, _ = decode_dna_with_mapping(cdna, rule, meta)
        assert recovered == reference, rule

    print("PASS: exported FASTQ -> upload parser -> SI grouping -> Q-score consensus -> rule decode")


if __name__ == "__main__":
    main()
