# DNA Data Storage System — whole-file pipeline

Active 6-step UI:

1. **Input** — upload the complete file/container (max 200 MB) and preview it.
2. **Compression** — exact no-compression bytes or the existing compression benchmark. The selected output becomes `stored_bytes`.
3. **Encoding** — map `stored_bytes` with SM / R0 / R1 / R2 / R∞. DNA structural-randomness analysis is available here.
4. **Strand Design** — FBR + SI + Payload + Filler + RBR, followed by sequencing-read simulation, explicit FASTQ generation/export, and SI-based consensus reconstruction. ECC is an independent layer and is currently `None`.
5. **Decoding** — use a known rule or auto-test SM/R0/R1/R2/R∞, rule-decode to `stored_bytes`, undo lossless wrappers where applicable, and recover the application file.
6. **Summarization** — original/stored/decoded comparison plus compression, DNA, randomness, strand, sequencing, consensus and recovery summaries.

## Core invariants

- No pixel representation is used by the DNA storage core.
- Compression remains in `compression_pipeline.py`.
- Rule mapping remains in `rulebase_codec_v4.py` / `dna_mapping.py`.
- Only `stored_bytes` are sent to DNA encoding.
- Randomness analysis evaluates the mapped DNA payload, excluding primer/SI/filler/read noise.
- FASTQ is generated from simulated sequencing reads with Phred+33 qualities. Bernoulli mode uses the Q equivalent of its fixed substitution probability; Q-score-driven mode uses Normal, linear-decay, or fixed-Q profiles and derives substitution probability from each Q score.
- With ECC=None, strand dropout/incomplete consensus cannot be repaired.

## Sequencing Q-score models and FASTQ

Step 4 supports two substitution probability models:

- `Independent Bernoulli per nucleotide`: one configured substitution probability is applied independently to every payload nucleotide. FASTQ quality uses the equivalent fixed Phred score.
- `Phred Q-score driven`: each payload base receives a simulated Q score and substitution probability is `p = 10^(-Q/10)`.

Q-score-driven profiles:

- `Normal Q distribution`
- `Linear quality decay`
- `Fixed Q`

The generated FASTQ uses Phred+33 quality strings whose length exactly matches each simulated read sequence. Insertions and deletions remain separate independent Bernoulli processes. FBR, SI, filler, and RBR remain protected in this simulator; whole-strand dropout is modeled separately.

## Decoding when the rule is unknown

Step 5 has two modes:

- `Known rule`: the user selects SM/R0/R1/R2/R∞.
- `Unknown rule — auto-detect`: all five decoding rules are attempted independently. The encoded-rule identity stored by the UI is not supplied to the candidate decoders.

For the in-app simulation, `stored_bytes` are retained as ground truth only for validation. A rule is identified when it reconstructs the exact stored byte stream. Container recognition/validation is also reported for every candidate. A truly external blind decode without any reference bytes would require a self-describing framing/checksum layer for robust rule identification of arbitrary raw binary files.
