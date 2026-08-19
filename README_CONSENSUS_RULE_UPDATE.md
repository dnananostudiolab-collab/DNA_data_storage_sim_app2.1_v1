# Consensus and rule-selection update

The active six-step whole-file DNA-storage UI keeps compression and SM/R0/R1/R2/R∞ unchanged.

## Step 4: sequencing / consensus

- Error probability model and Q-score generation remain part of sequencing-read/FASTQ simulation.
- FASTQ qualities are Phred+33.
- Consensus method now has two independent options:
  1. **Majority voting** — each read has equal weight at a nucleotide position.
  2. **Q-score weighted consensus** — a maximum-likelihood call uses each base's Phred Q score; high-Q calls contribute more evidence than low-Q calls.
- The selected consensus method is recorded in the final summary.

## Step 5: decoding

- **Known rule**: the user selects SM, R0, R1, R2, or R∞.
- **Unknown rule — auto-detect**: the app automatically tries all five decoders and selects the rule that reconstructs the stored file/container exactly in the current simulation workflow.
- Detailed candidate diagnostics are optional and hidden in an expander.
