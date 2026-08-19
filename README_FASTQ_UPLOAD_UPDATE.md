# FASTQ upload in Decoding

Step 5 Decoding now supports three DNA sources:

1. Consensus DNA generated in Strand Design / sequencing simulation.
2. Upload FASTQ.
3. Designed strands (clean baseline).

For Upload FASTQ, the app parses single-end FASTQ records, identifies strands from the current SI/strand design, extracts payload and Phred+33 qualities, performs either Majority voting or Q-score weighted consensus, reassembles the mapped DNA payload, and then applies Known-rule or Unknown-rule auto-detection decoding.

The current upload parser expects FASTQ reads corresponding to the current Strand Design. ECC remains None.
