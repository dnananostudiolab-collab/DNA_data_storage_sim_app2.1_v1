# FASTQ decode fix

This build fixes the Streamlit state handling for uploaded FASTQ reconstruction.

- Uploaded FASTQ consensus state is preserved when the app reruns for Decode.
- Changing FASTQ consensus method invalidates only the old FASTQ consensus/decode result.
- Known-rule decoding is displayed as `Selected rule`; only blind mode is displayed as `Auto-detected rule`.
- Uploaded FASTQ byte-recovery mismatch now shows an explicit diagnostic warning.
- SI carried by an uploaded read is checked against the current strand row before the read is accepted.

The DNA mapping/rulebase and compression cores are unchanged.

- Newly exported FASTQ headers include a short design fingerprint; uploading FASTQ from a different encoded payload is rejected before consensus/decode.
