# No-compression version

Flow: Input → Encoding → Strand Design / Sequencing / FASTQ / Consensus → Decoding → Summarization.

Run:
```bash
streamlit run app.py
```

There is no Compression step and no intermediate button. As soon as a file is uploaded, `stored_bytes` is set to the exact uploaded file byte stream and Encoding becomes available.
