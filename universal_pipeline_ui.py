from __future__ import annotations

import bz2
import gzip
import hashlib
import io
import json
import lzma
import os
import zipfile
from pathlib import Path
from typing import Any, Dict, List, Tuple

import pandas as pd
import streamlit as st

from config import DNA_PREVIEW_HEIGHT, MAPPING_OPTIONS, WORK_ROOT
from dna_codec import gc_content, homopolymer_stats
from dna_mapping import decode_dna_with_mapping, encode_bytes_to_dna, validate_container
from fragments import clean_dna, choose_auto_strand_design, prepare_dna_strands, strand_rows_to_csv
from randomness_ui import render_randomness_analysis
from restore_analysis import write_restored_file
from sequencing_simulator_v2 import consensus_by_si, generate_fastq, parse_fastq_reads, simulate_reads
from ui_helpers import fmt_bytes, magic_dict, preview_file, save_upload
from utils_bits_v2 import bytes_to_bitstring, detect_magic, sha256_bytes


FBR_DEFAULT = "ACACGACGCTCTTCCGATCT"
RBR_DEFAULT = "AGATCGGAAGAGCACACGTCT"


def _display_rule(name: str) -> str:
    return {
        "Simple Mapping": "SM",
        "R0_B9": "R0",
        "R1_B12": "R1",
        "R2_B15": "R2",
        "RINF_B16": "R∞",
    }.get(str(name), str(name))


def _preview_seq(seq: str, n: int = 600) -> str:
    seq = clean_dna(seq)
    return seq if len(seq) <= n else seq[:n] + f"\n... ({len(seq):,} nt total)"


def _clear(keys: List[str]) -> None:
    for key in keys:
        st.session_state.pop(key, None)


def _clear_after_input() -> None:
    _clear([
        "stored_bytes", "stored_file_path", "storage_method", "storage_meta", "compression_candidates",
        "dna", "bits", "codec_meta", "encoding_mapping", "strand_rows", "strand_config",
        "read_rows", "read_stats", "fastq_text", "fastq_stats",
        "consensus_rows", "consensus_dna", "consensus_stats", "consensus_evidence",
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", 
        "wholefile_randomness_report", "wholefile_randomness_hash", "wholefile_randomness_error",
    ])


def _clear_after_storage() -> None:
    _clear([
        "dna", "bits", "codec_meta", "encoding_mapping", "strand_rows", "strand_config",
        "read_rows", "read_stats", "fastq_text", "fastq_stats",
        "consensus_rows", "consensus_dna", "consensus_stats", "consensus_evidence",
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", 
        "wholefile_randomness_report", "wholefile_randomness_hash", "wholefile_randomness_error",
    ])


def _clear_after_encoding() -> None:
    _clear([
        "strand_rows", "strand_config", "read_rows", "read_stats", "fastq_text", "fastq_stats",
        "consensus_rows", "consensus_dna", "consensus_stats", "consensus_evidence",
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", 
        "wholefile_randomness_report", "wholefile_randomness_hash", "wholefile_randomness_error",
    ])


def _clear_after_strands() -> None:
    _clear([
        "read_rows", "read_stats", "fastq_text", "fastq_stats",
        "consensus_rows", "consensus_dna", "consensus_stats", "consensus_evidence",
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", 
    ])


def _clear_after_reads() -> None:
    _clear([
        "fastq_text", "fastq_stats", "consensus_rows", "consensus_dna", "consensus_stats", "consensus_evidence",
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", 
    ])


def _clear_after_consensus() -> None:
    _clear([
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule", "auto_decode_table",
        "uploaded_fastq_hash", "uploaded_fastq_name", "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats", "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna", "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence",
        "uploaded_fastq_consensus_method", "rule_resolution",
    ])


def _clear_after_uploaded_fastq_consensus() -> None:
    """Clear only downstream decode results; preserve uploaded FASTQ reconstruction state."""
    _clear([
        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
        "recovered_file_name", "restored_info", "decode_source_label", "decode_rule_mode", "detected_rule",
        "auto_decode_table", "rule_resolution",
    ])


def _byte_accuracy(reference: bytes, recovered: bytes) -> float:
    a = bytes(reference or b"")
    b = bytes(recovered or b"")
    denom = max(len(a), len(b), 1)
    matches = sum(x == y for x, y in zip(a, b))
    return matches / denom


def _recover_application_file(stored: bytes, method: str, input_name: str) -> Tuple[bytes, str, str]:
    method = str(method or "")
    data = bytes(stored or b"")
    name = os.path.basename(input_name or "recovered.bin")
    try:
        if method.startswith("gzip_lvl"):
            return gzip.decompress(data), name, "GZIP wrapper decompressed"
        if method.startswith("bz2_lvl"):
            return bz2.decompress(data), name, "BZ2 wrapper decompressed"
        if method.startswith("xz_p"):
            return lzma.decompress(data, format=lzma.FORMAT_XZ), name, "XZ wrapper decompressed"
        if method.startswith("zip_store") or method.startswith("zip_deflate"):
            with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
                members = [m for m in zf.namelist() if not m.endswith("/")]
                if not members:
                    raise ValueError("ZIP wrapper contains no file")
                member = members[0]
                return zf.read(member), os.path.basename(member) or name, "ZIP wrapper extracted"
    except Exception as exc:
        return data, name, f"Wrapper recovery failed; returning stored container ({exc})"

    m = detect_magic(data)
    if method == "No compression" or method == "keep_original":
        return data, name, "Original file bytes recovered"
    ext = m.ext if m else Path(name).suffix or ".bin"
    stem = Path(name).stem or "recovered"
    return data, f"{stem}{ext}", "Stored representation recovered"


def _render_file_preview(data: bytes, name: str, key: str) -> None:
    if not data:
        return
    out = WORK_ROOT / "ui_preview"
    out.mkdir(parents=True, exist_ok=True)
    safe_name = os.path.basename(name or "file.bin")
    path = out / f"{key}_{safe_name}"
    try:
        path.write_bytes(data)
        preview_file(str(path), name or "Preview")
    except Exception:
        st.caption("Preview unavailable for this file type.")


def _render_segmented_strand(row: Dict[str, Any]) -> None:
    st.markdown("**Designed strand**")
    pieces = [
        ("FBR", clean_dna(row.get("FBR", "")) or "—"),
        ("SI", clean_dna(row.get("Strand index", "")) or "—"),
        ("Payload", clean_dna(row.get("Payload", "")) or "—"),
        ("Filler", clean_dna(row.get("Filler", "")) or "—"),
        ("RBR", clean_dna(row.get("RBR", "")) or "—"),
    ]
    st.markdown("  ".join(f"**{label}**: `{value}`" for label, value in pieces))


def render_step_1_input() -> None:
    with st.container(border=True):
        st.markdown("## 1 Input")
        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("#### 📁 Input")
            uploaded = st.file_uploader(
                "Input file",
                type=None,
                key="whole_file_upload",
                label_visibility="collapsed",
            )
            st.caption("200MB per file")
            if uploaded is not None:
                raw = uploaded.getvalue()
                sig = f"{uploaded.name}|{len(raw)}|{hashlib.sha256(raw).hexdigest()}"
                if st.session_state.get("upload_signature") != sig:
                    path, data = save_upload(uploaded)
                    _clear_after_input()
                    st.session_state.update({
                        "upload_signature": sig,
                        "input_path": path,
                        "input_name": uploaded.name,
                        "input_bytes": data,
                        # No-compression route: the exact whole-file bytes are immediately ready for DNA encoding.
                        "stored_bytes": data,
                        "stored_file_path": path,
                        "storage_method": "No compression",
                        "storage_meta": {
                            "kind": "file_bytes",
                            "lossy": False,
                            "input_name": uploaded.name,
                        },
                        "compression_candidates": [],
                    })

            raw = st.session_state.get("input_bytes", b"") or b""
            path = st.session_state.get("input_path", "")
            if raw:
                m = detect_magic(raw)
                st.markdown("##### 📄 File properties")
                st.dataframe(pd.DataFrame([
                    {"Property": "File name", "Value": st.session_state.get("input_name", os.path.basename(path))},
                    {"Property": "File size", "Value": fmt_bytes(len(raw))},
                    {"Property": "Detected container", "Value": m.kind if m else "unknown / raw binary"},
                    {"Property": "Storage route", "Value": "No compression — exact file bytes"},
                    {"Property": "SHA-256", "Value": sha256_bytes(raw)},
                ]), use_container_width=True, hide_index=True)
                with st.expander("Binary preview", expanded=False):
                    st.text_area("Whole-file bytes → bits", bytes_to_bitstring(raw[:256]), height=120, disabled=True)

        with right:
            st.markdown("#### 🖼️ Preview")
            raw = st.session_state.get("input_bytes", b"") or b""
            if not raw:
                st.info("Upload a file to start.")
                return
            _render_file_preview(raw, st.session_state.get("input_name", "Input preview"), "input")

def render_step_3_encoding() -> None:
    with st.container(border=True):
        st.markdown("## 2 Encoding")
        stored = st.session_state.get("stored_bytes", b"") or b""
        if not stored:
            st.info("Upload a file first.")
            return

        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("#### 🧬 Design")
            previous = st.session_state.get("encoding_mapping", "Simple Mapping")
            if previous not in MAPPING_OPTIONS:
                previous = "Simple Mapping"
            mapping = st.selectbox(
                "Rule Selection",
                MAPPING_OPTIONS,
                index=MAPPING_OPTIONS.index(previous),
                format_func=_display_rule,
                key="encoding_mapping_widget",
            )
            if st.button("Run DNA Encoding", type="primary", use_container_width=True, key="run_dna_encoding_v5"):
                dna, bits, meta = encode_bytes_to_dna(stored, mapping)
                _clear_after_encoding()
                st.session_state.update({
                    "encoding_mapping": mapping,
                    "dna": dna,
                    "bits": bits,
                    "codec_meta": meta,
                })

            meta = st.session_state.get("codec_meta", {}) or {}
            st.markdown("##### 📄 Encoded data properties")
            rows = [
                {"Property": "Encoded data", "Value": fmt_bytes(len(stored))},
                {"Property": "Input bits", "Value": f"{len(stored) * 8:,} bits"},
                {"Property": "Rule", "Value": _display_rule(st.session_state.get("encoding_mapping", mapping))},
            ]
            if meta.get("dna_design_mode"):
                rows.append({"Property": "Rule architecture", "Value": meta.get("dna_design_mode")})
            if meta.get("rule_block_bytes"):
                rows.append({"Property": "Rule block", "Value": f"{meta.get('rule_block_bytes')} bytes"})
            st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)
            st.download_button("Download encoded binary", stored, "encoded_data.bin", use_container_width=True)

        with right:
            st.markdown("#### 🧬 DNA output")
            dna = clean_dna(st.session_state.get("dna", ""))
            if not dna:
                st.info("Run DNA Encoding first.")
            else:
                hstats = homopolymer_stats(dna)
                st.dataframe(pd.DataFrame([
                    {"Property": "DNA design", "Value": _display_rule(st.session_state.get("encoding_mapping", "—"))},
                    {"Property": "DNA length", "Value": f"{len(dna):,} nt"},
                    {"Property": "GC content", "Value": f"{gc_content(dna):.4f}"},
                    {"Property": "Longest homopolymer", "Value": hstats.get("longest", 0)},
                    {"Property": "Homopolymer segments (2+ bases)", "Value": hstats.get("segments_ge_2", hstats.get("count", "—"))},
                    {"Property": "DNA expansion", "Value": f"{len(dna)/max(1,len(stored)):.3f} nt/byte"},
                ]), use_container_width=True, hide_index=True)
                st.text_area("DNA payload preview", _preview_seq(dna), height=DNA_PREVIEW_HEIGHT, disabled=True)
                st.download_button("Download DNA payload", dna.encode(), "dna_payload.txt", "text/plain", use_container_width=True)

        dna = clean_dna(st.session_state.get("dna", ""))
        if dna:
            render_randomness_analysis(
                dna,
                state_prefix="wholefile",
                original_bit_count=len(stored) * 8,
                default_expanded=False,
            )


def render_step_4_strand_design() -> None:
    with st.container(border=True):
        st.markdown("## 3 Strand Design")
        dna = clean_dna(st.session_state.get("dna", ""))
        if not dna:
            st.info("Run DNA Encoding first.")
            return

        st.markdown("#### 🧵 Strand Design")
        with st.expander("Strand design", expanded=not bool(st.session_state.get("strand_rows"))):
            c1, c2, c3 = st.columns(3)
            target_len = c1.number_input("Total strand length", min_value=80, max_value=250, value=125, step=1, key="strand_total_len_v5")
            index_len = c2.number_input("SI length", min_value=2, max_value=24, value=8, step=1, key="strand_si_len_v5")
            ecc = c3.selectbox("ECC", ["None"], key="ecc_mode_v5")
            fbr = st.text_input("FBR", FBR_DEFAULT, key="strand_fbr_v5")
            rbr = st.text_input("RBR", RBR_DEFAULT, key="strand_rbr_v5")
            if st.button("Run Strand Design", type="primary", use_container_width=True, key="design_strands_v5"):
                cfg = choose_auto_strand_design(
                    len(dna), len(clean_dna(fbr)), len(clean_dna(rbr)), int(index_len),
                    min_total_len=int(target_len), max_total_len=int(target_len),
                )
                rows = prepare_dna_strands(
                    dna,
                    fbr=clean_dna(fbr),
                    rbr=clean_dna(rbr),
                    index_len=int(index_len),
                    target_total_len=int(cfg["target_total_len"]),
                    add_filler=True,
                )
                _clear_after_strands()
                st.session_state.update({
                    "strand_rows": rows,
                    "strand_config": {
                        **cfg,
                        "ecc": ecc,
                        "fbr": clean_dna(fbr),
                        "rbr": clean_dna(rbr),
                        "index_len": int(index_len),
                    },
                })

        rows = st.session_state.get("strand_rows", []) or []
        if not rows:
            st.info("Run Strand Design first.")
            return

        total_nt = sum(len(clean_dna(r.get("Full strand", ""))) for r in rows)
        total_filler = sum(int(r.get("Filler length", 0) or 0) for r in rows)
        st.markdown("##### 📄 Strand statistics")
        s1, s2, s3, s4, s5 = st.columns(5)
        s1.metric("Strands", f"{len(rows):,}")
        s2.metric("Total synthesized DNA", f"{total_nt:,} nt")
        s3.metric("Filler", f"{total_filler:,} nt")
        s4.metric("SI", f"{int(st.session_state.get('strand_config', {}).get('index_len', 0))} nt")
        s5.metric("ECC", st.session_state.get("strand_config", {}).get("ecc", "None"))
        st.caption("Architecture: FBR + SI + Payload + Filler + RBR. ECC is an independent layer and is disabled in this build.")

        selected_index = int(st.number_input(
            "Strand ID",
            min_value=1,
            max_value=max(1, len(rows)),
            value=1,
            step=1,
            key="inspect_prepared_strand_v5",
        ))
        _render_segmented_strand(rows[selected_index - 1])

        st.markdown("##### 🧾 Strand table")
        table_rows = []
        for r in rows[:500]:
            table_rows.append({
                "No.": r.get("No.", "—"),
                "SI": r.get("Strand index", "—"),
                "Payload length": r.get("Payload length", "—"),
                "Filler length": r.get("Filler length", "—"),
                "Total length": r.get("Total length", "—"),
                "GC content": r.get("GC content", "—"),
                "Longest homopolymer": r.get("Longest homopolymer", "—"),
            })
        st.dataframe(pd.DataFrame(table_rows), use_container_width=True, hide_index=True, height=280)
        if len(rows) > 500:
            st.caption(f"Showing first 500 of {len(rows):,} strands.")
        st.download_button("Download designed strands", strand_rows_to_csv(rows), "designed_strands.csv", "text/csv", use_container_width=True)

        st.markdown("---")
        st.markdown("#### 🧪 Sequencing Read Simulation")
        st.caption("Coverage creates independent copies of every prepared strand. Errors are applied at the nucleotide level, and each simulated read can be exported as FASTQ.")

        q1, q2, q3 = st.columns(3)
        coverage = q1.number_input("Coverage / noisy read sets", min_value=1, max_value=200, value=10, step=1, key="seq_coverage_v6")
        error_model_display = q2.selectbox(
            "Error probability model",
            ["Error rate", "Q-score"],
            key="seq_model_v7",
        )

        # Keep user-facing labels simple while preserving the simulator's
        # canonical internal model names.
        if error_model_display == "Error rate":
            error_model = "Independent Bernoulli per nucleotide"
        else:
            error_model = "Phred Q-score driven"

        seed = q3.number_input("Random seed", min_value=0, value=11, step=1, key="seq_seed_v6")

        qscore_profile = "Fixed Q"
        q_mean, q_std, q_start, q_end, q_fixed = 30.0, 4.0, 35.0, 20.0, 30.0
        if error_model_display == "Q-score":
            qscore_profile = st.selectbox(
                "Q-score profile",
                ["Normal Q distribution", "Linear quality decay", "Fixed Q"],
                key="seq_q_profile_v6",
            )
            if qscore_profile == "Normal Q distribution":
                a, b = st.columns(2)
                q_mean = a.number_input("Mean Q", min_value=2.0, max_value=40.0, value=30.0, step=1.0, key="seq_q_mean_v6")
                q_std = b.number_input("Q standard deviation", min_value=0.0, max_value=15.0, value=4.0, step=0.5, key="seq_q_std_v6")
                st.caption("Each base receives an independently sampled Phred Q score. Substitution probability is p = 10^(-Q/10).")
            elif qscore_profile == "Linear quality decay":
                a, b = st.columns(2)
                q_start = a.number_input("Start Q", min_value=2.0, max_value=40.0, value=35.0, step=1.0, key="seq_q_start_v6")
                q_end = b.number_input("End Q", min_value=2.0, max_value=40.0, value=20.0, step=1.0, key="seq_q_end_v6")
                st.caption("Quality decreases linearly from the beginning to the end of each payload; each Q defines p = 10^(-Q/10).")
            else:
                q_fixed = st.number_input("Fixed Q", min_value=2.0, max_value=40.0, value=30.0, step=1.0, key="seq_q_fixed_v6")
                st.caption("All payload bases use the same Phred Q score; substitution probability is p = 10^(-Q/10).")

        e1, e2, e3, e4 = st.columns(4)
        sub_rate = e1.number_input(
            "Substitution rate",
            min_value=0.0,
            max_value=0.50,
            value=0.0200,
            step=0.001,
            format="%.4f",
            key="seq_sub_v6",
            disabled=(error_model_display == "Q-score"),
            help="Used only by the independent Bernoulli model. Under the Q-score-driven model, substitution probability is derived from Q at each nucleotide.",
        )
        ins_rate = e2.number_input("Insertion rate", min_value=0.0, max_value=0.20, value=0.0000, step=0.0005, format="%.4f", key="seq_ins_v6")
        del_rate = e3.number_input("Deletion rate", min_value=0.0, max_value=0.20, value=0.0000, step=0.0005, format="%.4f", key="seq_del_v6")
        dropout_rate = e4.number_input("Strand dropout rate", min_value=0.0, max_value=1.0, value=0.0000, step=0.001, format="%.4f", key="seq_dropout_v6")
        st.caption("Insertion/deletion remain independent Bernoulli events. Errors are applied to payload nucleotides; strand dropout removes a complete prepared strand before coverage copies are generated.")

        if st.button("Generate Sequencing Reads", type="primary", use_container_width=True, key="generate_reads_v6"):
            read_rows, stats = simulate_reads(
                rows,
                coverage=int(coverage),
                substitution_rate=float(sub_rate),
                insertion_rate=float(ins_rate),
                deletion_rate=float(del_rate),
                strand_dropout_rate=float(dropout_rate),
                seed=int(seed),
                error_probability_model=error_model,
                qscore_profile=qscore_profile,
                q_mean=float(q_mean),
                q_std=float(q_std),
                q_start=float(q_start),
                q_end=float(q_end),
                q_fixed=float(q_fixed),
            )
            _clear_after_reads()
            st.session_state.update({"read_rows": read_rows, "read_stats": stats})

        reads = st.session_state.get("read_rows", []) or []
        stats = st.session_state.get("read_stats", {}) or {}
        if reads:
            r1, r2, r3, r4, r5 = st.columns(5)
            r1.metric("Sequencing reads", f"{len(reads):,}")
            r2.metric("Substitution events", f"{int(stats.get('substitution_events', 0)):,}")
            r3.metric("Insertion events", f"{int(stats.get('insertion_events', 0)):,}")
            r4.metric("Deletion events", f"{int(stats.get('deletion_events', 0)):,}")
            r5.metric("Dropped strands", f"{int(stats.get('dropped_strand_count', 0)):,}")
            st.caption(
                f"Model: {stats.get('error_probability_model', '—')} · "
                f"Q profile: {stats.get('qscore_profile', '—')} · "
                f"Mean payload Q: {float(stats.get('mean_payload_q', 0.0)):.2f} · "
                f"Observed substitution {100*float(stats.get('observed_substitution_rate',0)):.4f}% · "
                f"insertion {100*float(stats.get('observed_insertion_rate',0)):.4f}% · "
                f"deletion {100*float(stats.get('observed_deletion_rate',0)):.4f}%"
            )

            max_set = max(int(r.get("Copy", 1)) for r in reads)
            inspect_set = st.selectbox("Inspect read set", list(range(1, max_set + 1)), key="inspect_read_set_v5")
            subset = [r for r in reads if int(r.get("Copy", 0)) == int(inspect_set)]
            st.markdown(f"##### Simulated read set {inspect_set}")
            st.dataframe(pd.DataFrame(subset[:500]), use_container_width=True, hide_index=True, height=260)
            if len(subset) > 500:
                st.caption(f"Showing first 500 of {len(subset):,} reads in this set.")

            st.markdown("##### 🧾 FASTQ")
            st.caption("FASTQ uses the simulated per-base Phred+33 quality string. Under the Q-score-driven model, the same Q values also determine nucleotide substitution probabilities.")
            if st.button("Generate FASTQ", type="primary", use_container_width=True, key="generate_fastq_v6"):
                fq, fq_stats = generate_fastq(reads)
                st.session_state["fastq_text"] = fq
                st.session_state["fastq_stats"] = fq_stats

            fastq_text = st.session_state.get("fastq_text", "")
            fastq_stats = st.session_state.get("fastq_stats", {}) or {}
            if fastq_text:
                f1, f2, f3, f4 = st.columns(4)
                f1.metric("FASTQ records", f"{int(fastq_stats.get('records', 0)):,}")
                f2.metric("FASTQ size", fmt_bytes(int(fastq_stats.get("bytes", 0))))
                f3.metric("Mean payload Q", f"{float(fastq_stats.get('mean_payload_q', 0.0)):.2f}")
                f4.metric("Quality encoding", "Phred+33")
                preview_lines = fastq_text.splitlines()[:24]
                st.text_area("FASTQ preview", "\n".join(preview_lines), height=220, disabled=True)
                st.download_button("Download FASTQ", fastq_text.encode("utf-8"), "simulated_reads.fastq", "text/plain", use_container_width=True)

            st.markdown("---")
            st.markdown("#### 🗳️ Consensus Reconstruction")
            consensus_method = st.selectbox(
                "Consensus method",
                ["Majority voting", "Q-score weighted consensus"],
                key="consensus_method_v7",
            )
            if consensus_method == "Majority voting":
                st.caption("At each nucleotide position, the base observed most often across reads is selected. All reads have equal voting weight.")
            else:
                st.caption("At each nucleotide position, Phred Q scores weight the evidence. High-Q base calls contribute more than low-Q calls.")
            st.caption("Reads are grouped by the designed SI. Short read-to-reference alignment is used when insertions/deletions are present.")
            if st.button("Run Consensus Reconstruction", type="primary", use_container_width=True, key="run_consensus_v7"):
                crows, cdna, cstats, evidence = consensus_by_si(rows, reads, len(dna), method=consensus_method)
                _clear_after_consensus()
                st.session_state.update({
                    "consensus_rows": crows,
                    "consensus_dna": cdna,
                    "consensus_stats": cstats,
                    "consensus_evidence": evidence,
                })

        crows = st.session_state.get("consensus_rows", []) or []
        cstats = st.session_state.get("consensus_stats", {}) or {}
        if crows:
            c1, c2, c3, c4 = st.columns(4)
            acc = cstats.get("consensus_accuracy")
            c1.metric("Consensus accuracy", f"{100*acc:.4f}%" if acc is not None else "Incomplete")
            c2.metric("Mean confidence", f"{100*float(cstats.get('mean_confidence',0)):.2f}%")
            c3.metric("Low-confidence positions", int(cstats.get("low_confidence_positions", 0)))
            c4.metric("Final DNA mismatches", cstats.get("final_dna_mismatches") if cstats.get("final_dna_mismatches") is not None else "—")
            st.markdown("##### Consensus result by strand")
            st.dataframe(pd.DataFrame(crows[:500]), use_container_width=True, hide_index=True, height=280)
            if len(crows) > 500:
                st.caption(f"Showing first 500 of {len(crows):,} consensus strands.")
            if not st.session_state.get("consensus_dna"):
                st.warning("Consensus is incomplete. With ECC=None, missing/incomplete strands cannot be reconstructed.")
            else:
                st.download_button("Download consensus DNA", st.session_state["consensus_dna"].encode(), "consensus_dna.txt", "text/plain", use_container_width=True)

            evidence = st.session_state.get("consensus_evidence", []) or []
            if evidence:
                with st.expander("Position-level nucleotide evidence", expanded=False):
                    sid = int(st.number_input("Inspect position-level evidence for strand", min_value=1, max_value=len(rows), value=1, step=1, key="evidence_sid_v5"))
                    subset = [r for r in evidence if int(r.get("Strand ID", -1)) == sid]
                    st.dataframe(pd.DataFrame(subset), use_container_width=True, hide_index=True, height=300)


def _clean_reassembled_payload() -> str:
    rows = st.session_state.get("strand_rows", []) or []
    dna = clean_dna(st.session_state.get("dna", ""))
    if not rows:
        return dna
    return clean_dna("".join(clean_dna(r.get("Payload", "")) for r in rows))[:len(dna)]


def _generic_decode_meta(mapping: str, expected_bytes: int) -> Dict[str, Any]:
    """Minimal framing metadata shared across all candidate rules in blind simulation.

    The candidate rule itself is not taken from the encoder session. The current
    simulation does know the stored byte length, which is used only to trim the
    final padded rule block and to evaluate exact recovery.
    """
    meta: Dict[str, Any] = {"bytes_len": int(expected_bytes)}
    if mapping == "Simple Mapping":
        meta["bits_len"] = int(expected_bytes) * 8
    return meta


def _try_all_rules(dna_input: str) -> Tuple[Dict[str, Any] | None, List[Dict[str, Any]]]:
    """Try SM/R0/R1/R2/R∞ without using the session's encoded-rule identity."""
    stored_ref = bytes(st.session_state.get("stored_bytes", b"") or b"")
    input_name = st.session_state.get("input_name", "recovered.bin")
    method = st.session_state.get("storage_method", "No compression")
    candidates: List[Dict[str, Any]] = []
    best: Dict[str, Any] | None = None

    for rule in MAPPING_OPTIONS:
        row: Dict[str, Any] = {
            "Rule": _display_rule(rule),
            "Decode": "Failed",
            "Exact stored bytes": "—",
            "Container": "—",
            "Bytes": 0,
            "Rule erasures": "—",
            "Result": "Rejected",
        }
        candidate: Dict[str, Any] | None = None
        try:
            meta_hint = _generic_decode_meta(rule, len(stored_ref))
            decoded, bits, meta = decode_dna_with_mapping(dna_input, rule, meta_hint)
            exact = bool(stored_ref) and decoded == stored_ref
            magic = detect_magic(decoded)
            container_ok = False
            container_note = "No recognizable file/container signature"
            if magic is not None:
                container_ok, container_note = validate_container(decoded, magic.kind)

            recovered_file, recovered_name, recovery_note = _recover_application_file(decoded, method, input_name)
            score = 0
            if exact:
                score += 1000
            if container_ok:
                score += 100
            if magic is not None:
                score += int(10 * float(magic.confidence))
            erasures = int(meta.get("corrupted_block_count", 0) or 0)
            if erasures == 0:
                score += 1

            row.update({
                "Decode": "Success",
                "Exact stored bytes": "PASS" if exact else "FAIL",
                "Container": (f"Valid {magic.kind}" if container_ok and magic else (f"Detected {magic.kind}, invalid" if magic else "Not recognized")),
                "Bytes": len(decoded),
                "Rule erasures": erasures,
                "Result": "Identified" if exact else ("Candidate" if container_ok else "Rejected"),
            })
            candidate = {
                "mapping": rule,
                "data": decoded,
                "bits": bits,
                "meta": meta,
                "exact": exact,
                "container_ok": container_ok,
                "container_note": container_note,
                "score": score,
                "recovered_file": recovered_file,
                "recovered_name": recovered_name,
                "recovery_note": recovery_note,
            }
        except Exception as exc:
            row["Note"] = str(exc)[:140]

        candidates.append(row)
        if candidate is not None:
            if best is None or candidate["score"] > best["score"]:
                best = candidate

    # With a simulation reference, require exact stored-byte recovery to call the
    # rule identified. Without a reference, accept a uniquely valid container.
    if stored_ref:
        exact_candidates = [c for c in [best] if c is not None and c.get("exact")]
        if not exact_candidates:
            # Search all exact rows defensively by re-running only if the highest score
            # was not exact (normally impossible because exact gets +1000).
            return None, candidates
        return best, candidates

    return (best if best is not None and best.get("container_ok") else None), candidates


def render_step_5_decoding() -> None:
    with st.container(border=True):
        st.markdown("## 4 Decoding")
        encoded_mapping = st.session_state.get("encoding_mapping")
        if not encoded_mapping:
            st.info("Run DNA Encoding first.")
            return

        strand_rows = st.session_state.get("strand_rows", []) or []
        clean_payload = _clean_reassembled_payload()
        consensus = clean_dna(st.session_state.get("consensus_dna", ""))

        sources: List[str] = []
        if consensus:
            sources.append("Consensus DNA")
        if strand_rows:
            sources.append("Upload FASTQ")
        if clean_payload:
            sources.append("Designed strands (clean baseline)")
        if not sources:
            st.info("Run Strand Design first.")
            return

        left, right = st.columns(2, gap="large")
        with left:
            st.markdown("#### 📁 DNA input")
            source = st.radio("DNA source", sources, horizontal=True, key="decode_source_v8")

            dna_input = ""
            source_label = source

            if source == "Consensus DNA":
                dna_input = consensus
            elif source == "Designed strands (clean baseline)":
                dna_input = clean_payload
            else:
                uploaded_fastq = st.file_uploader(
                    "Upload FASTQ",
                    type=["fastq", "fq", "txt"],
                    key="decode_fastq_upload_v1",
                    help="Upload a single-end FASTQ read set generated from the current Strand Design.",
                )
                fastq_method = st.selectbox(
                    "FASTQ consensus method",
                    ["Majority voting", "Q-score weighted consensus"],
                    key="decode_fastq_consensus_method_v1",
                )
                if fastq_method == "Majority voting":
                    st.caption("Reads are grouped by SI and the most frequent nucleotide is selected at each payload position.")
                else:
                    st.caption("Reads are grouped by SI and uploaded Phred+33 Q-scores weight the nucleotide evidence at each payload position.")

                # A consensus reconstructed with one method must not silently be reused
                # after the user changes the method. Keep parsed reads, clear only the
                # consensus/decode products and require reconstruction again.
                previous_fastq_method = st.session_state.get("uploaded_fastq_consensus_method")
                if previous_fastq_method is not None and previous_fastq_method != fastq_method:
                    _clear([
                        "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna",
                        "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence",
                        "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error",
                        "recovered_file_bytes", "recovered_file_name", "restored_info",
                        "decode_source_label", "detected_rule", "auto_decode_table", "rule_resolution",
                    ])
                st.session_state["uploaded_fastq_consensus_method"] = fastq_method

                if uploaded_fastq is not None:
                    raw_fastq = uploaded_fastq.getvalue()
                    fq_hash = hashlib.sha256(raw_fastq).hexdigest()
                    if st.session_state.get("uploaded_fastq_hash") != fq_hash:
                        _clear([
                            "uploaded_fastq_read_rows", "uploaded_fastq_parse_stats",
                            "uploaded_fastq_consensus_rows", "uploaded_fastq_consensus_dna",
                            "uploaded_fastq_consensus_stats", "uploaded_fastq_consensus_evidence", "uploaded_fastq_consensus_method",
                            "decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error",
                            "recovered_file_bytes", "recovered_file_name", "restored_info",
                            "decode_source_label", "detected_rule", "auto_decode_table",
                        ])
                        st.session_state["uploaded_fastq_hash"] = fq_hash
                        st.session_state["uploaded_fastq_name"] = uploaded_fastq.name

                    if st.button(
                        "Reconstruct DNA from FASTQ",
                        type="primary",
                        use_container_width=True,
                        key="reconstruct_uploaded_fastq_v1",
                    ):
                        try:
                            fastq_text = raw_fastq.decode("utf-8-sig")
                        except UnicodeDecodeError as exc:
                            st.session_state["decode_error"] = f"FASTQ must be a text file encoded as UTF-8/ASCII: {exc}"
                        else:
                            try:
                                uploaded_reads, parse_stats = parse_fastq_reads(fastq_text, strand_rows)
                                if not uploaded_reads:
                                    raise ValueError("No FASTQ reads could be matched to the current Strand Design/SI table.")
                                original_dna_len = len(clean_dna(st.session_state.get("dna", "")))
                                crows, cdna, cstats, evidence = consensus_by_si(
                                    strand_rows,
                                    uploaded_reads,
                                    original_dna_len,
                                    method=fastq_method,
                                )
                                _clear_after_uploaded_fastq_consensus()
                                st.session_state.update({
                                    "uploaded_fastq_read_rows": uploaded_reads,
                                    "uploaded_fastq_parse_stats": parse_stats,
                                    "uploaded_fastq_consensus_rows": crows,
                                    "uploaded_fastq_consensus_dna": cdna,
                                    "uploaded_fastq_consensus_stats": cstats,
                                    "uploaded_fastq_consensus_evidence": evidence,
                                    "uploaded_fastq_consensus_method": fastq_method,
                                    "decode_error": "",
                                })
                            except Exception as exc:
                                st.session_state["decode_error"] = str(exc)

                parse_stats = st.session_state.get("uploaded_fastq_parse_stats", {}) or {}
                fastq_cstats = st.session_state.get("uploaded_fastq_consensus_stats", {}) or {}
                dna_input = clean_dna(st.session_state.get("uploaded_fastq_consensus_dna", ""))
                source_label = f"Uploaded FASTQ · {fastq_method}"

                if parse_stats:
                    f1, f2, f3, f4 = st.columns(4)
                    f1.metric("FASTQ records", f"{int(parse_stats.get('records', 0)):,}")
                    f2.metric("Matched reads", f"{int(parse_stats.get('matched_reads', 0)):,}")
                    f3.metric("Recognized strands", f"{int(parse_stats.get('recognized_strands', 0)):,}")
                    f4.metric("Mean Q", f"{float(parse_stats.get('mean_q', 0.0)):.2f}")
                    fastq_design_id = str(parse_stats.get("fastq_design_id", ""))
                    if fastq_design_id and fastq_design_id != "not embedded":
                        st.caption(f"FASTQ design fingerprint: `{fastq_design_id}` · matches current Strand Design")
                    elif fastq_design_id == "not embedded":
                        st.caption("FASTQ has no embedded design fingerprint (legacy/external FASTQ); compatibility is checked from SI/strand structure only.")
                    if int(parse_stats.get("unmatched_reads", 0)) > 0:
                        st.warning(f"{int(parse_stats.get('unmatched_reads', 0)):,} FASTQ reads could not be matched to the current SI/strand design.")
                    if int(parse_stats.get("structural_mismatch_reads", 0)) > 0:
                        st.warning(
                            f"{int(parse_stats.get('structural_mismatch_reads', 0)):,} matched reads contain a mismatch in FBR/SI/filler/RBR. "
                            "Payload extraction still uses the current strand boundaries."
                        )

                if fastq_cstats:
                    acc = fastq_cstats.get("consensus_accuracy")
                    c1, c2, c3 = st.columns(3)
                    c1.metric("FASTQ consensus", f"{100*acc:.4f}%" if acc is not None else "Incomplete")
                    c2.metric("Missing strands", int(fastq_cstats.get("missing_strands", 0)))
                    c3.metric("Low-confidence positions", int(fastq_cstats.get("low_confidence_positions", 0)))
                    if not dna_input:
                        st.warning("FASTQ consensus is incomplete. With ECC=None, missing/incomplete strands cannot be recovered.")

            rule_mode = st.radio(
                "Rule knowledge",
                ["Known rule", "Unknown rule — auto-detect"],
                horizontal=True,
                key="decode_rule_mode_v8",
            )

            selected_mapping = encoded_mapping
            if rule_mode == "Known rule":
                selected_mapping = st.selectbox(
                    "Rule",
                    MAPPING_OPTIONS,
                    index=MAPPING_OPTIONS.index(encoded_mapping) if encoded_mapping in MAPPING_OPTIONS else 0,
                    format_func=_display_rule,
                    key="known_decode_rule_v8",
                )
                st.caption("Select the rule used for decoding.")
            else:
                selected_mapping = None
                st.caption("The app automatically tests SM, R0, R1, R2 and R∞, then selects the rule that successfully reconstructs the stored file/container.")

            st.dataframe(pd.DataFrame([
                {"Property": "Rule knowledge", "Value": rule_mode},
                {"Property": "Rule", "Value": _display_rule(selected_mapping) if rule_mode == "Known rule" else "Auto-detect"},
                {"Property": "Input DNA", "Value": source},
                {"Property": "Input length", "Value": f"{len(dna_input):,} nt" if dna_input else "Not reconstructed"},
                {"Property": "ECC correction", "Value": "None"},
            ]), use_container_width=True, hide_index=True)
            st.text_area(
                "Input DNA payload preview",
                _preview_seq(dna_input) if dna_input else "Upload/process FASTQ to reconstruct the DNA payload.",
                height=120,
                disabled=True,
            )

            button_label = "Decode Recovered File" if rule_mode == "Known rule" else "Auto-detect Rule & Decode"
            if st.button(
                button_label,
                type="primary",
                use_container_width=True,
                key="decode_recovered_v8",
                disabled=not bool(dna_input),
            ):
                _clear(["decoded_stored_bytes", "decoded_bits", "decoded_meta", "decode_error", "recovered_file_bytes",
                        "recovered_file_name", "restored_info", "decode_source_label", "detected_rule", "auto_decode_table"])
                try:
                    if rule_mode == "Unknown rule — auto-detect":
                        candidate, candidate_rows = _try_all_rules(dna_input)
                        st.session_state["auto_decode_table"] = candidate_rows
                        if candidate is None:
                            raise ValueError("Auto-detection failed: none of the five rules reconstructed the stored file exactly.")
                        decoded = candidate["data"]
                        bits = candidate["bits"]
                        meta = candidate["meta"]
                        selected_mapping = candidate["mapping"]
                        recovered_file = candidate["recovered_file"]
                        recovered_name = candidate["recovered_name"]
                        recovery_note = candidate["recovery_note"]
                        st.session_state["detected_rule"] = selected_mapping
                        st.session_state["rule_resolution"] = "Auto-detected"
                    else:
                        decode_meta = (
                            st.session_state.get("codec_meta", {}) or {}
                            if selected_mapping == encoded_mapping
                            else _generic_decode_meta(selected_mapping, len(st.session_state.get("stored_bytes", b"") or b""))
                        )
                        decoded, bits, meta = decode_dna_with_mapping(dna_input, selected_mapping, decode_meta)
                        input_name = st.session_state.get("input_name", "recovered.bin")
                        method = st.session_state.get("storage_method", "No compression")
                        recovered_file, recovered_name, recovery_note = _recover_application_file(decoded, method, input_name)
                        st.session_state["detected_rule"] = selected_mapping
                        st.session_state["rule_resolution"] = "Selected"

                    out_dir = WORK_ROOT / "decoded_output"
                    out_dir.mkdir(parents=True, exist_ok=True)
                    recovered_path = out_dir / os.path.basename(recovered_name)
                    recovered_path.write_bytes(recovered_file)
                    restored_info = write_restored_file(decoded, str(out_dir), preferred_name="decoded_stored")
                    restored_info.update({
                        "application_file_path": str(recovered_path),
                        "application_file_name": recovered_name,
                        "recovery_note": recovery_note,
                    })
                    st.session_state.update({
                        "decoded_stored_bytes": decoded,
                        "decoded_bits": bits,
                        "decoded_meta": meta,
                        "recovered_file_bytes": recovered_file,
                        "recovered_file_name": recovered_name,
                        "restored_info": restored_info,
                        "decode_source_label": source_label,
                        "decode_error": "",
                    })
                except Exception as exc:
                    st.session_state["decode_error"] = str(exc)

            auto_rows = st.session_state.get("auto_decode_table", []) or []
            if rule_mode == "Unknown rule — auto-detect" and auto_rows:
                detected = st.session_state.get("detected_rule")
                if detected:
                    st.success(f"Auto-detected rule: {_display_rule(detected)}")
                with st.expander("Show auto-detection diagnostics"):
                    st.dataframe(pd.DataFrame(auto_rows), use_container_width=True, hide_index=True)

            if st.session_state.get("decode_error"):
                st.error(st.session_state["decode_error"])

        with right:
            st.markdown("#### 🖼️ Decoded output")
            decoded = st.session_state.get("decoded_stored_bytes")
            if decoded is None:
                st.info("Run Decode first.")
                return
            stored = st.session_state.get("stored_bytes", b"") or b""
            recovered_file = st.session_state.get("recovered_file_bytes", b"") or b""
            recovered_name = st.session_state.get("recovered_file_name", "recovered.bin")
            meta = st.session_state.get("decoded_meta", {}) or {}
            acc = _byte_accuracy(stored, decoded)
            exact = stored == decoded
            detected_rule = st.session_state.get("detected_rule")
            rule_resolution = st.session_state.get("rule_resolution", "Selected")
            c1, c2, c3, c4 = st.columns(4)
            c1.metric("Stored-byte accuracy", f"{100*acc:.6f}%")
            c2.metric("Stored SHA-256", "PASS" if exact else "FAIL")
            c3.metric("Auto-detected rule" if rule_resolution == "Auto-detected" else "Selected rule", _display_rule(detected_rule) if detected_rule else "—")
            c4.metric("Rule erasures", int(meta.get("corrupted_block_count", 0)))
            if str(st.session_state.get("decode_source_label", "")).startswith("Uploaded FASTQ") and not exact:
                st.warning(
                    "The uploaded FASTQ did not reconstruct the stored byte stream exactly. "
                    "Check FASTQ consensus accuracy, recognized/missing strands, the current Strand Design, "
                    "and whether this FASTQ was generated from the current encoded payload."
                )

            m = detect_magic(decoded)
            valid_note = "No recognizable container signature"
            if m:
                valid, valid_note = validate_container(decoded, m.kind)
                valid_note = ("Valid: " if valid else "Invalid: ") + valid_note
            st.caption(valid_note)
            st.dataframe(pd.DataFrame([
                {"Property": "Decoded stored bytes", "Value": fmt_bytes(len(decoded))},
                {"Property": "Recovered file", "Value": recovered_name},
                {"Property": "Recovered file size", "Value": fmt_bytes(len(recovered_file))},
                {"Property": "Storage method", "Value": st.session_state.get("storage_method", "—")},
                {"Property": "DNA source", "Value": st.session_state.get("decode_source_label", "—")},
                {"Property": "Decoded rule", "Value": _display_rule(detected_rule) if detected_rule else "—"},
            ]), use_container_width=True, hide_index=True)
            _render_file_preview(recovered_file, recovered_name, "decoded")
            d1, d2 = st.columns(2)
            d1.download_button("Download recovered file", recovered_file, file_name=recovered_name, use_container_width=True)
            d2.download_button("Download decoded stored container", decoded, file_name=f"decoded_stored{magic_dict(decoded).get('ext', '.bin')}", use_container_width=True)

def _summary_payload() -> Dict[str, Any]:
    raw = st.session_state.get("input_bytes", b"") or b""
    stored = st.session_state.get("stored_bytes", b"") or b""
    decoded = st.session_state.get("decoded_stored_bytes", b"") or b""
    recovered = st.session_state.get("recovered_file_bytes", b"") or b""
    dna = clean_dna(st.session_state.get("dna", ""))
    strands = st.session_state.get("strand_rows", []) or []
    read_stats = st.session_state.get("read_stats", {}) or {}
    consensus_stats = st.session_state.get("consensus_stats", {}) or {}
    rand = st.session_state.get("wholefile_randomness_report", {}) or {}
    quality = rand.get("sequence_quality", {}) or {}
    ising = rand.get("inverse_ising", {}) or {}
    logic = rand.get("logic_3input_1output", {}) or {}
    return {
        "input": {
            "name": st.session_state.get("input_name"),
            "bytes": len(raw),
            "sha256": sha256_bytes(raw) if raw else None,
        },
        "compression": {
            "method": st.session_state.get("storage_method"),
            "stored_bytes": len(stored),
            "compression_ratio": (len(raw) / max(1, len(stored))) if raw else None,
            "lossy": bool((st.session_state.get("storage_meta", {}) or {}).get("lossy")),
        },
        "encoding": {
            "rule": _display_rule(st.session_state.get("encoding_mapping", "—")),
            "dna_length_nt": len(dna),
            "gc_content": gc_content(dna) if dna else None,
            "longest_homopolymer": homopolymer_stats(dna).get("longest", 0) if dna else None,
        },
        "randomness": {
            "base_entropy_normalized": quality.get("base_entropy_normalized"),
            "gc_deviation_from_0_5": quality.get("gc_deviation_from_0_5"),
            "abs_polarization_gamma": ising.get("abs_polarization_gamma"),
            "alpha_rmse_from_0_5": logic.get("alpha_rmse_from_0_5"),
            "beta_rmse_from_0_125": logic.get("beta_rmse_from_0_125"),
        },
        "strand_design": {
            "strand_count": len(strands),
            "total_synthesized_nt": sum(len(clean_dna(r.get("Full strand", ""))) for r in strands),
            "ecc": (st.session_state.get("strand_config", {}) or {}).get("ecc", "None"),
        },
        "sequencing": read_stats,
        "consensus": consensus_stats,
        "recovery": {
            "decode_source": st.session_state.get("decode_source_label"),
            "detected_rule": _display_rule(st.session_state.get("detected_rule", "—")),
            "decoded_stored_bytes": len(decoded),
            "stored_byte_accuracy": _byte_accuracy(stored, decoded) if decoded or stored else None,
            "stored_sha256_match": bool(stored == decoded) if decoded or stored else None,
            "recovered_file_bytes": len(recovered),
            "whole_file_match": bool(raw == recovered) if raw or recovered else None,
        },
    }


def render_step_6_summary() -> None:
    with st.container(border=True):
        st.markdown("## 5 Summarization")
        decoded = st.session_state.get("decoded_stored_bytes")
        if decoded is None:
            st.info("Run Decode first.")
            return

        raw = st.session_state.get("input_bytes", b"") or b""
        stored = st.session_state.get("stored_bytes", b"") or b""
        recovered = st.session_state.get("recovered_file_bytes", b"") or b""
        recovered_name = st.session_state.get("recovered_file_name", "recovered.bin")

        st.markdown("#### 📊 Summary")
        c1, c2, c3 = st.columns(3, gap="large")
        with c1:
            st.markdown("##### Original")
            _render_file_preview(raw, st.session_state.get("input_name", "Original"), "summary_original")
            st.dataframe(pd.DataFrame([
                {"Property": "File", "Value": st.session_state.get("input_name", "—")},
                {"Property": "Size", "Value": fmt_bytes(len(raw))},
                {"Property": "SHA-256", "Value": sha256_bytes(raw) if raw else "—"},
            ]), use_container_width=True, hide_index=True)
        with c2:
            st.markdown("##### Stored (No compression)")
            _render_file_preview(stored, f"stored{magic_dict(stored).get('ext', '.bin')}", "summary_stored")
            st.dataframe(pd.DataFrame([
                {"Property": "Method", "Value": st.session_state.get("storage_method", "—")},
                {"Property": "Size", "Value": fmt_bytes(len(stored))},
                {"Property": "Ratio", "Value": f"{len(raw)/max(1,len(stored)):.3f}×" if raw else "—"},
            ]), use_container_width=True, hide_index=True)
        with c3:
            st.markdown("##### Decoded")
            _render_file_preview(recovered, recovered_name, "summary_decoded")
            st.dataframe(pd.DataFrame([
                {"Property": "File", "Value": recovered_name},
                {"Property": "Size", "Value": fmt_bytes(len(recovered))},
                {"Property": "Whole-file match", "Value": "PASS" if recovered == raw else "DIFF"},
            ]), use_container_width=True, hide_index=True)

        summary = _summary_payload()
        st.markdown("##### Storage analysis")
        st.dataframe(pd.DataFrame([
            {"Metric": "Original size", "Value": fmt_bytes(summary["input"]["bytes"])},
            {"Metric": "Stored size", "Value": fmt_bytes(summary["compression"]["stored_bytes"])},
            {"Metric": "Storage ratio (input/stored)", "Value": f"{summary['compression']['compression_ratio']:.3f}×" if summary["compression"]["compression_ratio"] is not None else "—"},
            {"Metric": "Method", "Value": summary["compression"]["method"]},
            {"Metric": "Lossy", "Value": "Yes" if summary["compression"]["lossy"] else "No"},
        ]), use_container_width=True, hide_index=True)

        st.markdown("##### Encode–strand analysis")
        st.dataframe(pd.DataFrame([
            {"Metric": "Mapping rule", "Value": summary["encoding"]["rule"]},
            {"Metric": "DNA length", "Value": f"{summary['encoding']['dna_length_nt']:,} nt"},
            {"Metric": "GC content", "Value": f"{summary['encoding']['gc_content']:.4f}" if summary["encoding"]["gc_content"] is not None else "—"},
            {"Metric": "Longest homopolymer", "Value": summary["encoding"]["longest_homopolymer"]},
            {"Metric": "Strand count", "Value": f"{summary['strand_design']['strand_count']:,}"},
            {"Metric": "Total synthesized DNA", "Value": f"{summary['strand_design']['total_synthesized_nt']:,} nt"},
            {"Metric": "ECC", "Value": summary["strand_design"]["ecc"]},
        ]), use_container_width=True, hide_index=True)

        rstats = summary.get("sequencing", {}) or {}
        cstats = summary.get("consensus", {}) or {}
        st.markdown("##### Sequencing / consensus analysis")
        st.dataframe(pd.DataFrame([
            {"Metric": "Coverage", "Value": rstats.get("coverage", "—")},
            {"Metric": "Reads", "Value": rstats.get("read_count", "—")},
            {"Metric": "Error probability model", "Value": ("Error rate" if rstats.get("error_probability_model") == "Independent Bernoulli per nucleotide" else "Q-score" if rstats.get("error_probability_model") == "Phred Q-score driven" else rstats.get("error_probability_model", "—"))},
            {"Metric": "Q-score profile", "Value": rstats.get("qscore_profile", "—")},
            {"Metric": "Mean payload Q", "Value": f"{float(rstats.get('mean_payload_q', 0.0)):.2f}" if rstats else "—"},
            {"Metric": "Substitution events", "Value": rstats.get("substitution_events", "—")},
            {"Metric": "Insertion events", "Value": rstats.get("insertion_events", "—")},
            {"Metric": "Deletion events", "Value": rstats.get("deletion_events", "—")},
            {"Metric": "Dropped strands", "Value": rstats.get("dropped_strand_count", "—")},
            {"Metric": "Consensus method", "Value": cstats.get("method", "—")},
            {"Metric": "Consensus accuracy", "Value": f"{100*float(cstats.get('consensus_accuracy')):.4f}%" if cstats.get("consensus_accuracy") is not None else "—"},
            {"Metric": "Final DNA mismatches", "Value": cstats.get("final_dna_mismatches", "—")},
        ]), use_container_width=True, hide_index=True)

        rand = summary.get("randomness", {}) or {}
        if any(v is not None for v in rand.values()):
            st.markdown("##### 🎲 Randomness summary")
            st.dataframe(pd.DataFrame([
                {"Metric": "Base entropy normalized", "Value": rand.get("base_entropy_normalized")},
                {"Metric": "GC deviation from 0.5", "Value": rand.get("gc_deviation_from_0_5")},
                {"Metric": "|Polarization γ|", "Value": rand.get("abs_polarization_gamma")},
                {"Metric": "Alpha RMSE from 0.5", "Value": rand.get("alpha_rmse_from_0_5")},
                {"Metric": "Beta RMSE from 0.125", "Value": rand.get("beta_rmse_from_0_125")},
            ]), use_container_width=True, hide_index=True)
        else:
            st.caption("Randomness analysis was not run in Step 3.")

        st.markdown("##### Recovery Quality Report")
        rec = summary["recovery"]
        st.dataframe(pd.DataFrame([
            {"Metric": "Decode source", "Value": rec.get("decode_source")},
            {"Metric": "Detected / selected rule", "Value": rec.get("detected_rule")},
            {"Metric": "Stored-byte accuracy", "Value": f"{100*float(rec.get('stored_byte_accuracy',0)):.6f}%"},
            {"Metric": "Stored SHA-256", "Value": "PASS" if rec.get("stored_sha256_match") else "FAIL"},
            {"Metric": "Recovered whole-file equality", "Value": "PASS" if rec.get("whole_file_match") else "DIFF"},
        ]), use_container_width=True, hide_index=True)

        st.download_button(
            "Download pipeline summary (JSON)",
            json.dumps(summary, indent=2, ensure_ascii=False, default=str).encode("utf-8"),
            "dna_storage_summary.json",
            "application/json",
            use_container_width=True,
        )


def _apply_style() -> None:
    st.markdown(
        """
<style>
.block-container {max-width: 1320px; padding-top: 1.2rem; padding-bottom: 3rem;}
h1, h2, h3, h4 {letter-spacing: -0.02em;}
div[data-testid="stMetric"] {border: 1px solid rgba(120,120,120,.16); border-radius: 14px; padding: .55rem .65rem;}
div[data-testid="stVerticalBlockBorderWrapper"] {border-radius: 16px;}
textarea, code, pre {font-family: ui-monospace, SFMono-Regular, Menlo, Consolas, monospace !important;}
</style>
""",
        unsafe_allow_html=True,
    )


def render_app() -> None:
    st.set_page_config(page_title="DNA Data Storage System — No Compression", page_icon="🧬", layout="wide")
    _apply_style()
    st.title("DNA Data Storage System")
    st.caption("Whole file/container → rule-based DNA → strand design → sequencing/FASTQ/consensus → recovered file")
    render_step_1_input()
    render_step_3_encoding()
    render_step_4_strand_design()
    render_step_5_decoding()
    render_step_6_summary()
