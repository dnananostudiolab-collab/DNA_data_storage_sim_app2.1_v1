from __future__ import annotations

import gzip
import io
import lzma
import bz2
import zipfile
from typing import Any, Dict, List, Optional, Tuple

import pandas as pd

try:
    from PIL import Image
except Exception:
    Image = None

import dna_codec
from utils_bits_v2 import bytes_to_bitstring, bitstring_to_bytes, detect_magic
from config import MAPPING_OPTIONS, IMAGE_KINDS


_RULE_INIT_DIMER = "TA"


def mapping_to_config(mapping_name: str) -> Dict[str, Any]:
    if mapping_name == "Simple Mapping":
        return {
            "mode": "SIMPLE",
            "scheme_name": "RINF_B16",
            "init_dimer": "TA",
            "whiten": False,
        }
    return {
        "mode": "TABLE",
        "scheme_name": mapping_name,
        "init_dimer": "TA",
        "whiten": False,
    }

def encode_bytes_to_dna(data: bytes, mapping_name: str) -> Tuple[str, str, Dict[str, Any]]:
    """Encode the exact stored byte stream into DNA.

    SM remains the direct 2-bit mapping baseline. R0/R1/R2/R∞ use the
    canonical rulebase-v4 implementation in rulebase_codec_v4.py.
    Compression is intentionally outside this function.
    """
    bits = bytes_to_bitstring(data)
    if mapping_name == "Simple Mapping":
        payload_bits = bits if bits else "0"
        dna, digits = dna_codec.simple_encode_bits_to_dna(payload_bits)
        meta = {
            "mapping": mapping_name,
            "dna_design_name": mapping_name,
            "dna_design_mode": "simple",
            "bits_len": len(bits),
            "bytes_len": len(data),
            "digits_len": len(digits),
            "dna_design_codec_version": "simple-v1",
        }
        return dna, bits, meta

    from rulebase_codec_v4 import encode_rulebase
    dna, meta = encode_rulebase(data, mapping_name)
    meta = dict(meta)
    meta.update({
        "mapping": mapping_name,
        "dna_design_name": mapping_name,
        "bits_len": len(bits),
        "bytes_len": len(data),
    })
    return dna, bits, meta


def decode_dna_with_mapping(
    dna: str,
    mapping_name: str,
    codec_meta: Optional[Dict[str, Any]] = None,
) -> Tuple[bytes, str, Dict[str, Any]]:
    """Decode DNA back to the stored byte stream using the selected design."""
    codec_meta = dict(codec_meta or {})

    if mapping_name == "Simple Mapping":
        decoded_bits, digits = dna_codec.simple_decode_dna_to_bits(dna)
        expected_bits = int(codec_meta.get("bits_len", len(decoded_bits)))
        decoded_bits = decoded_bits[:expected_bits] if expected_bits >= 0 else decoded_bits
        data, pad_bits = bitstring_to_bytes(decoded_bits, pad_to_byte=True)
        expected_bytes = codec_meta.get("bytes_len")
        if expected_bytes is not None:
            data = data[:int(expected_bytes)]
        meta = {
            "mapping": mapping_name,
            "dna_design_name": mapping_name,
            "dna_design_mode": "simple",
            "bits_len": len(decoded_bits),
            "bytes_len": len(data),
            "pad_bits_to_byte": pad_bits,
            "digits_len": len(digits),
            "dna_design_codec_version": "simple-v1",
        }
        return data, decoded_bits, meta

    if mapping_name in {"R0_B9", "R1_B12", "R2_B15", "RINF_B16"}:
        from rulebase_codec_v4 import decode_rulebase
        data, meta = decode_rulebase(dna, mapping_name, codec_meta)
        bits = bytes_to_bitstring(data)
        meta = dict(meta)
        meta.update({
            "mapping": mapping_name,
            "dna_design_name": mapping_name,
            "bits_len": len(bits),
            "bytes_len": len(data),
        })
        return data, bits, meta

    raise ValueError(f"Unsupported DNA design: {mapping_name}")


def decode_dna_with_design(
    dna: str,
    dna_design_name: str,
    dna_design_meta: Optional[Dict[str, Any]] = None,
) -> Tuple[bytes, str, Dict[str, Any]]:
    """Canonical alias used by rulebase regression tests and newer UI code."""
    return decode_dna_with_mapping(dna, dna_design_name, dna_design_meta)


def validate_container(data: bytes, magic_kind: str) -> Tuple[bool, str]:
    """Lightweight validation beyond magic signature."""
    try:
        if magic_kind == "zip" or magic_kind in {"docx", "pptx", "xlsx", "epub"}:
            with zipfile.ZipFile(io.BytesIO(data), "r") as zf:
                bad = zf.testzip()
                if bad is not None:
                    return False, f"ZIP test failed at {bad}"
            return True, "ZIP container opened successfully"
        if magic_kind == "gzip":
            gzip.decompress(data)
            return True, "GZIP decompressed successfully"
        if magic_kind == "xz":
            lzma.decompress(data, format=lzma.FORMAT_XZ)
            return True, "XZ decompressed successfully"
        if magic_kind == "bz2":
            bz2.decompress(data)
            return True, "BZ2 decompressed successfully"
        if magic_kind in IMAGE_KINDS and Image is not None:
            img = Image.open(io.BytesIO(data))
            img.verify()
            return True, "Image verified successfully"
        return True, "Magic signature accepted"
    except Exception as e:
        return False, str(e)

def blind_decode_dna(dna_text: str) -> Tuple[Dict[str, Any], pd.DataFrame]:
    """
    Decode DNA by trying available mappings.

    Important guard: New Design decoding can be intentionally expensive because it
    enumerates block candidates for repair. If a non-New-Design mapping already
    produces a verified self-describing file/container, we stop before trying New
    Design. This avoids treating ordinary mapped DNA as a repair-coded New Design
    stream and prevents long/hanging auto-detection.
    """
    dna = dna_codec.clean_dna_text(dna_text)
    rows: List[Dict[str, Any]] = []
    best: Optional[Dict[str, Any]] = None

    for mapping in MAPPING_OPTIONS:
        if mapping == "Toolkit RS Baseline":
            rows.append({
                "Mapping": mapping,
                "Status": "Manual only",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": 0.0,
                "Note": "Toolkit RS Baseline needs codec metadata such as original file size and data/parity columns; use Manual mapping after encoding in the same session.",
            })
            continue

        # Fast-stop before the expensive New Design candidate search if a normal
        # mapping has already produced a valid file/container.
        if mapping == "New Design" and best is not None and best.get("row", {}).get("Status") == "Valid":
            rows.append({
                "Mapping": mapping,
                "Status": "Skipped",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": 0.0,
                "Note": "Skipped because an earlier mapping already produced a valid self-describing file.",
            })
            break

        row: Dict[str, Any] = {"Mapping": mapping}
        try:
            data, bits, meta = decode_dna_with_mapping(dna, mapping)
            m = detect_magic(data)
            score = 0.0
            valid = False
            if data:
                score += 1.0
            if m:
                score += 10.0 * float(m.confidence)
                ok, note = validate_container(data, m.kind)
                valid = bool(ok)
                if ok:
                    score += 5.0
                else:
                    score -= 2.0
                row.update({
                    "Status": "Valid" if ok else "Weak",
                    "Magic": m.kind,
                    "Ext": m.ext,
                    "Confidence": m.confidence,
                    "Bytes": len(data),
                    "Score": score,
                    "Note": note,
                })
            else:
                row.update({
                    "Status": "No magic",
                    "Magic": "—",
                    "Ext": "—",
                    "Confidence": 0.0,
                    "Bytes": len(data),
                    "Score": score,
                    "Note": "No recognizable file/container signature",
                })

            candidate = {
                "mapping": mapping,
                "data": data,
                "bits": bits,
                "meta": meta,
                "magic": m,
                "score": score,
                "row": row,
            }
            if m is not None and (best is None or candidate["score"] > best["score"]):
                best = candidate


        except Exception as e:
            row.update({
                "Status": "Failed",
                "Magic": "—",
                "Ext": "—",
                "Confidence": 0.0,
                "Bytes": 0,
                "Score": -1.0,
                "Note": str(e)[:160],
            })
        rows.append(row)

    df = pd.DataFrame(rows)
    if best is None:
        raise ValueError("Auto-detection failed: no mapping produced a recognizable self-describing byte stream.")
    return best, df
