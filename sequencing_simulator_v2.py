from __future__ import annotations

import math
import random
import hashlib
from collections import Counter, defaultdict
from difflib import SequenceMatcher
from typing import Any, Dict, List, Sequence, Tuple

from fragments import clean_dna

BASES = "ACGT"
Q_MIN = 2
Q_MAX = 40
STRUCTURAL_Q = 40


def _strand_design_fingerprint(strand_rows: Sequence[Dict[str, Any]]) -> str:
    """Stable short fingerprint for the exact prepared-strand design/payload."""
    parts: List[str] = []
    for fallback, row in enumerate(strand_rows, start=1):
        parts.append("|".join([
            str(int(row.get("No.", fallback))),
            clean_dna(row.get("FBR", "")),
            clean_dna(row.get("Strand index", "")),
            clean_dna(row.get("Payload", "")),
            clean_dna(row.get("Filler", "")),
            clean_dna(row.get("RBR", "")),
        ]))
    return hashlib.sha256("\n".join(parts).encode("ascii")).hexdigest()[:16]


def _align_read_to_reference(reference: str, read: str) -> List[str | None]:
    """Map a short noisy read onto reference coordinates.

    Inserted read bases are ignored for position-wise consensus; deleted reference
    positions become None. This is intended for the simulator's short payload
    segments and low/moderate indel rates, not as a production aligner.
    """
    reference = clean_dna(reference)
    read = clean_dna(read)
    mapped: List[str | None] = [None] * len(reference)
    sm = SequenceMatcher(None, reference, read, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag == "equal":
            for k in range(i2 - i1):
                mapped[i1 + k] = read[j1 + k]
        elif tag == "replace":
            n = min(i2 - i1, j2 - j1)
            for k in range(n):
                mapped[i1 + k] = read[j1 + k]
        elif tag == "delete":
            pass
        elif tag == "insert":
            pass
    return mapped



def _align_read_and_q_to_reference(
    reference: str,
    read: str,
    quality_string: str,
) -> Tuple[List[str | None], List[int | None]]:
    """Map a noisy payload and its Phred+33 qualities to reference coordinates."""
    reference = clean_dna(reference)
    read = clean_dna(read)
    if len(quality_string) != len(read):
        quality_string = "I" * len(read)
    qscores = [max(0, ord(ch) - 33) for ch in quality_string]

    mapped_base: List[str | None] = [None] * len(reference)
    mapped_q: List[int | None] = [None] * len(reference)
    sm = SequenceMatcher(None, reference, read, autojunk=False)
    for tag, i1, i2, j1, j2 in sm.get_opcodes():
        if tag in {"equal", "replace"}:
            n = min(i2 - i1, j2 - j1)
            for k in range(n):
                mapped_base[i1 + k] = read[j1 + k]
                mapped_q[i1 + k] = qscores[j1 + k]
        # Insertions have no reference coordinate; deletions have no observed call.
    return mapped_base, mapped_q


def _q_weighted_base_call(observations: Sequence[Tuple[str, int]]) -> Tuple[str, float, Dict[str, float]]:
    """Maximum-likelihood consensus call using per-base Phred Q scores.

    For a candidate true base b, an observed matching call contributes log(1-p),
    while a mismatching call contributes log(p/3), where p=10^(-Q/10).
    """
    if not observations:
        return "N", 0.0, {b: float("-inf") for b in BASES}

    log_scores: Dict[str, float] = {b: 0.0 for b in BASES}
    for called, q in observations:
        p_err = min(0.999999, max(1e-12, _q_to_error_probability(q)))
        for candidate in BASES:
            prob = (1.0 - p_err) if called == candidate else (p_err / 3.0)
            log_scores[candidate] += math.log(max(prob, 1e-300))

    best_score = max(log_scores.values())
    winners = [b for b in BASES if abs(log_scores[b] - best_score) < 1e-12]
    chosen = winners[0]

    # Convert log-likelihoods to a normalized confidence for display.
    exp_scores = {b: math.exp(log_scores[b] - best_score) for b in BASES}
    denom = sum(exp_scores.values())
    confidence = exp_scores[chosen] / denom if denom else 0.0
    return chosen, confidence, log_scores


def _clamp_q(q: float) -> int:
    return int(max(Q_MIN, min(Q_MAX, round(float(q)))))


def _q_to_error_probability(q: int | float) -> float:
    return 10.0 ** (-float(q) / 10.0)


def _p_to_q(p: float) -> int:
    p = float(p)
    if p <= 0:
        return Q_MAX
    return _clamp_q(-10.0 * math.log10(min(1.0, p)))


def _quality_string(qscores: Sequence[int]) -> str:
    # Sanger FASTQ / Phred+33. Q is clamped to 2..40 in this simulator.
    return "".join(chr(int(q) + 33) for q in qscores)


def _make_qscores(
    length: int,
    rng: random.Random,
    profile: str,
    q_mean: float,
    q_std: float,
    q_start: float,
    q_end: float,
    q_fixed: float,
) -> List[int]:
    n = max(0, int(length))
    profile = str(profile or "Fixed Q")
    if n == 0:
        return []
    if profile == "Normal Q distribution":
        return [_clamp_q(rng.gauss(float(q_mean), max(0.0, float(q_std)))) for _ in range(n)]
    if profile == "Linear quality decay":
        if n == 1:
            return [_clamp_q(q_start)]
        return [
            _clamp_q(float(q_start) + (float(q_end) - float(q_start)) * (i / (n - 1)))
            for i in range(n)
        ]
    return [_clamp_q(q_fixed)] * n


def _mutate_payload_with_qscores(
    payload: str,
    rng: random.Random,
    substitution_model: str,
    configured_substitution_rate: float,
    qscores: Sequence[int],
    insertion_rate: float,
    deletion_rate: float,
) -> Tuple[str, List[int], Dict[str, int]]:
    """Mutate one payload and keep a quality score aligned to every emitted base."""
    seq = clean_dna(payload)
    out_bases: List[str] = []
    out_q: List[int] = []
    sub = ins = dele = 0

    for i, base in enumerate(seq):
        q = int(qscores[i]) if i < len(qscores) else Q_MAX

        # Deletion removes both the base call and its quality character.
        if rng.random() < float(deletion_rate):
            dele += 1
        else:
            p_sub = (
                _q_to_error_probability(q)
                if substitution_model == "Phred Q-score driven"
                else float(configured_substitution_rate)
            )
            called = base
            if rng.random() < p_sub:
                called = rng.choice([b for b in BASES if b != base])
                sub += 1
            out_bases.append(called)
            out_q.append(q)

        # Insertions are independent of the substitution/Q model in this build.
        if rng.random() < float(insertion_rate):
            out_bases.append(rng.choice(BASES))
            out_q.append(q)
            ins += 1

    return "".join(out_bases), out_q, {
        "Substitute count": sub,
        "Insertion count": ins,
        "Deletion count": dele,
    }


def simulate_reads(
    strand_rows: Sequence[Dict[str, Any]],
    coverage: int = 10,
    substitution_rate: float = 0.0,
    insertion_rate: float = 0.0,
    deletion_rate: float = 0.0,
    strand_dropout_rate: float = 0.0,
    seed: int = 11,
    error_probability_model: str = "Independent Bernoulli per nucleotide",
    qscore_profile: str = "Fixed Q",
    q_mean: float = 30.0,
    q_std: float = 4.0,
    q_start: float = 35.0,
    q_end: float = 20.0,
    q_fixed: float = 30.0,
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Create independent noisy copies of every non-dropped prepared strand.

    Error models
    ------------
    Independent Bernoulli per nucleotide
        Every payload nucleotide has the same configured substitution probability.
        FASTQ Q values are the Phred score implied by that probability.

    Phred Q-score driven
        A per-base Q profile is generated for each read. The substitution
        probability at a position is p=10^(-Q/10). Supported profiles are
        Normal Q distribution, Linear quality decay, and Fixed Q.

    Insertions/deletions remain independent Bernoulli events and are intentionally
    separate from the Q-driven substitution model. FBR/SI/filler/RBR are protected
    in this simulator so read grouping remains explicit; strand dropout removes a
    whole strand before coverage copies are generated.
    """
    coverage = max(1, int(coverage))
    rng = random.Random(int(seed))
    model = str(error_probability_model or "Independent Bernoulli per nucleotide")
    if model not in {"Independent Bernoulli per nucleotide", "Phred Q-score driven"}:
        raise ValueError(f"Unsupported error probability model: {model}")

    reads: List[Dict[str, Any]] = []
    dropped: List[int] = []
    sub_total = ins_total = del_total = 0
    payload_nt_total = 0
    all_payload_q: List[int] = []
    design_id = _strand_design_fingerprint(strand_rows)

    for fallback, row in enumerate(strand_rows, start=1):
        sid = int(row.get("No.", fallback))
        si = clean_dna(row.get("Strand index", ""))
        payload = clean_dna(row.get("Payload", ""))
        fbr = clean_dna(row.get("FBR", ""))
        filler = clean_dna(row.get("Filler", ""))
        rbr = clean_dna(row.get("RBR", ""))

        if rng.random() < float(strand_dropout_rate):
            dropped.append(sid)
            continue

        for copy_no in range(1, coverage + 1):
            # A deterministic per-read RNG keeps results stable for the same seed.
            read_rng = random.Random(f"{seed}|strand={sid}|copy={copy_no}")
            if model == "Phred Q-score driven":
                payload_q = _make_qscores(
                    len(payload), read_rng, qscore_profile,
                    q_mean, q_std, q_start, q_end, q_fixed,
                )
            else:
                payload_q = [_p_to_q(float(substitution_rate))] * len(payload)

            noisy_payload, noisy_payload_q, counts = _mutate_payload_with_qscores(
                payload,
                read_rng,
                substitution_model=model,
                configured_substitution_rate=float(substitution_rate),
                qscores=payload_q,
                insertion_rate=float(insertion_rate),
                deletion_rate=float(deletion_rate),
            )
            sub = int(counts.get("Substitute count", 0))
            ins = int(counts.get("Insertion count", 0))
            dele = int(counts.get("Deletion count", 0))
            sub_total += sub
            ins_total += ins
            del_total += dele
            payload_nt_total += len(payload)
            all_payload_q.extend(payload_q)

            full_read = fbr + si + noisy_payload + filler + rbr
            full_q = (
                [STRUCTURAL_Q] * len(fbr)
                + [STRUCTURAL_Q] * len(si)
                + noisy_payload_q
                + [STRUCTURAL_Q] * len(filler)
                + [STRUCTURAL_Q] * len(rbr)
            )
            if len(full_q) != len(full_read):
                raise RuntimeError("FASTQ quality length does not match simulated read length")

            reads.append({
                "Read ID": f"s{sid}_c{copy_no}",
                "Design ID": design_id,
                "Strand ID": sid,
                "SI": si,
                "Copy": copy_no,
                "Reference payload": payload,
                "Read payload": noisy_payload,
                "Read sequence": full_read,
                "Quality string": _quality_string(full_q),
                "Payload quality string": _quality_string(noisy_payload_q),
                "Mean payload Q": (sum(payload_q) / len(payload_q)) if payload_q else 0.0,
                "Substitutions": sub,
                "Insertions": ins,
                "Deletions": dele,
            })

    total_events = sub_total + ins_total + del_total
    return reads, {
        "coverage": coverage,
        "design_id": design_id,
        "read_count": len(reads),
        "strand_count": len(strand_rows),
        "dropped_strand_count": len(dropped),
        "dropped_strands": dropped,
        "substitution_events": sub_total,
        "insertion_events": ins_total,
        "deletion_events": del_total,
        "total_error_events": total_events,
        "reference_payload_nt_observed": payload_nt_total,
        "observed_substitution_rate": sub_total / max(1, payload_nt_total),
        "observed_insertion_rate": ins_total / max(1, payload_nt_total),
        "observed_deletion_rate": del_total / max(1, payload_nt_total),
        "configured_substitution_rate": float(substitution_rate),
        "configured_insertion_rate": float(insertion_rate),
        "configured_deletion_rate": float(deletion_rate),
        "configured_dropout_rate": float(strand_dropout_rate),
        "error_probability_model": model,
        "qscore_profile": qscore_profile if model == "Phred Q-score driven" else "Fixed from Bernoulli p",
        "mean_payload_q": (sum(all_payload_q) / len(all_payload_q)) if all_payload_q else 0.0,
        "min_payload_q": min(all_payload_q) if all_payload_q else 0,
        "max_payload_q": max(all_payload_q) if all_payload_q else 0,
        "q_mean": float(q_mean),
        "q_std": float(q_std),
        "q_start": float(q_start),
        "q_end": float(q_end),
        "q_fixed": float(q_fixed),
        "error_scope": "Payload only; SI/primers/filler protected in this simulator",
        "seed": int(seed),
    }


def consensus_by_si(
    strand_rows: Sequence[Dict[str, Any]],
    read_rows: Sequence[Dict[str, Any]],
    original_dna_length: int,
    method: str = "Majority voting",
) -> Tuple[List[Dict[str, Any]], str, Dict[str, Any], List[Dict[str, Any]]]:
    """Group reads by SI and reconstruct one payload consensus per strand.

    Supported methods
    -----------------
    Majority voting
        At each reference position, choose the nucleotide observed most often.

    Q-score weighted consensus
        Use each observed base's Phred Q score in a maximum-likelihood base call,
        so high-quality reads contribute more evidence than low-quality reads.
    """
    method = str(method or "Majority voting")
    if method not in {"Majority voting", "Q-score weighted consensus"}:
        raise ValueError(f"Unsupported consensus method: {method}")

    by_si: Dict[str, List[Dict[str, Any]]] = defaultdict(list)
    for read in read_rows:
        by_si[clean_dna(read.get("SI", ""))].append(dict(read))

    consensus_rows: List[Dict[str, Any]] = []
    evidence_rows: List[Dict[str, Any]] = []
    payloads: List[str] = []
    missing = incomplete = mismatches = low_conf = ties = 0
    confidence_values: List[float] = []

    for fallback, row in enumerate(strand_rows, start=1):
        sid = int(row.get("No.", fallback))
        si = clean_dna(row.get("Strand index", ""))
        reference = clean_dna(row.get("Payload", ""))
        reads = by_si.get(si, [])
        if not reads:
            missing += 1
            consensus_rows.append({
                "Strand ID": sid,
                "SI": si,
                "Coverage": 0,
                "Reference length": len(reference),
                "Consensus length": 0,
                "Consensus payload": "",
                "Mismatches": "—",
                "Mean confidence": 0.0,
                "Status": "Missing",
            })
            continue

        aligned_pairs: List[Tuple[List[str | None], List[int | None]]] = []
        for r in reads:
            read_payload = clean_dna(r.get("Read payload", ""))
            payload_qual = str(r.get("Payload quality string", ""))
            aligned_pairs.append(_align_read_and_q_to_reference(reference, read_payload, payload_qual))

        consensus_chars: List[str] = []
        strand_conf: List[float] = []
        strand_mismatch = 0
        complete = True

        for pos in range(len(reference)):
            observations: List[Tuple[str, int]] = []
            for bases, qs in aligned_pairs:
                base = bases[pos] if pos < len(bases) else None
                q = qs[pos] if pos < len(qs) else None
                if isinstance(base, str) and base in BASES and q is not None:
                    observations.append((str(base), int(q)))

            observed_bases = [b for b, _ in observations]
            counts = Counter(observed_bases)
            q_sums = {b: sum(q for called, q in observations if called == b) for b in BASES}

            if not observations:
                complete = False
                chosen = "N"
                conf = 0.0
                winners: List[str] = []
            elif method == "Majority voting":
                best = max(counts.values())
                winners = [b for b in BASES if counts.get(b, 0) == best]
                chosen = winners[0]
                conf = best / len(observations)
                ties += int(len(winners) > 1)
            else:
                chosen, conf, log_scores = _q_weighted_base_call(observations)
                best_score = max(log_scores.values())
                winners = [b for b in BASES if abs(log_scores[b] - best_score) < 1e-12]
                ties += int(len(winners) > 1)

            if conf < 0.70:
                low_conf += 1
            consensus_chars.append(chosen)
            strand_conf.append(conf)
            confidence_values.append(conf)
            if chosen != reference[pos]:
                strand_mismatch += 1
                mismatches += 1

            evidence_rows.append({
                "Strand ID": sid,
                "SI": si,
                "Position": pos + 1,
                "Coverage": len(observations),
                "A": counts.get("A", 0),
                "C": counts.get("C", 0),
                "G": counts.get("G", 0),
                "T": counts.get("T", 0),
                "A Q-sum": q_sums["A"],
                "C Q-sum": q_sums["C"],
                "G Q-sum": q_sums["G"],
                "T Q-sum": q_sums["T"],
                "Consensus": chosen,
                "Reference": reference[pos],
                "Confidence": conf,
                "Correct": chosen == reference[pos],
            })

        consensus = "".join(consensus_chars)
        if "N" in consensus:
            complete = False
        if not complete:
            incomplete += 1
        else:
            payloads.append(consensus)

        consensus_rows.append({
            "Strand ID": sid,
            "SI": si,
            "Coverage": len(reads),
            "Reference length": len(reference),
            "Consensus length": len(consensus.replace("N", "")) if not complete else len(consensus),
            "Consensus payload": consensus,
            "Mismatches": strand_mismatch,
            "Mean confidence": sum(strand_conf) / max(1, len(strand_conf)),
            "Status": "Recovered" if complete else "Incomplete",
        })

    all_complete = (missing == 0 and incomplete == 0 and len(consensus_rows) == len(strand_rows))
    consensus_dna = ""
    if all_complete:
        consensus_dna = clean_dna("".join(r["Consensus payload"] for r in consensus_rows))[:int(original_dna_length)]

    reference_dna = clean_dna("".join(clean_dna(r.get("Payload", "")) for r in strand_rows))[:int(original_dna_length)]
    final_mismatches = None
    accuracy = None
    if consensus_dna and len(consensus_dna) == len(reference_dna):
        final_mismatches = sum(a != b for a, b in zip(consensus_dna, reference_dna))
        accuracy = 1.0 - final_mismatches / max(1, len(reference_dna))

    return consensus_rows, consensus_dna, {
        "method": method,
        "all_strands_complete": all_complete,
        "missing_strands": missing,
        "incomplete_strands": incomplete,
        "final_dna_mismatches": final_mismatches,
        "consensus_accuracy": accuracy,
        "mean_confidence": sum(confidence_values) / max(1, len(confidence_values)),
        "low_confidence_positions": low_conf,
        "ties": ties,
        "consensus_dna_length": len(consensus_dna),
    }, evidence_rows



def parse_fastq_reads(
    fastq_text: str,
    strand_rows: Sequence[Dict[str, Any]],
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """Parse uploaded single-end FASTQ reads into the internal read-row format.

    The parser uses the current Strand Design to recover the SI and payload region
    from each read. It is therefore intended for FASTQ data generated from the
    prepared strands of the current project/session. Phred+33 quality strings are
    preserved so Q-score weighted consensus can use the uploaded base qualities.
    """
    text = str(fastq_text or "")
    raw_lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    # Ignore empty separator lines, but preserve all non-empty FASTQ fields.
    lines = [line.strip() for line in raw_lines if line.strip() != ""]
    if not lines:
        raise ValueError("FASTQ file is empty.")
    if len(lines) % 4 != 0:
        raise ValueError(
            f"Invalid FASTQ structure: expected groups of 4 lines, found {len(lines)} non-empty lines."
        )
    if not strand_rows:
        raise ValueError("Current Strand Design is required before FASTQ decoding.")

    rows_by_sid: Dict[int, Dict[str, Any]] = {}
    rows_by_si: Dict[str, Dict[str, Any]] = {}
    for fallback, row in enumerate(strand_rows, start=1):
        sid = int(row.get("No.", fallback))
        si = clean_dna(row.get("Strand index", ""))
        rows_by_sid[sid] = dict(row)
        if si:
            rows_by_si[si] = dict(row)

    first = dict(strand_rows[0])
    common_fbr = clean_dna(first.get("FBR", ""))
    index_len = len(clean_dna(first.get("Strand index", "")))

    parsed: List[Dict[str, Any]] = []
    unmatched = 0
    structural_mismatch = 0
    current_design_id = _strand_design_fingerprint(strand_rows)
    header_design_ids = set()
    q_values_all: List[int] = []
    recognized_sids = set()

    import re

    for rec_idx in range(0, len(lines), 4):
        header, seq_raw, plus, qual = lines[rec_idx:rec_idx + 4]
        record_no = rec_idx // 4 + 1
        if not header.startswith("@"):
            raise ValueError(f"Invalid FASTQ record {record_no}: header must start with '@'.")
        if not plus.startswith("+"):
            raise ValueError(f"Invalid FASTQ record {record_no}: third line must start with '+'.")

        design_match = re.search(r"(?:^|[|;])design=([0-9a-fA-F]{16})(?:$|[|;])", header[1:])
        if design_match:
            header_design_id = design_match.group(1).lower()
            header_design_ids.add(header_design_id)
            if header_design_id != current_design_id:
                raise ValueError(
                    "FASTQ design fingerprint does not match the current Strand Design. "
                    "Upload the FASTQ generated from the current encoded payload, or recreate the matching encoding/strand design first."
                )

        seq = str(seq_raw).upper().replace(" ", "")
        if any(base not in BASES for base in seq):
            raise ValueError(
                f"FASTQ record {record_no} contains unsupported nucleotide(s). "
                "This build accepts A/C/G/T reads."
            )
        if len(seq) != len(qual):
            raise ValueError(
                f"FASTQ record {record_no} has sequence length {len(seq)} but quality length {len(qual)}."
            )

        # First try the app-generated header (@s<strand>_c<copy>), then fall back
        # to SI extraction from the protected structural prefix.
        row = None
        sid = None
        copy_no = 1
        m = re.search(r"(?:^|[^A-Za-z0-9])s(\d+)_c(\d+)(?:$|[^0-9])", header[1:])
        if m:
            sid = int(m.group(1))
            copy_no = int(m.group(2))
            row = rows_by_sid.get(sid)

        if row is None and index_len > 0 and len(seq) >= len(common_fbr) + index_len:
            si_guess = clean_dna(seq[len(common_fbr):len(common_fbr) + index_len])
            row = rows_by_si.get(si_guess)
            if row is not None:
                sid = int(row.get("No.", 0) or 0)

        if row is None:
            unmatched += 1
            continue

        si = clean_dna(row.get("Strand index", ""))
        fbr = clean_dna(row.get("FBR", ""))
        filler = clean_dna(row.get("Filler", ""))
        rbr = clean_dna(row.get("RBR", ""))

        # Even when an app-generated header names a strand (e.g. s12_c3),
        # verify that the SI carried by the sequence agrees with the CURRENT
        # Strand Design. This prevents a FASTQ from an older/different encoding
        # from being silently assigned by strand number and decoded as garbage.
        if si and len(seq) >= len(fbr) + len(si):
            observed_si = clean_dna(seq[len(fbr):len(fbr) + len(si)])
            if observed_si != si:
                structural_mismatch += 1
                unmatched += 1
                continue

        prefix_len = len(fbr) + len(si)
        suffix_len = len(filler) + len(rbr)
        if len(seq) < prefix_len + suffix_len:
            unmatched += 1
            continue

        expected_prefix = fbr + si
        expected_suffix = filler + rbr
        if seq[:prefix_len] != expected_prefix or (suffix_len and seq[-suffix_len:] != expected_suffix):
            structural_mismatch += 1

        end = len(seq) - suffix_len if suffix_len else len(seq)
        payload = seq[prefix_len:end]
        payload_qual = qual[prefix_len:end]
        qs = [max(0, ord(ch) - 33) for ch in payload_qual]
        q_values_all.extend(qs)
        recognized_sids.add(int(sid or row.get("No.", 0) or 0))

        parsed.append({
            "Read ID": header[1:] or f"uploaded_{record_no}",
            "Design ID": next(iter(header_design_ids)) if header_design_ids else "",
            "Strand ID": int(sid or row.get("No.", 0) or 0),
            "SI": si,
            "Copy": int(copy_no),
            "Reference payload": clean_dna(row.get("Payload", "")),
            "Read payload": payload,
            "Read sequence": seq,
            "Quality string": qual,
            "Payload quality string": payload_qual,
            "Mean payload Q": (sum(qs) / len(qs)) if qs else 0.0,
            "Substitutions": None,
            "Insertions": None,
            "Deletions": None,
            "Source": "Uploaded FASTQ",
        })

    total_records = len(lines) // 4
    return parsed, {
        "records": total_records,
        "matched_reads": len(parsed),
        "unmatched_reads": unmatched,
        "recognized_strands": len(recognized_sids),
        "structural_mismatch_reads": structural_mismatch,
        "mean_q": (sum(q_values_all) / len(q_values_all)) if q_values_all else 0.0,
        "quality_encoding": "Phred+33",
        "current_design_id": current_design_id,
        "fastq_design_id": next(iter(header_design_ids)) if len(header_design_ids) == 1 else ("multiple" if header_design_ids else "not embedded"),
        "design_fingerprint_match": (len(header_design_ids) == 1 and next(iter(header_design_ids)) == current_design_id) if header_design_ids else None,
    }

def reads_to_fastq(read_rows: Sequence[Dict[str, Any]]) -> str:
    lines: List[str] = []
    for r in read_rows:
        seq = clean_dna(r.get("Read sequence", ""))
        if not seq:
            continue
        qual = str(r.get("Quality string", ""))
        if len(qual) != len(seq):
            # Backward compatibility for read sets created by older sessions.
            qual = "I" * len(seq)
        read_id = str(r.get("Read ID", "read"))
        design_id = str(r.get("Design ID", "")).strip()
        header = f"@{read_id}" + (f"|design={design_id}" if design_id and "design=" not in read_id else "")
        lines.extend([header, seq, "+", qual])
    return "\n".join(lines) + ("\n" if lines else "")


def generate_fastq(read_rows: Sequence[Dict[str, Any]]) -> Tuple[str, Dict[str, Any]]:
    """Build an exportable FASTQ read set using the simulated per-base Q scores."""
    text = reads_to_fastq(read_rows)
    lengths = [len(clean_dna(r.get("Read sequence", ""))) for r in read_rows if clean_dna(r.get("Read sequence", ""))]
    mean_q_values = [float(r.get("Mean payload Q", 0.0)) for r in read_rows if r.get("Mean payload Q") is not None]
    return text, {
        "records": len(lengths),
        "bytes": len(text.encode("utf-8")),
        "quality_model": "Simulated Phred+33 per-base Q scores",
        "min_read_length": min(lengths) if lengths else 0,
        "max_read_length": max(lengths) if lengths else 0,
        "mean_read_length": (sum(lengths) / len(lengths)) if lengths else 0.0,
        "mean_payload_q": (sum(mean_q_values) / len(mean_q_values)) if mean_q_values else 0.0,
    }
