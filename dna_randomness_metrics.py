"""
dna_randomness_metrics.py
=========================

Metrics for evaluating structural randomness of DNA-storage sequences.

Implemented from the uploaded manuscript where equations/mappings are explicit:
1. Sequence quality: GC, longest homopolymer, information density.
2. Translational active-particle model:
   A=(1,0), C=(0,1), G=(-1,0), T=(0,-1)
   MSD(n) = 4 D n + V^2 n^2
3. 3-input/1-output model:
   alpha = N0/(N0+N1), beta = (N0+N1)/total.
4. Polarization gamma = mean(spin), with 0 -> -1 and 1 -> +1.

Reproducibility note
--------------------
The uploaded main manuscript refers to Supporting Information S4 for the exact
closed-form definitions of inverse-Ising lambda and h, but S4 is absent.
Therefore lambda and h are estimated by standard pseudo-likelihood:
    logit P(s_i=+1 | neighbours) = 2h + 2J*sum(neighbour spins)
The returned ising_lambda is J. It is scientifically meaningful, but must not
be claimed to be numerically identical to the missing S4 formula.

The manuscript also omits the full cube vertex-transition table for the
rotational model. rotational_cycle_metrics() is a transparent paper-inspired
reconstruction based on the six cyclic permutations described in the text.

Dependency: NumPy only.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from math import pi, sqrt
from typing import Any, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

import numpy as np


DNA_ALPHABET = frozenset("ACGT")
DNA_TO_BITS: Mapping[str, str] = {"A": "00", "C": "01", "G": "10", "T": "11"}
TRANSLATIONAL_VECTORS: Mapping[str, Tuple[float, float]] = {
    "A": (1.0, 0.0),
    "C": (0.0, 1.0),
    "G": (-1.0, 0.0),
    "T": (0.0, -1.0),
}
LOGIC_CONTEXTS: Tuple[str, ...] = tuple(f"{i:03b}" for i in range(8))


class DNASequenceError(ValueError):
    """Raised when a DNA sequence cannot be analysed safely."""


def clean_dna(sequence: str, *, convert_u_to_t: bool = False) -> str:
    """Normalize a DNA sequence and reject non-ACGT symbols."""
    if not isinstance(sequence, str):
        raise TypeError("sequence must be a string")
    seq = "".join(sequence.upper().split())
    if convert_u_to_t:
        seq = seq.replace("U", "T")
    if not seq:
        raise DNASequenceError("DNA sequence is empty")
    invalid = sorted(set(seq) - DNA_ALPHABET)
    if invalid:
        raise DNASequenceError(
            "DNA sequence contains unsupported symbols: "
            + ", ".join(invalid)
            + ". Only A, C, G and T are accepted."
        )
    return seq


def dna_to_bits(sequence: str) -> str:
    """Convert DNA to bits using A=00, C=01, G=10, T=11."""
    seq = clean_dna(sequence)
    return "".join(DNA_TO_BITS[base] for base in seq)


def _longest_homopolymer(seq: str) -> int:
    longest = current = 1
    for previous, current_base in zip(seq, seq[1:]):
        if current_base == previous:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
    return longest


def _homopolymer_segment_count(seq: str, minimum_length: int = 2) -> int:
    if minimum_length < 1:
        raise ValueError("minimum_length must be >= 1")
    count = 0
    run = 1
    for previous, current in zip(seq, seq[1:]):
        if current == previous:
            run += 1
        else:
            if run >= minimum_length:
                count += 1
            run = 1
    if run >= minimum_length:
        count += 1
    return count


def _normalized_entropy_from_counts(counts: Sequence[int]) -> float:
    values = np.asarray(counts, dtype=float)
    total = float(values.sum())
    if total <= 0:
        return float("nan")
    probabilities = values[values > 0] / total
    entropy = -float(np.sum(probabilities * np.log2(probabilities)))
    max_entropy = np.log2(len(values))
    return entropy / max_entropy if max_entropy > 0 else 0.0


def sequence_quality_metrics(
    sequence: str,
    *,
    original_bit_count: Optional[int] = None,
) -> Dict[str, float]:
    """Compute basic DNA-storage sequence metrics."""
    seq = clean_dna(sequence)
    n = len(seq)
    counts = {base: seq.count(base) for base in "ACGT"}
    gc_fraction = (counts["G"] + counts["C"]) / n
    result: Dict[str, float] = {
        "length_nt": float(n),
        "a_fraction": counts["A"] / n,
        "c_fraction": counts["C"] / n,
        "g_fraction": counts["G"] / n,
        "t_fraction": counts["T"] / n,
        "gc_fraction": gc_fraction,
        "gc_deviation_from_0_5": abs(gc_fraction - 0.5),
        "longest_homopolymer": float(_longest_homopolymer(seq)),
        "homopolymer_segments_ge_2": float(_homopolymer_segment_count(seq, 2)),
        "base_entropy_normalized": _normalized_entropy_from_counts(
            [counts["A"], counts["C"], counts["G"], counts["T"]]
        ),
    }
    if original_bit_count is not None:
        if original_bit_count < 0:
            raise ValueError("original_bit_count must be >= 0")
        result["information_density_bits_per_nt"] = original_bit_count / n
    return result


def _fit_nonnegative_two_term_model(
    n_values: np.ndarray,
    y_values: np.ndarray,
    *,
    linear_multiplier: float,
) -> Tuple[float, float, float, np.ndarray]:
    """Fit y = linear_multiplier*diffusion*n + speed^2*n^2."""
    n_values = np.asarray(n_values, dtype=float)
    y_values = np.asarray(y_values, dtype=float)
    if n_values.ndim != 1 or y_values.ndim != 1:
        raise ValueError("n_values and y_values must be one-dimensional")
    if len(n_values) != len(y_values) or len(n_values) < 2:
        raise ValueError("At least two matched observations are required")

    design = np.column_stack([linear_multiplier * n_values, n_values**2])
    candidates: List[np.ndarray] = []
    coef, *_ = np.linalg.lstsq(design, y_values, rcond=None)
    if np.all(coef >= 0):
        candidates.append(coef)

    quadratic = design[:, 1]
    b = max(0.0, float(np.dot(quadratic, y_values) / np.dot(quadratic, quadratic)))
    candidates.append(np.array([0.0, b]))

    linear = design[:, 0]
    a = max(0.0, float(np.dot(linear, y_values) / np.dot(linear, linear)))
    candidates.append(np.array([a, 0.0]))
    candidates.append(np.array([0.0, 0.0]))

    best_coef = min(candidates, key=lambda c: float(np.sum((y_values - design @ c) ** 2)))
    fitted = design @ best_coef
    residual_ss = float(np.sum((y_values - fitted) ** 2))
    total_ss = float(np.sum((y_values - np.mean(y_values)) ** 2))
    r_squared = 1.0 - residual_ss / total_ss if total_ss > 0 else 1.0
    diffusion = float(best_coef[0])
    speed = sqrt(max(0.0, float(best_coef[1])))
    return diffusion, speed, r_squared, fitted


def translational_active_particle_metrics(
    sequence: str,
    *,
    max_lag: int = 20,
    step_length: float = 1.0,
) -> Dict[str, Any]:
    """Paper-faithful 2D active-particle MSD analysis."""
    seq = clean_dna(sequence)
    if max_lag < 2:
        raise ValueError("max_lag must be >= 2")
    if step_length <= 0:
        raise ValueError("step_length must be > 0")
    if len(seq) <= max_lag:
        raise DNASequenceError(
            f"Sequence length ({len(seq)}) must be greater than max_lag ({max_lag})."
        )

    steps = np.asarray([TRANSLATIONAL_VECTORS[b] for b in seq], dtype=float) * step_length
    positions = np.vstack([np.zeros((1, 2)), np.cumsum(steps, axis=0)])
    lags = np.arange(1, max_lag + 1, dtype=int)
    starts = np.arange(0, len(seq) - max_lag + 1)
    msd = np.empty(max_lag, dtype=float)
    for index, lag in enumerate(lags):
        displacement = positions[starts + lag] - positions[starts]
        msd[index] = float(np.mean(np.sum(displacement**2, axis=1)))

    diffusion, velocity, r_squared, fitted = _fit_nonnegative_two_term_model(
        lags.astype(float), msd, linear_multiplier=4.0
    )
    return {
        "velocity_V_L_per_nt": velocity,
        "diffusion_D_L2_per_nt": diffusion,
        "fit_r_squared": r_squared,
        "max_lag_nt": int(max_lag),
        "window_count": int(len(starts)),
        "lags_nt": lags.tolist(),
        "msd": msd.tolist(),
        "msd_fitted": fitted.tolist(),
    }


def bits_to_matrix(bits: str, *, shape: Optional[Tuple[int, int]] = None) -> np.ndarray:
    """Reshape bits to a strict supplied shape or the largest near-square matrix."""
    if not bits or any(bit not in "01" for bit in bits):
        raise ValueError("bits must be a non-empty string containing only 0 and 1")
    values = np.fromiter((int(bit) for bit in bits), dtype=np.int8)
    if shape is not None:
        rows, cols = shape
        if rows < 2 or cols < 3:
            raise ValueError("shape must have at least 2 rows and 3 columns")
        required = rows * cols
        if len(values) < required:
            raise ValueError(f"Need {required} bits for shape {shape}, but only {len(values)} exist")
        return values[:required].reshape(rows, cols)

    rows = int(np.floor(np.sqrt(len(values))))
    while rows >= 2:
        cols = len(values) // rows
        if cols >= 3:
            return values[: rows * cols].reshape(rows, cols)
        rows -= 1
    raise ValueError("Not enough bits to construct a matrix of at least 2 x 3")


def dna_to_binary_matrix(
    sequence: str, *, shape: Optional[Tuple[int, int]] = None
) -> np.ndarray:
    return bits_to_matrix(dna_to_bits(sequence), shape=shape)


def _fit_inverse_ising_pseudolikelihood(
    spins: np.ndarray,
    *,
    max_iter: int = 100,
    tolerance: float = 1e-9,
    ridge: float = 1e-8,
) -> Tuple[float, float, bool, int]:
    """Estimate uniform Ising coupling J and field h by pseudo-likelihood."""
    neighbour_sum = np.zeros_like(spins, dtype=float)
    degree = np.zeros_like(spins, dtype=float)
    neighbour_sum[1:, :] += spins[:-1, :]
    degree[1:, :] += 1
    neighbour_sum[:-1, :] += spins[1:, :]
    degree[:-1, :] += 1
    neighbour_sum[:, 1:] += spins[:, :-1]
    degree[:, 1:] += 1
    neighbour_sum[:, :-1] += spins[:, 1:]
    degree[:, :-1] += 1

    valid = degree > 0
    x = neighbour_sum[valid].astype(float)
    y = ((spins[valid] + 1) / 2).astype(float)
    design = np.column_stack([np.ones_like(x), x])
    beta = np.zeros(2, dtype=float)  # [2h, 2J]
    converged = False

    for iteration in range(1, max_iter + 1):
        linear = np.clip(design @ beta, -40.0, 40.0)
        probability = 1.0 / (1.0 + np.exp(-linear))
        weights = np.maximum(probability * (1.0 - probability), 1e-9)
        gradient = design.T @ (y - probability) - ridge * beta
        hessian_positive = design.T @ (weights[:, None] * design) + ridge * np.eye(2)
        try:
            step = np.linalg.solve(hessian_positive, gradient)
        except np.linalg.LinAlgError:
            step = np.linalg.pinv(hessian_positive) @ gradient
        beta_new = beta + step
        if np.linalg.norm(beta_new - beta) < tolerance:
            beta = beta_new
            converged = True
            break
        beta = beta_new

    return float(beta[1] / 2), float(beta[0] / 2), converged, iteration


def inverse_ising_metrics(
    sequence: str, *, matrix_shape: Optional[Tuple[int, int]] = None
) -> Dict[str, Any]:
    """Compute gamma and standard pseudo-likelihood lambda/h estimates."""
    matrix = dna_to_binary_matrix(sequence, shape=matrix_shape)
    spins = 2.0 * matrix.astype(float) - 1.0
    gamma = float(np.mean(spins))
    horizontal = spins[:, :-1] * spins[:, 1:]
    vertical = spins[:-1, :] * spins[1:, :]
    nearest_neighbour_correlation = float(
        (horizontal.sum() + vertical.sum()) / (horizontal.size + vertical.size)
    )
    interaction, field, converged, iterations = _fit_inverse_ising_pseudolikelihood(spins)
    return {
        "polarization_gamma": gamma,
        "abs_polarization_gamma": abs(gamma),
        "ising_lambda": interaction,
        "abs_ising_lambda": abs(interaction),
        "ising_h": field,
        "abs_ising_h": abs(field),
        "nearest_neighbour_correlation": nearest_neighbour_correlation,
        "matrix_rows": int(matrix.shape[0]),
        "matrix_cols": int(matrix.shape[1]),
        "estimation_method": "standard_pseudolikelihood_not_missing_paper_S4_formula",
        "optimizer_converged": bool(converged),
        "optimizer_iterations": int(iterations),
    }


def logic_3input_1output_metrics(
    sequence: str,
    *,
    matrix_shape: Optional[Tuple[int, int]] = None,
    layout: str = "paper_grid",
) -> Dict[str, Any]:
    """
    Compute alpha and beta for all eight 3-bit inputs.

    paper_grid inputs for output (row,col), row>=1,col>=2:
        [left-two, left-one, above] -> current output.
    This yields (rows-1)*(cols-2), matching the manuscript's sample count.
    """
    bits = dna_to_bits(sequence)
    counts_zero = {context: 0 for context in LOGIC_CONTEXTS}
    counts_one = {context: 0 for context in LOGIC_CONTEXTS}

    if layout == "paper_grid":
        matrix = bits_to_matrix(bits, shape=matrix_shape)
        rows, cols = matrix.shape
        for row in range(1, rows):
            for col in range(2, cols):
                context = f"{matrix[row, col-2]}{matrix[row, col-1]}{matrix[row-1, col]}"
                output = int(matrix[row, col])
                (counts_zero if output == 0 else counts_one)[context] += 1
        expected_sample_count = (rows - 1) * (cols - 2)
        matrix_rows, matrix_cols = rows, cols
    elif layout == "linear":
        for index in range(len(bits) - 3):
            context = bits[index:index+3]
            output = int(bits[index+3])
            (counts_zero if output == 0 else counts_one)[context] += 1
        expected_sample_count = max(0, len(bits) - 3)
        matrix_rows = matrix_cols = None
    else:
        raise ValueError("layout must be 'paper_grid' or 'linear'")

    totals = {c: counts_zero[c] + counts_one[c] for c in LOGIC_CONTEXTS}
    total_samples = sum(totals.values())
    if total_samples == 0:
        raise DNASequenceError("Not enough data for 3-input/1-output analysis")

    alpha: Dict[str, Optional[float]] = {}
    beta: Dict[str, float] = {}
    for context in LOGIC_CONTEXTS:
        total = totals[context]
        alpha[context] = counts_zero[context] / total if total > 0 else None
        beta[context] = total / total_samples

    observed_alpha = np.asarray([v for v in alpha.values() if v is not None], dtype=float)
    all_beta = np.asarray(list(beta.values()), dtype=float)
    alpha_rmse = float(np.sqrt(np.mean((observed_alpha - 0.5) ** 2)))
    beta_rmse = float(np.sqrt(np.mean((all_beta - 0.125) ** 2)))

    return {
        "alpha": alpha,
        "beta": beta,
        "alpha_rmse_from_0_5": alpha_rmse,
        "beta_rmse_from_0_125": beta_rmse,
        "alpha_rms_raw": float(np.sqrt(np.mean(observed_alpha**2))),
        "beta_rms_raw": float(np.sqrt(np.mean(all_beta**2))),
        "missing_alpha_contexts": [c for c, v in alpha.items() if v is None],
        "counts_zero": counts_zero,
        "counts_one": counts_one,
        "sample_count": int(total_samples),
        "expected_sample_count": int(expected_sample_count),
        "layout": layout,
        "matrix_rows": matrix_rows,
        "matrix_cols": matrix_cols,
    }


_ROTATIONAL_CYCLES: Mapping[str, Tuple[float, float, float]] = {
    "ACGT": (0.0, 0.0, 1.0),
    "ATGC": (0.0, 0.0, -1.0),
    "ATCG": (1.0, 0.0, 0.0),
    "AGCT": (-1.0, 0.0, 0.0),
    "ACTG": (0.0, 1.0, 0.0),
    "AGTC": (0.0, -1.0, 0.0),
}


def _all_rotations(text: str) -> Iterable[str]:
    for shift in range(len(text)):
        yield text[shift:] + text[:shift]


_ROTATIONAL_LOOKUP: Dict[str, Tuple[float, float, float]] = {}
for canonical_cycle, direction in _ROTATIONAL_CYCLES.items():
    for rotation in _all_rotations(canonical_cycle):
        _ROTATIONAL_LOOKUP[rotation] = direction


def rotational_cycle_metrics(
    sequence: str,
    *,
    max_lag: int = 20,
    angular_step: float = pi / 2.0,
) -> Dict[str, Any]:
    """Paper-inspired periodic-cycle MSAD approximation."""
    seq = clean_dna(sequence)
    if max_lag < 2:
        raise ValueError("max_lag must be >= 2")
    if len(seq) < max_lag + 3:
        raise DNASequenceError(
            f"Sequence length ({len(seq)}) must be at least {max_lag + 3}."
        )

    increments: List[Tuple[float, float, float]] = []
    classified = 0
    for start in range(len(seq) - 3):
        window = seq[start:start+4]
        direction = _ROTATIONAL_LOOKUP.get(window)
        if direction is None:
            increments.append((0.0, 0.0, 0.0))
        else:
            increments.append(tuple(angular_step * value for value in direction))
            classified += 1

    angular_steps = np.asarray(increments, dtype=float)
    positions = np.vstack([np.zeros((1, 3)), np.cumsum(angular_steps, axis=0)])
    lags = np.arange(1, max_lag + 1, dtype=int)
    starts = np.arange(0, len(angular_steps) - max_lag + 1)
    msad = np.empty(max_lag, dtype=float)
    for index, lag in enumerate(lags):
        displacement = positions[starts + lag] - positions[starts]
        msad[index] = float(np.mean(np.sum(displacement**2, axis=1)))

    rotational_diffusion, angular_velocity, r_squared, fitted = _fit_nonnegative_two_term_model(
        lags.astype(float), msad, linear_multiplier=2.0
    )
    return {
        "angular_velocity_omega_rad_per_step": angular_velocity,
        "rotational_diffusion_DR_rad2_per_step": rotational_diffusion,
        "fit_r_squared": r_squared,
        "max_lag_steps": int(max_lag),
        "window_count": int(len(starts)),
        "classified_fourmer_fraction": classified / len(increments),
        "lags": lags.tolist(),
        "msad": msad.tolist(),
        "msad_fitted": fitted.tolist(),
        "method": "paper_inspired_cycle_reconstruction_not_exact_missing_transition_table",
    }


@dataclass
class DNASequenceRandomnessReport:
    sequence_quality: Dict[str, Any]
    translational_active_particle: Dict[str, Any]
    inverse_ising: Dict[str, Any]
    logic_3input_1output: Dict[str, Any]
    rotational_active_particle: Optional[Dict[str, Any]]
    metadata: Dict[str, Any]

    def to_dict(self) -> Dict[str, Any]:
        return asdict(self)

    def app_summary(self) -> Dict[str, Any]:
        q = self.sequence_quality
        t = self.translational_active_particle
        i = self.inverse_ising
        l = self.logic_3input_1output
        summary: Dict[str, Any] = {
            "DNA length (nt)": int(q["length_nt"]),
            "GC fraction": q["gc_fraction"],
            "GC deviation from 0.5": q["gc_deviation_from_0_5"],
            "Longest homopolymer": int(q["longest_homopolymer"]),
            "Homopolymer segments (2+ bases)": int(q["homopolymer_segments_ge_2"]),
            "Translational velocity V": t["velocity_V_L_per_nt"],
            "Translational diffusion D": t["diffusion_D_L2_per_nt"],
            "|Polarization gamma|": i["abs_polarization_gamma"],
            "|Inverse-Ising lambda|": i["abs_ising_lambda"],
            "|Inverse-Ising h|": i["abs_ising_h"],
            "Alpha RMSE from 0.5": l["alpha_rmse_from_0_5"],
            "Beta RMSE from 0.125": l["beta_rmse_from_0_125"],
        }
        if "information_density_bits_per_nt" in q:
            summary["Information density (bits/nt)"] = q["information_density_bits_per_nt"]
        if self.rotational_active_particle is not None:
            r = self.rotational_active_particle
            summary["Angular velocity omega"] = r["angular_velocity_omega_rad_per_step"]
            summary["Rotational diffusion DR"] = r["rotational_diffusion_DR_rad2_per_step"]
        return summary


def analyze_dna_randomness(
    sequence: str,
    *,
    original_bit_count: Optional[int] = None,
    matrix_shape: Optional[Tuple[int, int]] = None,
    max_lag: int = 20,
    logic_layout: str = "paper_grid",
    include_rotational_approximation: bool = True,
) -> DNASequenceRandomnessReport:
    """Run the complete metric set for one DNA sequence."""
    seq = clean_dna(sequence)
    quality = sequence_quality_metrics(seq, original_bit_count=original_bit_count)
    translational = translational_active_particle_metrics(seq, max_lag=max_lag)
    ising = inverse_ising_metrics(seq, matrix_shape=matrix_shape)
    logic = logic_3input_1output_metrics(
        seq, matrix_shape=matrix_shape, layout=logic_layout
    )
    rotational = (
        rotational_cycle_metrics(seq, max_lag=max_lag)
        if include_rotational_approximation
        else None
    )
    return DNASequenceRandomnessReport(
        sequence_quality=quality,
        translational_active_particle=translational,
        inverse_ising=ising,
        logic_3input_1output=logic,
        rotational_active_particle=rotational,
        metadata={
            "dna_to_bits_mapping": dict(DNA_TO_BITS),
            "max_lag": max_lag,
            "matrix_shape_requested": matrix_shape,
            "logic_layout": logic_layout,
            "lower_is_better_metrics": [
                "gc_deviation_from_0_5",
                "longest_homopolymer",
                "velocity_V_L_per_nt",
                "diffusion_D_L2_per_nt",
                "abs_polarization_gamma",
                "abs_ising_lambda",
                "abs_ising_h",
                "alpha_rmse_from_0_5",
                "beta_rmse_from_0_125",
                "angular_velocity_omega_rad_per_step",
            ],
            "caveats": [
                "Inverse-Ising lambda/h use pseudo-likelihood because manuscript SI S4 is missing.",
                "Rotational omega/DR use a paper-inspired cycle reconstruction because the full cube transition table is missing.",
                "Do not compare runs with different matrix shapes, max_lag values, or logic layouts.",
            ],
        },
    )


def compare_sequences(
    sequences: Mapping[str, str],
    *,
    original_bit_counts: Optional[Mapping[str, int]] = None,
    matrix_shape: Optional[Tuple[int, int]] = None,
    max_lag: int = 20,
    logic_layout: str = "paper_grid",
    include_rotational_approximation: bool = True,
) -> List[Dict[str, Any]]:
    """Analyse named sequences and return rows suitable for pandas/Streamlit."""
    rows: List[Dict[str, Any]] = []
    for name, sequence in sequences.items():
        original_bits = original_bit_counts.get(name) if original_bit_counts else None
        report = analyze_dna_randomness(
            sequence,
            original_bit_count=original_bits,
            matrix_shape=matrix_shape,
            max_lag=max_lag,
            logic_layout=logic_layout,
            include_rotational_approximation=include_rotational_approximation,
        )
        row = {"Sequence": name}
        row.update(report.app_summary())
        rows.append(row)
    return rows
