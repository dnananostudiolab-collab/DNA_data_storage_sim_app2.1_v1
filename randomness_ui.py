from __future__ import annotations

import hashlib
from typing import Any, Dict, Optional

import pandas as pd
import plotly.graph_objects as go
import streamlit as st

from dna_randomness_metrics import DNASequenceError, analyze_dna_randomness


def _safe_float(value: Any, default: float = 0.0) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _metric_text(value: Any, digits: int = 4) -> str:
    try:
        return f"{float(value):.{digits}f}"
    except (TypeError, ValueError):
        return "—"


def _report_hash(dna: str, original_bit_count: Optional[int], max_lag: int, include_rotational: bool) -> str:
    payload = f"{dna}|{original_bit_count}|{max_lag}|{int(include_rotational)}".encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _reference_line(fig: go.Figure, value: float, label: str) -> None:
    fig.add_hline(y=value, line_dash="dash", annotation_text=label, annotation_position="top left")


def render_randomness_analysis(
    dna: str,
    *,
    state_prefix: str,
    original_bit_count: Optional[int] = None,
    default_expanded: bool = False,
) -> None:
    """Render and cache structural-randomness metrics for one encoded DNA payload."""
    if not dna:
        return

    report_key = f"{state_prefix}_randomness_report"
    hash_key = f"{state_prefix}_randomness_hash"
    error_key = f"{state_prefix}_randomness_error"

    with st.expander("🎲 DNA Structural Randomness", expanded=default_expanded):
        st.caption(
            "Evaluates the encoded DNA payload only. Primers, strand IDs, filler, sequencing noise and consensus DNA are excluded."
        )

        c1, c2 = st.columns(2)
        max_allowed = max(2, min(100, len(dna) - 1))
        default_lag = min(20, max_allowed)
        max_lag = c1.number_input(
            "Maximum lag",
            min_value=2,
            max_value=max_allowed,
            value=default_lag,
            step=1,
            key=f"{state_prefix}_randomness_max_lag",
        )
        include_rotational = c2.checkbox(
            "Include experimental rotational model",
            value=False,
            key=f"{state_prefix}_randomness_rotational",
            help="Uses a paper-inspired reconstruction because the complete transition table is unavailable.",
        )

        current_hash = _report_hash(dna, original_bit_count, int(max_lag), bool(include_rotational))
        stale = st.session_state.get(hash_key) != current_hash
        button_label = "Run Randomization Analysis" if report_key not in st.session_state else "Recalculate Randomization Analysis"

        if st.button(
            button_label,
            key=f"{state_prefix}_run_randomness_analysis",
            type="primary",
            use_container_width=True,
        ):
            try:
                report = analyze_dna_randomness(
                    dna,
                    original_bit_count=original_bit_count,
                    max_lag=int(max_lag),
                    logic_layout="paper_grid",
                    include_rotational_approximation=bool(include_rotational),
                ).to_dict()
                st.session_state[report_key] = report
                st.session_state[hash_key] = current_hash
                st.session_state.pop(error_key, None)
            except (DNASequenceError, ValueError, TypeError, RuntimeError) as exc:
                st.session_state.pop(report_key, None)
                st.session_state.pop(hash_key, None)
                st.session_state[error_key] = str(exc)

        if stale and report_key in st.session_state:
            st.info("DNA or analysis settings changed. Recalculate to update the results.")
        if st.session_state.get(error_key):
            st.warning(f"Randomization analysis is unavailable: {st.session_state[error_key]}")

        report: Optional[Dict[str, Any]] = st.session_state.get(report_key)
        if not report:
            st.caption("Run the analysis to generate summary metrics and plots.")
            return

        quality = report.get("sequence_quality", {})
        translational = report.get("translational_active_particle", {})
        ising = report.get("inverse_ising", {})
        logic = report.get("logic_3input_1output", {})
        rotational = report.get("rotational_active_particle")

        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric("Base entropy", _metric_text(quality.get("base_entropy_normalized")))
        m2.metric("GC deviation", _metric_text(quality.get("gc_deviation_from_0_5")))
        m3.metric("|Polarization γ|", _metric_text(ising.get("abs_polarization_gamma")))
        m4.metric("Alpha RMSE", _metric_text(logic.get("alpha_rmse_from_0_5")))
        m5.metric("Beta RMSE", _metric_text(logic.get("beta_rmse_from_0_125")))
        st.caption("For entropy, values closer to 1 are better. For deviation, |γ| and RMSE metrics, values closer to 0 are better.")

        tab1, tab2, tab3, tab4, tab5 = st.tabs([
            "Base composition", "Translational MSD", "Inverse Ising", "Logic α / β", "Summary"
        ])

        with tab1:
            base_df = pd.DataFrame({
                "Base": ["A", "C", "G", "T"],
                "Fraction": [
                    _safe_float(quality.get("a_fraction")),
                    _safe_float(quality.get("c_fraction")),
                    _safe_float(quality.get("g_fraction")),
                    _safe_float(quality.get("t_fraction")),
                ],
            })
            fig = go.Figure(go.Bar(x=base_df["Base"], y=base_df["Fraction"], text=base_df["Fraction"].map(lambda x: f"{x:.3f}"), textposition="auto"))
            _reference_line(fig, 0.25, "Ideal = 0.25")
            fig.update_layout(title="Nucleotide composition", xaxis_title="Nucleotide", yaxis_title="Fraction", yaxis_range=[0, max(0.35, float(base_df["Fraction"].max()) * 1.2)])
            st.plotly_chart(fig, use_container_width=True, key=f"{state_prefix}_base_composition_chart")

        with tab2:
            lags = translational.get("lags_nt", [])
            msd = translational.get("msd", [])
            fitted = translational.get("msd_fitted", [])
            fig = go.Figure()
            fig.add_trace(go.Scatter(x=lags, y=msd, mode="lines+markers", name="Observed MSD"))
            fig.add_trace(go.Scatter(x=lags, y=fitted, mode="lines", name="Fitted MSD"))
            fig.update_layout(title="Translational active-particle MSD", xaxis_title="Lag (nt)", yaxis_title="Mean squared displacement")
            st.plotly_chart(fig, use_container_width=True, key=f"{state_prefix}_translation_msd_chart")
            t1, t2, t3 = st.columns(3)
            t1.metric("Velocity V", _metric_text(translational.get("velocity_V_L_per_nt"), 6))
            t2.metric("Diffusion D", _metric_text(translational.get("diffusion_D_L2_per_nt"), 6))
            t3.metric("Fit R²", _metric_text(translational.get("fit_r_squared"), 4))

        with tab3:
            labels = ["|γ|", "|λ|", "|h|", "|Neighbour corr.|"]
            values = [
                abs(_safe_float(ising.get("abs_polarization_gamma"))),
                abs(_safe_float(ising.get("abs_ising_lambda"))),
                abs(_safe_float(ising.get("abs_ising_h"))),
                abs(_safe_float(ising.get("nearest_neighbour_correlation"))),
            ]
            fig = go.Figure(go.Bar(x=labels, y=values, text=[f"{v:.4f}" for v in values], textposition="auto"))
            fig.update_layout(title="Local binary structure", xaxis_title="Metric", yaxis_title="Absolute value")
            st.plotly_chart(fig, use_container_width=True, key=f"{state_prefix}_ising_chart")
            st.caption("λ and h are standard pseudo-likelihood estimates, not the unavailable Supporting Information S4 closed-form formula.")

        with tab4:
            alpha = logic.get("alpha", {}) or {}
            beta = logic.get("beta", {}) or {}
            contexts = [f"{i:03b}" for i in range(8)]
            alpha_values = [alpha.get(ctx) for ctx in contexts]
            beta_values = [_safe_float(beta.get(ctx)) for ctx in contexts]

            fig_a = go.Figure(go.Bar(x=contexts, y=alpha_values, name="Alpha"))
            _reference_line(fig_a, 0.5, "Ideal = 0.5")
            fig_a.update_layout(title="Conditional output balance α", xaxis_title="3-bit context", yaxis_title="P(output = 0 | context)", yaxis_range=[0, 1])
            st.plotly_chart(fig_a, use_container_width=True, key=f"{state_prefix}_alpha_chart")

            fig_b = go.Figure(go.Bar(x=contexts, y=beta_values, name="Beta"))
            _reference_line(fig_b, 0.125, "Ideal = 0.125")
            fig_b.update_layout(title="Context occurrence balance β", xaxis_title="3-bit context", yaxis_title="Context fraction")
            st.plotly_chart(fig_b, use_container_width=True, key=f"{state_prefix}_beta_chart")

        with tab5:
            rows = [
                ("DNA length", int(_safe_float(quality.get("length_nt"))), "nt"),
                ("Base entropy normalized", _safe_float(quality.get("base_entropy_normalized")), ""),
                ("GC fraction", _safe_float(quality.get("gc_fraction")), ""),
                ("GC deviation from 0.5", _safe_float(quality.get("gc_deviation_from_0_5")), ""),
                ("Longest homopolymer", int(_safe_float(quality.get("longest_homopolymer"))), "nt"),
                ("Homopolymer segments ≥2", int(_safe_float(quality.get("homopolymer_segments_ge_2"))), "segments"),
                ("Translational velocity V", _safe_float(translational.get("velocity_V_L_per_nt")), "L/nt"),
                ("Translational diffusion D", _safe_float(translational.get("diffusion_D_L2_per_nt")), "L²/nt"),
                ("|Polarization γ|", _safe_float(ising.get("abs_polarization_gamma")), ""),
                ("|Inverse-Ising λ|", _safe_float(ising.get("abs_ising_lambda")), ""),
                ("|Inverse-Ising h|", _safe_float(ising.get("abs_ising_h")), ""),
                ("Alpha RMSE from 0.5", _safe_float(logic.get("alpha_rmse_from_0_5")), ""),
                ("Beta RMSE from 0.125", _safe_float(logic.get("beta_rmse_from_0_125")), ""),
            ]
            if "information_density_bits_per_nt" in quality:
                rows.insert(1, ("Information density", _safe_float(quality.get("information_density_bits_per_nt")), "bits/nt"))
            if rotational:
                rows.extend([
                    ("Angular velocity ω", _safe_float(rotational.get("angular_velocity_omega_rad_per_step")), "rad/step"),
                    ("Rotational diffusion DR", _safe_float(rotational.get("rotational_diffusion_DR_rad2_per_step")), "rad²/step"),
                    ("Classified four-mer fraction", _safe_float(rotational.get("classified_fourmer_fraction")), ""),
                ])
            st.dataframe(pd.DataFrame(rows, columns=["Metric", "Value", "Unit"]), use_container_width=True, hide_index=True)
            st.download_button(
                "Download randomness report (JSON)",
                data=pd.Series(report).to_json(force_ascii=False, indent=2).encode("utf-8"),
                file_name=f"{state_prefix}_dna_randomness_report.json",
                mime="application/json",
                key=f"{state_prefix}_download_randomness_report",
            )
            if rotational:
                st.warning("The rotational model is a paper-inspired approximation because the complete cube transition table is unavailable.")
