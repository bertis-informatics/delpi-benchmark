import polars as pl
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from benchmark.dataset import datasets
from benchmark.result_reader import ResultReader
from benchmark.constant import tool_color_map, tool_display_name_map
from delpi.utils.fdr import calculate_q_value
from benchmark import PROJECT_DIR


def _draw_and_save_venn(
    set_list,
    set_labels,
    set_colors,
    out_path,
    subset_fontsize=16,
    label_fontsize=20,
    fontweight="bold",
):
    plt.figure()
    if len(set_list) == 3:
        v = venn3(set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8)
        subset_ids = ("100", "010", "110", "001", "101", "011", "111")
    elif len(set_list) == 2:
        v = venn2(set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8)
        subset_ids = ("10", "01", "11")
    else:
        raise ValueError(
            f"Expected 2 or 3 tool result sets, got {len(set_list)}: {set_labels}"
        )

    # Style the numeric subset labels (inside the venn regions).
    for sid in subset_ids:
        t = v.get_label_by_id(sid)
        if t is not None:
            t.set_fontsize(subset_fontsize)
            t.set_fontweight(fontweight)

    # Style the set labels (outside, tool names).
    for t in v.set_labels or []:
        if t is not None:
            t.set_fontsize(label_fontsize)
            t.set_fontweight(fontweight)

    plt.savefig(out_path, bbox_inches="tight")
    plt.close()


def _plain_pep_score_col(tool):
    """Score column used for plain-peptide aggregation per tool."""
    if tool.startswith("diann"):
        return "global_precursor_q_value", True
    if tool == "msgf":
        return "spec_e_value", True
    return "posterior_error", True


def _venn_plain_pep_by_fdr(results_dict, fdr_threshold=0.01):
    """Original FDR-based plain-peptide sets per tool."""
    out = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        pep_df = df.with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"\([^)]*\)", "")
            .alias("plain_peptide")
        )
        if tool == "diann-1.8.1":
            pep_df = pep_df.filter(
                (pl.col("is_decoy") == False)
                & (pl.col("global_precursor_q_value") <= fdr_threshold)
            )
        else:
            score_col = "spec_e_value" if tool == "msgf" else "posterior_error"
            pep_df = pep_df.group_by("plain_peptide").agg(
                pl.col(score_col, "is_decoy").sort_by(score_col).first()
            )
            pep_df = calculate_q_value(
                pep_df,
                out_column="global_pep_q_value",
                score_column=score_col,
                score_sort_descending=False,
            )
            pep_df = pep_df.filter(
                (pl.col("is_decoy") == False)
                & (pl.col("global_pep_q_value") <= fdr_threshold)
            )
        out.append((tool, set(pep_df["plain_peptide"].unique())))
    return out


def _venn_plain_pep_by_fdp(reader, results_dict, fdp_threshold=0.01):
    """FDP-based plain-peptide sets per tool."""
    out = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        score_col, score_ascending = _plain_pep_score_col(tool)
        cutoff = reader.plain_peptide_fdp_cutoff(
            df,
            score_col,
            score_ascending=score_ascending,
            fdp_threshold=fdp_threshold,
        )
        pep_df = df.filter(pl.col("is_decoy") == False).with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"\([^)]*\)", "")
            .alias("plain_peptide")
        )
        pep_df = (
            pep_df.sort(score_col, descending=not score_ascending)
            .group_by("plain_peptide")
            .agg(pl.col(score_col).first().alias("score"))
        )
        if cutoff is None:
            pep_df = pep_df.head(0)
        elif score_ascending:
            pep_df = pep_df.filter(pl.col("score") <= cutoff)
        else:
            pep_df = pep_df.filter(pl.col("score") >= cutoff)
        out.append((tool, set(pep_df["plain_peptide"].unique())))
    return out


def _venn_precursor_by_fdr(results_dict, fdr_threshold=0.01):
    """Original FDR-based modified-sequence sets per tool."""
    out = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        q_val_col = (
            "global_precursor_q_value"
            if tool.startswith("diann")
            else "global_peptide_q_value"
        )
        tmp_df = df.filter(
            (pl.col("is_decoy") == False) & (pl.col(q_val_col) <= fdr_threshold)
        )
        out.append((tool, set(tmp_df["modified_sequence"].unique())))
    return out


def _venn_precursor_by_fdp(reader, results_dict, fdp_threshold=0.01):
    """FDP-based modified-sequence sets per tool."""
    out = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        score_col, ascending, cutoff = reader.precursor_fdp_cutoff(
            df, tool, fdp_threshold=fdp_threshold
        )
        if cutoff is None:
            modified_sequences = set()
        else:
            pass_filter = (
                (pl.col(score_col) <= cutoff)
                if ascending
                else (pl.col(score_col) >= cutoff)
            )
            tmp_df = df.filter((pl.col("is_decoy") == False) & pass_filter)
            modified_sequences = set(tmp_df["modified_sequence"].unique())
        out.append((tool, modified_sequences))
    return out


def _emit_venn(pairs, save_dir, filename):
    set_list, set_labels, set_colors = [], [], []
    for tool, s in pairs:
        set_list.append(s)
        set_labels.append(tool_display_name_map[tool])
        set_colors.append(tool_color_map[tool])
    _draw_and_save_venn(set_list, set_labels, set_colors, save_dir / filename)


def generate_venn_diagrams_by_plain_peptide(
    save_dir, ds_name, fdr_threshold=0.01, fdp_threshold=0.01
):
    """Generate both FDR- and FDP-based plain-peptide Venn diagrams."""
    reader = ResultReader(ds_name)
    results_dict = reader.load()

    fdr_pairs = _venn_plain_pep_by_fdr(results_dict, fdr_threshold)
    _emit_venn(fdr_pairs, save_dir, f"{ds_name}_venn_diagram_plain_pep.png")

    fdp_pairs = _venn_plain_pep_by_fdp(reader, results_dict, fdp_threshold)
    _emit_venn(fdp_pairs, save_dir, f"{ds_name}_venn_diagram_plain_pep_by_fdp.png")


def generate_venn_diagrams(save_dir, ds_name, fdr_threshold=0.01, fdp_threshold=0.01):
    """Generate both FDR- and FDP-based precursor (modified-sequence) Venn diagrams."""
    reader = ResultReader(ds_name)
    results_dict = reader.load()

    fdr_pairs = _venn_precursor_by_fdr(results_dict, fdr_threshold)
    _emit_venn(fdr_pairs, save_dir, f"{ds_name}_venn_diagram.png")

    fdp_pairs = _venn_precursor_by_fdp(reader, results_dict, fdp_threshold)
    _emit_venn(fdp_pairs, save_dir, f"{ds_name}_venn_diagram_by_fdp.png")


def venn_diagrams_all():
    save_dir = PROJECT_DIR / "reports/figures"

    save_dir.mkdir(exist_ok=True, parents=True)

    for ds_name in list(datasets):
        print(ds_name)
        generate_venn_diagrams(save_dir, ds_name)

    for ds_name in list(datasets):
        print(ds_name)
        generate_venn_diagrams_by_plain_peptide(save_dir, ds_name)
