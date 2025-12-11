import polars as pl
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from benchmark.dataset import datasets
from benchmark.result_reader import tool_display_name_map, tool_color_map, ResultReader
from delpi.utils.fdr import calculate_q_value


def gen_venn_diagrams_by_plain_peptide(ds_name):
    # venn diagram
    # ds_name = "2023-LFQ-single-run"
    # ds_name = "CPTAC_LUAD_Phospho"
    ds_name = "2024-Mouse-Phospho"
    reader = ResultReader(ds_name)
    results_dict = reader.load()

    set_list = []
    set_labels = []
    set_colors = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        # q_val_col = (
        #     "global_precursor_q_value"
        #     if tool.startswith("diann")
        #     else "global_peptide_q_value"
        # )
        pep_df = df.with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"\([^)]*\)", "")
            .alias("plain_peptide")
        )

        if tool == "diann-1.8.1":
            pep_df = pep_df.filter(
                (pl.col("is_decoy") == False)
                & (pl.col("global_precursor_q_value") <= 0.01)
            )
        else:
            if tool == "msgf":
                score_col = "spec_e_value"
            else:
                score_col = "posterior_error"

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
                (pl.col("is_decoy") == False) & (pl.col("global_pep_q_value") <= 0.01)
            )
        modified_sequences = set(pep_df["plain_peptide"].unique())

        set_list.append(modified_sequences)
        set_labels.append(tool_display_name_map[tool])
        set_colors.append(tool_color_map[tool])

    plt.figure()
    # Draw venn2 or venn3 depending on number of tools included
    if len(set_list) == 3:
        venn_obj = venn3(
            set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
        )
    elif len(set_list) == 2:
        venn_obj = venn2(
            set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
        )
    else:
        raise ValueError(
            f"Expected 2 or 3 tool result sets, got {len(set_list)}: {set_labels}"
        )
    plt.savefig(rf"./benchmark/figures/{ds_name}_venn_diagram_plain_pep.png")


def gen_venn_diagrams(ds_name):
    # venn diagram
    # ds_name = "2023-LFQ-single-run"
    ds_name = "CPTAC_LUAD_Phospho"
    reader = ResultReader(ds_name)
    results_dict = reader.load()

    set_list = []
    set_labels = []
    set_colors = []
    for tool, df in results_dict.items():
        if tool == "diabert":
            continue
        q_val_col = (
            "global_precursor_q_value"
            if tool.startswith("diann")
            else "global_peptide_q_value"
        )

        tmp_df = df.filter((pl.col("is_decoy") == False) & (pl.col(q_val_col) <= 0.01))
        modified_sequences = set(tmp_df["modified_sequence"].unique())
        set_list.append(modified_sequences)
        set_labels.append(tool_display_name_map[tool])
        set_colors.append(tool_color_map[tool])

    plt.figure()
    # Draw venn2 or venn3 depending on number of tools included
    if len(set_list) == 3:
        venn_obj = venn3(
            set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
        )
    elif len(set_list) == 2:
        venn_obj = venn2(
            set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
        )
    else:
        raise ValueError(
            f"Expected 2 or 3 tool result sets, got {len(set_list)}: {set_labels}"
        )
    plt.savefig(rf"./benchmark/figures/{ds_name}_venn_diagram.png")


def venn_diagrams_all():
    for ds_name in list(datasets):
        print(ds_name)
        gen_venn_diagrams(ds_name)
