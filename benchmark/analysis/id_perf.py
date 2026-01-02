import polars as pl
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from benchmark.dataset import datasets
from benchmark.result_reader import ResultReader
from benchmark.constant import tool_color_map, tool_display_name_map


def generate_id_performance_comparison_report():
    for ds_name in list(datasets):
        print(ds_name)
        reader = ResultReader(ds_name)
        results_dict = reader.save()
        # results_dict = reader.load()
        summary_df = reader.generate_summary(results_dict)
        summary_df.write_csv(r"./reports/" + f"{ds_name}_result_summary.csv")


# def generate_venn_diagrams(ds_name):
#     # venn diagram
#     # ds_name = "2023-LFQ-single-run"
#     reader = ResultReader(ds_name)
#     results_dict = reader.load()

#     set_list = []
#     set_labels = []
#     set_colors = []
#     for tool, df in results_dict.items():
#         if tool == "diabert":
#             continue
#         q_val_col = (
#             "global_precursor_q_value"
#             if tool.startswith("diann")
#             else "global_peptide_q_value"
#         )

#         tmp_df = df.filter((pl.col("is_decoy") == False) & (pl.col(q_val_col) <= 0.01))
#         modified_sequences = set(tmp_df["modified_sequence"].unique())
#         set_list.append(modified_sequences)
#         set_labels.append(tool_display_name_map[tool])
#         set_colors.append(tool_color_map[tool])

#     plt.figure()
#     # Draw venn2 or venn3 depending on number of tools included
#     if len(set_list) == 3:
#         venn_obj = venn3(
#             set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
#         )
#     elif len(set_list) == 2:
#         venn_obj = venn2(
#             set_list, set_labels=set_labels, set_colors=set_colors, alpha=0.8
#         )
#     else:
#         raise ValueError(
#             f"Expected 2 or 3 tool result sets, got {len(set_list)}: {set_labels}"
#         )
#     plt.savefig(rf"./benchmark/figures/{ds_name}_venn_diagram.png")


# def venn_diagrams_all():
#     for ds_name in list(datasets):
#         print(ds_name)
#         gen_venn_diagrams(ds_name)
