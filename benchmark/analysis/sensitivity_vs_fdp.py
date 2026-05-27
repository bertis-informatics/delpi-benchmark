import polars as pl
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from benchmark.dataset import datasets
from benchmark.result_reader import ResultReader
from benchmark.constant import tool_color_map, tool_display_name_map


#    for ds_name in list(datasets):
ds_name = "2023-LFQ-single-run"
print(ds_name)
reader = ResultReader(ds_name)

results_dict = reader.load()
# results_dict = reader.save()
# results_dict = reader.load()
# summary_df = reader.generate_summary(results_dict)
# summary_df.write_csv(r"./reports/" + f"{ds_name}_result_summary.csv")
self = reader

results_dict["delpi"]

entrap_db = self.dataset_config["entrapment"]
target_db = self.dataset_config["target"]
target_db = [target_db] if isinstance(target_db, str) else target_db
prot_r, pep_r = self.estimate_entrapment_ratio()

exclude_filter = pl.col("fasta_id").str.contains(entrap_db) & (
    pl.any_horizontal([pl.col("fasta_id").str.contains(p) for p in target_db])
)

filtered_results_dict = dict()

for tool, df in results_dict.items():
    tmp_df = df.filter(
        (pl.col("is_decoy") == False) & (pl.col("global_precursor_q_value") <= 0.1)
    )
    tmp_df = tmp_df.filter(~exclude_filter).with_columns(
        is_fp=pl.col("fasta_id").str.contains(entrap_db)
    )
    filtered_results_dict[tool] = tmp_df

filtered_results_dict["delpi"]

# filtered_results_dict["alphadia"].columns


df = pl.read_csv(
    "/data1/benchmark/DIA/2023-LFQ-single-run/diann-1.8.1/report.tsv", separator="\t"
)
df.select(["PEP", "CScore", "Q.Value", "Lib.Q.Value"]).sort(
    "PEP", descending=False
).head(20)

summary_dict = dict()
for tool, df in filtered_results_dict.items():
    counts = dict()
    ################ precursor counting ####################
    tmp_df = df.filter(
        (pl.col("is_decoy") == False) & (pl.col("global_precursor_q_value") <= 0.01)
    )
    num_precursors = tmp_df.n_unique(["modified_sequence", "precursor_charge"])

    tmp_df = tmp_df.filter(~exclude_filter).with_columns(
        is_fp=pl.col("fasta_id").str.contains(entrap_db)
    )
    tmp_df = tmp_df.unique(
        [
            "modified_sequence",
            "precursor_charge",
            "is_fp",
            "global_precursor_q_value",
        ]
    )

    T = tmp_df.filter(~pl.col("is_fp")).shape[0]
    E = tmp_df.filter(pl.col("is_fp")).shape[0]
    fdp_lb = E / (T + E)
    fdp_comb = (E * (1 + 1 / pep_r)) / (T + E)

    counts["precursors"] = num_precursors
    counts["FDP_lb"] = fdp_lb * 100
    counts["FDP_comb"] = fdp_comb * 100

    ################ protein group counting ####################
    tmp_df = df.filter(
        (pl.col("is_decoy") == False)
        & (pl.col("is_decoy_protein_group") == False)
        & (pl.col("global_protein_group_q_value") <= 0.01)
        & (pl.col("protein_group").is_not_null())
    )

    num_protein_groups = tmp_df["protein_group"].n_unique()

    # If a protein group included ≥2 proteins and at least
    # one of them was from the original target database,
    # it was taken as an original target protein group.
    # https://www.nature.com/articles/s41592-025-02719-x
    tmp_df = tmp_df.with_columns(
        is_fp=pl.col("protein_group")
        .str.split(";")
        .list.eval(pl.element().str.contains(entrap_db))
        .list.all()
    )
    tmp_df = tmp_df.unique(["protein_group", "is_fp"])
    T = tmp_df.filter(~pl.col("is_fp")).shape[0]
    E = tmp_df.filter(pl.col("is_fp")).shape[0]

    fdp_lb = E / (T + E)
    fdp_comb = (E * (1 + 1 / prot_r)) / (T + E)
    assert num_protein_groups == T + E

    counts["protein groups"] = num_protein_groups
    counts["FDP_lb_pg"] = fdp_lb * 100
    counts["FDP_comb_pg"] = fdp_comb * 100
    summary_dict[tool] = counts

summary_df = (
    pd.DataFrame.from_dict(summary_dict, orient="index")
    .reset_index()
    .rename(columns={"index": "tool"})
)
