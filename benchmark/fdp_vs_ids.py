import polars as pl
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt


from delpi.database.fasta_parser import FastaParser
from delpi.search.config import SearchConfig
from benchmark.dataset import tools
from benchmark.tools.msgf import MSGFReader
from benchmark.tools.sage import SageReader
from benchmark.tools.diann import DIANNReader
from benchmark.tools.alphadia import AlphaDIAReader
from benchmark.tools.diabert import DIABertReader
from benchmark.tools.delpi import DelPiReader
from benchmark.dataset import datasets
from benchmark.result_reader import ResultReader
from delpi.utils.fdr import calculate_q_value


def plot_fdp_vs_ids():

    ds_name = "2023-LFQ"

    reader = ResultReader(ds_name)
    # results_dict = reader.save()
    results_dict = reader.load()

    # pass
    self = reader

    entrap_db = self.dataset_config["entrapment"]
    target_db = self.dataset_config["target"]
    target_db = [target_db] if isinstance(target_db, str) else target_db
    prot_r, pep_r = self.estimate_entrapment_ratio()

    exclude_filter = pl.col("fasta_id").str.contains(entrap_db) & (
        pl.any_horizontal([pl.col("fasta_id").str.contains(p) for p in target_db])
    )

    summary_dict = dict()
    for tool, df in results_dict.items():

        ################ precursor level ####################
        tmp_df = df.filter(
            (pl.col("is_decoy") == False)
            & (pl.col("global_precursor_q_value") <= 0.01)
        )

        tmp_df = tmp_df.filter(~exclude_filter).with_columns(
            is_fp=pl.col("fasta_id").str.contains(entrap_db)
        )
        
        tmp_df = tmp_df.unique(
            ["modified_sequence", "precursor_charge", "is_fp", "global_precursor_q_value"]
        )

        # if tool == "delpi":
        #     precursor_fdp_df = tmp_df.sort("score").with_columns(
        #         n_fp=pl.col("is_fp").cum_sum(),
        #         n_total=pl.int_range(1, pl.len() + 1, dtype=pl.UInt32)
        #     ).with_columns(
        #         fdp_lb_precursor=(pl.col("n_fp") / pl.col("n_total")).cast(pl.Float32),
        #         fdp_comb_precursor=((pl.col("n_fp") * (1 + 1 / pep_r)) / pl.col("n_total")).cast(pl.Float32)
        #     )
        #     precursor_fdp_df.unique("fdp_comb_precursor", maintain_order=True, keep="last")

        precursor_fdps = []
        for qval_cutoff in np.arange(0.001, 0.011, 0.001):
            T = tmp_df.filter(pl.col("global_precursor_q_value") <= qval_cutoff).filter(~pl.col("is_fp")).shape[0]
            E = tmp_df.filter(pl.col("global_precursor_q_value") <= qval_cutoff).filter(pl.col("is_fp")).shape[0]
            fdp_lb = E / (T + E)
            fdp_comb = (E * (1 + 1 / pep_r)) / (T + E)
            precursor_fdps.append((qval_cutoff, T + E, fdp_lb, fdp_comb))
        precursor_fdp_df = pl.DataFrame(precursor_fdps, schema=["q_value_cutoff", "n_total", "fdp_lb", "fdp_comb_precursor"])
        # Calculate cumulative FDP for each row at precursor level
        

        
        ################ protein group level ####################
        tmp_df = df.filter(
            (pl.col("is_decoy") == False)
            & (pl.col("is_decoy_protein_group") == False)
            & (pl.col("global_protein_group_q_value") <= 0.01)
            & (pl.col("protein_group").is_not_null())
        )

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

        # Calculate cumulative FDP for each row at protein group level
        protein_group_fdps = []
        for qval_cutoff in np.arange(0.001, 0.011, 0.001):
            T = tmp_df.filter(pl.col("global_protein_group_q_value") <= qval_cutoff).filter(~pl.col("is_fp")).shape[0]
            E = tmp_df.filter(pl.col("global_protein_group_q_value") <= qval_cutoff).filter(pl.col("is_fp")).shape[0]
            fdp_lb = E / (T + E)
            fdp_comb = (E * (1 + 1 / prot_r)) / (T + E)
            protein_group_fdps.append((qval_cutoff, T + E, fdp_lb, fdp_comb))
        protein_group_fdp_df = pl.DataFrame(protein_group_fdps, schema=["q_value_cutoff", "n_total", "fdp_lb", "fdp_comb_pg"])
        # protein_group_fdp_df = tmp_df.with_columns(
        #     n_fp=pl.col("is_fp").cum_sum(),
        #     n_total=pl.int_range(1, pl.len() + 1, dtype=pl.UInt32)
        # ).with_columns(
        #     fdp_lb_pg=(pl.col("n_fp") / pl.col("n_total")).cast(pl.Float32),
        #     fdp_comb_pg=((pl.col("n_fp") * (1 + 1 / prot_r)) / pl.col("n_total")).cast(pl.Float32)
        # )

        summary_dict[tool] = {
            "precursor_fdp_df": precursor_fdp_df, 
            "protein_group_fdp_df": protein_group_fdp_df
        }

    # Plot number of precursors vs combined FDP
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
    
    # Precursor level plot
    for tool, data in summary_dict.items():
        precursor_df = data["precursor_fdp_df"]
        ax1.plot(
            precursor_df["fdp_comb_precursor"].to_numpy(),
            precursor_df["n_total"].to_numpy(),
            label=tool,
            marker='o',
            markersize=4,
            linewidth=1.5
        )
    
    ax1.set_xlabel("Combined FDP", fontsize=12)
    ax1.set_ylabel("Number of Precursors", fontsize=12)
    ax1.set_title("Precursor Level: Number of IDs vs Combined FDP", fontsize=14)
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Protein group level plot
    for tool, data in summary_dict.items():
        protein_df = data["protein_group_fdp_df"]
        ax2.plot(
            protein_df["fdp_comb_pg"].to_numpy(),
            protein_df["n_total"].to_numpy(),
            label=tool,
            marker='s',
            markersize=4,
            linewidth=1.5
        )
    
    ax2.set_xlabel("Combined FDP", fontsize=12)
    ax2.set_ylabel("Number of Protein Groups", fontsize=12)
    ax2.set_title("Protein Group Level: Number of IDs vs Combined FDP", fontsize=14)
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig("benchmark/figures/ids_vs_fdp_comb.png")
    # plt.show()

        