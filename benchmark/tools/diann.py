import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from benchmark.tools.base import BaseToolReader


class DIANNReader(BaseToolReader):

    def read(self) -> pl.DataFrame:

        seq_df = self.seq_df

        # file_path = r"/data1/benchmark/DIA/2025-SCP/diann-1.8.1/report.tsv"
        file_path = self.output_dir / "report.tsv"
        num_runs = pl.read_csv(file_path, separator="\t", columns=["Run"])[
            "Run"
        ].n_unique()

        col_map = {
            "Run": "run",
            "Modified.Sequence": "modified_sequence",
            "Precursor.Charge": "precursor_charge",
            "MS2.Scan": "frame_num",
            "Protein.Group": "protein_group",
            "Protein.Ids": "fasta_id",
        }

        if num_runs == 1:
            col_map.update(
                {
                    "Q.Value": "global_precursor_q_value",
                    "PG.Q.Value": "global_protein_group_q_value",
                }
            )
        else:
            col_map.update(
                {
                    "Lib.Q.Value": "global_precursor_q_value",
                    "Lib.PG.Q.Value": "global_protein_group_q_value",
                }
            )

        df = (
            pl.read_csv(file_path, separator="\t", columns=list(col_map.keys()))
            .rename(col_map)
            .with_columns(pl.lit(False).alias("is_decoy"))
        )

        prot_id_to_fasta_id = {
            row[0]: row[1]
            for row in seq_df.select(pl.col("protein_id", "fasta_id")).iter_rows()
        }

        def convert_prot_id(prot_id: str) -> str:
            ids = prot_id.split(";")
            fasta_ids = [
                prot_id_to_fasta_id[id_] if id_ in prot_id_to_fasta_id else id_
                for id_ in ids
            ]
            return ";".join(fasta_ids)

        df = df.with_columns(
            protein_group=pl.col("protein_group").map_elements(convert_prot_id)
        ).with_columns(fasta_id=pl.col("fasta_id").map_elements(convert_prot_id))

        self.df = df
        return df

    def run_protein_inference(self):
        return self.df.with_columns(is_decoy_protein_group=pl.lit(False))

    def read_diann_pg_matrix(self):

        seq_df = self.seq_df
        file_path = self.output_dir / "report.pg_matrix.tsv"
        q_df = pl.read_csv(file_path, separator="\t")

        col_map = {"Protein.Group": "protein_group"}
        col_map.update(
            {
                c: Path(c).stem
                for c in q_df.columns
                if c.startswith("/data1/benchmark/DIA/2023-LFQ/mzML")
            }
        )

        q_df = q_df.rename(col_map).select(sorted(list(col_map.values())))
        prot_id_to_fasta_id = {
            row[0]: row[1]
            for row in seq_df.select(pl.col("protein_id", "fasta_id")).iter_rows()
        }

        def convert_prot_id(prot_id: str) -> str:
            ids = prot_id.split(";")
            fasta_ids = [
                prot_id_to_fasta_id[id_] if id_ in prot_id_to_fasta_id else id_
                for id_ in ids
            ]
            return ";".join(fasta_ids)

        q_df = q_df.with_columns(
            protein_group=pl.col("protein_group").map_elements(convert_prot_id)
        )

        return q_df


def test():
    from delpi.database.fasta_parser import FastaParser

    fasta_parser = FastaParser(
        r"/data1/benchmark/DIA/2023-LFQ/FASTA/Fasta-Merged_MusScereArabIRT-Review_10090-Total_559292-Review_3702-iRT_Fusion-20220302.fasta",
        # r"/data1/FASTA/2022-02-18-reviewed-UP000005640-UP000000625.fas",
    )
    seq_df = fasta_parser.parse()
    seq_df = seq_df.with_columns(
        pl.col("fasta_id").str.split("|").list.get(0).alias("db"),
        pl.col("fasta_id").str.split("|").list.get(1).alias("protein_id"),
        pl.col("fasta_id").str.split("|").list.get(2).alias("protein_name"),
    )

    # reader = DIANNReader(r"/data1/benchmark/DIA/2025-SCP/diann-1.8.1", seq_df)
    # reader = DIANNReader(r"/data1/benchmark/DIA/2022-NeatPlasma/diann-1.8.1", seq_df)
    reader = DIANNReader(r"/data1/benchmark/DIA/2023-LFQ/diann-1.8.1", seq_df)

    df = reader.read()
    reader.count_protein_groups()
    reader.count_precursors()
