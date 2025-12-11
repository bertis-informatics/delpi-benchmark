import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from delpi.chem.modification import Modification
from delpi.utils.fdr import calculate_q_value
from benchmark.tools.base import BaseToolReader
from delpi import PROJECT_DIR


class DIABertReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/diabert_params.yaml"

    def read(self) -> pl.DataFrame:

        seq_df = self.seq_df
        file_path = list(self.output_dir.rglob("**/*precursor.csv"))[0]
        run_name = file_path.parent.stem

        col_map = {
            "ProteinID": "protein_id",
            "RT": "rt_in_seconds",
        }

        # # df1 = pl.read_csv(r'/data1/benchmark/DIA/2023-LFQ/diabert/HF_20210517_M200-Y800_DIA_2/fdr_stats/HF_20210517_M200-Y800_DIA_2.mzML_score_0.2_precursor.csv')
        # # df1.filter(pl.col("q_value") <= 0.01).filter(pl.col("label") == 1).shape
        # df.join(
        #     df1.select(pl.col("transition_group_id", "label", "q_value")),
        #     left_on="PrecursorID",
        #     right_on="transition_group_id",
        #     how="left",
        # )

        df = pl.read_csv(file_path)
        df = (
            df.with_columns(
                precursor_charge=pl.col("PrecursorID").str.tail(1).cast(pl.Int32),
                modified_sequence=pl.col("PrecursorID").str.head(-1),
                run=pl.lit(run_name, dtype=pl.String),
                precursor_q_value=pl.lit(0, dtype=pl.Float32),
                protein_group_q_value=pl.lit(0, dtype=pl.Float32),
                global_precursor_q_value=pl.lit(0, dtype=pl.Float32),
                global_protein_group_q_value=pl.lit(0, dtype=pl.Float32),
                is_decoy=pl.lit(False, dtype=pl.Boolean),
            )
            .rename(col_map)
            .drop(
                [
                    "iRT",
                    "ProteinName",
                    "PrecursorQuant",
                    "PeptideSequence",
                    "FileName",
                    "PrecursorID",
                ]
            )
        )

        df = df.join(
            seq_df.select(pl.col("protein_id", "fasta_id")), on="protein_id", how="left"
        ).drop("protein_id")

        self.df = df

        return df

    def run_protein_inference(self):
        df = self.df.with_columns(
            protein_group=pl.col("fasta_id"), is_decoy_protein_group=pl.lit(False)
        )
        self.df = df
        return df


def test():
    from delpi.database.fasta_parser import FastaParser

    fasta_parser = FastaParser(
        # r"/data1/FASTA/Fasta-Merged_MusScereIRT-Review_10090-Total_559292-iRT_Fusion-20211103.fasta"
        r"/data1/benchmark/DIA/2023-LFQ/FASTA/Fasta-Merged_MusScereArabIRT-Review_10090-Total_559292-Review_3702-iRT_Fusion-20220302.fasta"
    )
    seq_df = fasta_parser.parse()
    seq_df = seq_df.with_columns(
        pl.col("fasta_id").str.split("|").list.get(0).alias("db"),
        pl.col("fasta_id").str.split("|").list.get(1).alias("protein_id"),
        pl.col("fasta_id").str.split("|").list.get(2).alias("protein_name"),
    )

    reader = DIABertReader(r"/data1/benchmark/DIA/2023-LFQ/diabert", seq_df)
    df = reader.read()
    df = reader.run_protein_inference()
    reader.count_protein_groups()
    reader.count_precursors()
