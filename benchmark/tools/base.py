import polars as pl
from pathlib import Path
from typing import Dict, Any, Union
from abc import ABC, abstractmethod

from pyproteininference.pipeline import ProteinInferencePipeline
from delpi.utils.yaml_file import load_yaml


class BaseToolReader(ABC):
    """Base class for tool-specific result readers"""

    prot_infer_params = None
    # prot_infer_score_column = "score"

    def __init__(self, output_dir: Union[str, Path], seq_df: pl.DataFrame):
        self.output_dir = Path(output_dir)
        self.seq_df = seq_df
        self.df = None

    @property
    def prot_infer_score_column(self) -> str:
        params = load_yaml(self.prot_infer_params)
        return params["parameters"]["score"]["psm_score"]

    @property
    def decoy_symbol(self) -> str:
        params = load_yaml(self.prot_infer_params)
        return params["parameters"]["identifiers"]["decoy_symbol"]

    def read_and_infer(self) -> pl.DataFrame:
        df = self.read()
        df = self.run_protein_inference()
        return df

    @abstractmethod
    def read(self) -> pl.DataFrame:
        """Read results from the tool's output file"""
        pass

    def run_protein_inference(self) -> pl.DataFrame:

        df = self.df
        output_dir = self.output_dir

        input_file = output_dir / "protein_inference_input.tsv"
        output_file = output_dir / "protein_inference_output.csv"

        if not output_file.exists():
            # make proteinIds with unique proteins
            interim_df = (
                df.filter(pl.col("global_precursor_q_value") <= 0.01)
                .select(
                    pl.col("modified_sequence").alias("peptide"),
                    pl.col(self.prot_infer_score_column),
                    pl.col("fasta_id").str.replace_all(";", " ").alias("proteinIds"),
                )
                .with_row_index("PSMId")
            )

            interim_df.write_csv(input_file, separator="\t")

            pipeline = ProteinInferencePipeline(
                parameter_file=self.prot_infer_params,
                combined_files=[input_file],
                output_filename=output_file,
            )
            pipeline.execute()

        df = self.merge_protein_infer_results(output_file)
        self.df = df

        return df

    def merge_protein_infer_results(self, output_file: Union[str, Path]) -> None:

        out_df = pl.read_csv(output_file)
        decoy_symbol = self.decoy_symbol

        pg_df = out_df.select(
            pl.col("Protein").str.split(" ").alias("protein_group"),
            pl.col("Q_Value").alias("global_protein_group_q_value"),
            pl.col("Peptides").str.split(" ").alias("modified_sequence"),
        )

        ## find protein groups consisting of only decoy proteins
        pg_df = pg_df.with_columns(
            pl.col("protein_group")
            .list.eval(pl.element().str.contains(decoy_symbol))
            .list.all()
            .alias("is_decoy_protein_group"),
        )

        pg_df = pg_df.with_columns(
            pl.col("protein_group").list.sort().list.join(";"),
        ).explode("modified_sequence")

        df = self.df.join(pg_df, on="modified_sequence", how="left")

        return df

    def count_precursors(self):
        df = self.df.filter(pl.col("is_decoy") == False)
        return (
            df.filter(pl.col("global_precursor_q_value") <= 0.01)
            .unique(["modified_sequence", "precursor_charge"])
            .shape[0]
        )

    def count_protein_groups(self):
        df = self.df.filter(pl.col("is_decoy") == False).filter(
            pl.col("is_decoy_protein_group") == False
        )

        return df.filter(pl.col("global_protein_group_q_value") <= 0.01)[
            "protein_group"
        ].n_unique()

    def protein_groups(self) -> pl.DataFrame:
        df = self.df.filter(pl.col("is_decoy") == False).filter(
            pl.col("is_decoy_protein_group") == False
        )

        return (
            df.filter(pl.col("global_protein_group_q_value") <= 0.01)["protein_group"]
            .unique()
            .to_list()
        )
