import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from benchmark.tools.base import BaseToolReader
from benchmark import PROJECT_DIR
from delpi.utils.fdr import calculate_q_value


class SageReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/sage_params.yaml"

    def read(self) -> pl.DataFrame:

        col_map = {
            "scannr": "frame_num",
            "charge": "precursor_charge",
            "peptide": "modified_sequence",
            "proteins": "fasta_id",
            "spectrum_q": "precursor_q_value",
            "peptide_q": "peptide_q_value",
            "is_decoy": "is_decoy",
            "posterior_error": "posterior_error",
            "sage_discriminant_score": "score",
        }
        file_path = self.output_dir / "results.sage.parquet"
        df = (
            pl.read_parquet(file_path, columns=list(col_map))
            .rename(col_map)
            .with_columns(pl.col("posterior_error").exp())
        )

        pre_df = df.group_by(["modified_sequence", "precursor_charge"]).agg(
            pl.col("score", "is_decoy").sort_by("score").last()
        )
        pre_df = calculate_q_value(pre_df, out_column="global_precursor_q_value")

        pep_df = df.group_by(["modified_sequence"]).agg(
            pl.col("score", "is_decoy").sort_by("score").last()
        )
        pep_df = calculate_q_value(pep_df, out_column="global_peptide_q_value")

        df = df.join(
            pre_df.select(
                pl.col(
                    "modified_sequence", "precursor_charge", "global_precursor_q_value"
                )
            ),
            on=["modified_sequence", "precursor_charge"],
            how="left",
        ).join(
            pep_df.select(pl.col("modified_sequence", "global_peptide_q_value")),
            on=["modified_sequence"],
            how="left",
        )

        df = (
            df.with_columns(
                pl.col("modified_sequence")
                .str.replace_all(r"[+15.9949]", "(UniMod:35)", literal=True)
                .str.replace_all(r"[+42.0106]", "(UniMod:1)", literal=True)
                .str.replace_all(r"[+79.9663]", "(UniMod:21)", literal=True)
                .str.replace_all(r"[+57.0215]", "(UniMod:4)", literal=True)
                .str.replace_all(r"[+229.1629]", "(UniMod:737)", literal=True)
                .str.replace_all(r"-", "", literal=True)
            )
            # .with_columns(modified_sequence=("_" + pl.col("modified_sequence") + "."))
            # .with_columns(proteotypic=~pl.col("fasta_id").str.contains(";"))
            # .select(pl.exclude("is_decoy"))
        )
        self.df = df

        return df


def test():
    from delpi.database.fasta_parser import FastaParser

    fasta_parser = FastaParser(
        r"/data1/FASTA/2022-02-18-reviewed-UP000005640-UP000000625.fas"
    )
    seq_df = fasta_parser.parse()
    seq_df = seq_df.with_columns(
        pl.col("fasta_id").str.split("|").list.get(0).alias("db"),
        pl.col("fasta_id").str.split("|").list.get(1).alias("protein_id"),
        pl.col("fasta_id").str.split("|").list.get(2).alias("protein_name"),
    )

    reader = SageReader(r"/data1/benchmark/DDA/CPTAC_COAD/sage", seq_df)
    df = reader.read()
    df = reader.run_protein_inference()
    reader.count_precursors()
    reader.count_protein_groups()
