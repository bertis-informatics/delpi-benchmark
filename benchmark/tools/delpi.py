import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from delpi.chem.modification import Modification
from delpi.database.utils import get_modified_sequence
from delpi.search.config import SearchConfig
from delpi.search.result_manager import ResultManager
from benchmark.tools.base import BaseToolReader
from delpi import PROJECT_DIR


class DelPiReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/delpi_params.yaml"

    def read(self) -> pl.DataFrame:

        search_config = SearchConfig(self.output_dir / "params.yaml")
        seq_df = pl.read_parquet(search_config.db_dir / "sequence_df.parquet")

        fasta_id_list = seq_df["fasta_id"].to_list()

        file_path = self.output_dir / "pmsm.bak.parquet"
        df = pl.read_parquet(file_path).filter(pl.col("precursor_q_value") <= 0.01)

        df = df.with_columns(
            pl.col("peptide_index", "peptidoform_index"),
            pl.when(pl.col("mod_ids").is_null())
            .then(pl.col("peptide"))
            .otherwise(
                pl.struct(["peptide", "mod_ids", "mod_sites"]).map_elements(
                    lambda x: get_modified_sequence(
                        x["peptide"],
                        x["mod_ids"],
                        x["mod_sites"],
                        use_unimod_id=True,
                    ),
                    return_dtype=pl.String,
                ),
            )
            .alias("modified_sequence"),
        )
        df = df.with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"_", "")
            .str.replace_all(r"\.", ""),
        ).with_columns(posterior_error=1 - (1 / (1 + (-pl.col("score")).exp())))

        ## fasta_id mapping from protien_index
        def map_fasta_id(protein_index_list, is_decoy) -> str:
            if is_decoy:
                fasta_ids = [f"rev_{fasta_id_list[pid]}" for pid in protein_index_list]
            else:
                fasta_ids = [fasta_id_list[pid] for pid in protein_index_list]
            return ";".join(fasta_ids)

        df = df.with_columns(
            fasta_id=pl.struct(["protein_index", "is_decoy"]).map_elements(
                lambda x: map_fasta_id(
                    x["protein_index"],
                    x["is_decoy"],
                ),
                return_dtype=pl.String,
            ),
        )
        df = df.drop(
            [
                "protein_group",
                "master_protein",
                "protein_group_q_value",
                "global_protein_group_q_value",
                "protein_index",
            ]
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

    reader = DelPiReader(r"/data1/benchmark/DIA/2025-SCP/delpi", seq_df)
    df = reader.read()
    df = reader.run_protein_inference()

    reader.count_protein_groups()
    reader.count_precursors()
