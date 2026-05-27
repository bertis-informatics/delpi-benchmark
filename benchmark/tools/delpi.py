import polars as pl

from delpi.database.utils import get_modified_sequence
from delpi.search.config import SearchConfig
from benchmark.tools.base import BaseToolReader
from benchmark import PROJECT_DIR

# def map_fasta_id(protein_index_list, is_decoy) -> str:
#     if is_decoy:
#         fasta_ids = [f"rev_{fasta_id_list[pid]}" for pid in protein_index_list]
#     else:
#         fasta_ids = [fasta_id_list[pid] for pid in protein_index_list]
#     return ";".join(fasta_ids)


class DelPiReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/delpi_params.yaml"

    def read(self) -> pl.DataFrame:

        search_config = SearchConfig(self.output_dir / "params.yaml")
        results_file = search_config.output_dir / "pmsm_results.parquet"

        if results_file.exists():
            df = pl.read_parquet(results_file)
        else:
            df = pl.read_csv(
                search_config.output_dir / "pmsm_results.tsv", separator="\t"
            )

        # add prefix ("rev_") to decoys
        df = df.with_columns(
            pl.when(pl.col("is_decoy"))
            .then(
                pl.col("fasta_id")
                .str.split(";")
                .list.eval("rev_" + pl.element())
                .list.join(";")
            )
            .otherwise(pl.col("fasta_id"))
            .alias("fasta_id"),
        )
        df = df.with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"_", "")
            .str.replace_all(r"\.", ""),
        )

        df = df.drop(
            [
                "protein_group",
                "master_protein",
                "protein_group_q_value",
                "global_protein_group_q_value",
                "protein_index",
            ],
            strict=False,
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
