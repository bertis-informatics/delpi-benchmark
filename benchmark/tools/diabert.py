import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from delpi.chem.modification import Modification
from delpi.utils.fdr import calculate_q_value
from benchmark.tools.base import BaseToolReader
from benchmark import PROJECT_DIR


class DIABertReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/diabert_params.yaml"

    def read(self) -> pl.DataFrame:

        seq_df = self.seq_df

        # interim per-precursor scoring file (target + decoy, all scores)
        eval_path = list(self.output_dir.rglob("**/finetune/output/result_*_eval.csv"))[
            0
        ]
        # filename: result_<run>.mzML_eval.csv
        run_name = eval_path.stem[len("result_") : -len(".mzML_eval")]

        df = pl.read_csv(eval_path)
        # decoy transition_group_id is prefixed with "DECOY_"; strip before
        # joining against the lib whose decoy peptides have no such prefix
        df = df.with_columns(
            modified_sequence=pl.col("transition_group_id")
            .str.head(-1)
            .str.strip_prefix("DECOY_"),
            precursor_charge=pl.col("transition_group_id").str.tail(1).cast(pl.Int32),
            is_decoy=(pl.col("label") == 0.0),
            run=pl.lit(run_name, dtype=pl.String),
            rt_in_seconds=pl.col("RT"),
        )

        # peptide -> protein mapping from spectral library (1:1 by design)
        lib_df = self._load_lib_mapping()
        df = df.join(
            lib_df,
            on=["modified_sequence", "precursor_charge"],
            how="left",
        )

        # precursor-level q-value from score (higher = better)
        df = calculate_q_value(
            df,
            score_column="score",
            score_sort_descending=True,
            out_column="precursor_q_value",
        )

        # join seq_df to get fasta_id (null for decoys: lib protein_id is "rev_sp|...")
        df = df.join(
            seq_df.select(pl.col("protein_id", "fasta_id")),
            on="protein_id",
            how="left",
        )

        # single-run: global q-value equals local
        df = df.with_columns(
            global_precursor_q_value=pl.col("precursor_q_value"),
        ).drop(["transition_group_id", "label", "file_name", "iRT", "RT"])

        self.df = df
        return df

    def _load_lib_mapping(self) -> pl.DataFrame:
        """Deduped (modified_sequence, precursor_charge) -> protein_id mapping
        from the DIA-BERT spectral library. Cached as parquet next to the tsv.
        """
        lib_tsv = self.output_dir / "diann-speclib.filtered.tsv"
        lib_pq = lib_tsv.with_suffix(".peptide_protein.parquet")
        if lib_pq.exists():
            return pl.read_parquet(lib_pq)

        lib_df = (
            pl.scan_csv(lib_tsv, separator="\t")
            .select(
                pl.col("FullUniModPeptideName").alias("modified_sequence"),
                pl.col("PrecursorCharge").cast(pl.Int32).alias("precursor_charge"),
                pl.col("ProteinID").alias("protein_id"),
            )
            .unique()
            .collect()
        )
        lib_df.write_parquet(lib_pq)
        return lib_df

    def run_protein_inference(self):
        # DIA-BERT speclib is built with unique peptide->protein mapping, so we
        # skip protein inference: every precursor's protein_id is its protein
        # group directly. Decoys (lib id like "rev_sp|...") don't match seq_df,
        # so we fall back to the raw protein_id for grouping.
        df = self.df.with_columns(
            protein_group=pl.coalesce(pl.col("fasta_id"), pl.col("protein_id")),
            is_decoy_protein_group=pl.col("is_decoy"),
        )

        # protein-level q-value: best precursor score per protein_group
        pg_df = df.group_by("protein_group").agg(
            pl.col("score").max().alias("score"),
            pl.col("is_decoy").any().alias("is_decoy"),
        )
        pg_df = calculate_q_value(
            pg_df,
            score_column="score",
            score_sort_descending=True,
            out_column="protein_group_q_value",
        ).select("protein_group", "protein_group_q_value")

        df = (
            df.join(pg_df, on="protein_group", how="left")
            .with_columns(
                global_protein_group_q_value=pl.col("protein_group_q_value"),
            )
            .drop("protein_id")
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
