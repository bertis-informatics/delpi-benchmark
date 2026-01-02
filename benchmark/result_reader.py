import polars as pl
import pandas as pd
from pathlib import Path

from benchmark.dataset import tools
from benchmark.dataset import datasets
from benchmark.constant import tool_to_reader_class
from delpi.database.fasta_parser import FastaParser
from delpi.search.config import SearchConfig


class ResultReader:
    """Unified result reader for mass spectrometry search tools"""

    def __init__(self, dataset_name: str):
        """
        Initialize the result reader with dataset configuration

        Args:
            dataset_name: Name of the dataset (key from datasets dict)
            dataset_config: Dataset configuration (value from datasets dict)
            base_data_dir: Base directory where tool results are stored
        """
        self.dataset_name = dataset_name
        self.dataset_config = datasets[dataset_name]
        self.acq_method = self.dataset_config["acq_method"]
        self.base_data_dir = Path(rf"/data1/benchmark/{self.acq_method}")
        ds_dir = self.base_data_dir / self.dataset_name
        self.search_config = SearchConfig(ds_dir / "delpi/params.yaml")
        fasta_file = self.search_config.fasta_file
        fasta_parser = FastaParser(fasta_file)

        seq_df = fasta_parser.parse()
        self.seq_df = seq_df.with_columns(
            pl.col("fasta_id").str.split("|").list.get(0).alias("db"),
            pl.col("fasta_id").str.split("|").list.get(1).alias("protein_id"),
            pl.col("fasta_id").str.split("|").list.get(2).alias("protein_name"),
        )

        tool_readers = dict()
        for tool in tools[self.acq_method]:
            if (ds_dir / tool).exists():
                tool_readers[tool] = tool_to_reader_class[tool](
                    ds_dir / tool, self.seq_df
                )

        self._tool_readers = tool_readers

    def load(self):
        results_dict = {}
        ds_dir = self.base_data_dir / self.dataset_name / "search_results"
        for tool, tool_reader in self._tool_readers.items():
            df = pl.read_parquet(ds_dir / f"{tool}.parquet")
            results_dict[tool] = df
        return results_dict

    def save(self):
        results_dict = {}
        ds_dir = self.base_data_dir / self.dataset_name / "search_results"
        ds_dir.mkdir(exist_ok=True, parents=True)
        for tool, tool_reader in self._tool_readers.items():
            df = tool_reader.read_and_infer()
            df.write_parquet(ds_dir / f"{tool}.parquet")
            results_dict[tool] = df

        return results_dict

    def estimate_entrapment_ratio(self):
        search_config = self.search_config

        entrap_db = self.dataset_config["entrapment"]
        # target_db = self.dataset_config["target"]
        # target_db = [target_db] if isinstance(target_db, str) else target_db

        # protein level entrapment ratio
        seq_df = (
            pl.scan_parquet(search_config.db_dir / "sequence_df.parquet")
            .select(pl.col("protein_index", "fasta_id"))
            .with_columns(
                is_entrapment_protein=pl.col("fasta_id").str.contains(entrap_db)
            )
            .collect()
        )

        N = seq_df.shape[0]
        E = seq_df.filter(pl.col("is_entrapment_protein")).shape[0]
        T = N - E
        seq_r = E / T

        # peptide/precursor level entrapment ratio
        pep_df = (
            pl.scan_parquet(search_config.db_dir / "peptide_df.parquet")
            .filter(pl.col("is_decoy") == False)
            .select(pl.col("peptide_index", "protein_index"))
            .collect()
        )

        pep_df = (
            pep_df.explode("protein_index")
            .join(seq_df, on="protein_index", how="left")
            .group_by("peptide_index")
            .agg(
                pl.col("fasta_id"),
                (pl.col("is_entrapment_protein").n_unique() > 1).alias(
                    "target_entrapment_overlap"
                ),
            )
            .with_columns(pl.col("fasta_id").list.join(";"))
            .filter(~pl.col("target_entrapment_overlap"))
        )

        # peptide or precursor level entrapment ratio
        N = pep_df.shape[0]
        E = pep_df.filter(pl.col("fasta_id").str.contains(entrap_db)).shape[0]
        T = N - E
        pep_r = E / T

        return seq_r, pep_r

    def generate_summary(self, results_dict):

        entrap_db = self.dataset_config["entrapment"]
        target_db = self.dataset_config["target"]
        target_db = [target_db] if isinstance(target_db, str) else target_db
        prot_r, pep_r = self.estimate_entrapment_ratio()

        exclude_filter = pl.col("fasta_id").str.contains(entrap_db) & (
            pl.any_horizontal([pl.col("fasta_id").str.contains(p) for p in target_db])
        )

        summary_dict = dict()
        for tool, df in results_dict.items():
            counts = dict()
            ################ precursor counting ####################
            tmp_df = df.filter(
                (pl.col("is_decoy") == False)
                & (pl.col("global_precursor_q_value") <= 0.01)
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
        return pl.from_pandas(summary_df)
