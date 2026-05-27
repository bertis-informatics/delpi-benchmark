import polars as pl
import pandas as pd
from pathlib import Path

from benchmark.dataset import tools
from benchmark.dataset import datasets
from benchmark.constant import tool_to_reader_class
from delpi.database.fasta_parser import FastaParser
from delpi.search.config import SearchConfig


def score_spec_for_tool(tool: str, level: str = "precursor"):
    """Return ``(score_col, ascending)`` for the FDP-based cutoff.

    `ascending=True` means lower score is better (q-value, e-value, PEP).
    DIA-NN does not expose a per-PSM raw score, so we fall back to the
    tool-reported global q-values for both precursor and protein levels.
    """
    if tool.startswith("diann"):
        col = (
            "global_precursor_q_value"
            if level == "precursor"
            else "global_protein_group_q_value"
        )
        return col, True
    if tool == "msgf":
        return "spec_e_value", True
    if tool == "diabert":
        return "score", False
    # delpi, sage, alphadia
    return "posterior_error", True


def cutoff_at_fdp(
    eval_df: pl.DataFrame,
    ratio: float,
    fdp_threshold: float = 0.01,
    ascending: bool = True,
):
    """Find the largest q where running combined entrapment-FDP <= fdp_threshold.

    `eval_df` must have columns 'q' and 'is_fp' (bool). `ascending=True`
    means lower 'q' is better (q-values, e-values, PEP); `False` means
    higher is better (raw classifier scores).
    `ratio` is the entrapment/target ratio at the appropriate level.
    Returns a dict with keys 'q_cutoff', 'T', 'E', 'fdp_lb', 'fdp_comb' or
    None if no cutoff satisfies the threshold.

    Tie-aware: rows sharing the same score are grouped into a single bin
    and either all included or all excluded, so the result is independent
    of intra-bin ordering (important when many rows share the boundary,
    e.g. q=0).
    """
    if eval_df.height == 0:
        return None
    binned = (
        eval_df.group_by("q")
        .agg(
            T_inc=(~pl.col("is_fp")).cast(pl.Int64).sum(),
            E_inc=pl.col("is_fp").cast(pl.Int64).sum(),
        )
        .sort("q", descending=not ascending)
        .with_columns(
            T=pl.col("T_inc").cum_sum(),
            E=pl.col("E_inc").cum_sum(),
        )
        .with_columns(
            fdp_comb=(pl.col("E") * (1.0 + 1.0 / ratio)) / (pl.col("T") + pl.col("E"))
        )
    )
    valid = binned.filter(pl.col("fdp_comb") <= fdp_threshold)
    if valid.height == 0:
        return None
    # Cutoff = the boundary bin (largest q for ascending, smallest for descending).
    last = valid.sort("q", descending=not ascending).tail(1)
    T = int(last["T"][0])
    E = int(last["E"][0])
    return {
        "q_cutoff": last["q"][0],
        "T": T,
        "E": E,
        "fdp_lb": E / (T + E) if (T + E) > 0 else 0.0,
        "fdp_comb": float(last["fdp_comb"][0]),
    }


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

    def generate_summary_by_fdp(self, results_dict, fdp_threshold=0.01):
        """Count IDs at a target *combined* entrapment-FDP threshold.

        Unlike `generate_summary`, which trusts each tool's reported
        q-values at 1% FDR, this method re-ranks PSMs/protein groups by
        the tool-reported score (q-value) and finds the largest cutoff at
        which the entrapment-based combined FDP stays under
        `fdp_threshold` (default 1%). Counts are reported at that cutoff.
        """

        entrap_db = self.dataset_config["entrapment"]
        target_db = self.dataset_config["target"]
        target_db = [target_db] if isinstance(target_db, str) else target_db
        prot_r, pep_r = self.estimate_entrapment_ratio()

        exclude_filter = pl.col("fasta_id").str.contains(entrap_db) & (
            pl.any_horizontal([pl.col("fasta_id").str.contains(p) for p in target_db])
        )

        summary_dict = dict()
        for tool, df in results_dict.items():
            # tool = 'delpi'
            # df = results_dict[tool]
            counts = dict()

            ################ precursor counting ####################
            prec_score_col, prec_asc = score_spec_for_tool(tool, "precursor")
            prec_df = df.filter(pl.col("is_decoy") == False)
            # Collapse to one row per precursor with its best raw score.
            # Pick the fasta_id from the best-scoring row to keep the
            # entrapment / target classification consistent with that PSM.
            prec_df = (
                prec_df.sort(prec_score_col, descending=not prec_asc)
                .group_by(["modified_sequence", "precursor_charge"])
                .agg(
                    pl.col(prec_score_col).first().alias("q"),
                    pl.col("fasta_id").first(),
                )
            )

            prec_eval = prec_df.filter(~exclude_filter).with_columns(
                is_fp=pl.col("fasta_id").str.contains(entrap_db)
            )

            res = cutoff_at_fdp(
                prec_eval.select("q", "is_fp"),
                pep_r,
                fdp_threshold,
                ascending=prec_asc,
            )
            if res is None:
                counts["precursors"] = 0
                counts["FDP_lb"] = 0.0
                counts["FDP_comb"] = 0.0
                counts["q_cutoff"] = None
            else:
                # Total precursors at this cutoff includes the shared
                # target+entrapment precursors that were excluded from
                # the FDP evaluation.
                if prec_asc:
                    pass_filter = pl.col("q") <= res["q_cutoff"]
                else:
                    pass_filter = pl.col("q") >= res["q_cutoff"]
                num_precursors = prec_df.filter(pass_filter).height
                counts["precursors"] = num_precursors
                counts["FDP_lb"] = res["fdp_lb"] * 100
                counts["FDP_comb"] = res["fdp_comb"] * 100
                counts["q_cutoff"] = res["q_cutoff"]
            counts["score_col"] = prec_score_col

            ################ protein group counting ####################
            pg_score_col, pg_asc = score_spec_for_tool(tool, "protein_group")
            pg_df = df.filter(
                (pl.col("is_decoy") == False)
                & (pl.col("is_decoy_protein_group") == False)
                & (pl.col("protein_group").is_not_null())
            )
            # One row per protein group with its best raw score.
            pg_score_agg = (
                pl.col(pg_score_col).min() if pg_asc else pl.col(pg_score_col).max()
            )
            pg_df = (
                pg_df.group_by("protein_group")
                .agg(pg_score_agg.alias("q"))
                .with_columns(
                    is_fp=pl.col("protein_group")
                    .str.split(";")
                    .list.eval(pl.element().str.contains(entrap_db))
                    .list.all()
                )
            )

            res = cutoff_at_fdp(
                pg_df.select("q", "is_fp"),
                prot_r,
                fdp_threshold,
                ascending=pg_asc,
            )
            if res is None:
                counts["protein groups"] = 0
                counts["FDP_lb_pg"] = 0.0
                counts["FDP_comb_pg"] = 0.0
                counts["q_cutoff_pg"] = None
            else:
                counts["protein groups"] = res["T"] + res["E"]
                counts["FDP_lb_pg"] = res["fdp_lb"] * 100
                counts["FDP_comb_pg"] = res["fdp_comb"] * 100
                counts["q_cutoff_pg"] = res["q_cutoff"]
            counts["score_col_pg"] = pg_score_col

            summary_dict[tool] = counts

        summary_df = (
            pd.DataFrame.from_dict(summary_dict, orient="index")
            .reset_index()
            .rename(columns={"index": "tool"})
        )
        return pl.from_pandas(summary_df)

    def _entrap_filters(self):
        """Return (entrap_db, target_db, exclude_filter) reused across helpers."""
        entrap_db = self.dataset_config["entrapment"]
        target_db = self.dataset_config["target"]
        target_db = [target_db] if isinstance(target_db, str) else target_db
        exclude_filter = pl.col("fasta_id").str.contains(entrap_db) & (
            pl.any_horizontal([pl.col("fasta_id").str.contains(p) for p in target_db])
        )
        return entrap_db, target_db, exclude_filter

    def precursor_fdp_cutoff(
        self,
        df: pl.DataFrame,
        tool: str,
        fdp_threshold: float = 0.01,
    ):
        """Return ``(score_col, ascending, cutoff)`` at the precursor level
        where the combined entrapment FDP first reaches ``fdp_threshold``
        for the given tool df. ``cutoff`` is in the original score units
        (q-value, e-value, PEP, or raw classifier score) and is ``None``
        if no cutoff satisfies the threshold. Callers should filter with
        ``col <= cutoff`` when ``ascending`` is True, ``col >= cutoff``
        otherwise.
        """
        entrap_db, _, exclude_filter = self._entrap_filters()
        _, pep_r = self.estimate_entrapment_ratio()
        score_col, ascending = score_spec_for_tool(tool, "precursor")

        prec_df = (
            df.filter(pl.col("is_decoy") == False)
            .sort(score_col, descending=not ascending)
            .group_by(["modified_sequence", "precursor_charge"])
            .agg(
                pl.col(score_col).first().alias("q"),
                pl.col("fasta_id").first(),
            )
        )
        prec_eval = prec_df.filter(~exclude_filter).with_columns(
            is_fp=pl.col("fasta_id").str.contains(entrap_db)
        )
        res = cutoff_at_fdp(
            prec_eval.select("q", "is_fp"),
            pep_r,
            fdp_threshold,
            ascending=ascending,
        )
        cutoff = None if res is None else res["q_cutoff"]
        return score_col, ascending, cutoff

    def plain_peptide_fdp_cutoff(
        self,
        df: pl.DataFrame,
        score_col: str,
        score_ascending: bool = True,
        fdp_threshold: float = 0.01,
    ):
        """Return the plain-peptide-level score cutoff (in original score
        units) at which the combined entrapment FDP first reaches
        `fdp_threshold`. The best score per plain peptide is used.

        `score_ascending=True` means lower is better (q-value, e-value,
        PEP). Returns the cutoff value or None.
        """
        entrap_db, _, exclude_filter = self._entrap_filters()
        _, pep_r = self.estimate_entrapment_ratio()

        pep_df = df.filter(pl.col("is_decoy") == False).with_columns(
            pl.col("modified_sequence")
            .str.replace_all(r"\([^)]*\)", "")
            .alias("plain_peptide")
        )
        pep_df = (
            pep_df.sort(score_col, descending=not score_ascending)
            .group_by("plain_peptide")
            .agg(
                pl.col(score_col).first().alias("score"),
                pl.col("fasta_id").first(),
            )
        )
        pep_df = pep_df.with_columns(pl.col("score").alias("q"))
        pep_eval = pep_df.filter(~exclude_filter).with_columns(
            is_fp=pl.col("fasta_id").str.contains(entrap_db)
        )
        res = cutoff_at_fdp(
            pep_eval.select("q", "is_fp"),
            pep_r,
            fdp_threshold,
            ascending=score_ascending,
        )
        return None if res is None else res["q_cutoff"]
