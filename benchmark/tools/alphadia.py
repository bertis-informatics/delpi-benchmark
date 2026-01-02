import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from delpi.chem.modification import Modification
from delpi.utils.fdr import calculate_q_value
from benchmark import PROJECT_DIR
from benchmark.tools.base import BaseToolReader


class AlphaDIAReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/alphadia_params.yaml"

    @staticmethod
    def get_modified_sequence(seq: str, mod_ids: str, mod_sites: str) -> str:
        if mod_ids is None or len(mod_ids) < 1:
            return seq

        seq = list(seq)
        mod_ids = mod_ids.split(";")
        mod_sites = map(int, mod_sites.split(";"))
        for mod_id, pos in zip(mod_ids, mod_sites):
            mod_name = mod_id.split("@")[0]
            unimod_id = Modification.get(mod_name).accession_num
            pos = len(seq) - 1 if pos == -1 else pos - 1
            seq[pos] += f"(UniMod:{unimod_id})"
        return "".join(seq)

    def read(self) -> pl.DataFrame:

        seq_df = self.seq_df
        output_dir = self.output_dir
        # file_path = r"/data1/benchmark/DIA/2022-NeatPlasma/alphadia"
        col_map = {
            "run": "run",
            "frame_center": "frame_num",
            "charge": "precursor_charge",
            "modified_sequence": "modified_sequence",
            "proteins": "protein_id",
            "decoy": "is_decoy",
            "precursor_idx": "precursor_index",
            "proba": "posterior_error",
            "score": "score",
        }

        dfs = [
            pl.read_parquet(psm_path)
            for psm_path in output_dir.glob("quant/**/psm.parquet")
        ]
        df = pl.concat(dfs)
        df = df.with_columns(
            pl.when(pl.col("mods").is_null())
            .then(pl.col("sequence"))
            .otherwise(
                pl.struct(["sequence", "mods", "mod_sites"]).map_elements(
                    lambda x: self.get_modified_sequence(
                        x["sequence"],
                        x["mods"],
                        x["mod_sites"],
                    ),
                    return_dtype=pl.String,
                )
            )
            .alias("modified_sequence")
        )

        df = df.select(list(col_map)).rename(col_map).cast({"is_decoy": pl.Boolean})

        ## Calculate global precursor q-value
        pre_df = (
            df.select(pl.col("precursor_index", "posterior_error", "is_decoy"))
            .group_by(["precursor_index"])
            .agg(pl.all().sort_by("posterior_error").first())
        )

        pre_df = calculate_q_value(
            pre_df,
            out_column="global_precursor_q_value",
            score_column="posterior_error",
            score_sort_descending=False,
        )

        pep_df = df.group_by(["modified_sequence"]).agg(
            pl.col("posterior_error", "is_decoy").sort_by("posterior_error").first()
        )
        pep_df = calculate_q_value(
            pep_df,
            out_column="global_peptide_q_value",
            score_column="posterior_error",
            score_sort_descending=False,
        )

        df = df.join(
            pre_df.select(pl.col("precursor_index", "global_precursor_q_value")),
            on=["precursor_index"],
            how="left",
        ).join(
            pep_df.select(pl.col("modified_sequence", "global_peptide_q_value")),
            on=["modified_sequence"],
            how="left",
        )

        ## replace protein_id with fasta_id
        prot_id_to_fasta_id = {
            row[0]: row[1]
            for row in seq_df.select(pl.col("protein_id", "fasta_id")).iter_rows()
        }

        def map_fasta_id(protein_id_str, is_decoy) -> str:
            protein_index_list = protein_id_str.split(";")
            if is_decoy:
                fasta_ids = [
                    f"rev_{prot_id_to_fasta_id[pid]}" for pid in protein_index_list
                ]
            else:
                fasta_ids = [prot_id_to_fasta_id[pid] for pid in protein_index_list]
            return ";".join(fasta_ids)

        df = df.with_columns(
            fasta_id=pl.struct(["protein_id", "is_decoy"]).map_elements(
                lambda x: map_fasta_id(
                    x["protein_id"],
                    x["is_decoy"],
                ),
                return_dtype=pl.String,
            ),
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

    reader = AlphaDIAReader(r"/data1/benchmark/DIA/2022-NeatPlasma/alphadia", seq_df)
    df = reader.read()
    df = reader.run_protein_inference()
    # reader.df = df
    reader.count_protein_groups()

    # df.filter(pl.col("global_precursor_q_value") <= 0.01).unique(["modified_sequence", "precursor_charge"]).shape
