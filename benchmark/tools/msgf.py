import polars as pl
from pathlib import Path
from typing import Dict, Any, Union

from benchmark.tools.base import BaseToolReader
from delpi.utils.fdr import calculate_q_value
from delpi import PROJECT_DIR


class MSGFReader(BaseToolReader):

    prot_infer_params = PROJECT_DIR / "benchmark/proteininfer/msgf_params.yaml"

    def get_output_files(self):
        return list(
            [
                f
                for f in Path(self.output_dir).glob("*.tsv")
                if not f.stem.startswith("protein_inference")
            ]
        )

    def convert_mzid_to_tsv(self) -> None:
        import subprocess

        # for run_name, input_file in zip(search_config.run_names, search_config.input_files):
        for f in Path(self.output_dir).glob("*.mzid"):
            run_name = f.stem
            input_file = self.output_dir / f"{run_name}.mzid"
            output_file = self.output_dir / f"{run_name}.tsv"
            cmd = [
                "java",
                "-Xmx3500M",
                "-cp",
                "/data/software/msgf/MSGFPlus.jar",
                "edu.ucsd.msjava.ui.MzIDToTsv",
                "-showDecoy",
                "1",
                "-i",
                f"{input_file}",
                "-o",
                f"{output_file}",
            ]
            result = subprocess.run(cmd, capture_output=True, text=True)
            print("STDOUT:", result.stdout)

    def read(self) -> pl.DataFrame:
        # file_path = r"/data1/benchmark/DDA/CPTAC_COAD/msgf"
        col_map = {
            "ScanNum": "frame_num",
            "Charge": "precursor_charge",
            "Peptide": "modified_sequence",
            "Protein": "fasta_id",
            "QValue": "precursor_q_value",
            "PepQValue": "peptide_q_value",
            "SpecEValue": "spec_e_value",
        }

        dfs = []
        for tsv_file in self.get_output_files():
            dfs.append(
                pl.read_csv(tsv_file, separator="\t", columns=list(col_map.keys()))
                .rename(col_map)
                .with_columns(pl.lit(tsv_file.stem).alias("run"))
            )

        df = pl.concat(dfs)

        df = df.with_columns(
            pl.col("fasta_id").str.contains("XXX_sp").alias("is_decoy")
        )

        # Remove cleavage info in fasta_id column (e.g. (pre=K,post=M))
        # Make sure unique proteins in fasta_id
        df = df.with_columns(
            pl.col("fasta_id").str.replace_all(r"\([^)]*\)", "")
        ).with_columns(pl.col("fasta_id").str.split(";").list.unique().list.join(";"))

        pre_df = df.group_by(["modified_sequence", "precursor_charge"]).agg(
            pl.col("spec_e_value", "is_decoy").sort_by("spec_e_value").first()
        )
        pre_df = calculate_q_value(
            pre_df,
            out_column="global_precursor_q_value",
            score_column="spec_e_value",
            score_sort_descending=False,
        )

        pep_df = df.group_by(["modified_sequence"]).agg(
            pl.col("spec_e_value", "is_decoy").sort_by("spec_e_value").first()
        )
        pep_df = calculate_q_value(
            pep_df,
            out_column="global_peptide_q_value",
            score_column="spec_e_value",
            score_sort_descending=False,
        )

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
                .str.replace_all(r"+15.995", "(UniMod:35)", literal=True)
                .str.replace_all(r"+42.011", "(UniMod:1)", literal=True)
                .str.replace_all(r"+79.966", "(UniMod:21)", literal=True)
                .str.replace_all(r"+57.021", "(UniMod:4)", literal=True)
                .str.replace_all(r"+229.163", "(UniMod:737)", literal=True)
            )
            # .with_columns(modified_sequence=("_" + pl.col("modified_sequence") + "."))
            # .with_columns(proteotypic=~pl.col("fasta_id").str.contains(";"))
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

    reader = MSGFReader(r"/data1/benchmark/DDA/CPTAC_COAD/msgf", seq_df)
    df = reader.read()
    df = reader.run_protein_inference()
    reader.count_precursors()
    reader.count_protein_groups()
