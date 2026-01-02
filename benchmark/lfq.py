from pathlib import Path
import polars as pl

from benchmark.result_reader import ResultReader
from benchmark import PROJECT_DIR
from delpi.search.dia.max_lfq import maxlfq


def read_lfq_experiment_design_df():
    exp_df = pl.read_csv(PROJECT_DIR / r"data/lfq_files.tsv", separator="\t")
    exp_df = exp_df.select(
        pl.col("File name").str.replace_all(".raw", "").alias("run_name"),
        pl.col("Experiment design").alias("experiment_design"),
    )
    return exp_df


def read_diann_pg_matrix():

    exp_df = read_lfq_experiment_design_df()
    run_names = exp_df["run_name"].to_list()

    reader = ResultReader("2023-LFQ")
    diann_reader = reader._tool_readers["diann-1.8.1"]
    pg_quant_df = diann_reader.read_diann_pg_matrix()
    pg_quant_df = pg_quant_df.filter(
        pl.col("protein_group")
        .str.split(";")
        .list.eval(pl.element().str.contains("MOUSE"))
        .list.all()
    )

    # run_cols = [c for c in pg_quant_df.columns if c != "protein_group"]
    pg_quant_df = pg_quant_df.unpivot(
        index=["protein_group"],
        on=run_names,
        variable_name="run_name",
        value_name="abundance",
    )

    return pg_quant_df


# q_df.shape
def read_delpi_pg_matrix():

    pmsm_df = (
        pl.read_parquet("/data1/benchmark/DIA/2023-LFQ/search_results/delpi.parquet")
        .filter(
            (pl.col("is_decoy") == False) & (pl.col("is_decoy_protein_group") == False)
        )
        .filter(
            pl.col("protein_group")
            .str.split(";")
            .list.eval(pl.element().str.contains("MOUSE"))
            .list.all()
        )
    )

    # Merge MS2 area info from final PMSM results
    final_pmsm_df = pl.read_csv(
        "/data1/benchmark/DIA/2023-LFQ/delpi/pmsm_results.tsv", separator="\t"
    )

    pmsm_df = pmsm_df.join(
        final_pmsm_df.select(pl.col("run_index", "run_name", "pmsm_index", "ms2_area")),
        on=["run_index", "pmsm_index"],
        how="left",
    )

    ## collapse_precursors_to_peptides
    df = (
        pmsm_df.filter(pl.col("ms2_area").is_not_null() & (pl.col("ms2_area") > 0))
        .group_by(["protein_group", "peptide_index", "run_name"])
        .agg(
            pl.col("ms2_area").median().alias("peptide_abundance"),
            pl.len().alias("n_precursors_used"),
        )
    )

    pg_quant_df = maxlfq(df, run_col="run_name", min_peptides_per_protein=1)
    # pg_quant_df = maxlfq_like(df, intensity_col="peptide_abundance")

    return pg_quant_df


def compute_lfq_stats(
    pg_quant_df: pl.DataFrame,
    min_reps: int = 3,
    intensity_min: float = 0.0,
) -> tuple[pl.DataFrame, pl.DataFrame, pl.DataFrame]:
    """Compute LFQ stats per ratio group.

    Parameters
    - pg_quant_df: wide-format protein group matrix from DIA-NN reader.
    - exp_df: mapping DataFrame with columns [run_name, experiment_design]. If None, it will be loaded via read_lfq_experiment_design_df().
    - min_reps: minimum number of quantified replicates to include a protein group (default 3).
    - intensity_min: values <= this threshold are considered non-quantified.

    Returns
    - counts_df: per ratio group count of protein groups with >= min_reps quantifications.
      columns: [experiment_design, n_pg_ge_min_reps]
    - box_stats_df: per ratio group CV distribution summary for proteins with n>=min_reps.
      columns: [experiment_design, median, q1, q3, iqr, min, max]
    - cv_df: per-protein CV values used to generate box stats (for plotting if desired)
      columns: [experiment_design, protein_group, n_reps, mean, sd, cv]
    """

    exp_df = read_lfq_experiment_design_df()

    # 1) Long format and annotate with experiment design (ratio)
    long_df = pg_quant_df.join(exp_df, on="run_name", how="inner")

    # Keep only quantified observations (non-null and > intensity_min)
    q_long = long_df.filter(
        pl.col("abundance").is_not_null() & (pl.col("abundance") > intensity_min)
    )

    # 2) Count proteins with at least min_reps quantifications per ratio
    per_ratio_pg_counts = q_long.group_by(["experiment_design", "protein_group"]).agg(
        pl.len().alias("n_quants")
    )

    counts_df = (
        per_ratio_pg_counts.filter(pl.col("n_quants") >= min_reps)
        .group_by("experiment_design")
        .agg(pl.len().alias("quantifiable_proteins"))
        .sort("experiment_design")
    )

    # 3) Compute per-protein CV within each ratio (only proteins with >= min_reps)
    per_ratio_pg_stats = q_long.group_by(["experiment_design", "protein_group"]).agg(
        pl.len().alias("n_reps"),
        pl.col("abundance").mean().alias("mean"),
        pl.col("abundance").std(ddof=1).alias("sd"),
    )

    cv_df = (
        per_ratio_pg_stats.filter((pl.col("n_reps") >= min_reps) & (pl.col("mean") > 0))
        .with_columns((pl.col("sd") / pl.col("mean")).alias("cv"))
        .select("experiment_design", "protein_group", "n_reps", "mean", "sd", "cv")
    )

    # 4) Summarize CV distribution for boxplots per ratio group
    box_stats_df = (
        cv_df.group_by("experiment_design")
        .agg(
            pl.col("cv").median().alias("median"),
            pl.col("cv").quantile(0.25).alias("q1"),
            pl.col("cv").quantile(0.75).alias("q3"),
            (pl.col("cv").quantile(0.75) - pl.col("cv").quantile(0.25)).alias("iqr"),
            pl.col("cv").min().alias("min"),
            pl.col("cv").max().alias("max"),
            pl.len().alias("quantifiable_proteins"),
        )
        .sort("experiment_design")
    )

    return counts_df, cv_df, box_stats_df


def estimate_lfq_performance():

    save_dir = PROJECT_DIR / "reports"
    exp_order = pl.DataFrame(
        [
            ("Ref", 0),
            ("1:4", 1),
            ("1:2", 2),
            ("2:3", 3),
            ("1:1", 4),
            ("3:2", 5),
            ("2:1", 6),
        ],
        orient="row",
        schema={"experiment_design": pl.String, "order": pl.Int32},
    )

    pg_quant_df = read_diann_pg_matrix()
    counts_df1, cv_df1, box_stats_df1 = compute_lfq_stats(pg_quant_df)

    pg_quant_df = read_delpi_pg_matrix()
    counts_df2, cv_df2, box_stats_df2 = compute_lfq_stats(pg_quant_df)

    cv_df1.select(pl.col("experiment_design", "protein_group", "cv")).rename(
        {"cv": "cv_diann"}
    ).join(
        cv_df2.select(pl.col("experiment_design", "protein_group", "cv")).rename(
            {"cv": "cv_delpi"}
        ),
        on=["experiment_design", "protein_group"],
        how="full",
        coalesce=True,
    ).join(
        exp_order, on="experiment_design", how="left"
    ).sort(
        "order"
    ).filter(
        pl.col("experiment_design") != "Ref"
    ).select(
        pl.col("experiment_design", "protein_group"),
        pl.col("cv_diann") * 100,
        pl.col("cv_delpi") * 100,
    ).write_csv(
        # r"/data1/benchmark/DIA/2023-LFQ/search_results/cv.tsv",
        save_dir / "lfq_cv.tsv",
        separator="\t",
    )

    counts_df = (
        counts_df1.rename({"quantifiable_proteins": "DIA-NN"})
        .join(
            counts_df2.rename({"quantifiable_proteins": "DelPi"}),
            on="experiment_design",
        )
        .join(exp_order, on="experiment_design", how="left")
        .sort("order")
    )

    counts_df.write_csv(
        save_dir / "lfq_quantifiable_proteins.tsv",
        separator="\t",
    )


def estimate_absolute_error():

    save_dir = PROJECT_DIR / "reports"

    exp_order = pl.DataFrame(
        [
            ("Ref", 0),
            ("1:4", 1),
            ("1:2", 2),
            ("2:3", 3),
            ("1:1", 4),
            ("3:2", 5),
            ("2:1", 6),
        ],
        orient="row",
        schema={"experiment_design": pl.String, "order": pl.Int32},
    )

    exp_df = read_lfq_experiment_design_df()
    exp_df = exp_df.join(
        pl.DataFrame(
            [
                ("1:4", 0.25),
                ("1:2", 0.5),
                ("2:3", 2 / 3),
                ("1:1", 1.0),
                ("Ref", 1.0),
                ("3:2", 3 / 2),
                ("2:1", 2.0),
            ],
            orient="row",
            schema={"experiment_design": pl.String, "true_ratio": pl.Float64},
        ),
        on="experiment_design",
        how="left",
    )

    df1 = read_diann_pg_matrix()
    df2 = read_delpi_pg_matrix()

    df1 = df1.filter(pl.col("abundance").is_not_null()).join(
        exp_df, on="run_name", how="left"
    )
    df2 = df2.filter(pl.col("abundance").is_not_null()).join(
        exp_df, on="run_name", how="left"
    )

    df1 = df1.group_by(["experiment_design", "protein_group"]).agg(
        pl.col("abundance").median().alias("median_ab"), pl.col("true_ratio").first()
    )

    df2 = df2.group_by(["experiment_design", "protein_group"]).agg(
        pl.col("abundance").median().alias("median_ab"), pl.col("true_ratio").first()
    )

    ref_df1 = (
        df1.filter(pl.col("experiment_design") == "Ref")
        .rename({"median_ab": "ref_ab"})
        .select(["protein_group", "ref_ab"])
    )
    ref_df2 = (
        df2.filter(pl.col("experiment_design") == "Ref")
        .rename({"median_ab": "ref_ab"})
        .select(["protein_group", "ref_ab"])
    )

    error_df1 = (
        df1.filter(pl.col("experiment_design") != "Ref")
        .join(ref_df1, on="protein_group", how="inner")
        .with_columns(
            (pl.col("median_ab") / pl.col("ref_ab")).alias("observed_ratio"),
        )
        .select(
            pl.col("experiment_design", "protein_group"),
            (
                (pl.col("observed_ratio") - pl.col("true_ratio")).abs()
                / pl.col("true_ratio")
            ).alias("relative_error"),
        )
        .join(exp_order, on="experiment_design", how="left")
        .sort("order")
    )

    error_df2 = (
        df2.filter(pl.col("experiment_design") != "Ref")
        .join(ref_df2, on="protein_group", how="inner")
        .with_columns(
            (pl.col("median_ab") / pl.col("ref_ab")).alias("observed_ratio"),
        )
        .select(
            pl.col("experiment_design", "protein_group"),
            (
                (pl.col("observed_ratio") - pl.col("true_ratio")).abs()
                / pl.col("true_ratio")
            ).alias("relative_error"),
        )
        .join(exp_order, on="experiment_design", how="left")
        .sort("order")
    )

    err_df = (
        error_df1.join(
            error_df2,
            on=["experiment_design", "protein_group", "order"],
            how="full",
            suffix="_delpi",
            coalesce=True,
        )
        .join(exp_order, on="experiment_design", how="left")
        .sort("order")
    )

    (
        error_df2.group_by("experiment_design")
        .agg(
            pl.col("relative_error").median().alias("med_relative_error"),
            pl.col("order").first(),
        )
        .sort("order")["med_relative_error"]
        - error_df1.group_by("experiment_design")
        .agg(
            pl.col("relative_error").median().alias("med_relative_error"),
            pl.col("order").first(),
        )
        .sort("order")["med_relative_error"]
    )

    err_df.with_columns(
        (pl.col("relative_error") * 100).alias("relative_error_diann"),
        (pl.col("relative_error_delpi") * 100).alias("relative_error_delpi"),
    ).select(
        pl.col(
            "experiment_design",
            "protein_group",
            "relative_error_diann",
            "relative_error_delpi",
        )
    ).write_csv(
        save_dir / "lfq_absolute_error.tsv",
        separator="\t",
    )
