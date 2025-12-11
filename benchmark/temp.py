import polars as pl


df0 = pl.read_csv(
    r"/data1/benchmark/DIA/2023-LFQ/delpi.v2/pmsm_results.tsv", separator="\t"
)

df1 = pl.read_csv(
    r"/data1/benchmark/DIA/2023-LFQ/delpi/pmsm_results.tsv", separator="\t"
)

df2 = pl.read_csv(
    r"/data1/benchmark/DIA/2023-LFQ/diann-1.8.1/report.tsv", separator="\t"
)


df2.filter(
    (pl.col("Lib.Q.Value") <= 0.01)
    # (pl.col("Q.Value") <= 0.01)
).group_by("Run").len().sort("Run")


df0.filter((pl.col("global_precursor_q_value") <= 0.01))["precursor_index"].n_unique()
df1.filter((pl.col("global_precursor_q_value") <= 0.01))["precursor_index"].n_unique()

df0.filter((pl.col("precursor_q_value") <= 0.01))["precursor_index"].n_unique()
df1.filter((pl.col("precursor_q_value") <= 0.01))["precursor_index"].n_unique()

df2.filter((pl.col("Lib.Q.Value") <= 0.01))["Precursor.Id"].n_unique()

df2.filter((pl.col("Q.Value") <= 0.01))["Precursor.Id"].n_unique()
# df2.filter((pl.col("Global.Q.Value") <= 0.01))["Precursor.Id"].n_unique()


df2.filter(
    # (pl.col("Global.Q.Value") < 0.01) &
    (pl.col("Lib.Q.Value") < 0.01)
).select(pl.col("Q.Value", "Lib.Q.Value", "Global.Q.Value", "PEP"))


df2.filter(
    (pl.col("PG.Q.Value") <= 0.01) & (pl.col("Lib.PG.Q.Value") <= 0.01)
).group_by("Run").agg(pl.col("Protein.Group").n_unique()).sort("Run")


df2.filter((pl.col("PG.Q.Value") < 0.01) & (pl.col("Lib.PG.Q.Value") > 0.01))
