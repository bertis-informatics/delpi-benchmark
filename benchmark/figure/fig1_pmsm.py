import polars as pl

pep_df = pl.read_parquet(
    r"/data1/DelPi_Data/speclib/2021-03-16-reviewed-contam-UP000002311/peptide_df.parquet"
)

pep_index = pep_df.filter(pl.col("sequence_length") == 8).sample(1).item(0, "peptide_index")

peptidoform_index = pl.scan_parquet(
    r"/data1/DelPi_Data/speclib/2021-03-16-reviewed-contam-UP000002311/modification_df.parquet"
).filter(pl.col("peptide_index") == pep_index).collect().item(0, "peptidoform_index")


precursor_index = pl.scan_parquet(
    r"/data1/DelPi_Data/speclib/2021-03-16-reviewed-contam-UP000002311/precursor_df.parquet"
).filter(pl.col("peptidoform_index") == peptidoform_index).collect().item(0, "precursor_index")


speclib_df = pl.scan_parquet(
    r"/data1/DelPi_Data/speclib/2021-03-16-reviewed-contam-UP000002311/speclib_df.parquet"
).filter(pl.col("precursor_index") == precursor_index).collect()




# if self.is_prefix:
#             return self.cleavage_index + 1
#         return len(self.sequence) - self.cleavage_index - 1


tmp_df = (
    speclib_df
    .with_columns(
        fragment=pl.when(pl.col("is_prefix"))
        .then(pl.col("cleavage_index") + 1)
        .otherwise(10 - pl.col("cleavage_index") - 1)
    )
    .select(pl.col("mz", "fragment", "is_prefix", "charge", "predicted_intensity"))
    .with_columns(
        pl.when(pl.col("is_prefix"))
        .then("b" + pl.col("fragment").cast(str))
        .otherwise("y" + pl.col("fragment").cast(str))
        .alias("ion_type")
    )
    .sort("mz")
)

# Plot the MS/MS spectrum
import matplotlib.pyplot as plt
import numpy as np

# Convert to pandas for easier plotting
df_pandas = tmp_df.to_pandas()

# Create the plot
fig, ax = plt.subplots(figsize=(12, 8))

# Create stem plot (mirror plot style)
b_ions = df_pandas[df_pandas["is_prefix"] == True]
y_ions = df_pandas[df_pandas["is_prefix"] == False]

# Plot b-ions (positive intensity)
if not b_ions.empty:
    ax.stem(
        b_ions["mz"],
        b_ions["predicted_intensity"],
        linefmt="b-",
        markerfmt="bo",
        basefmt=" ",
        label="b-ions",
    )

    # Add ion labels with charge information
    for _, row in b_ions.iterrows():
        charge_suffix = "+" * int(row["charge"])
        label = f"{row['ion_type']}{charge_suffix}"
        ax.annotate(
            label,
            (row["mz"], row["predicted_intensity"]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            color="blue",
        )

# Plot y-ions (positive intensity as well)
if not y_ions.empty:
    ax.stem(
        y_ions["mz"],
        y_ions["predicted_intensity"],
        linefmt="r-",
        markerfmt="ro",
        basefmt=" ",
        label="y-ions",
    )

    # Add ion labels with charge information
    for _, row in y_ions.iterrows():
        charge_suffix = "+" * int(row["charge"])
        label = f"{row['ion_type']}{charge_suffix}"
        ax.annotate(
            label,
            (row["mz"], row["predicted_intensity"]),
            xytext=(0, 10),
            textcoords="offset points",
            ha="center",
            va="bottom",
            fontsize=10,
            color="red",
        )

# Customize the plot
ax.set_xlabel("m/z", fontsize=12)
ax.set_ylabel("Predicted Intensity", fontsize=12)
ax.set_title("MS/MS Spectrum - Fragment Ion Intensities", fontsize=14)
ax.grid(True, alpha=0.3)
ax.legend()

# Set y-axis limits with some padding (only positive values)
max_intensity = df_pandas["predicted_intensity"].max()
ax.set_ylim(0, max_intensity * 1.2)

plt.tight_layout()
plt.savefig("./benchmark/temp_figures/msms_spectrum_fragment_ions.png")

# Print the data for reference
print("Fragment Ion Data:")
print(tmp_df)
