import polars as pl
import pandas as pd
import numpy as np
from pathlib import Path
from matplotlib import pyplot as plt
from matplotlib_venn import venn2, venn3

from benchmark.dataset import datasets
from benchmark.result_reader import ResultReader


def generate_id_performance_comparison_report(save_dir):

    save_dir = Path(r"/home/jungkap/workspace/delpi_root/delpi-benchmark/reports")
    # ds_name = '2024-Mouse-Phospho'
    # ds_name = "2023-LFQ-single-run"
    for ds_name in list(datasets):
        print(ds_name)
        fdp_threshold = 0.05 if ds_name == "2024-Mouse-Phospho" else 0.01
        reader = ResultReader(ds_name)
        results_dict = reader.save()
        # results_dict = reader.load()
        summary_df = reader.generate_summary(results_dict)
        summary_df.write_csv(save_dir / f"{ds_name}_result_summary.csv")

        # FDP-based summary: cut at 1% combined entrapment FDP per tool.
        summary_fdp_df = reader.generate_summary_by_fdp(
            results_dict, fdp_threshold=fdp_threshold
        )
        summary_fdp_df.write_csv(save_dir / f"{ds_name}_result_summary_by_fdp.csv")
