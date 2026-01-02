from benchmark.analysis.id_perf import generate_id_performance_comparison_report
from benchmark.analysis.id_perf import (
    generate_venn_diagrams,
    generate_venn_diagrams_by_plain_peptide,
)
from benchmark import PROJECT_DIR
from benchmark.dataset import datasets


if __name__ == "__main__":
    save_dir = PROJECT_DIR / "reports/figures"

    save_dir.mkdir(exist_ok=True, parents=True)

    for ds_name in list(datasets):
        print(ds_name)
        generate_venn_diagrams(save_dir, ds_name)

    for ds_name in list(datasets):
        print(ds_name)
        generate_venn_diagrams_by_plain_peptide(save_dir, ds_name)
