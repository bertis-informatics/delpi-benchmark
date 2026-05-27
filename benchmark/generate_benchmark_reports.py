from benchmark.analysis.id_perf import generate_id_performance_comparison_report
from benchmark.analysis.venn_diagram import (
    generate_venn_diagrams,
    generate_venn_diagrams_by_plain_peptide,
)
from benchmark import PROJECT_DIR
from benchmark.dataset import datasets

if __name__ == "__main__":

    ## ID performance comparison report
    save_dir = PROJECT_DIR / "reports"
    save_dir.mkdir(exist_ok=True, parents=True)
    generate_id_performance_comparison_report(save_dir)

    fig_save_dir = save_dir / "figures"
    fig_save_dir.mkdir(exist_ok=True, parents=True)
    for ds_name in list(datasets):
        print(ds_name)
        fdp_threshold = 0.05 if ds_name == "2024-Mouse-Phospho" else 0.01
        generate_venn_diagrams(fig_save_dir, ds_name, fdp_threshold=fdp_threshold)

    for ds_name in list(datasets):
        print(ds_name)
        fdp_threshold = 0.05 if ds_name == "2024-Mouse-Phospho" else 0.01
        generate_venn_diagrams_by_plain_peptide(
            fig_save_dir, ds_name, fdp_threshold=fdp_threshold
        )
