# DelPi Benchmark

Reproducible benchmark analysis for the DelPi (Deep Learning-based Peptide Identification) study. This repository contains all the code and analysis scripts used to evaluate DelPi's performance against other peptide identification tools.


### Tools Compared

| Tool Name | Version | Search Mode | Repository |
|-----------|---------|-------------|------------|
| DIA-NN | v1.8.1 | DIA | [github.com/vdemichev/DiaNN](https://github.com/vdemichev/DiaNN) |
| AlphaDIA | v1.12.1 | DIA | [github.com/MannLabs/alphadia](https://github.com/MannLabs/alphadia) |
| DIA-BERT | v1.1 | DIA | [github.com/guomics-lab/DIA-BERT](https://github.com/guomics-lab/DIA-BERT) |
| MS-GF+ | v2024.03.26 | DDA | [github.com/MSGFPlus/msgfplus](https://github.com/MSGFPlus/msgfplus) |
| Sage | v0.14.7 | DDA | [github.com/lazear/sage](https://github.com/lazear/sage) |
| DelPi | - | DIA & DDA | [github.com/bertis-informatics/delpi](https://github.com/bertis-informatics/delpi) |

## Datasets

| Name | Acquisition | Species | Instrument | Repository | Reference |
|------|-------------|---------|------------|------------|-----------|
| Dilution series | DIA | Mouse, Yeast | Q Exactive HF | [PXD034709](https://www.ebi.ac.uk/pride/archive/projects/PXD034709) | [Lou et al., Nat. Commun., 2023](https://doi.org/10.1038/s41467-023-39828-0) |
| Single cell | DIA | Human | Astral | [PXD049412](https://www.ebi.ac.uk/pride/archive/projects/PXD049412) | [Bubis et al., Nat. Methods, 2025](https://doi.org/10.1038/s41592-024-02544-w) |
| Neat Plasma | DIA | Human | Q Exactive HF-X | [PXD038394](https://www.ebi.ac.uk/pride/archive/projects/PXD038394) | [Scott et al., Nat. Commun., 2023](https://doi.org/10.1038/s41467-023-38488-4) |
| Phosphoproteomics | DIA | Mouse | Astral | [MSV000093613](https://massive.ucsd.edu/ProteoSAFe/dataset.jsp?task=MSV000093613) | [Lancaster et al., Nat. Commun., 2024](https://doi.org/10.1038/s41467-024-53985-y) |
| Global proteomics | DDA | Human | Q-Exactive Plus | [PDC000109](https://pdc.cancer.gov/pdc/study/PDC000109) | [Vasaikar et al., Cell, 2019](https://doi.org/10.1016/j.cell.2019.03.030) |
| Phosphoproteomics (TMT-labeled) | DDA | Human | Fusion Lumos | [PDC000149](https://pdc.cancer.gov/pdc/study/PDC000149) | [Gillette et al., Cell, 2020](https://doi.org/10.1016/j.cell.2020.06.052) |


## Setup with Conda

**Prerequisites:** [Anaconda](https://www.anaconda.com/download) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html) must be installed.

```bash
# Clone the repositories
git clone https://github.com/bertis-informatics/delpi-benchmark.git
git clone https://github.com/bertis-informatics/delpi.git

# Create conda environment with Python 3.12
conda create -n delpi-benchmark python=3.12 -y

# Activate the conda environment
conda activate delpi-benchmark

# Install DelPi in development mode
cd delpi
pip install -e .
cd ..

# Install benchmark package in development mode
cd delpi-benchmark
pip install -e .
```


## Running Benchmark Analyses


Each benchmark script is designed to be run independently. Make sure you have the required data files in the `data/` directory before running analyses.

```bash
# Example: Run label-free quantification benchmark
python benchmark/generate_benchmark_reports.py
```


## Authors
* Jungkap Park <jungkap.park@bertis.com>
