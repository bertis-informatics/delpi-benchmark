tools = {
    "DIA": ["alphadia", "diabert", "diann-1.8.1", "delpi"],
    "DDA": ["msgf", "sage", "delpi"],
}

datasets = {
    "2023-LFQ": {
        "acq_method": "DIA",
        "entrapment": "ARATH",  # Arabidopsis
        "target": ["MOUSE", "YEAST"],
        # "fasta": "Fasta-Merged_MusScereArabIRT-Review_10090-Total_559292-Review_3702-iRT_Fusion-20220302.fasta",
    },
    "2023-LFQ-single-run": {
        "acq_method": "DIA",
        "entrapment": "ARATH",  # Arabidopsis
        "target": ["MOUSE", "YEAST"],
        # "fasta": "Fasta-Merged_MusScereArabIRT-Review_10090-Total_559292-Review_3702-iRT_Fusion-20220302.fasta",
    },
    "2022-NeatPlasma": {
        "acq_method": "DIA",
        "entrapment": "ECOLI",
        "target": "HUMAN",
        # "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
    },
    "2025-SCP": {
        "acq_method": "DIA",
        "entrapment": "ECOLI",
        "target": "HUMAN",
        # "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
    },
    "2024-Mouse-Phospho": {
        "acq_method": "DIA",
        "entrapment": "YEAST",
        "target": "MOUSE",
        # "fasta": "Fasta-Merged_MusScereIRT-Review_10090-Total_559292-iRT_Fusion-20211103.fasta",
    },
    "CPTAC_COAD": {
        "acq_method": "DDA",
        "entrapment": "ECOLI",
        "target": "HUMAN",
        # "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
    },
    "CPTAC_LUAD_Phospho": {
        "acq_method": "DDA",
        "entrapment": "ECOLI",
        "target": "HUMAN",
        # "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
    },
}

dda_datasets = {k: v for k, v in datasets.items() if v["acq_method"] == "DDA"}
dia_datasets = {k: v for k, v in datasets.items() if v["acq_method"] == "DIA"}
