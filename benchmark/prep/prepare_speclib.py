from delpi.data_prep.env_var import DELPI_DATA_DIR, FASTA_DATA_DIR
from delpi.database import PeptideDatabase


speclib_settings = [
    {
        "fasta": "Fasta-Merged_MusScereIRT-Review_10090-Total_559292-iRT_Fusion-20211103.fasta",
        "tmt": False,
        "phospho": False,
    },
    {
        "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
        "tmt": False,
        "phospho": False,
    },
    {
        "fasta": "Fasta-Merged_MusScereIRT-Review_10090-Total_559292-iRT_Fusion-20211103.fasta",
        "tmt": False,
        "phospho": True,
    },

    {
        "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
        "tmt": True,
        "phospho": False,
    },
    {
        "fasta": "2022-02-18-reviewed-UP000005640-UP000000625.fas",
        "tmt": True,
        "phospho": True,
    },
]
global_mods = [
    ("Carbamidomethyl", "C", "anywhere", True),
    ("Acetyl", "*", "protein_n_term", False),
    ("Oxidation", "M", "anywhere", False),
]
phospho_mods = [
    ("Carbamidomethyl", "C", "anywhere", True),
    ("Phospho", "S", "anywhere", False),
    ("Phospho", "T", "anywhere", False),
    ("Phospho", "Y", "anywhere", False),
]

tmt_global_mods = [
    ("Carbamidomethyl", "C", "anywhere", True),
    ("TMT6plex", "*", "peptide_n_term", True),
    ("TMT6plex", "K", "anywhere", True),
    ("Acetyl", "*", "protein_n_term", False),
    ("Oxidation", "M", "anywhere", False),
]

tmt_phospho_mods = [
    ("Carbamidomethyl", "C", "anywhere", True),
    ("TMT6plex", "*", "peptide_n_term", True),
    ("TMT6plex", "K", "anywhere", True),
    ("Phospho", "S", "anywhere", False),
    ("Phospho", "T", "anywhere", False),
    ("Phospho", "Y", "anywhere", False),
]



def generate_protein_db(speclib_setting, device_idx=0):

    fasta_file = speclib_setting["fasta"]
    fasta_path = FASTA_DATA_DIR / fasta_file

    is_tmt = speclib_setting["tmt"]
    is_phospho = speclib_setting["phospho"]

    dir_name = fasta_path.stem
    if is_tmt:
        dir_name += "_tmt"

    if is_phospho:
        dir_name += "_phospho"

    save_dir = DELPI_DATA_DIR / "speclib" / dir_name
    if save_dir.exists():
        print(f"Directory {save_dir} already exists. Skipping...")
        return

    if is_tmt:
        mod_params = tmt_phospho_mods if is_phospho else tmt_global_mods
    else:
        mod_params = phospho_mods if is_phospho else global_mods
    max_mods = 3 if is_phospho else 2

    db = PeptideDatabase()
    db = db.build(
        fasta_file=fasta_path,
        min_len=7,
        max_len=30,
        max_missed_cleavages=1,
        mod_param_set=mod_params,
        max_mods=max_mods,
        decoy="pseudo_reverse",
        device=f"cuda:{device_idx}"
    )
    db.save(save_dir, overwrite=True)
    print(f"Saved to {save_dir}")


if __name__ == "__main__":
    for speclib_setting in speclib_settings[3:]:
        print(speclib_setting)
        generate_protein_db(speclib_setting, device_idx=1)
