from pathlib import Path
import pandas as pd
import yaml
from deposition.make_cifs import assemble_group_changed_state_cifs

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

cfg = config["groupdepflow"]
GROUP_DEP_DIR = cfg["groupdep_directory"]
REFINEMENT_CSV = cfg["refinement_csv"]
EVENT_CSV = cfg["event_csv"]
TEMPLATE_CIF = cfg["template_cif"]
LIGAND_CSV = cfg["ligand_csv"]
CHUNK_SIZE = 5


def run_assemble_group_changed_state_cifs():
    refinement_df = pd.read_csv(REFINEMENT_CSV)
    event_df = pd.read_csv(EVENT_CSV)
    ligand_df = pd.read_csv(LIGAND_CSV)
    group_dep_dir = GROUP_DEP_DIR
    assemble_group_changed_state_cifs(
        refinement_df,
        event_df,
        group_dep_dir,
        TEMPLATE_CIF,
        ligand_df)
    # check that files are generated for each xtal_id in refinement_df
    for xtal_id in refinement_df['xtal_id']:
        sf_cif_path = f"{group_dep_dir}/{xtal_id}-sf.cif"
        assert Path(sf_cif_path).exists(), f"missing {sf_cif_path}"
        cif_path = f"{group_dep_dir}/{xtal_id}.cif"
        assert Path(cif_path).exists(), f"missing {cif_path}"

if __name__ == "__main__":
    run_assemble_group_changed_state_cifs()
