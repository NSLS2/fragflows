from deposition.make_cifs import assemble_group_changed_state_cifs
import pandas as pd
from pathlib import Path

group_dep_dir = ""
REFINEMENT_CSV = ""
EVENT_CSV = ""
LIGAND_CSV = ""

def test_assemble_group_changed_state_cifs():
    refinement_df = pd.read_csv(REFINEMENT_CSV)
    event_df = pd.read_csv(EVENT_CSV)
    ligand_df = pd.read_csv(LIGAND_CSV)
    assemble_group_changed_state_cifs(
        refinement_df,
        event_df,
        group_dep_dir,
        'data_template_Xray_GroupDep.20260423.cif',
        ligand_df)
    # check that files are generated for each xtal_id in refinement_df
    for xtal_id in refinement_df['xtal_id']:
        sf_cif_path = f"{group_dep_dir}/{xtal_id}-sf.cif"
        assert Path(sf_cif_path).exists(), f"missing {sf_cif_path}"
        cif_path = f"{group_dep_dir}/{xtal_id}.cif"
        assert Path(cif_path).exists(), f"missing {cif_path}"