from deposition.make_cifs import assemble_group_changed_state_cifs
import pandas as pd
from pathlib import Path

group_dep_dir = "/nsls2/data/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/test_group_deposition_20260424"

def test_assemble_group_changed_state_cifs():
    refinement_df = pd.read_csv("cypd.20260430.refinement.csv")
    event_df = pd.read_csv("cypd.event_table.20260422.csv")
    ligand_df = pd.read_csv("cypd.20251001.ligands_fixed.csv")
    group_dep_dir = "../test_group_deposition_20260430"
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