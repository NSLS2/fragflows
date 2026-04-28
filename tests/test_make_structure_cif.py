from pathlib import Path
from deposition.make_cifs import make_changed_state_cif
from deposition.load import ispyb_xml_to_cif_block
import pandas as pd

def test_load_ispyb_xml_to_cif_block():
    xml_path = "/nsls2/data4/amx/proposals/2025-3/pass-319624/319624-20250809-dtime/mx319624-1/cypd-1/1/FGZ-007_1/autoProcOutput/autoPROC.xml"
    assert Path(xml_path).exists()
    block = ispyb_xml_to_cif_block(xml_path)
    assert block is not None

def test_make_changed_state_cif():
    TEST_STRUCTURE = "cypd-1"
    refinement_df = pd.read_csv("cypd.20260424.refinement.csv")
    ligand_df = pd.read_csv("cypd.20251001.ligands_fixed.csv")
    doc = make_changed_state_cif(
        refinement_df,
        TEST_STRUCTURE,
        "data_template_Xray_GroupDep.20260423.cif",
        ligand_df,
    )

    doc.write_file(f"{TEST_STRUCTURE}.mmcif")

    assert doc is not None
    assert Path(f"{TEST_STRUCTURE}.mmcif").exists()