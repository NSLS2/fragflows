from pathlib import Path
from deposition.make_cifs import make_changed_state_cif
from deposition.load import ispyb_xml_to_cif_block
import pandas as pd

def test_load_ispyb_xml_to_cif_block():
    xml_path = ""
    assert Path(xml_path).exists()
    block = ispyb_xml_to_cif_block(xml_path)
    assert block is not None

def test_make_changed_state_cif():
    TEST_STRUCTURE = ""
    refinement_df = pd.read_csv("")
    ligand_df = pd.read_csv("")
    doc = make_changed_state_cif(
        refinement_df,
        TEST_STRUCTURE,
        "data_template_Xray_GroupDep.20260423.cif",
        ligand_df,
    )

    doc.write_file(f"{TEST_STRUCTURE}.mmcif")

    assert doc is not None
    assert Path(f"{TEST_STRUCTURE}.mmcif").exists()