from pathlib import Path
from deposition.make_cifs import make_changed_state_sf_cif
import pandas as pd


def test_make_changed_state_sf_cif():
    TEST_STRUCTURE = "cypd-1"
    refinement_df = pd.read_csv("cypd.20260424.refinement.csv")
    event_df = pd.read_csv("cypd.event_table.20260422.csv")
    doc = make_changed_state_sf_cif(
        refinement_df,
        event_df,
        TEST_STRUCTURE,
    )
    doc.write_file(f"{TEST_STRUCTURE}-sf.mmcif")
    assert doc is not None
    assert Path(f"{TEST_STRUCTURE}-sf.mmcif").exists()