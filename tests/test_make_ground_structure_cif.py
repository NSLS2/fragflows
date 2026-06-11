import os

from deposition.legacy import (
    build_ground_structure_cif_metadata_loops,
    get_mtz_input_from_dimple_log,
    build_ground_structure_factor_cif,
    build_ground_structure_cif,
    _collect_dimple_metadata,
)
from deposition.cif_blocks import dimple_mtz_to_cif_block
import pandas as pd
from pathlib import Path
import gemmi

group_dep_dir = ""
group_dep_set = [
    k[:-4]
    for k in os.listdir(
        group_dep_dir
    )
    if k.endswith(".cif") and "-sf.cif" not in k
]
print(group_dep_set)
PANDDA_ANALYSE_ALL_INFO = ""
all_data_combined_df = pd.read_csv(
    PANDDA_ANALYSE_ALL_INFO
)
unmodelled_xtal_ids = [
    d for d in all_data_combined_df["dtag"] if d not in group_dep_set
]
unmodelled_xtal_ids = sorted(unmodelled_xtal_ids, key=lambda x: int(x.split("-")[1]))
print(unmodelled_xtal_ids)

GROUP_DEP_CSV = ""
LIGAND_CSV = ""
TEMPLATE_CIF_PATH = ""
test_id1 = ""
test_id2 = ""
test_id3 = ""
test_target_name = ""
test_dimple_log = ""
test_dimple_mtz = ""
test_spacegroup = ""
dimple_dir = ""


def test_legacy_build_ground_structure_cif_loops():
    df = pd.read_csv(GROUP_DEP_CSV)
    unmodelled_xtal_ids = [test_id1, test_id2]
    new_block = build_ground_structure_cif_metadata_loops(df, unmodelled_xtal_ids)
    assert new_block.find_loop_item("_exptl_crystal.id").loop.width() == 3
    assert new_block is not None
    assert new_block.name == "XXXX"
    d = gemmi.cif.Document()
    d.add_copied_block(new_block)
    # d.write_file("test_ground_structure_cif.cif")
    # assert Path("test_ground_structure_cif.cif").exists()


def test_get_mtz_input_from_dimple_log():
    mtz_path = get_mtz_input_from_dimple_log(
        test_dimple_log
    )
    assert Path(mtz_path).exists()


def test_dimple_mtz_to_cif_block():
    mtz_path = test_dimple_mtz
    block = dimple_mtz_to_cif_block(
        mtz_path, gemmi.SpaceGroup(test_spacegroup), xtal_id=test_id3
    )
    assert block is not None
    # assert block.name == "cypd-116"


def test_legacy_build_ground_structure_factor_cif():
    df = pd.read_csv(GROUP_DEP_CSV)
    dimple_dir = ""
    ligand_csv = ""
    d = build_ground_structure_factor_cif(
        df,
        pd.read_csv(ligand_csv),
        dimple_dir,
        unmodelled_xtal_ids[:100],
    )
    assert d is not None
    d.write_file("ground_batch1-sf.cif")
    assert Path("ground_batch1-sf.cif").exists()

    d = build_ground_structure_factor_cif(
        df,
        pd.read_csv(ligand_csv),
        dimple_dir,
        unmodelled_xtal_ids[100:],
    )
    assert d is not None
    d.write_file("ground_batch2-sf.cif")
    assert Path("ground_batch2-sf.cif").exists()


def test_legacy_build_ground_structure_cif():
    df = pd.read_csv(GROUP_DEP_CSV)

    import os

    # unmodelled_xtal_ids = ["cypd-1", "cypd-119"]
    template_cif_path = TEMPLATE_CIF_PATH
    dimple_metadata = _collect_dimple_metadata(dimple_dir, df, unmodelled_xtal_ids)
    doc = build_ground_structure_cif(
        df,
        unmodelled_xtal_ids[:100],
        template_cif_path,
        dimple_metadata,
        target_name=test_target_name,
    )
    assert doc is not None
    doc.write_file("ground_batch1.cif")
    assert Path("ground_batch1.cif").exists()

    doc = build_ground_structure_cif(
        df,
        unmodelled_xtal_ids[100:],
        template_cif_path,
        dimple_metadata,
        target_name=test_target_name,
    )
    assert doc is not None
    doc.write_file("ground_batch2.cif")
    assert Path("ground_batch2.cif").exists()
