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


group_dep_set = [
    k[:-4]
    for k in os.listdir(
        "/nsls2/data/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/fragflows_group_deposition_20260505"
    )
    if k.endswith(".cif") and "-sf.cif" not in k
]
print(group_dep_set)
all_data_combined_df = pd.read_csv(
    "/nsls2/data/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/pandda_analyse_ccp4-7.0_20251016-1/analyses/all_datasets_info.csv"
)
unmodelled_xtal_ids = [
    d for d in all_data_combined_df["dtag"] if d not in group_dep_set
]
unmodelled_xtal_ids = sorted(unmodelled_xtal_ids, key=lambda x: int(x.split("-")[1]))
print(unmodelled_xtal_ids)


def test_legacy_build_ground_structure_cif_loops():
    df = pd.read_csv("cypd.20260422.group_dep.csv")
    unmodelled_xtal_ids = ["cypd-1", "cypd-87"]
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
        "/nsls2/data4/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/models_20251016-1/cypd-116/dimple.log"
    )
    assert Path(mtz_path).exists()


def test_dimple_mtz_to_cif_block():
    mtz_path = "/nsls2/data4/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/models_20251016-1/cypd-116/cypd-116.dimple.mtz"
    block = dimple_mtz_to_cif_block(
        mtz_path, gemmi.SpaceGroup("C2"), xtal_id="cypd-116"
    )
    assert block is not None
    # assert block.name == "cypd-116"


def test_legacy_build_ground_structure_factor_cif():
    df = pd.read_csv("cypd.20260422.group_dep.csv")
    dimple_dir = "/nsls2/data4/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/models_20251016-1"

    d = build_ground_structure_factor_cif(
        df,
        pd.read_csv("cypd.20251001.ligands_fixed.csv"),
        dimple_dir,
        unmodelled_xtal_ids[:100],
    )
    assert d is not None
    d.write_file("cypd_ground_batch1-sf.cif")
    assert Path("cypd_ground_batch1-sf.cif").exists()

    d = build_ground_structure_factor_cif(
        df,
        pd.read_csv("cypd.20251001.ligands_fixed.csv"),
        dimple_dir,
        unmodelled_xtal_ids[100:],
    )
    assert d is not None
    d.write_file("cypd_ground_batch2-sf.cif")
    assert Path("cypd_ground_batch2-sf.cif").exists()


def test_legacy_build_ground_structure_cif():
    df = pd.read_csv("cypd.20260422.group_dep.csv")
    dimple_dir = "/nsls2/data4/amx/proposals/2025-3/pass-319624/319624-cypd2-processing/models_20251016-1"

    import os

    # unmodelled_xtal_ids = ["cypd-1", "cypd-119"]
    template_cif_path = "data_template_Xray_GroupDep.20260423.cif"
    dimple_metadata = _collect_dimple_metadata(dimple_dir, df, unmodelled_xtal_ids)
    doc = build_ground_structure_cif(
        df,
        unmodelled_xtal_ids[:100],
        template_cif_path,
        dimple_metadata,
        target_name="Cyclophilin D",
    )
    assert doc is not None
    doc.write_file("cypd_ground_batch1.cif")
    assert Path("cypd_ground_batch1.cif").exists()

    doc = build_ground_structure_cif(
        df,
        unmodelled_xtal_ids[100:],
        template_cif_path,
        dimple_metadata,
        target_name="Cyclophilin D",
    )
    assert doc is not None
    doc.write_file("cypd_ground_batch2.cif")
    assert Path("cypd_ground_batch2.cif").exists()
