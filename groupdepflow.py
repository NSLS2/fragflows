import argparse
from pathlib import Path
import pandas as pd
import yaml
from deposition.make_cifs import assemble_group_changed_state_cifs
from deposition.cif_blocks import dimple_mtz_to_cif_block
from validation.cif_checks import (
    diffrn_ids_disjoint_check,
    ground_mtz_to_cif_check,
    has_ligand,
    count_ligands,
    count_event_blocks,
    cif_loop_check,
    diffrn_id_check,
    map_ligands_to_events,
    soaked_compound_check,
)
from collections import defaultdict
import os

import os

from deposition.legacy import (
    build_ground_structure_factor_cif,
    build_ground_structure_cif,
    _collect_dimple_metadata,
)


# Load yaml config parameters

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

cfg = config["groupdepflow"]
GROUP_DEP_DIR = cfg["groupdep_directory"]
REFINEMENT_CSV = cfg["refinement_csv"]
EVENT_CSV = cfg["event_csv"]
TEMPLATE_CIF = cfg["template_cif"]
LIGAND_CSV = cfg["ligand_csv"]
PANDDA_INPUT_DIR = cfg["pandda_input_dir"]
PANDDA_ALL_DATASETS_CSV = cfg["pandda_all_datasets_csv"]
GROUP_DEP_CSV = cfg["groupdep_csv"]
TARGET_NAME = cfg["target_name"] #fancy target name for the CIF metadata, e.g. "Cyclophilin D"
CHUNK_SIZE = 5

# load dataframes
refinement_df = pd.read_csv(REFINEMENT_CSV)
group_dep_df = pd.read_csv(GROUP_DEP_CSV)
event_df = pd.read_csv(EVENT_CSV)
ligand_df = pd.read_csv(LIGAND_CSV)
group_dep_dir = GROUP_DEP_DIR

def run_assemble_group_changed_state_cifs(only_validate: bool = False):

    if not only_validate:
        assemble_group_changed_state_cifs(
            refinement_df, event_df, group_dep_dir, TEMPLATE_CIF, ligand_df
        )
    # check that files are generated for each xtal_id in refinement_df
    for xtal_id in refinement_df["xtal_id"]:
        sf_cif_path = f"{group_dep_dir}/{xtal_id}-sf.cif"
        assert Path(sf_cif_path).exists(), f"missing {sf_cif_path}"
        cif_path = f"{group_dep_dir}/{xtal_id}.cif"
        assert Path(cif_path).exists(), f"missing {cif_path}"

        # validation checks
        assert has_ligand(cif_path, "UNL"), f"ligand UNL not found in {cif_path}"
        ligand_count = count_ligands(cif_path, "UNL")
        event_count = count_event_blocks(sf_cif_path)
        assert (
            event_count <= ligand_count
        ), f"more event blocks ({event_count}) than ligands ({ligand_count}) in {xtal_id}"
        cif_loop_check(cif_path)
        diffrn_id_check(sf_cif_path, xtal_id)
        diffrn_id_check(cif_path, xtal_id)
        map_ligands_to_events(cif_path, sf_cif_path)
        soaked_compound_check(sf_cif_path, ligand_df)


def create_ground_state_cifs():
    group_dep_set = [
    k[:-4]
    for k in os.listdir(
        GROUP_DEP_DIR
    )
    if k.endswith(".cif") and "-sf.cif" not in k
    ]
    print(group_dep_set)
    all_data_combined_df = pd.read_csv(
        PANDDA_ALL_DATASETS_CSV
    )
    unmodelled_xtal_ids = [
        d for d in all_data_combined_df["dtag"] if d not in group_dep_set
    ]
    unmodelled_xtal_ids = sorted(unmodelled_xtal_ids, key=lambda x: int(x.split("-")[1]))
    print(unmodelled_xtal_ids)

    # ground state identifier
    ground_state_identifiers = set([s.split("-")[0] for s in unmodelled_xtal_ids])

    # structure factor cifs, multi-block single file
    if len(ground_state_identifiers) != 1:
        raise ValueError(f"expected exactly one ground state identifier, found: {ground_state_identifiers}")
    ground_state_id = ground_state_identifiers.pop()

    d = build_ground_structure_factor_cif(
        group_dep_df,
        ligand_df,
        PANDDA_INPUT_DIR,
        unmodelled_xtal_ids,
    )
    assert d is not None
    d.write_file(f"{GROUP_DEP_DIR}/{ground_state_id}_ground-sf.cif")
    assert Path(f"{GROUP_DEP_DIR}/{ground_state_id}_ground-sf.cif").exists()

    # now the structure cif, which requires collecting metadata from the dimple logs and PDBs
    # there will be one block and one structure for all
    dimple_metadata = _collect_dimple_metadata(PANDDA_INPUT_DIR, group_dep_df, unmodelled_xtal_ids)

    ground_mtz_to_cif_check(
        dimple_metadata,
        f"{GROUP_DEP_DIR}/{ground_state_id}_ground-sf.cif"
    )

    doc = build_ground_structure_cif(
        group_dep_df,
        unmodelled_xtal_ids,
        TEMPLATE_CIF,
        dimple_metadata,
        target_name=TARGET_NAME,
    )
    assert doc is not None
    doc.write_file(f"{GROUP_DEP_DIR}/{ground_state_id}_ground.cif")
    assert Path(f"{GROUP_DEP_DIR}/{ground_state_id}_ground.cif").exists()


def create_group_dep_index(group_dep_dir: str):
    f_index = defaultdict(dict)
    group_dep_path = Path(group_dep_dir)

    for f in sorted(group_dep_path.glob("*.cif")):
        name = f.name
        if name.endswith("-sf.cif"):
            k = name[:-7]
            f_index[k]["sf"] = name
        else:
            k = name[:-4]
            f_index[k]["model"] = name

    missing = [
        k for k, v in f_index.items() if "model" not in v or "sf" not in v
    ]
    if missing:
        raise ValueError(f"missing model/sf CIF pair(s) for: {', '.join(sorted(missing))}")

    index_path = group_dep_path / "index.txt"
    with open(index_path, "w") as i:
        for k in sorted(f_index):
            v = f_index[k]
            i.writelines(f"label: {k}\nmodel: {v['model']}\nsf: {v['sf']}\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Only run validation checks on existing CIF files; skip CIF assembly.",
    )
    args = parser.parse_args()
    run_assemble_group_changed_state_cifs(only_validate=args.validate)
    create_ground_state_cifs()
    diffrn_ids_disjoint_check(GROUP_DEP_DIR)
    create_group_dep_index(GROUP_DEP_DIR)
