import argparse
from pathlib import Path
import pandas as pd
import yaml
from deposition.make_cifs import assemble_group_changed_state_cifs
from validation.cif_checks import (
    has_ligand,
    count_ligands,
    count_event_blocks,
    cif_loop_check,
    diffrn_id_check,
    map_ligands_to_events,
    soaked_compound_check,
)

with open("config.yaml", "r") as yaml_file:
    config = yaml.safe_load(yaml_file)

cfg = config["groupdepflow"]
GROUP_DEP_DIR = cfg["groupdep_directory"]
REFINEMENT_CSV = cfg["refinement_csv"]
EVENT_CSV = cfg["event_csv"]
TEMPLATE_CIF = cfg["template_cif"]
LIGAND_CSV = cfg["ligand_csv"]
CHUNK_SIZE = 5


def run_assemble_group_changed_state_cifs(only_validate: bool = False):
    refinement_df = pd.read_csv(REFINEMENT_CSV)
    event_df = pd.read_csv(EVENT_CSV)
    ligand_df = pd.read_csv(LIGAND_CSV)
    group_dep_dir = GROUP_DEP_DIR
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--validate",
        action="store_true",
        help="Only run validation checks on existing CIF files; skip CIF assembly.",
    )
    args = parser.parse_args()
    run_assemble_group_changed_state_cifs(only_validate=args.validate)
