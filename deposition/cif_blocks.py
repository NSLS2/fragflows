import gemmi
import numpy as np

# map mtz column labels from various programs to pdbx refln compliant names
# https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/


def validate_mtz_columns(mtz: gemmi.Mtz, mtz_columns: list):
    # Ensure that we are adding correct columns
    for column in mtz_columns:
        if column not in mtz.column_labels():
            raise ValueError(f"Mtz file: {mtz_path} is missing {column} label")


def populate_minimal_block_pairs(
    block: gemmi.cif.Block,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    wavelength: str,
    details: str,
    cell_precision: int = 4,
):

    # create block and populate header pairs
    block.set_pair("_cell.length_a", str(np.round(cell.a, cell_precision)))
    block.set_pair("_cell.length_b", str(np.round(cell.b, cell_precision)))
    block.set_pair("_cell.length_c", str(np.round(cell.c, cell_precision)))
    block.set_pair("_cell.angle_alpha", str(np.round(cell.alpha, cell_precision)))
    block.set_pair("_cell.angle_beta", str(np.round(cell.beta, cell_precision)))
    block.set_pair("_cell.angle_gamma", str(np.round(cell.gamma, cell_precision)))

    block.set_pair("_diffrn.id", "1")
    block.set_pair("_diffrn.crystal_id", "1")
    block.set_pair("_diffrn.details", f'"{details}"')

    block.set_pair("_diffrn_radiation_wavelength.id", "1")
    block.set_pair("_diffrn_radiation_wavelength.wavelength", wavelength)

    block.set_pair("_reflns_crystal.id", "1")
    block.set_pair("_reflns_scale.group_code", "1")

    block.set_pair("_symmetry.space_group_name_H-M", f'"{spacegroup.hm}"')
    block.set_pair("_symmetry.Int_Tables_number", f"{spacegroup.number}")

    return


def insert_pair_into_cif_block(block: gemmi.cif.Block, pair_key: str, *args):
    """
    Inserts new key-value pairs into a CIF block after a specified category.

    Traverses the CIF block and inserts one or more key-value
    pairs immediately after the category specified by
    `pair_key`. It preserves the order of the original items and avoids
    overwriting existing keys. Outputs a copy of the CIF block with inserted values.

    Parameters:
        block (gemmi.cif.Block): The original CIF block to modify.
        pair_key (str): The category name (e.g., "struct") after which the new pairs will be inserted.
        *args (tuple): One or more key-value tuples representing the new pairs to insert,
                       e.g., ("_struct.title", "My Title").

    Returns:
        gemmi.cif.Block: A new CIF block with the inserted key-value pairs.

    Raises:
        ValueError: If any of the new keys already exist in the block, or
                    if the specified `pair_key` category is not found, or
                    if CIF formatting errors (e.g., multiple categories in a loop) are detected.
    """
    existing_tags = {
        tag
        for category in block.get_mmcif_category_names()
        for tag in block.find_mmcif_category(category).tags
    }
    incoming_keys = {arg[0] for arg in args}
    if incoming_keys & existing_tags:
        raise ValueError("Attempting to overwrite existing key-value pair.")
    if pair_key not in [s.split(".")[0] for s in list(existing_tags)]:
        raise ValueError("Attempting to insert after non-existent key")

    new_block = gemmi.cif.Block(block.name)

    def get_item_category(item):
        if item.pair:
            return item.pair[0].split(".")[0]
        elif item.loop:
            tags = set(s.split(".")[0] for s in item.loop.tags)
            if len(tags) != 1:
                raise ValueError(
                    "Error in CIF loop table formatting: found multiple category keys."
                )
            return tags.pop()
        else:
            raise ValueError("Found CIF item that is neither pair nor loop")

    previous_category = get_item_category(block[0])
    new_block.add_item(block[0])

    for b in list(block)[1:]:
        current_category = get_item_category(b)

        if current_category != previous_category:
            if previous_category == pair_key:
                for arg in args:
                    new_block.set_pair(*arg)
            previous_category = current_category

        new_block.add_item(b)

    return new_block


def refinement_cif_to_cif_block(
    refinement_cif_path: str,
    block_name: str = "xxxxsf",
    details: str = "data from final ensemble refinement with ligand",
) -> gemmi.cif.Block:

    try:
        refinement_cif = gemmi.cif.read_file(refinement_cif_path)
        block = refinement_cif.sole_block()
    except RuntimeError as e:
        print(f"Caught RuntimeError: {e}")
    block.name = block_name
    refinement_block = insert_pair_into_cif_block(
        block,
        "_cell",
        ("_diffrn.id", "1"),
        ("_diffrn.crystal_id", "1"),
        ("_diffrn.details", f'"{details}"'),
    )
    return refinement_block


def event_map_to_cif_block(
    map_path: str,
    high_res: float,
    block_name: str = "xxxxBsf",
    details: str = "test event",
    wavelength: str = "0.9201",
) -> gemmi.cif.Block:

    # read map
    ccp4map = gemmi.read_ccp4_map(map_path, setup=True)
    ccp4map.grid.spacegroup = gemmi.SpaceGroup("P 1")
    ccp4map.update_ccp4_header()

    # inverse FFT
    fphi = gemmi.transform_map_to_f_phi(ccp4map.grid)
    sf = fphi.prepare_asu_data(dmin=high_res)

    # create block and populate header pairs
    block = gemmi.cif.Block(block_name)
    populate_minimal_block_pairs(
        block, sf.unit_cell, ccp4map.grid.spacegroup, wavelength, details
    )

    block.init_mmcif_loop(
        "_refln",
        [
            "crystal_id",
            "wavelength_id",
            "scale_group_code",
            "index_h",
            "index_k",
            "index_l",
            "pdbx_FWT",
            "pdbx_PHWT",
        ],
    )

    loop = block.find_loop("_refln.crystal_id").get_loop()

    for k in range(0, len(sf)):
        row = [
            "1",
            "1",
            "1",
        ]

        row += [str(j) for j in sf.miller_array[k, :]]

        row += [
            str(round(np.abs(sf.value_array[k]), 2)),
            str(round(np.degrees(np.angle(sf.value_array[k])) % 360, 2)),
        ]

        loop.add_row(row)

    return block


def original_mtz_to_cif_block(
    mtz_path: str,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    block_name: str = "xxxxAsf",
    details: str = "data from original reflections",
    mtz_columns: list = ["IMEAN", "SIGIMEAN", "F", "SIGF"],
) -> gemmi.cif.Block:

    mtz = gemmi.read_mtz_file(mtz_path)

    try:
        validate_mtz_columns(mtz, mtz_columns)
    except ValueError as e:
        print(f"Caught ValueError during validation: {e}")
        raise e

    # verify that cell dimensions are consistent between original mtz and refined cif
    if not np.all(
        np.isclose(
            [cell.a, cell.b, cell.c], [mtz.cell.a, mtz.cell.b, mtz.cell.c], atol=0.02
        )
    ):
        raise ValueError(
            "Inconsistent cell dimensions: "
            f"{[cell.a, cell.b, cell.c]} vs. "
            f"{[mtz.cell.a, mtz.cell.b, mtz.cell.c]}"
        )

    # translational symmetry operations may not be defined in original intensity file
    if mtz.spacegroup != spacegroup:
        print(
            "Refinement cif and original intensities have different spacegroups,"
            "but same point groups."
        )

    mtz_to_cif = gemmi.MtzToCif()
    mtz_to_cif.wavelength = mtz.dataset(1).wavelength
    mtz_to_cif.spec_lines = [
        "H H index_h",
        "K H index_k",
        "L H index_l",
        "? IMEAN J intensity_meas",
        "& SIGIMEAN Q intensity_sigma_meas",
        "? F F F_meas_au",
        "& SIGF Q F_meas_sigma_au",
    ]

    block = gemmi.cif.read_string(mtz_to_cif.write_cif_to_string(mtz)).sole_block()
    block.name = block_name

    # need to add new block pairs and then move them towards top of file
    block.set_pairs(
        "_diffrn.", {"id": "1", "crystal_id": "1", "details": f'"{details}"'}, raw=True
    )
    idx = [
        block.get_index(k)
        for k in ["_diffrn.id", "_diffrn.crystal_id", "_diffrn.details"]
    ]
    [block.move_item(j, k) for j, k in zip(idx, [10, 11, 12])]

    return block


def dimple_mtz_to_cif_block(
    mtz_path: str,
    spacegroup: gemmi.SpaceGroup,
    block_name: str = "xxxxsf",
    details: str = "ground state model xxxx",
    mtz_columns: list = ["F", "SIGF", "FC", "PHIC", "FWT", "PHWT", "FOM", "FreeR_flag"],
    default_rfree_value: np.int32 = 0,
) -> gemmi.cif.Block:

    mtz = gemmi.read_mtz_file(mtz_path)

    try:
        validate_mtz_columns(mtz, mtz_columns)
    except ValueError as e:
        print(f"Caught ValueError during validation: {e}")
        raise e

    if spacegroup != mtz.spacegroup:
        raise ValueError(f"Unexpected spacegroup {mtz.spacegroup} in {mtz_path}")

    mtz_to_cif = gemmi.MtzToCif()
    mtz_to_cif.spec_lines = [
        "H H index_h",
        "K H index_k",
        "L H index_l",
        "? F F F_meas_au",
        "& SIGF Q F_meas_sigma_au",
        "? FreeR_flag I pdbx_r_free_flag",
        "& FreeR_flag I status S",
        "? FC F F_calc_au",
        "& PHIC P phase_calc",
        "? FWT F pdbx_FWT",
        "& PHWT P pdbx_PHWT",
        "? DELFWT F pdbx_DELFWT",
        "& PHDELWT P pdbx_PHDELWT",
        "FOM W fom",
    ]

    mtz_to_cif.free_flag_value = default_rfree_value

    block = gemmi.cif.read_string(mtz_to_cif.write_cif_to_string(mtz)).sole_block()
    block.name = block_name

    # need to add new block pairs and then move them towards top of file
    block.set_pairs(
        "_diffrn.", {"id": "1", "crystal_id": "1", "details": f'"{details}"'}, raw=True
    )
    idx = [
        block.get_index(k)
        for k in ["_diffrn.id", "_diffrn.crystal_id", "_diffrn.details"]
    ]
    [block.move_item(j, k) for j, k in zip(idx, [10, 11, 12])]

    return block
