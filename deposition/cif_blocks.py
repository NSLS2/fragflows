import gemmi
import numpy as np
from collections import defaultdict

# map mtz column labels from various programs to pdbx refln compliant names
# https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/


def validate_mtz_columns(mtz: gemmi.Mtz, mtz_columns: list):
    # Ensure that we are adding correct columns
    for column in mtz_columns:
        if column not in mtz.column_labels():
            raise ValueError(f"Mtz file: {mtz.title} is missing {column} label")


def populate_minimal_block_pairs(
    block: gemmi.cif.Block,
    cell: gemmi.UnitCell,
    spacegroup: gemmi.SpaceGroup,
    wavelength: str,
    details: str,
    cell_precision: int = 4,
    xtal_id: str = "1",
):
    diffrn_id = xtal_id

    # create block and populate header pairs
    block.set_pair("_cell.length_a", str(np.round(cell.a, cell_precision)))
    block.set_pair("_cell.length_b", str(np.round(cell.b, cell_precision)))
    block.set_pair("_cell.length_c", str(np.round(cell.c, cell_precision)))
    block.set_pair("_cell.angle_alpha", str(np.round(cell.alpha, cell_precision)))
    block.set_pair("_cell.angle_beta", str(np.round(cell.beta, cell_precision)))
    block.set_pair("_cell.angle_gamma", str(np.round(cell.gamma, cell_precision)))

    block.set_pair("_diffrn.id", diffrn_id)
    block.set_pair("_diffrn.crystal_id", xtal_id)
    block.set_pair("_diffrn.details", f'"{details}"')

    block.set_pair("_diffrn_radiation_wavelength.id", "1")
    block.set_pair("_diffrn_radiation_wavelength.wavelength", wavelength)

    block.set_pair("_symmetry.space_group_name_H-M", f'"{spacegroup.hm}"')
    block.set_pair("_symmetry.Int_Tables_number", f"{spacegroup.number}")

    return


def insert_pair_into_cif_block(block: gemmi.cif.Block, pair_key: str, *args)->gemmi.cif.Block:
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
    block_name: str = "XXXXsf",
    xtal_id: str = "1",
    crystal_treatment: str = "?",
    details: str = "data from final ensemble refinement",
) -> gemmi.cif.Block:
    
    diffrn_id = xtal_id

    try:
        refinement_cif = gemmi.cif.read_file(refinement_cif_path)
        block = refinement_cif.sole_block()
    except RuntimeError as e:
        print(f"Caught RuntimeError: {e}")
    block.name = block_name
    refinement_block = insert_pair_into_cif_block(
        block,
        "_cell",
        ("_diffrn.id", diffrn_id),
        ("_diffrn.crystal_id", xtal_id),
        ("_diffrn.crystal_treatment", f'"{crystal_treatment}"'),
        ("_diffrn.details", f'"{details}"'),
    )

    for item in refinement_block:
        if item.pair and '.crystal_id' in item.pair[0]:
            if item.pair[1] != xtal_id:
                raise ValueError(
                    f"Found inconsistent crystal_id values: {item.pair[1]} vs. {xtal_id}"
                )
        if item.loop and {'_refln'} == set([s.split('.')[0] for s in item.loop.tags]):
            if "_refln.crystal_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.crystal_id"], xtal_id)
            if "_refln.wavelength_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.wavelength_id"], "1")


    return refinement_block


def event_map_to_cif_block(
    map_path: str,
    high_res: float,
    wavelength: str,
    block_name: str = "XXXXBsf",
    details: str = "test event",
    xtal_id: str = "1",
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
        block, sf.unit_cell, ccp4map.grid.spacegroup, wavelength, details, xtal_id=xtal_id
    )

    block.init_mmcif_loop(
        "_refln",
        [
            "crystal_id",
            "wavelength_id",
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
            xtal_id,
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
    block_name: str = "XXXXAsf",
    details: str = "data from original reflections",
    mtz_columns: list = ["IMEAN", "SIGIMEAN", "F", "SIGF"],
    xtal_id: str = "1",
) -> gemmi.cif.Block:
    
    diffrn_id = xtal_id

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
        "& SIGIMEAN Q intensity_sigma",
        "? F F F_meas_au",
        "& SIGF Q F_meas_sigma_au",
    ]

    block = gemmi.cif.read_string(mtz_to_cif.write_cif_to_string(mtz)).sole_block()
    block.name = block_name

    # need to add new block pairs and then move them towards top of file
    block.set_pairs(
        "_diffrn.", {"id": diffrn_id, "crystal_id": xtal_id, "details": f'"{details}"'}, raw=True
    )
    idx = [
        block.get_index(k)
        for k in ["_diffrn.id", "_diffrn.crystal_id", "_diffrn.details"]
    ]
    [block.move_item(j, k) for j, k in zip(idx, [10, 11, 12])]

    for item in block:
        if item.loop and {'_refln'} == set([s.split('.')[0] for s in item.loop.tags]):
            if "_refln.crystal_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.crystal_id"], xtal_id)
            if "_refln.wavelength_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.wavelength_id"], "1")

    return block


def dimple_mtz_to_cif_block(
    mtz_path: str,
    spacegroup: gemmi.SpaceGroup | None=None,
    block_name: str = "XXXXsf",
    crystal_treatment: str = "?",
    diffrn_details: str = "not refined to convergence",
    mtz_columns: list = ["F", "SIGF", "FC", "PHIC", "FWT", "PHWT", "FOM", "FreeR_flag"],
    default_rfree_value: np.int32 = 0,
    xtal_id: str = "1",
) -> gemmi.cif.Block:

    mtz = gemmi.read_mtz_file(mtz_path)

    diffrn_id = xtal_id

    try:
        validate_mtz_columns(mtz, mtz_columns)
    except ValueError as e:
        print(f"Caught ValueError during validation: {e} for {xtal_id}")
        raise e

    if spacegroup.point_group_hm() != mtz.spacegroup.point_group_hm():
        raise ValueError(f"point groups don't match for {spacegroup} and {mtz.spacegroup} in {mtz_path}")

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
        "_diffrn.", {"id": diffrn_id, "crystal_id": xtal_id, "crystal_treatment": f'"{crystal_treatment}"', "details": f'"{diffrn_details}"'}, raw=True
    )
    idx = [
        block.get_index(k)
        for k in ["_diffrn.id", "_diffrn.crystal_id", "_diffrn.crystal_treatment", "_diffrn.details"]
    ]
    [block.move_item(j, k) for j, k in zip(idx, [10, 11, 12, 13])]

    # add crystal_id and wavelength_id columns to refln loop if not already present
    for item in block:
        if item.loop and {'_refln'} == set([s.split('.')[0] for s in item.loop.tags]):
            if "_refln.crystal_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.crystal_id"], xtal_id)
            if "_refln.wavelength_id" not in item.loop.tags:
                item.loop.add_columns(["_refln.wavelength_id"], "1")

    return block

def deduplicate_cif_loops(
        source_block: gemmi.cif.Block,
        target_block: gemmi.cif.Block,
        reference_tag: str='_entity.id',
        only_drop: bool = False,
    ) -> tuple[gemmi.cif.Block, gemmi.cif.Block]:

    """take two cif blocks, a source and target, and merge loops giving preference to source
    in case of duplicate keys."""

    source_item = source_block.find_loop_item(reference_tag)
    target_item = target_block.find_loop_item(reference_tag)

    if not only_drop:
        if source_item is None or target_item is None:
            raise ValueError("One or both blocks are missing the required loop with reference_tag: " + reference_tag)
        
        if not set(source_item.loop.tags).issubset(set(target_item.loop.tags)):
            raise ValueError("Source loop contains tags that are not present in target loop, cannot merge.")
        
        source_item.loop.add_columns([t for t in target_item.loop.tags if t not in source_item.loop.tags], "?")

        if source_item.loop.tags != target_item.loop.tags:
            print("source loop tags: ",source_item.loop.tags)
            print("target loop tags: ",target_item.loop.tags)
            raise ValueError("After adding missing columns, source and target loops have different tags, cannot merge.")
        idx_to_move = [
            row_index
            for row_index, entity_id in enumerate(source_block.find_loop(reference_tag))
            if entity_id not in target_block.find_loop(reference_tag)
        ]

        for row_index in idx_to_move:
            row = [source_item.loop[row_index,c] for c in range(source_item.loop.width())]
            target_item.loop.add_row(row)

    # source block without the redundant loop
    fresh_block = gemmi.cif.Block(source_block.name)
    for item in source_block:
        if item.loop and reference_tag in item.loop.tags:
            continue
        else:
            fresh_block.add_item(item)

    return fresh_block, target_block

def prune_empty_loops(block: gemmi.cif.Block) -> gemmi.cif.Block:
    """remove loops from a cif block that have no rows"""
    new_block = gemmi.cif.Block(block.name)
    for item in block:
        if item.pair:
            new_block.add_item(item)                     
        if item.loop and item.loop.length() > 0:
            new_block.add_item(item)
    return new_block

def filter_mmcif_categories(block: gemmi.cif.Block, allowed_categories: list[str]) -> gemmi.cif.Block:
    """remove items from a cif block that do not belong to allowed categories"""
    new_block = gemmi.cif.Block(block.name)
    for item in block:
        if item.pair:
            category = item.pair[0].split(".")[0] + '.'
            if category in allowed_categories:
                new_block.add_item(item)
        elif item.loop:
            tags = set(s.split(".")[0] + '.' for s in item.loop.tags)
            if len(tags) != 1:
                raise ValueError(
                    "Error in CIF loop table formatting: found multiple category keys."
                )
            category = tags.pop()
            if category in allowed_categories:
                new_block.add_item(item)
        else:
            raise ValueError("Found CIF item that is neither pair nor loop")
    return new_block

def prepare_cif_block_for_merging(block: gemmi.cif.Block, allowed_categories: list[str]) -> gemmi.cif.Block:
    """prepare a cif block for merging by pruning empty loops and filtering categories"""
    pruned_block = prune_empty_loops(block)
    filtered_block = filter_mmcif_categories(pruned_block, allowed_categories)
    return filtered_block


def update_entity_id_loops(block: gemmi.cif.Block, exclude_polymer_entity_ids: bool=True):

    entity_id_to_asym = defaultdict(list)
    for entity_id, asym_id in zip(block.find_loop('_struct_asym.entity_id'), block.find_loop('_struct_asym.id')):
        entity_id_to_asym[entity_id].append(asym_id)

    asym_id_to_resname = defaultdict(list)
    for asym_id, resname in zip(block.find_loop('_atom_site.label_asym_id'), block.find_loop('_atom_site.label_comp_id')):
        asym_id_to_resname[asym_id].append(resname)

    entity_id_to_exclude = [entity_id for entity_id, entity_type in zip(block.find_loop('_entity.id'), block.find_loop('_entity.type')) if entity_type == 'polymer'] if exclude_polymer_entity_ids else []

    for entity_id in entity_id_to_exclude:
        entity_id_to_asym.pop(entity_id, None)
    
    entity_id_to_resname = {}
    for entity_id, asym_ids in entity_id_to_asym.items():
        resname = set().union(*[r for r in [asym_id_to_resname[asym_id] for asym_id in asym_ids]])
        if len(resname) > 1:
            # could be a more helpful exception message
            raise ValueError(f"for {block} found multiple labels for asym, should be unique")
        entity_id_to_resname[entity_id] = resname.pop()


    # will need updating for each project
    resname_to_description = {
        'HOH': {'src_method': 'nat', 'pdbx_description': "water"},
        'CA': {'src_method': 'syn', 'pdbx_description': "\"CALCIUM ION\""},
        'CL': {'src_method': 'syn', 'pdbx_description': "\"CHLORIDE ION\""},
        'DMS': {'src_method': 'syn', 'pdbx_description': "\"DIMETHYL SULFOXIDE\""},
        'NA': {'src_method': 'syn', 'pdbx_description': "\"SODIUM ION\""},
        'UNL': {'src_method': 'syn', 'pdbx_description': "LIGAND"},
        'MG': {'src_method': 'syn', 'pdbx_description': "\"MAGNESIUM ION\""},
    }

    #update the entity_id table src_method and pdbx_description for non-polymer
    loop = block.find_loop_item('_entity.id').loop
    src_method_idx = loop.tags.index('_entity.src_method')
    pdbx_description_idx = loop.tags.index('_entity.pdbx_description')
    type_idx = loop.tags.index('_entity.type')
    entity_id_idx = loop.tags.index('_entity.id')

    for i in range(loop.length()): #cannot assume sequential entity ids
        entity_type = loop[i, type_idx]
        if entity_type != 'polymer':
            entity_id = loop[i, entity_id_idx]
            loop[i, src_method_idx] = resname_to_description[entity_id_to_resname[entity_id]]['src_method']
            loop[i, pdbx_description_idx] = resname_to_description[entity_id_to_resname[entity_id]]['pdbx_description']


def convert_cif_pairs_to_loop(block: gemmi.cif.Block, category: str) -> gemmi.cif.Block:
    """convert pairs in a cif block that belong to a specified category into a loop"""

    existing_loop_tags = set().union(*[set(item.loop.tags) for item in block if item.loop])
    if any([category==tag.split('.')[0] for tag in existing_loop_tags]):
        raise ValueError(f"Block already contains loop items from category {category}")

    new_block = gemmi.cif.Block(block.name)
    loop_tags = []
    loop_values = []
    
    for item in block:
        if item.pair and item.pair[0].split('.')[0] == category:
            loop_tags.append(item.pair[0].split('.')[1])
            loop_values.append(item.pair[1])
        else:
            new_block.add_item(item)
    if loop_tags:
        new_block.init_mmcif_loop(category, loop_tags)
        loop = new_block.find_loop(f"{category}.{loop_tags[0]}").get_loop()
        loop.add_row([str(v) for v in loop_values])
    return new_block