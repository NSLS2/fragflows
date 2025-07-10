import gemmi
import pandas as pd
import string
from .cif_blocks import (
    refinement_cif_to_cif_block,
    original_mtz_to_cif_block,
    event_map_to_cif_block,
)
from .structure import remove_ground_state
from .utils import letter_generator

def make_changed_state_sf_cif(
    table1: pd.DataFrame,
    table2: pd.DataFrame,
    xtal_id: str,
) -> gemmi.cif.Document:

    doc = gemmi.cif.Document()
    # everything after primary refinement block will be named xxxxAsf, xxxxBsf, etc.
    letter_gen = letter_generator()

    # read data from disk to cif blocks

    # main refinement block
    try:
        REFINEMENT_SF_CIF = table1.loc[
            table1["xtal_id"] == xtal_id, "refined_reflection_file"
        ].iloc[0]
    except KeyError as e:
        print(f"caught {e}: missing refinement reflection cif for {xtal_id}")
    refinement_block = refinement_cif_to_cif_block(REFINEMENT_SF_CIF)
    rblock = gemmi.as_refln_blocks(gemmi.cif.read_file(REFINEMENT_SF_CIF))[0]

    # original intensity measurements
    REFLECTION_DATA = table1.loc[
        table1["xtal_id"] == xtal_id, "reflection_data_file"
    ].iloc[0]
    original_mtz_block = original_mtz_to_cif_block(
        REFLECTION_DATA,
        rblock.cell,
        rblock.spacegroup,
        block_name=f"xxxx{next(letter_gen)}sf",
    )

    # wavelength set by original reflection file
    wl = original_mtz_block.find_value("_diffrn_radiation_wavelength.wavelength")

    # update blocks and append to doc

    # refinement_block
    refinement_block.set_pair("_diffrn_radiation_wavelength.id", "1")
    refinement_block.set_pair("_diffrn_radiation_wavelength.wavelength", f"{wl}")
    refinement_block.move_item(-3, -1)
    doc.add_copied_block(refinement_block)

    # intensity block
    doc.add_copied_block(original_mtz_block)

    # event map block, resolution determined by measured intensities
    high_res = gemmi.read_mtz_file(REFLECTION_DATA).resolution_high()
    for idx, row in table2[table2["xtal_id"] == xtal_id].iterrows():
        EVENT_MAP = row["event_map_file"]
        # xtal_uid = row['uid']
        # catalog_id = table1.loc[table1['xtal_id'] == xtal_id, 'catalog_id'].iloc[0]
        x, y, z = row["x"], row["y"], row["z"]
        event_background_density_correction = row["1-BDC"]

        diffrn_details = f"ligand evidence PanDDA event map;1-BDC=({event_background_density_correction});event_site=({x},{y},{z})"
        event_map_block = event_map_to_cif_block(
            EVENT_MAP,
            high_res=high_res,
            block_name=f"xxxx{next(letter_gen)}sf",
            wavelength=wl,
            details=diffrn_details,
        )
        doc.add_copied_block(event_map_block)

    return doc


def make_changed_state_cif(
    table1: pd.DataFrame,
    xtal_id: str,
    template_path: str,
    ligand_csv: str,
    block_append_identifier: str = "comp_",
    block_name: str = "xxxx",
):
    """This function will need information from the template cif file.

    Convert ensemble cif to structure for ground state removal, then convert structure to cif
    block for combining with template block.

    The ensemble cif contains multiple blocks corresponding to ligand and/or ion definitions
    after the atomic model. These definitions are included in the final cif document.
    """
    doc = gemmi.cif.Document()
    print(xtal_id)
    structure_cif = table1.loc[
        table1["xtal_id"] == xtal_id, "refined_structure_file"
    ].iloc[0]
    ensemble_structure = gemmi.read_structure(structure_cif)
    changed_state_structure = remove_ground_state(ensemble_structure)
    sblock = changed_state_structure.make_mmcif_block()

    # information from template cif
    template_block = gemmi.cif.read_file(template_path).sole_block()

    # inject structure specific information into template_block
    # this information will ultimately include: collection_date, wavelength, catalog_id

    # catalog_id
    ligand_df = pd.read_csv(ligand_csv)
    xtal_id_soak_entry = ligand_df[ligand_df['xtal_id'] == xtal_id]
    if not xtal_id_soak_entry.empty:
        catalog_id = xtal_id_soak_entry['catalog_id'].iloc[0]
    else:
        catalog_id = None
    
    try:
        s = template_block.find_pair('_struct.title')[-1]
        template_block.set_pair('_struct.title',s.replace('XXXX',catalog_id))
    except TypeError:
        print(f'Error updating template cif block ligand catalog_id for {xtal_id}')

    doc.add_copied_block(template_block, pos=-1)

    if "ensemble" not in sblock.name:
        raise ValueError("input does not appear to be pandda ensemble model")

    for x in sblock.get_mmcif_category_names():
        pair_key = f"{x}entry_id"
        pair = sblock.find_value(pair_key)
        if pair is not None:
            sblock.set_pair(pair_key, block_name)

    # update block name to reflect ground state removal, verify ensemble input
    # assign default entry_id to be updated by PDB
    sblock.name = block_name
    sblock.set_pair("_entry.id", block_name)

    # ground state removed, changed state model
    doc.add_copied_block(sblock, pos=-1)

    # Append additional cif blocks
    structure_cif_blocks = gemmi.cif.read_file(structure_cif)
    for block in structure_cif_blocks:
        if block_append_identifier in block.name:
            doc.add_copied_block(block, pos=-1)

    return doc

def assemble_group_changed_state_cifs(
        refinement_table: pd.DataFrame,
        event_table: pd.DataFrame,
        group_dep_dir: str,
        **kwargs
):
    xtal_ids = list(refinement_table['xtal_id'])
    if len(set(xtal_ids)) != len(xtal_ids):
        raise Exception('redundant xtal_ids found in refinement_table')
    
    for xtal_id in xtal_ids:
        print(f'generating files for {xtal_id}')
        changed_state_sf_doc = make_changed_state_sf_cif(refinement_table, event_table, xtal_id)
        changed_state_doc = make_changed_state_cif(refinement_table, xtal_id, **kwargs)
        changed_state_sf_doc.write_file(f'{group_dep_dir}/{xtal_id}-sf.cif')
        changed_state_doc.write_file(f'{group_dep_dir}/{xtal_id}.cif')

    return