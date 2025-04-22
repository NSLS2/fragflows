import gemmi
import pandas as pd
import string
from .cif_blocks import (
    refinement_cif_to_cif_block,
    original_mtz_to_cif_block,
    event_map_to_cif_block,
)
from .structure import remove_ground_state


def letter_generator():
    def get_letter(num):
        result = []
        while num > 0:
            num, remainder = divmod(num - 1, 26)
            result.append(string.ascii_uppercase[remainder])
            return "".join(reversed(result))

    for k in range(1, 15000):
        yield get_letter(k)


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
    # rblock = gemmi.as_refln_blocks(gemmi.cif.read_file(REFINEMENT_SF_CIF))[0]

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
    doc.add_copied_block(template_block, pos=-1)

    # update block name to reflect ground state removal, verify ensemble input
    # assign default entry_id to be updated by PDB
    sblock.name = block_name
    sblock.set_pair("_entry.id", block_name)

    for x in sblock.get_mmcif_category_names():
        pair_key = f"{x}entry_id"
        pair = sblock.find_value(pair_key)
        if pair is not None:
            if "ensemble" not in pair:
                raise ValueError("input does not appear to be pandda ensemble model")
            else:
                sblock.set_pair(pair_key, block_name)

    # ground state removed, changed state model
    doc.add_copied_block(sblock, pos=-1)

    # Append additional cif blocks
    structure_cif_blocks = gemmi.cif.read_file(structure_cif)
    for block in structure_cif_blocks:
        if block_append_identifier in block.name:
            doc.add_copied_block(block, pos=-1)

    return doc
