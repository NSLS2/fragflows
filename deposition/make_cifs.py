import gemmi
import pandas as pd
import numpy as np
from pathlib import Path
import json
from .cif_blocks import (
    refinement_cif_to_cif_block,
    original_mtz_to_cif_block,
    event_map_to_cif_block,
    deduplicate_cif_loops,
    prepare_cif_block_for_merging,
)
from .structure import remove_ground_state
from .utils import letter_generator
from urllib.parse import urlencode
from .load import ispyb_xml_to_cif_block
from .beamline_parameters import BEAMLINE_PARAMETERS
from .utils import convert_iso_date_to_ymd
from .structure import tweak_occupancy, resolve_entities

def make_changed_state_sf_cif(
    table1: pd.DataFrame,
    table2: pd.DataFrame,
    xtal_id: str,
    xtal_id_key: str = "xtal_id",
) -> gemmi.cif.Document:

    diffrn_id = xtal_id

    doc = gemmi.cif.Document()
    # everything after primary refinement block will be named xxxxAsf, xxxxBsf, etc.
    letter_gen = letter_generator()

    # read data from disk to cif blocks

    row = table1.loc[table1[xtal_id_key] == xtal_id].iloc[0]
    diffrn_crystal_treatment = {
            "method": "soak",
            "catalog_id": row.get("catalog_id", "unknown"),
            "smiles": row.get("smiles", "unknown"),
            "solvent": "DMSO",
        }
    
    diffrn_crystal_treatment_json = f"'{json.dumps(diffrn_crystal_treatment)}'"

    # main refinement block
    try:
        REFINEMENT_SF_CIF = table1.loc[
            table1[xtal_id_key] == xtal_id, "refined_reflection_file"
        ].iloc[0]
    except KeyError as e:
        print(f"caught {e}: missing refinement reflection cif for {xtal_id}")
    refinement_block = refinement_cif_to_cif_block(
                            REFINEMENT_SF_CIF,
                            xtal_id=xtal_id,
                            crystal_treatment=diffrn_crystal_treatment_json,
                        )
    
    rblock = gemmi.as_refln_blocks(gemmi.cif.read_file(REFINEMENT_SF_CIF))[0]

    # original intensity measurements
    REFLECTION_DATA = table1.loc[
        table1[xtal_id_key] == xtal_id, "reflection_data_file"
    ].iloc[0]
    original_mtz_block = original_mtz_to_cif_block(
        REFLECTION_DATA,
        rblock.cell,
        rblock.spacegroup,
        block_name=f"xxxx{next(letter_gen)}sf",
        xtal_id=xtal_id,
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
    for idx, row in table2[table2[xtal_id_key] == xtal_id].iterrows():
        EVENT_MAP = row["event_map_file"]
        # xtal_uid = row['uid']
        # catalog_id = table1.loc[table1[xtal_id_key] == xtal_id, 'catalog_id'].iloc[0]
        x, y, z = row["x"], row["y"], row["z"]
        event_background_density_correction = row["1-BDC"]

        diffrn_details = {
            "ligand_evidence": "PanDDA_event_map",
            "1-BDC": event_background_density_correction,
            "event_site_x": np.round(x,3),
            "event_site_y": np.round(y,3),
            "event_site_z": np.round(z,3),
        }
        diffrn_details_json = f"'{json.dumps(diffrn_details)}'"

        event_map_block = event_map_to_cif_block(
            EVENT_MAP,
            high_res=high_res,
            block_name=f"xxxx{next(letter_gen)}sf",
            wavelength=wl,
            details=diffrn_details_json,
            xtal_id=xtal_id,
        )
        doc.add_copied_block(event_map_block)

    

    return doc


def make_changed_state_cif(
    table1: pd.DataFrame,
    xtal_id: str,
    template_path: str,
    ligand_df: pd.DataFrame,
    block_append_identifier: str = "comp_",
    block_name: str = "xxxx",
    xtal_id_key: str = "xtal_id",
):
    """This function will need information from the template cif file.

    Convert ensemble cif to structure for ground state removal, then convert structure to cif
    block for combining with template block.

    The ensemble cif contains multiple blocks corresponding to ligand and/or ion definitions
    after the atomic model. These definitions are included in the final cif document.
    """
    diffrn_id = xtal_id

    print(xtal_id)
    row = table1.loc[table1[xtal_id_key] == xtal_id].iloc[0]
    structure_cif = row["refined_structure_file"]
    ensemble_structure = gemmi.read_structure(structure_cif)
    n_atoms = len([a for m in ensemble_structure for c in m for r in c for a in r])
    tweak_occupancy(ensemble_structure)
    resolve_entities(ensemble_structure)

    if len([a for m in ensemble_structure for c in m for r in c for a in r]) != n_atoms:
        raise ValueError("number of atoms changed, something went wrong")
    
    sblock = ensemble_structure.make_mmcif_block()

    sblock = prepare_cif_block_for_merging(
        sblock,
        allowed_categories=[
            '_entry.',
            '_cell.',
            '_symmetry.',
            '_entity.',
            '_chem_comp.',
            '_struct_asym.',
            '_entity_poly_seq.',
            '_atom_site.'
        ]
    )


    sblock_metadata_blocks = gemmi.cif.read_file(structure_cif)
    sblock_metadata_block = sblock_metadata_blocks[0]
    sblock_metadata_block = prepare_cif_block_for_merging(
        sblock_metadata_block,
        allowed_categories=[
            '_software.',
            '_refine.',
            '_refine_ls_shell.',
            '_refine_ls_restr.',
            '_ccp4_form_factor.',
            '_atom_type.',
            '_ccp4_reflns_twin.',
            '_ccp4_refine_ls.',
            '_ccp4_refine.',
        ]
    )

    # information from template cif
    template_block = gemmi.cif.read_file(template_path).sole_block()

    # update .crystal_id pairs from template block
    crystal_id = row.get(xtal_id_key, None)
    if pd.isna(crystal_id) or str(crystal_id).strip() == "":
        raise ValueError(
            f"Missing or empty '{xtal_id_key}' for xtal_id '{xtal_id}' in refinement table row"
        )
    template_block.set_pair("_exptl_crystal_grow.crystal_id", str(crystal_id))
    template_block.set_pair("_diffrn.crystal_id", str(crystal_id))

    # we need to resolve entity and entity_poly duplications
    # the template_block is modified in place, sblock returns a modified copy because we needed
    # to drop loops.

    sblock, template_block = deduplicate_cif_loops(sblock, template_block, reference_tag='_entity.id')
    sblock, template_block = deduplicate_cif_loops(
        sblock, template_block, reference_tag='_entity_poly.entity_id', only_drop=True
    )

    # inject structure specific information into template_block
    # this information will ultimately include:
    # reflection stats, collection_date, wavelength, catalog_id

    # collection_date and other beamline parameters
    if BEAMLINE_PARAMETERS.get(row["det_serial_no"]) is None:
        raise ValueError(f"missing beamline parameters for detector {row['det_serial_no']}")
    
    update_cif_block_from_dataframe_collection_info(
        template_block, row, BEAMLINE_PARAMETERS.get(row["det_serial_no"])
    )

    # load processing stats from ispyb xml
    rstats_xml = row["xml_path"]
    if not Path(rstats_xml).exists():
        raise ValueError(f"missing ispyb xml file at {rstats_xml} for {xtal_id}")
    
    rstats_block = ispyb_xml_to_cif_block(rstats_xml)

    for idx, item in enumerate(rstats_block):
        template_block.add_item(item)

    # catalog_id
    xtal_id_soak_entry = ligand_df[ligand_df[xtal_id_key] == xtal_id]
    if not xtal_id_soak_entry.empty:
        catalog_id = xtal_id_soak_entry['catalog_id'].iloc[0]
    else:
        catalog_id = None
    
    try:
        s = template_block.find_pair('_struct.title')[-1]
        if not catalog_id:
            print(f'Warning: missing catalog_id for {xtal_id}, cannot update template cif block title')
            catalog_id = 'DMSO'

        template_block.set_pair('_struct.title',s.replace('XXXX',catalog_id))
    except TypeError:
        print(f'Error updating template cif block ligand catalog_id for {xtal_id}')

    #if "ensemble" not in sblock.name:
    #    raise ValueError("input does not appear to be pandda ensemble model")

    for x in sblock.get_mmcif_category_names():
        pair_key = f"{x}entry_id"
        pair = sblock.find_value(pair_key)
        if pair is not None:
            sblock.set_pair(pair_key, block_name)

    # update block name
    sblock.name = block_name
    sblock.set_pair("_entry.id", block_name)



    # rather than append template block to structure block, loop through
    # template and add elements to structure block to avoid separate blocks
    # for idx, item in enumerate(template_block):
    #    sblock.add_item(item)
    #    sblock.move_item(-1,idx)


    #for item in sblock:
    #    template_block.add_item(item)

    doc = gemmi.cif.Document()
    #cif_text = template_block.as_string()
    #fixed_template_block = gemmi.cif.read_string(cif_text).sole_block()
    
    for item in sblock_metadata_block:
        template_block.add_item(item)
    for item in sblock:
        template_block.add_item(item)

    for item in template_block:
        if item.pair:
            if 'pdbx_diffrn_id' == item.pair[0].split('.')[-1] or 'diffrn_id' == item.pair[0].split('.')[-1]:
                template_block.set_pair(item.pair[0], diffrn_id)
            elif '_diffrn.id' == item.pair[0]:
                template_block.set_pair(item.pair[0], diffrn_id)
            elif '_exptl_crystal_grow.crystal_id' == item.pair[0]:
                template_block.set_pair(item.pair[0], crystal_id)
            elif '_diffrn.crystal_id' == item.pair[0]:
                template_block.set_pair(item.pair[0], crystal_id)

    template_block.name = block_name
    doc.add_copied_block(template_block)

    for block in sblock_metadata_blocks:
        if block_append_identifier in block.name:
            doc.add_copied_block(block)

    return doc

def assemble_group_changed_state_cifs(
        refinement_table: pd.DataFrame,
        event_table: pd.DataFrame,
        group_dep_dir: str,
        template_cif_path: str,
        ligand_table: pd.DataFrame,
        xtal_id_key: str = "xtal_id",
        block_append_identifier: str = "comp_",
        block_name: str = "xxxx",
):
    xtal_ids = list(refinement_table[xtal_id_key])
    if len(set(xtal_ids)) != len(xtal_ids):
        raise Exception('redundant xtal_ids found in refinement_table')
    
    for xtal_id in xtal_ids:
        print(f'generating files for {xtal_id}')
        changed_state_sf_doc = make_changed_state_sf_cif(
            refinement_table,
            event_table,
            xtal_id,
            xtal_id_key=xtal_id_key,
        )
        changed_state_doc = make_changed_state_cif(
            refinement_table,
            xtal_id,
            template_cif_path,
            ligand_table,
            block_append_identifier=block_append_identifier,
            block_name=block_name,
            xtal_id_key=xtal_id_key,
        )
        changed_state_sf_doc.write_file(f'{group_dep_dir}/{xtal_id}-sf.cif')
        changed_state_doc.write_file(f'{group_dep_dir}/{xtal_id}.cif')

    return

def update_cif_block_from_dataframe_collection_info(
    cif_block: gemmi.cif.Block,
    row: pd.Series,
    beamline_params: dict[str, str],
    required_fields: list[str] = ['data_collection_date', 'wavelength', 'det_serial_no'],
) -> dict[str, str]:
    """
    Update a CIF block with values from a dataframe row.
    
    Args:
        cif_block: The CIF block to update
        row: A dataframe row (pd.Series)
        beamline_params: Dict mapping {df_column: cif_key}
        required_fields: List of required dataframe columns that must exist

    Returns:
        dict with status info including missing_fields and errors
    
    Raises:
        ValueError: If required fields are missing
    """
    required_fields = required_fields or []
    
    missing_in_df = [field for field in required_fields if pd.isna(row.get(field))]
    if missing_in_df:
        raise ValueError(f"Missing required fields in dataframe: {missing_in_df}")
    
    updated_fields = []
    
    for mmcif_key, mmcif_value in beamline_params['cif_update'].items():
        # Check if value exists in dataframe
        try:
            cif_block.set_pair(mmcif_key, str(mmcif_value))
            updated_fields.append(mmcif_key)
        except Exception as e:
            raise ValueError(f"Error updating CIF block for key {mmcif_key} with value {mmcif_value}: {e}")
        
        # update wavelength in two places
    cif_block.set_pair("_diffrn_source.pdbx_wavelength_list", str(row.get("wavelength")))
    updated_fields.append("_diffrn_source.pdbx_wavelength_list")
    cif_block.set_pair("_diffrn_radiation_wavelength.wavelength", str(row.get("wavelength")))
    updated_fields.append("_diffrn_radiation_wavelength.wavelength")
    dt = convert_iso_date_to_ymd(row.get("data_collection_date"), beamline_params.get("timezone"))
    cif_block.set_pair("_diffrn_detector.pdbx_collection_date", dt)
    updated_fields.append("_diffrn_detector.pdbx_collection_date")
    
    return {
        "updated_fields": updated_fields,
    }