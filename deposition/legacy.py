import os
from pathlib import Path
import pandas as pd
from deposition.pipeline_programs import xml_path_to_pipeline_programs
from fragflows_db.crud import mtz_from_xml
from fragflows_db.utils import sha256sum
from deposition.beamline_parameters import BEAMLINE_PARAMETERS
from deposition.load import ispyb_xml_to_cif_block
from deposition.cif_blocks import convert_cif_pairs_to_loop, dimple_mtz_to_cif_block
import json
from deposition.utils import letter_generator, convert_iso_date_to_ymd
import gemmi
from deposition.structure import remap_entity_ids
from deposition.cif_blocks import (
        prepare_cif_block_for_merging,
        insert_pair_into_cif_block,
        deduplicate_cif_loops
)

# some helpful functions for combining dataframes from different
# fragflows versions

# 20260421
# helper function for extending the original csv/dataframe used for 
# dimpleflow with ispyb xml and hdf5 info needed for group deposition.

def merge_dfs_for_group_dep(
    df1: pd.DataFrame, 
    df2: pd.DataFrame, 
    columns_to_add: list = None
) -> pd.DataFrame:
    """
    Merge dataframes by matching MTZ file checksums.
    
    Computes SHA256 checksums for MTZ files in df1 and df2, then maps columns
    from df2 to df1 based on matching checksums.
    
    Args:
        df1: DataFrame with 'filepath' column
        df2: DataFrame with 'xml_path' column (contains MTZ reference)
        columns_to_add: Columns to map from df2 to df1. Defaults to ispyb metadata.
    
    Returns:
        df1 with new columns added via checksum mapping
    """
    if columns_to_add is None:
        columns_to_add = [
            "xml_path", 
            "data_collection_date", 
            "wavelength", 
            "det_description", 
            "det_serial_no"
        ]
    
    # Compute checksums for df1 files
    for idx, row in df1.iterrows():
        mtz_filepath = row["filepath"]
        if not Path(mtz_filepath).exists():
            raise FileNotFoundError(f"File {mtz_filepath} does not exist.")
        mtz_checksum = sha256sum(mtz_filepath)
        print(f"Computed checksum for {mtz_filepath}: {mtz_checksum}")
        df1.at[idx, "mtz_checksum"] = mtz_checksum

    # Compute checksums for df2 files (from XML)
    for idx, row in df2.iterrows():
        df2.at[idx, "mtz_checksum"] = sha256sum(mtz_from_xml(row["xml_path"]))

    # Map columns from df2 to df1 using checksums
    for col in columns_to_add:
        mapping = df2.dropna(subset=["mtz_checksum", col]).set_index("mtz_checksum")[col]
        df1[col] = df1["mtz_checksum"].map(mapping)
    
    return df1.copy()

# for now put this in legacy, eventually we will extend db to include batch jobs and pipeline info

def build_ground_structure_cif_metadata_loops(
        group_dep_df: pd.DataFrame,
        unmodelled_xtal_ids: list[str],
        block_name: str="XXXX",
        **kwargs) -> gemmi.cif.Block:
    """
    Build CIF block with experimental metadata for ground structure deposition.
    
    Creates loops for crystal, diffraction source, detector, radiation, and
    reflection statistics for all crystals in the group deposition.
    
    Args:
        group_dep_df: DataFrame with columns: xtal_id, wavelength, data_collection_date, 
                      det_serial_no, xml_path
        unmodelled_xtal_ids: List of crystal IDs to include
        block_name: Name for output CIF block
    Kwargs:
        optional dict with default tags for any loops
    
    Returns:
        gemmi.cif.Block with populated metadata loops
    
    Raises:
        ValueError: If xtal_id duplicates or missing values found
        FileNotFoundError: If referenced files don't exist
    """
    
    # Validation
    if len(group_dep_df['xtal_id'].unique()) < len(group_dep_df):
        raise ValueError("Duplicate xtal_id values found in input dataframe.")
    
    for xtal_id in unmodelled_xtal_ids:
        if xtal_id not in group_dep_df['xtal_id'].values:
            raise ValueError(f"xtal_id {xtal_id} not found in input dataframe.")

    # Initialize block and loops
    new_block = gemmi.cif.Block(block_name)

    new_block.init_mmcif_loop("_exptl_crystal", ["id","density_matthews","density_percent_sol"])
    new_block.init_mmcif_loop(
        "_exptl_crystal_grow", ["crystal_id","method", "pH", "temp", "pdbx_details", "pdbx_pH_range"])
    new_block.init_mmcif_loop("_diffrn", ["id", "crystal_id", "ambient_temp"])
    new_block.init_mmcif_loop(
        "_diffrn_source",
        ["diffrn_id", "source", "type", "pdbx_wavelength_list", "pdbx_synchrotron_beamline", "pdbx_synchrotron_site"]
    )
    new_block.init_mmcif_loop(
        "_diffrn_detector",
        ["diffrn_id", "details", "detector", "type", "pdbx_collection_date"]
    )
    new_block.init_mmcif_loop(
        "_diffrn_radiation",
        ["diffrn_id", "pdbx_scattering_type", "pdbx_diffrn_protocol", "wavelength_id", "pdbx_monochromatic_or_laue_m_l", "monochromator"]
    )
    new_block.init_mmcif_loop(
        "_diffrn_radiation_wavelength",
        ["id", "wavelength"]
    )

    # Populate loops for each crystal
    pdbx_ordinal = 1
    for i, xtal_id in enumerate(unmodelled_xtal_ids):
        row = group_dep_df[group_dep_df["xtal_id"] == xtal_id].iloc[0]
        df_wavelength = row["wavelength"]
        df_collection_date = row["data_collection_date"]
        serial_no = row["det_serial_no"]
        tz = BEAMLINE_PARAMETERS.get(serial_no, {}).get("timezone", "UTC")
        
        # Build beamline metadata dictionary
        df_collection_date_dict = (
            BEAMLINE_PARAMETERS.get(serial_no, {}).get("cif_update") 
            | {
                "_diffrn_source.pdbx_wavelength_list": str(df_wavelength),
                "_diffrn_source.diffrn_id": xtal_id,
                "_diffrn_detector.diffrn_id": xtal_id,
                "_diffrn_radiation.diffrn_id": xtal_id,
                "_diffrn_radiation.wavelength_id": xtal_id,
                "_diffrn_detector.pdbx_collection_date": convert_iso_date_to_ymd(df_collection_date, tz),
            }
            | kwargs
            | {"_exptl_crystal_grow.crystal_id": xtal_id}
        )

        # Add crystal metadata rows
        exptl_crystal_loop = new_block.find_loop_item("_exptl_crystal.id").loop
        exptl_crystal_loop.add_row([xtal_id, "?", "?"])
        
        exptl_crystal_grow_loop = new_block.find_loop_item("_exptl_crystal_grow.crystal_id").loop
        exptl_crystal_grow_loop.add_row([df_collection_date_dict.get(tag, "?") for tag in exptl_crystal_grow_loop.tags])
        
        # Add diffraction metadata rows
        diffrn_loop = new_block.find_loop_item("_diffrn.id").loop
        diffrn_loop.add_row([xtal_id, xtal_id, "100"])
        
        diffrn_source_loop = new_block.find_loop_item("_diffrn_source.diffrn_id").loop
        diffrn_source_loop.add_row([df_collection_date_dict.get(tag, "?") for tag in diffrn_source_loop.tags])
        
        diffrn_detector_loop = new_block.find_loop_item("_diffrn_detector.diffrn_id").loop
        diffrn_detector_loop.add_row([df_collection_date_dict.get(tag, "?") for tag in diffrn_detector_loop.tags])
        
        diffrn_radiation_loop = new_block.find_loop_item("_diffrn_radiation.diffrn_id").loop
        diffrn_radiation_loop.add_row([df_collection_date_dict.get(tag, "?") for tag in diffrn_radiation_loop.tags])

        if i == 0:
            diffrn_radiation_wavelength_loop = new_block.find_loop_item("_diffrn_radiation_wavelength.id").loop
            diffrn_radiation_wavelength_loop.add_row([xtal_id, str(df_wavelength)])

        # Add reflection statistics rows
        xml_path = row["xml_path"]
        if not Path(xml_path).exists():
            raise FileNotFoundError(f"File {xml_path} does not exist.")
        
        rstats_block = ispyb_xml_to_cif_block(xml_path)
        rstats_dict = {i.pair[0]: i.pair[1] for i in rstats_block}
        
        # Initialize rstats loops on first crystal
        if new_block.find_loop_item("_reflns.pdbx_diffrn_id") is None:
            print("Initializing refinement statistics loops...")
            reflns_tags = [k.split('.')[1] for k in rstats_dict.keys() if k.split('.')[0] == '_reflns']
            shell_tags = [k.split('.')[1] for k in rstats_dict.keys() if k.split('.')[0] == '_reflns_shell']
            new_block.init_mmcif_loop("_reflns", reflns_tags)
            new_block.init_mmcif_loop("_reflns_shell", shell_tags)
        
        # Add diffrn_id and ordinal to statistics
        rstats_dict.update({
            "_reflns.pdbx_diffrn_id": xtal_id,
            "_reflns_shell.pdbx_diffrn_id": xtal_id,
            "_reflns.pdbx_ordinal": str(pdbx_ordinal),
            "_reflns_shell.pdbx_ordinal": str(pdbx_ordinal),
        })
        
        reflns_loop = new_block.find_loop_item("_reflns.pdbx_diffrn_id").loop
        reflns_loop.add_row([rstats_dict.get(tag, "?") for tag in reflns_loop.tags])
        
        reflns_shell_loop = new_block.find_loop_item("_reflns_shell.pdbx_diffrn_id").loop
        reflns_shell_loop.add_row([rstats_dict.get(tag, "?") for tag in reflns_shell_loop.tags])
        
        pdbx_ordinal += 1

    return new_block

def get_mtz_input_from_dimple_log(dimple_log_path: str) -> str:
    """
    Extract MTZ input file path from dimple.log.
    
    Args:
        dimple_log_path: Path to dimple.log file
    
    Returns:
        Path to the MTZ input file
    
    Raises:
        FileNotFoundError: If MTZ file doesn't exist
        ValueError: If 'data_file:' entry not found in log
    """
    with open(dimple_log_path, "r") as f:
        for line in f:
            if (parts := line.split()) and 'data_file:' == parts[0]:
                mtz_path = parts[1]
                if Path(mtz_path).exists():
                    return mtz_path
                else:
                    raise FileNotFoundError(f"MTZ file {mtz_path} does not exist.")
    raise ValueError(f"No line starting with 'data_file:' found in {dimple_log_path}.")


def _collect_dimple_metadata(
    models_dir: str,
    group_dep_df: pd.DataFrame,
    unmodelled_xtal_ids: list[str],
    log_file_name: str = "dimple.log"
) -> dict:
    """
    Scan models directory and collect dimple output metadata.
    
    Args:
        models_dir: Root directory containing dimple outputs
        group_dep_df: DataFrame with checksums for MTZ validation
        unmodelled_xtal_ids: List of crystals to include
        log_file_name: Name of dimple log file
    
    Returns:
        Dictionary mapping xtal_id to dimple metadata
    
    Raises:
        ValueError: If xtal_id not found or checksum mismatch
        FileNotFoundError: If required files don't exist
    """
    # Collect xtal_ids from dimple log directories
    xtal_id_from_dimple_dirs = []
    for r, d, f in os.walk(models_dir):
        xtal_id_from_dimple_dirs.extend(
            [os.path.basename(r) for f_ in f if f_ == log_file_name]
        )
    
    # Validate all requested crystals have dimple outputs
    if any(xtal_id not in xtal_id_from_dimple_dirs for xtal_id in unmodelled_xtal_ids):
        raise ValueError(
            "Some xtal_id values from unmodelled list not in dimple log directories."
        )

    # Build metadata dictionary
    dimple_dict = {}
    for r, d, f in os.walk(models_dir):
        for f_ in f:
            if f_ == log_file_name:
                dimple_log_path = os.path.join(r, f_)
                xtal_id = os.path.basename(r)
                
                # Get input MTZ path from log
                mtz_path = get_mtz_input_from_dimple_log(dimple_log_path)
                
                if not Path(mtz_path).exists():
                    raise FileNotFoundError(f"MTZ file {mtz_path} does not exist.")

                # Validate MTZ via checksum
                dimple_log_mtz_checksum = sha256sum(mtz_path)
                if dimple_log_mtz_checksum not in group_dep_df["mtz_checksum"].values:
                    raise ValueError(
                        f"MTZ file {mtz_path} checksum not found in group_dep_df."
                    )
                
                # Verify xtal_id matches
                expected_xtal_id = group_dep_df[
                    group_dep_df["mtz_checksum"] == dimple_log_mtz_checksum
                ]["xtal_id"].values[0]
                if expected_xtal_id != xtal_id:
                    raise ValueError(
                        f"MTZ checksum mismatch: expected {expected_xtal_id}, got {xtal_id}"
                    )

                # Get dimple MTZ output
                dimple_mtz = os.path.join(r, f"{xtal_id}.dimple.mtz")
                if not Path(dimple_mtz).exists():
                    raise FileNotFoundError(
                        f"Dimple MTZ file {dimple_mtz} does not exist."
                    )
                
                # Get dimple PDB output
                dimple_pdb = os.path.join(r, f"{xtal_id}.dimple.pdb")
                if not Path(dimple_pdb).exists():
                    raise FileNotFoundError(
                        f"Dimple PDB file {dimple_pdb} does not exist."
                    )

                dimple_dict[xtal_id] = {
                    "xtal_id": xtal_id,
                    "dimple_log_path": dimple_log_path,
                    "dimple_mtz": dimple_mtz,
                    "dimple_pdb": dimple_pdb,
                }

    return dimple_dict


def build_ground_structure_factor_cif(
        group_dep_df: pd.DataFrame,
        ligand_df: pd.DataFrame,
        models_dir: str,
        unmodelled_xtal_ids: list[str],
        log_file_name: str = "dimple.log",
) -> gemmi.cif.Document:
    """
    Build CIF document with structure factor blocks from dimple outputs.
    
    Creates one structure factor block per crystal using dimple-refined MTZ files.
    Links crystal metadata from ligand_df.
    
    Args:
        group_dep_df: DataFrame with xtal_id, xml_path, and checksums
        ligand_df: DataFrame with xtal_id, catalog_id, smiles, and crystal treatment info
        models_dir: Root directory containing dimple model folders
        unmodelled_xtal_ids: List of crystal IDs to include
        log_file_name: Name of dimple log file
    
    Returns:
        gemmi.cif.Document with structure factor blocks
    
    Raises:
        FileNotFoundError: If required files don't exist
        ValueError: If data consistency issues found
    """
    unmodelled_xtal_ids = sorted(unmodelled_xtal_ids, key=lambda x: int(x.split('-')[-1]))

    # Compute checksums for all input MTZ files
    group_dep_df_copy = group_dep_df.copy()
    for idx, row in group_dep_df_copy.iterrows():
        xml_filepath = row["xml_path"]
        if not Path(xml_filepath).exists():
            raise FileNotFoundError(f"File {xml_filepath} does not exist.")
        
        mtz_from_xml_filepath = mtz_from_xml(xml_filepath)
        mtz_checksum = sha256sum(mtz_from_xml_filepath)
        print(f"Computed checksum for {mtz_from_xml_filepath}: {mtz_checksum}")
        group_dep_df_copy.at[idx, "mtz_checksum"] = mtz_checksum

    # Collect dimple metadata
    dimple_dict = _collect_dimple_metadata(
        models_dir, group_dep_df_copy, unmodelled_xtal_ids, log_file_name
    )

    # Build document with blocks sorted by xtal_id
    doc = gemmi.cif.Document()
    block_letter = ""
    lg = letter_generator()
    
    for xtal_id in unmodelled_xtal_ids:
        dimple_entry = dimple_dict[xtal_id]
        
        # Get ligand/crystal treatment info
        ligand_row = ligand_df.loc[ligand_df["xtal_id"] == xtal_id].iloc[0]
        diffrn_crystal_treatment = {
            "method": "soak",
            "catalog_id": ligand_row.get("catalog_id", "unknown"),
            "smiles": ligand_row.get("smiles", "unknown"),
            "solvent": "DMSO",
        }
        diffrn_crystal_treatment_json = f"'{json.dumps(diffrn_crystal_treatment)}'"

        # Convert dimple MTZ to CIF block
        block = dimple_mtz_to_cif_block(
            mtz_path=dimple_entry["dimple_mtz"],
            spacegroup=gemmi.SpaceGroup("C2"),
            crystal_treatment=diffrn_crystal_treatment_json,
            block_name=f"XXXX{block_letter}sf",
            xtal_id=xtal_id,
        )
        block_letter = next(lg)
        doc.add_copied_block(block)

        print(
            f"Added block for {xtal_id} with treatment: "
            f"{diffrn_crystal_treatment_json}"
        )
                
    return doc

def convert_pdb_to_cif_block(pdb_path: str) -> gemmi.cif.Block:
    import subprocess
    """
    Convert a PDB file to a CIF block using gemmi.
    
    Args:
        pdb_path: Path to the input PDB file
    
    Returns:
        gemmi.cif.Block containing the structure data
    
    Raises:
        FileNotFoundError: If the PDB file does not exist
        ValueError: If the PDB file cannot be parsed
    """
    if not Path(pdb_path).exists():
        raise FileNotFoundError(f"PDB file {pdb_path} does not exist.")
    
    try:
        cif_path = pdb_path.replace(".pdb", ".cif")
        subprocess.run(
            ["/nsls2/software/mx/ccp4-9/bin/gemmi","convert", pdb_path, cif_path],
            check=True, cwd=os.path.dirname(pdb_path)
        )
        return gemmi.cif.read_file(cif_path).sole_block()
    except Exception as e:
        raise ValueError(f"Error parsing PDB file {pdb_path}: {e}")

def build_ground_structure_cif(
        group_dep_df: pd.DataFrame,
        unmodelled_xtal_ids: list[str],
        template_cif_path: str,
        dimple_metadata: dict,
        target_name: str | None = None,
        ) -> gemmi.cif.Document:
    
    template_block = gemmi.cif.read_file(template_cif_path).sole_block()
    # get exptl_crystal_grow details from template
    exptl_crystal_grow = {
        "_exptl_crystal_grow.method": "\"VAPOR DIFFUSION, SITTING DROP\"",
        "_exptl_crystal_grow.pH": "6.7",
        "_exptl_crystal_grow.temp": "294",
        "_exptl_crystal_grow.pdbx_details": "\"2.2 M DL-malic acid\""
    }
    unmodelled_xtal_ids = sorted(unmodelled_xtal_ids, key=lambda x: int(x.split('-')[-1]))
    first_id = unmodelled_xtal_ids[0]
    representative_sblock = convert_pdb_to_cif_block(dimple_metadata[first_id]["dimple_pdb"])
    remap_entity_ids(representative_sblock)
    representative_mblock = build_ground_structure_cif_metadata_loops(
        group_dep_df, unmodelled_xtal_ids, **exptl_crystal_grow
        )

    # copy strategy used for changed state structure cif
    template_block = prepare_cif_block_for_merging(
        template_block,
        allowed_categories=[
            '_audit_conform.',
            '_pdbx_database_status.',
            '_pdbx_audit_support',
            '_audit_author.',
            '_pdbx_contact_author.',
            '_pdbx_deposit_group.',
            '_citation.',
            '_citation_author.',
            '_struct.',
            '_struct_keywords.',
            '_pdbx_struct_assembly_depositor_info.',
            '_pdbx_struct_assembly_gen_depositor_info.',
            '_pdbx_struct_assembly_auth_evidence_depositor_info.',
            '_pdbx_entry_details.',
            '_exptl.',
            '_entity.',
            '_entity_poly.',
            '_entity_src_gen.',
            '_pdbx_initial_refinement_model.',

        ]
    )
    template_block.set_pair("_pdbx_entry_details.has_ligand_of_interest", "N")
    template_block.set_pair("_exptl.crystals_number", str(len(unmodelled_xtal_ids)))

    template_block.set_pair("_struct.title", f"\"Crystal structure of {target_name} with no ligand modelled\"")

    representative_sblock = prepare_cif_block_for_merging(
        representative_sblock,
        allowed_categories=[
            '_entry.',
            '_cell.',
            '_symmetry.',
            '_entity.',
            '_chem_comp.',
            '_struct_asym.',
            '_entity_poly_seq.',
            '_atom_site.',
            '_software.',
            '_refine.',
            '_refine_ls_shell.',
            '_refine_ls_restr.',
            '_atom_type.',
        ],
        disallowed_tags = ["_refine_ls_shell.number_reflns_work"]
    )

    representative_sblock, template_block = deduplicate_cif_loops(
        representative_sblock, template_block, reference_tag='_entity.id'
    )
    representative_sblock, template_block = deduplicate_cif_loops(
        representative_sblock, template_block , reference_tag='_entity_poly.entity_id', only_drop=True
    )

    # update and fix to conform to PDBx/mmCIF dictionary
    # artifact of ccp4-9 provided gemmi convert function pdb -> cif

    st = gemmi.read_pdb(dimple_metadata[first_id]["dimple_pdb"])
    d_res_h = [s.split()[-1] for s in st.make_pdb_headers().split('\n') if 'BIN RESOLUTION RANGE HIGH' in s][0]
    d_res_l = [s.split()[-1] for s in st.make_pdb_headers().split('\n') if 'BIN RESOLUTION RANGE LOW' in s][0]
    refine_ls_rfree = [s.split()[-1] for s in st.make_pdb_headers().split('\n') if 'BIN FREE R VALUE SET COUNT' in s][0]
    refine_ls_reflection_number_all = [s.split()[-1] for s in st.make_pdb_headers().split('\n') if 'REFLECTION IN BIN' in s][0]
    representative_sblock.set_pair("_refine_ls_shell.d_res_high", d_res_h)
    representative_sblock.set_pair("_refine_ls_shell.d_res_low", d_res_l)
    representative_sblock.set_pair("_refine_ls_shell.number_reflns_R_free", refine_ls_rfree)
   
    representative_sblock = insert_pair_into_cif_block(
                                representative_sblock,
                                "_refine_ls_shell",
                                ("_refine_ls_shell.number_reflns_R_work", refine_ls_reflection_number_all),
                            )
    
    representative_sblock.set_pair("_refine_ls_shell.pdbx_refine_id", "\"X-RAY DIFFRACTION\"")

    # insert default cross-validation method to metadata block
    representative_sblock = insert_pair_into_cif_block(
                                representative_sblock,
                                "_refine",
                                ("_refine.pdbx_method_to_determine_struct", "\"MOLECULAR REPLACEMENT\""),
                            )


    # load processing stats from ispyb xml

    rstats_xml = [group_dep_df.loc[group_dep_df["xtal_id"] == xid, "xml_path"].values[0] for xid in unmodelled_xtal_ids]

    for xid, xml_path in zip(unmodelled_xtal_ids, rstats_xml):
        if not Path(xml_path).exists():
            raise ValueError(f"missing ispyb xml file at {xml_path} for {xid}")
        
    list_of_dicts = []
    for xml_path in rstats_xml:
        list_of_dicts.extend(xml_path_to_pipeline_programs(xml_path))
    unique_programs = [dict(t) for t in {tuple(sorted(d.items())) for d in list_of_dicts}]
    
    # software table updates
    # we need to populate this with some additional info about processing pipeline
    representative_sblock = convert_cif_pairs_to_loop(representative_sblock, '_software')

    # extract pipeline programs as list[dict] from ispyb xml and insert into template block
    software_item = representative_sblock.find_loop_item("_software.pdbx_ordinal")
    software_tags = software_item.loop.tags if software_item else []
    pdbx_ordinal = int(max(representative_sblock.find_loop("_software.pdbx_ordinal"))) + 1 if software_item else 1
    for pipeline_program in unique_programs:
        field_map = {f"{'_software'}.{k}": v for k, v in pipeline_program.items()} | {"_software.pdbx_ordinal": str(pdbx_ordinal)}
        values = []
        for tag in software_tags:
            if tag in field_map:
                values.append(field_map[tag])
            else:
                values.append("?")
        pdbx_ordinal += 1
        software_item.loop.add_row(values)


    doc = gemmi.cif.Document()
    new_block = gemmi.cif.Block("XXXX")
    for item in template_block:
        new_block.add_item(item)
    for item in representative_mblock:
        new_block.add_item(item)
    for item in representative_sblock:
        new_block.add_item(item)
    doc.add_copied_block(new_block)
    
    return doc
