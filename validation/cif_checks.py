import os

import gemmi
import json
import pandas as pd

def has_ligand(structure_cif: str, ligand_name: "UNL") -> bool:
    st = gemmi.read_structure(structure_cif)
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.name == ligand_name:
                    return True
    return False

def count_ligands(structure_cif: str, ligand_name: "UNL") -> int:
    st = gemmi.read_structure(structure_cif)
    count = 0
    for model in st:
        for chain in model:
            for residue in chain:
                if residue.name == ligand_name:
                    count += 1
    return count

def count_event_blocks(structure_factor_cif: str, start_idx: int = 2)->int:
    d = gemmi.cif.read_file(structure_factor_cif)
    return len(d[start_idx:])

def cif_loop_check(cif_file: str):
    d = gemmi.cif.read_file(cif_file)
    for block in d:
        for item in block:
            if item.loop:
                if len(set([s.split(".")[0] for s in item.loop.tags])) > 1:
                    raise ValueError(f"Loop item {item.tags} contains tags from multiple categories")
                
def diffrn_id_check(cif_file: str, expected_id: str):
    d = gemmi.cif.read_file(cif_file)
    diffrn_id_found = False
    for block in d:
        diffrn_id = block.find_pair('_diffrn.id')
        if diffrn_id:
            if diffrn_id[1] != expected_id:
                raise ValueError(f"_diffrn.id {diffrn_id[1]} does not match expected {expected_id}")
            diffrn_id_found = True
    
    if not diffrn_id_found:
        raise ValueError(f"_diffrn.id not found in {cif_file}")
                
def map_ligands_to_events(structure_cif: str, structure_factor_cif: str, ligand_names: list[str]=["UNL"], event_block_start_idx: int = 2) -> dict[str, list[int]]:
    st = gemmi.read_structure(structure_cif)
    d = gemmi.cif.read_file(structure_factor_cif)
    ligands = []
    for m in st:
        for c in m:
            for r in c:
                if r.name in ligand_names and len(r) > 0:
                    x = sum(a.pos.x for a in r) / len(r)
                    y = sum(a.pos.y for a in r) / len(r)
                    z = sum(a.pos.z for a in r) / len(r)
                    ligands.append((c.name, r.name, r.seqid.num, round(x, 2), round(y, 2), round(z, 2)))
    event_blocks = d[event_block_start_idx:]
    event_block_jsons = [gemmi.cif.as_string(b.find_value('_diffrn.details')).strip() for b in event_blocks]
    event_block_jsons = [s[1:-1] if s.startswith("'") and s.endswith("'") else s for s in event_block_jsons]
    event_block_data = [json.loads(j) for j in event_block_jsons]
    event_block_pos = [(ebd['event_site_x'], ebd['event_site_y'], ebd['event_site_z']) for ebd in event_block_data]
    for ligand in ligands:
        closest_event = min(event_block_pos, key=lambda ebp: ((ligand[3] - ebp[0])**2 + (ligand[4] - ebp[1])**2 + (ligand[5] - ebp[2])**2)**0.5)
        dist = ((ligand[3] - closest_event[0])**2 + (ligand[4] - closest_event[1])**2 + (ligand[5] - closest_event[2])**2)**0.5
        print(f"{structure_cif}, Ligand: {ligand}, Closest event: {closest_event}, Distance: {round(dist, 2)}")

def soaked_compound_check(structure_factor_cif: str, ligand_df: pd.DataFrame):
    d = gemmi.cif.read_file(structure_factor_cif)
    crystal_id = d[0].find_pair('_diffrn.crystal_id')[1]
    crystal_treatment = gemmi.cif.as_string(d[0].find_value('_diffrn.crystal_treatment')).strip()
    crystal_treatment = crystal_treatment[1:-1]
    crystal_treatment_json = json.loads(crystal_treatment)
    smiles_matches = ligand_df.loc[
        ligand_df['xtal_id'].astype(str).str.strip() == str(crystal_id).strip(),
        'smiles',
    ]
    if smiles_matches.empty:
        raise ValueError(f"Crystal {crystal_id} not found in ligand table")
    smiles_from_df = smiles_matches.iloc[0]
    smiles_from_cif = crystal_treatment_json.get('smiles')
    if smiles_from_cif != smiles_from_df:
        raise ValueError(f"Smiles from CIF {smiles_from_cif} does not match expected {smiles_from_df} for crystal {crystal_id}")
    print(f"Crystal {crystal_id} has matching smiles: {smiles_from_cif}")

def ground_mtz_to_cif_check(dimple_metadata: dict, ground_structure_factor_cif: str):
    print("Checking dimensions of structure factor amplitudes in ground state models...")
    doc = gemmi.cif.read_file(ground_structure_factor_cif)
    for block in doc:
        if block.find_pair('_diffrn.id'):
            diffrn_id = block.find_pair('_diffrn.id')[1]
            dimple_mtz = dimple_metadata[diffrn_id]['dimple_mtz']
            mtz = gemmi.read_mtz_file(dimple_mtz)
            f_mtz = mtz.column_with_label('F')
            f_cif = block.find_loop('_refln.F_meas_au')
            assert len(f_mtz) == len(f_cif), f"Number of reflections in MTZ ({len(f_mtz)}) does not match number of reflections in CIF ({len(f_cif)}) for crystal {diffrn_id}"

def diffrn_ids_disjoint_check(group_dep_dir: str):
    diffrn_ids_ground = set()
    for f in os.listdir(group_dep_dir):
        if f.endswith("ground-sf.cif"):
            d = gemmi.cif.read_file(f"{group_dep_dir}/{f}")
            for block in d:
                if block.find_pair('_diffrn.id'):
                    diffrn_id = block.find_pair('_diffrn.id')[1]
                    if diffrn_id in diffrn_ids_ground:
                        raise ValueError(f"Duplicate _diffrn.id {diffrn_id} found in {f}")
                    diffrn_ids_ground.add(diffrn_id)

    diffrn_ids_hits = set()
    for f in os.listdir(group_dep_dir):
        if f.endswith("-sf.cif") and not f.endswith("ground-sf.cif"):
            d = gemmi.cif.read_file(f"{group_dep_dir}/{f}") # only refinement block
            block = d[0]
            if block.find_pair('_diffrn.id'):
                diffrn_id = block.find_pair('_diffrn.id')[1]
                if diffrn_id in diffrn_ids_hits:
                    raise ValueError(f"Duplicate _diffrn.id {diffrn_id} found in {f}")
                if diffrn_id in diffrn_ids_ground:
                    raise ValueError(f"_diffrn.id {diffrn_id} in {f} also found in ground state CIFs")
                    diffrn_ids_hits.add(diffrn_id) 

    assert diffrn_ids_ground.isdisjoint(diffrn_ids_hits), f"_diffrn.id values are not disjoint between ground state and hit CIFs in {group_dep_dir}" 
    print(f"All _diffrn.id values are unique across CIF files in {group_dep_dir}")