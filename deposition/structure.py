import gemmi
from .utils import letter_generator
from collections import defaultdict


def fix_formal_charges(func):
    """
    Phenix refine (v1.20) appears to write formal charge sign
    at the end of the charge, e.g. 2+, 1- which causes an issue
    with gemmi type checking when creating structure object from cif
    block. Before handing Block to make_structure_from_block, fix the
    charges so that sign is in front.

    Args:
        gemmi.make_structure_from_block (callable): modify input Block

    Returns:
        callable: modified gemmi.make_structure_from_block

    """

    def wrapper(*args, **kwargs):
        column_name = "_atom_site.pdbx_formal_charge"
        column = args[0].find_loop(column_name)
        for idx, charge in enumerate(column):
            if charge != "?":
                if "+" in charge:
                    new_charge = charge.replace("+", "")
                    column[idx] = new_charge
                elif "-" in charge:
                    new_charge = "".join(["-", charge.replace("-", "")])
                    column[idx] = new_charge
        return gemmi.make_structure_from_block(args[0])

    return wrapper


@fix_formal_charges
def make_structure_from_block_wrapper(*args, **kwargs):
    return gemmi.make_structure_from_block(*args, **kwargs)


def merge_altlocs(st: gemmi.Structure):
    """
    Merge alternate locations (altlocs) in a gemmi Structure object.

    This function iterates through a gemmi Structure object, detects atoms with
    alternate locations (altlocs), and merges them if they are within the RMSD
    threshold.

    Parameters
    ----------
    st : gemmi.Structure
        The molecular structure to process.

    Returns
    -------
    None
        The function modifies the input `st` in place.
    """
    for model in st:
        for chain in model:
            for residue in chain:
                if any(atom.has_altloc() for atom in residue):
                    merge_residue_altlocs(st, residue, threshold=0.02)


def merge_residue_altlocs(st: gemmi.Structure, res: gemmi.Residue, **kwargs):
    """
    Merges alternate location (altloc) atoms at residue level
    within a given RMSD threshold.

    Parameters
    ----------
    st : gemmi.Structure
        The molecular structure being processed. Included to renumber models
        after atom removal.
    res : gemmi.Residue
        The residue containing alternate location atoms.
    **kwargs : dict
        Additional parameters for `residue_altloc_dist`.

    Returns
    -------
    None
        The function modifies `res` in place and updates the structure numbering.

    Notes
    -----
    - If only one altloc remains, its identifier is removed, i.e. set to NULL.
    - The RMSD threshold is set via `kwargs`.
    - Calls `remove_atoms_from_res()` and `residue_altloc_dist()` for processing.
    - bilinear interpolation to assign occupancy, b_iso of merged conformer
    - assumptions used for bilinear interp may not hold with aniso b values

    Examples
    --------
    >>> import gemmi
    >>> structure = gemmi.read_structure("example.pdb")
    >>> residue = structure[0]['A'][10]  # Model 0, Chain A, Residue 10
    >>> merge_residue_altlocs(structure, residue, threshold=0.02)
    """

    altlocs = set()
    [altlocs.add(a.altloc) for a in res]
    altlocs = sorted(altlocs)[::-1]
    for i, al_i in enumerate(altlocs):
        atoms_i = atoms_from_altloc(res, al_i)
        for al_j in altlocs[(i + 1) :]:
            atoms_j = atoms_from_altloc(res, al_j)
            if residue_altloc_dist(atoms_i, atoms_j, **kwargs):
                print(f'removing residue: {res.name} {res.seqid} in {st.name}')
                for a_i, a_j in zip(atoms_i, atoms_j):
                    print(a_i.name, a_i.occ)
                    print(a_j.name, a_j.occ)
                    a_ij_occ = a_j.occ + a_i.occ
                    a_j.b_iso = (a_i.occ / a_ij_occ) * a_i.b_iso + (
                        a_j.occ / a_ij_occ
                    ) * a_j.b_iso
                    a_j.occ = a_ij_occ
                remove_atoms_from_res(res, [al_i])

    # remove altloc identifier if only single altloc remains
    altlocs = set()
    [altlocs.add(a.altloc) for a in res]
    altlocs = sorted(altlocs)[::-1]
    if len(altlocs) == 1:
        if "A" not in altlocs:
            raise ValueError("expected altloc of A for single altloc")
        else:
            for a in res:
                a.altloc = "\0"  # null altloc ''

    st.renumber_models()

    return


def remove_atoms_from_res(res: gemmi.Residue, altlocs: list):
    """Delete atoms from residue object that references structure
    based on altloc identifier.

    Parameters
    ----------
    res : gemmi.Residue
        The residue identified for modification/atom removal.
    altlocs : list
        altloc identifiers to be removed (can be multiple)

    Returns
    -------
    None
        The function modifies the structure in place.

    Notes
    -----
    - Only increment if not deleting an atom, otherwise may cause
    strange side effects due to modifying in place
    """

    i = 0
    while i < len(res):
        if res[i].altloc in altlocs:
            del res[i]
        else:
            i += 1


def atoms_from_altloc(res: gemmi.Residue, altloc="A") -> list:
    """Return a list of atoms from a given residue"""

    atoms = [a for a in res if a.altloc == altloc]
    if len(atoms) == 0:
        raise ValueError("empty altloc atoms")
    else:
        return atoms


def residue_altloc_dist(grp1: list, grp2: list, threshold=0.05) -> bool:
    """Calculate Euclidean distance for all pairs of atoms corresponding
    to different conformers, i.e. altlocs. If any one distance for a given
    pair exceeds the threshold, it is assumed that the conformers are in
    fact different. Otherwise, conformers are assumed to have the same
    xyz coordinates. Summing all distances is avoided, numerical
    imprecision may cause false positives.

    Parameters
    ----------
    grp1 : list
        List of gemmi.Atom objects
    grp2 : list
        List of gemmi.Atom objects, must match exactly list 1 for correspondence.

    Returns
    -------
    bool
        True if same, False if different conformers.
    """

    if len(grp1) != len(grp2):
        raise ValueError("different sized groups of atoms")

    for i, j in zip(grp1, grp2):
        if i.pos.dist(j.pos) > threshold:
            if i.name != j.name:
                raise ValueError("atom name mismatch")
            return False
    else:
        return True


def remove_ground_state(
    st: gemmi.Structure, model_idx: int = 0, res_name: str = "UNL"
) -> gemmi.Structure:
    """Removes ground state conformation from pandda-style generated
    ensemble model.

    Parameters
    ----------
    st : gemmi.Structure
        Structure object to be split.
    model_idx : int
        Integer of model from gemmi structure object to be split.
    res_name : str
        Residue name of ligand, used to identify 'changed' state.

    Returns
    -------
    None
        Modifies gemmi structure in place.
    """

    letter_gen = letter_generator()
    altlocs = set()

    # survey altlocs
    for chain in st[model_idx]:
        for res in chain:
            if res.name == res_name:
                for atom in res:
                    altlocs.add(atom.altloc)

    # remove these atoms (ground state)
    n_atoms_before = st[model_idx].count_atom_sites()
    sel = gemmi.Selection(",".join([":"] + list(altlocs)))
    sel.remove_not_selected(st[model_idx])
    st.renumber_models()
    n_atoms_after = st[model_idx].count_atom_sites()
    if n_atoms_after >= n_atoms_before:
        raise ValueError

    # create name change map
    sorted_altlocs = sorted(list(altlocs))
    altlocs_map = {}
    letter_gen = letter_generator()
    for k in sorted(list(altlocs)):
        altlocs_map[k] = next(letter_gen)

    # update altloc names
    for chain in st[model_idx]:
        for res in chain:
            for atom in res:
                if atom.has_altloc():
                    atom.altloc = altlocs_map[atom.altloc]

    merge_altlocs(st)

    return st

def get_residue(st: gemmi.Structure, chain: gemmi.Chain, seqid: int):
    for m in st:
        for c in m:
            if c.name == chain.name:
                for r in c:
                    if r.seqid.num == seqid:
                        return r
    raise ValueError('residue not found')



def tweak_occupancy(
    st: gemmi.Structure,
    threshold: float = 1.0,
    occ_tweak: float = 2e-3,
    epsilon: float = 3e-7,
):
    """Tweak occupancy values for atoms in a structure to ensure that they sum to <1.0
    for any given set of alternate locations (altlocs). Numerical imprecision
    from refinement mmcif can cause occupancies to sum to slightly above 1.0 which
    can cause issues with PDB validation. This difference can be as small as machine
    epsilon for float64.

    Parameters
    ----------
    st : gemmi.Structure
        The molecular structure to process.
    threshold : float
        The occupancy sum threshold above which a warning is raised. Default is 1.0 to allow for some numerical imprecision.
    occ_tweak : float
        The amount by which to reduce the occupancy of the most occupied atom altloc when the sum exceeds the threshold. Default is 0.002.
    epsilon : float
        A small value to allow for numerical imprecision when comparing occupancy sums to the threshold. Default is 3e-7, which is near machine
        epsilon for float32.
    Returns
    -------
    None
        The function modifies the input `st` in place.
    """
    for model in st:
        for chain in model:
            for residue in chain:
                atoms_by_name = defaultdict(list)

                for atom in residue:
                    atoms_by_name[atom.name].append(atom)

                for atoms in atoms_by_name.values():
                    # Only adjust true alternate conformers (no blank-altloc atoms).
                    altloc_atoms = [atom for atom in atoms if atom.has_altloc()]
                    if len(altloc_atoms) <= 1:
                        continue

                    total_occ = sum(atom.occ for atom in altloc_atoms)
                    # Floating-point sums around 1.0 can land just below threshold.
                    if total_occ < (threshold - epsilon):
                        continue

                    atom_to_adjust = max(altloc_atoms, key=lambda atom: atom.occ)
                    # Ensure strict < threshold after adjustment.
                    required = total_occ - (threshold - epsilon)
                    decrement = max(occ_tweak, required)
                    atom_to_adjust.occ = max(atom_to_adjust.occ - decrement, 0.0)

                    adjusted_total = sum(atom.occ for atom in altloc_atoms)
                    if adjusted_total >= threshold:
                        raise ValueError(
                            f"Occupancy sum for altloc atoms {', '.join(atom.name for atom in altloc_atoms)} "
                            f"in residue {residue.name} {residue.seqid} is {adjusted_total:.6f} and still >= {threshold} "
                            "after tweaking."
                        )
                   

def resolve_entities(st: gemmi.Structure, one_polymer=True):
    """Normalize entity assignment and IDs for a Gemmi structure.

    This allows Gemmi deduplication to merge split polymer entities (for example,
    when chains differ only by modeled vs. unmodeled terminal residues). When the sequence
    is inferred from the structure by refmac5 when outputing mmcif during refinement,
    there may be discrepancies in entity assignment if there are multiple chain copies 
    in the asymmetric unit.

    The workflow performed is:
    1. Recompute residue entity types and subchain labels.
    2. Rebuild st.entities from subchains and deduplicate polymers with setup_entities().
    3. Renumber entity names as numeric strings and relink residue entity IDs.

    Args:
        st: Input structure modified in place.
        one_polymer: If True, enforce a single deduplicated polymer entity.

    Raises:
        NotImplementedError: If one_polymer is False.
    """
    if one_polymer:
        reference_sequence = max([entity.full_sequence for entity in st.entities], key=len)
        for entity in st.entities:
            if entity.entity_type == gemmi.EntityType.Polymer:
                entity.full_sequence = reference_sequence

        st.setup_entities()
        for i, e in enumerate(st.entities, start=1):
            e.name = str(i)
        st.add_entity_ids(overwrite=True)
    else:
        raise NotImplementedError("Multiple polymer entities not yet supported")