from itertools import combinations
from typing import Optional
import gemmi
import numpy as np
import string
from collections import OrderedDict
from sklearn.cluster import DBSCAN
from sklearn.metrics import pairwise_distances
from .sync_solvent import sync_solvent_labels


# map altlocs
# expand altlocs
# remove reference altlocs
# relabel donor structure altlocs
# merge structures
# merge altlocs if all match, relabel as reference altlocs
# apply occupancy policy

def atoms_from_altloc(res: gemmi.Residue, altloc: str="A") -> list[gemmi.Atom]:
    """Return a list of atoms from a given residue"""
    atoms = [a for a in res if a.altloc == altloc]
    if len(atoms) == 0:
        raise ValueError("empty altloc atoms")
    else:
        return atoms

def atom_from_altloc_and_name(res: gemmi.Residue, altloc: str=None, atom_name: str=None)->gemmi.Atom:
    atoms = [a for a in res if a.altloc == altloc and a.name == atom_name]
    if len(atoms) > 1:
        raise ValueError(f'multiple atoms found in {res} for {atom_name} and altloc {altloc}')
    if len(atoms) == 0:
        raise ValueError(f'no atoms found in {res} for {atom_name} and altloc {altloc}')
    return atoms[0]

def generate_pair_list(items: set) -> Optional[list[str]]:
    return [set(pair) if len(items) > 1 else {} for pair in list(combinations(sorted(items),2))]

def atom_altlocs(gemmi_residue: gemmi.Residue=None)->dict:
    atom_names = set()
    atom_altlocs = {}
    for a in gemmi_residue:
        atom_names.add(a.name)

    for an in atom_names:
       atom_altlocs[an] = set()

    for a in gemmi_residue:
        atom_altlocs[a.name].add(a.altloc)
    return atom_altlocs

def map_altlocs(st: gemmi.Structure, model_idx: int=0)->list[dict]:
    residues = []
    for c in st[model_idx]:
        for r in c:
            #chain/seqid/icode used for fast lookup
            residue = {
                'structure_ref': st,
                'resname': r.name,
                'model': model_idx,
                'chain': c.name,
                'seqid': r.seqid.num,
                'icode': r.seqid.icode,
                'atom_altlocs': {},
            }
            residue['atom_altlocs'] = atom_altlocs(r)
            residues.append(residue)
    return residues


def max_altlocs(residues: list[dict]) -> set:
    """
    Determine the maximal set of alternate-location identifiers (altlocs)
    present across a list of residue records.

    For each residue, the union of all altloc identifiers present in its
    atoms is computed. The function returns the altloc set with the
    greatest cardinality. If multiple residues share this maximal size,
    their altloc sets must agree exactly; otherwise an exception is
    raised, indicating inconsistent conformer labeling in the input.

    Parameters
    ----------
    residues : list of dict
        Residue dictionaries, each containing an ``"atom_altlocs"`` field
        mapping atom names to sets of altloc identifiers.

    Returns
    -------
    set of str
        The maximal altloc set across all residues.

    Raises
    ------
    ValueError
        If residues with the same maximal number of altlocs do not share
        an identical set of altloc identifiers.
    """
    
    # Collect each residue's union of altloc identifiers, flatten
    altloc_sets = [
        set().union(*res['atom_altlocs'].values())
        if res.get("atom_altlocs") else set()
        for res in residues
    ]

    if not altloc_sets:
        return set()

    # Identify the set with maximum size
    max_set = max(altloc_sets, key=len)
    same_size_sets = [s for s in altloc_sets if len(s) == len(max_set)]

    # Validate that all residues with maximal size agree on altloc IDs
    if len(same_size_sets) > 1:
        combined = set().union(*same_size_sets)
        if combined != max_set:
            raise ValueError(
                f"Inconsistent altloc IDs across residues (expected {max_set}, found {combined})."
            )

    return max_set


def get_gemmi_residue(residue: dict)->gemmi.Residue:
    st = residue['structure_ref']
    model = st[residue['model']]
    return model.sole_residue(residue['chain'], gemmi.SeqId(residue['seqid'],residue['icode']))

def copy_insert_reference_altloc(residue: dict, copied_altloc_id: str=None)->None:
    if copied_altloc_id is None:
        raise Exception('copied_conformer_id not defined')

    gemmi_res = get_gemmi_residue(residue)
    atoms_to_add = []
    for a in gemmi_res:
        if a.altloc == '\x00':
            if copied_altloc_id in residue['atom_altlocs'][a.name]:
                raise ValueError(f'altloc {copied_altloc_id} already exists for {a} in {residue}')
            a_to_add = a.clone()
            a_to_add.altloc = copied_altloc_id
            atoms_to_add.append(a_to_add)
            
    [gemmi_res.add_atom(a) for a in atoms_to_add] #avoid simultaneous adding/iterating
    
    # update atom altlocs after copying
    residue['atom_altlocs'] = atom_altlocs(gemmi_res)

def remove_altloc(residue: dict, remove_altloc_id: str=None)->None:
    if remove_altloc_id not in set().union(*residue.get('atom_altlocs',{}).values()):
        raise ValueError(f'altloc {remove_altloc_id} not in residue {residue}')
    gemmi_res = get_gemmi_residue(residue)
    for i in range(len(gemmi_res)-1,-1,-1): #iterate through backwards
        if gemmi_res[i].altloc == remove_altloc_id:
            del gemmi_res[i]
    residue['atom_altlocs'] = atom_altlocs(gemmi_res)

def remove_altlocs_by_atom(gemmi_res: gemmi.Residue, altlocs: list[str]=None, atom_name: str=None):
    for i in range(len(gemmi_res)-1,-1,-1):
        if gemmi_res[i].altloc in altlocs and gemmi_res[i].name == atom_name:
            del gemmi_res[i]

  
def relabel_altloc(residue: dict, from_label: str=None, to_label: str=None):
    gemmi_res = get_gemmi_residue(residue)
    if from_label == to_label:
        raise ValueError('cannot change label, labels are the same')

    else:
        atoms_to_relabel = atoms_from_altloc(gemmi_res, from_label)
        for a in atoms_to_relabel:
            if a.name == to_label:
                raise ValueError(f'cannot change {a.name} from: {from_label} to: {to_label} for {residue}, atoms in to_label exist')
            if a.altloc == from_label:
                a.altloc = to_label
    residue['atom_altlocs'] = atom_altlocs(gemmi_res)


def get_reference_atom(gemmi_residue: gemmi.Residue, atom_name: str)->gemmi.Atom:
    for a in gemmi_residue:
        if a.altloc == '\x00' and a.name == atom_name:
            return a
    raise ValueError(f'no reference altloc found for {gemmi_residue} and {atom_name}')

def fill_altlocs(residue: dict, altlocs_to_fill: set=None):

    gemmi_res = get_gemmi_residue(residue)
    all_atom_altlocs = set().union(*residue['atom_altlocs'].values())
    if not all_atom_altlocs.issubset(altlocs_to_fill|{'\x00'}):
        raise ValueError(f'current altlocs for: {residue} are not subset of max altlocs')

    atoms_to_add = []
    for a in gemmi_res:
        if not a.has_altloc():
            reference_atom = get_reference_atom(gemmi_res, a.name)
            altlocs_to_add = altlocs_to_fill - residue['atom_altlocs'][a.name]
            for altloc in sorted(altlocs_to_add):
                atom_to_add = reference_atom.clone()
                atom_to_add.altloc = altloc
                atoms_to_add.append(atom_to_add)

    [gemmi_res.add_atom(a) for a in atoms_to_add]
    residue['atom_altlocs'] = atom_altlocs(gemmi_res)

def fill_all_residue_altlocs(residues: list[dict])->None:
    altlocs_to_fill = max_altlocs(residues)
    for residue in residues:
        fill_altlocs(residue, altlocs_to_fill)
    
def remove_reference_altlocs(residues: list[dict])->None:
    for residue in residues:
        if '\x00' in set().union(*residue['atom_altlocs'].values()):
            remove_altloc(residue, '\x00')

def relabel_residue_altlocs(residue: dict, basis_label: str=None)->None:
    if len(basis_label) > 1:
        raise ValueError(f'basis_label {basis_label} can only be single char')
    if basis_label not in string.ascii_uppercase:
        raise ValueError(f'basis_label must be capital letter but is {basis_label}')
    if basis_label is None:
        raise ValueError('no basis_label provided')
    st = residue['structure_ref']
    gemmi_res = get_gemmi_residue(residue)
    _ra = sorted(set().union(*residue['atom_altlocs'].values()), reverse=True)
    
    # update structure
    update_map = OrderedDict(zip(_ra, map(lambda x: string.ascii_uppercase[ord(x) - 64 + ord(basis_label) - 64 - 1], _ra)))
    for a in gemmi_res:
        a.altloc = update_map[a.altloc]
    # update residue dict
    residue['atom_altlocs'] = atom_altlocs(gemmi_res)

def relabel_all_residue_altlocs(residues: list[dict], basis_label='A')->None:
    for residue in residues:
        relabel_residue_altlocs(residue, basis_label)

def generate_expanded_structure(residues: list[dict], basis_label: str=None)->None:
    fill_all_residue_altlocs(residues)
    remove_reference_altlocs(residues)
    if basis_label:
        relabel_all_residue_altlocs(residues, basis_label)
    return

def map_residues_to_merge(donor: list[dict], acceptor: list[dict], keys=("chain","seqid","icode"))->list[dict]:
    generate_expanded_structure(acceptor)
    generate_expanded_structure(donor, sorted(max_altlocs(acceptor))[-1])

    def key_tuple(d: dict)->tuple:
        return tuple(d[k] for k in keys)

    donor_index = {key_tuple(d): d for d in donor}
    acceptor_index = {key_tuple(d): d for d in acceptor}
    results = []

    # reference if no corresponding acceptor
    acceptor_structures = set()
    for residue in acceptor:
        acceptor_structures.add(residue['structure_ref'])
    if len(acceptor_structures) > 1:
        raise ValueError('multiple acceptor structures detected')
    else:
        acceptor_structure_ref = acceptor_structures.pop()
    
    for k in donor_index.keys() | acceptor_index.keys():
        d = donor_index.get(k)
        a = acceptor_index.get(k)
        entry = {
            "acceptor_structure_ref": acceptor_structure_ref,
            "acceptor": a,
            "donor": d,
            "complete": (d is not None and a is not None),
        }
        results.append(entry)

    return results

def get_chain_from_structure(st: gemmi.Structure, model_idx: int=0, chain: str=None)->gemmi.Chain:
    for c in st[model_idx]:
        if c.name == chain:
            return c
    raise ValueError(f'chain with name {chain} not found')

def merge_residue(residue_to_merge: dict)->None:
    donor_residue = residue_to_merge['donor']
    acceptor_residue = residue_to_merge['acceptor']

    # case residue placed in donor with no corresponding acceptor
    if acceptor_residue is None and donor_residue is not None:
        acceptor_structure_ref = residue_to_merge['acceptor_structure_ref']
        c = get_chain_from_structure(acceptor_structure_ref, 0, donor_residue['chain'])
        gemmi_donor_residue = get_gemmi_residue(donor_residue)
        c.add_residue(gemmi_donor_residue)

        # update
        residue = {
            'chain': c.name,
            'structure_ref': acceptor_structure_ref,
            'model': donor_residue['model'],
            'seqid': gemmi_donor_residue.seqid.num,
            'icode': gemmi_donor_residue.seqid.icode,
            'atom_altlocs': atom_altlocs(gemmi_donor_residue)
        }
        residue_to_merge['acceptor'] = residue

    elif acceptor_residue is not None and donor_residue is None:
        pass #we leave acceptor as is

    # this case should not happen
    elif acceptor_residue is None and donor_residue is None:
        raise ValueError('missing residues, both donor/acceptor are None')
        
    else:
        gemmi_donor_residue = get_gemmi_residue(donor_residue)
        gemmi_acceptor_residue = get_gemmi_residue(acceptor_residue)
        for a in gemmi_donor_residue:
            a_ = a.clone()
            gemmi_acceptor_residue.add_atom(a_)
        acceptor_residue['atom_altlocs'] = atom_altlocs(gemmi_acceptor_residue)

def check_pairs(residue: dict, threshold: float=0.01)->None:
    """Determine whether or not two altlocs have the same conformation.
    This information will be used to decide if altlocs will be merged or if refinement
    restraints need to be generated for a given pair. This check is performed per atom.

    Populates:
      residue['atom_altloc_pair_list'][atom_name] -> list of pairs of altloc IDs
      residue['atom_altloc_pair_match'][atom_name][frozenset({alt1, alt2})] -> bool
    """

    gemmi_residue = get_gemmi_residue(residue)
    
    residue['atom_altloc_pair_list']: dict[str, list] = {}
    residue['atom_altloc_pair_match']: dict[str, dict[frozenset[str], bool]]= {}
    
    for atom_name in residue['atom_altlocs'].keys():
        residue['atom_altloc_pair_list'][atom_name] = generate_pair_list(residue['atom_altlocs'][atom_name])
 
    for atom_name, altloc_pairs in residue['atom_altloc_pair_list'].items():
        residue['atom_altloc_pair_match'][atom_name] = {}
        for pair in altloc_pairs:
            if len(pair) != 2:
                raise ValueError(f'malformed pair for {residue} pair {pair}')
                
            altloc_1, altloc_2 = pair
            atom_1 = atom_from_altloc_and_name(gemmi_residue, altloc_1, atom_name)
            atom_2 = atom_from_altloc_and_name(gemmi_residue, altloc_2, atom_name)
            pair_match_key = frozenset({altloc_1, altloc_2})
            
            if atom_1.pos.dist(atom_2.pos) <= threshold:
                residue['atom_altloc_pair_match'][atom_name][pair_match_key] = True
            else:
                residue['atom_altloc_pair_match'][atom_name][pair_match_key] = False

def telescope_altlocs(residue: dict, altlocs: set):
    """Collapses multiple identical conformers into a single reference altloc
    when all altlocs for that atom match within the threshold. Only occurs
    if the residue has the maximum number of possible altlocs.
    After modification, the residue's altloc mapping is rebuilt and
    pairwise comparisons are recomputed via ``check_pairs``.

    Parameters
    ----------
    residue : dict
        A residue dictionary with keys such as ``'atom_altlocs'`` and
        ``'atom_altloc_pair_match'``. The residue must be compatible with
        ``get_gemmi_residue`` and contain up-to-date altloc information.
    altlocs : set of str
        The maximum/expected set of altloc identifiers for the current
        residue. An atom's altlocs must exactly match this set in order to
        be considered for collapsing into the reference altloc.

    Returns
    -------
    None
        The function modifies ``residue`` in place. Updated fields include:
        ``residue['atom_altlocs']`` and
        ``residue['atom_altloc_pair_match']``.

    Notes
    -----
    This operation is only applied when all conformers for a given atom
    are geometrically equivalent (i.e., within the positional threshold
    used when generating ``atom_altloc_pair_match``). If any altloc pair
    does not match, no collapsing occurs and the altlocs are retained.

    Examples
    --------
    Collapse altlocs for a residue after performing pairwise comparison::

        residue = {...}
        check_pairs(residue)
        telescope_altlocs(residue, altlocs={'A', 'B', 'C'})
    """
    
    gemmi_residue = get_gemmi_residue(residue)
    for atom_name, pair_match in residue['atom_altloc_pair_match'].items():
        if all(pair_match.values()) and residue['atom_altlocs'][atom_name] == altlocs:
            # remove all but one
            sorted_altlocs = sorted(residue['atom_altlocs'][atom_name])
            a = atom_from_altloc_and_name(gemmi_residue, sorted_altlocs[0], atom_name)
            a.altloc = '\x00'
            remove_altlocs_by_atom(gemmi_residue, sorted_altlocs[1:], atom_name)

        else:
            # make restraints / do nothing
            pass
    
    residue['atom_altlocs'] = atom_altlocs(gemmi_residue)
    check_pairs(residue)

def generate_refmac_dist_restraints(residue: dict, sigma: float=0.02, restraint_bond_type: int=1):
    gemmi_residue = get_gemmi_residue(residue)
    residue['atom_altloc_refmac_restraints'] = {}
    for atom_name, pairs in residue['atom_altloc_pair_match'].items():
        residue['atom_altloc_refmac_restraints'][atom_name] = []
        if not all(pairs.values()) or not set().union(*residue['atom_altlocs'].values()) == max_altlocs:            
            matched_pairs = [p for p,v in pairs.items() if v]
            for matched_pair in matched_pairs:
                altloc_1, altloc_2 = matched_pair
                restraint = (
                    f"exte dist first  chain {residue['chain']} resi {residue['seqid']:>4} "
                    f"alte {altloc_1} atom {atom_name:<4} "
                    f"second chain {residue['chain']} resi {residue['seqid']:>4} "
                    f"alte {altloc_2} atom {atom_name:<4} "
                    f"value 0.0 sigma {sigma:.2f} type {restraint_bond_type}"
                )
                residue['atom_altloc_refmac_restraints'][atom_name].append(restraint)


def update_atom_occupancy(residue: dict, atom_name: str, altloc: str, occ: np.float64):
    gemmi_residue = get_gemmi_residue(residue)
    atom_to_update = [a for a in gemmi_residue if a.name == atom_name and a.altloc == altloc]

    if len(atom_to_update) > 1:
        raise ValueError(f'multiple {atom_name} found in {residue}')
    elif len(atom_to_update) == 0:
        raise ValueError(f'no {atom_name} found in {residue}')
    else:
        atom_to_update[0].occ = occ

def apply_occupancy_policy(residue: dict, bdc: np.float64=0.8):
    """The acceptor structure has been updated to include all atoms in the ensemble.
    We can recover the original altlocs present in the acceptor or ground state
    by looking at {current acceptor} - {donor}. This information is needed to
    properly scale occupancy. Operations are done per atom.
    Acceptor atom occupancy is bdc*(1/n_a), and donor occ is (1-bdc)*(1/n_b), bdc
    is the background density correction and is an approximation of which is
    subsequently refined.
    
    There are three cases:
    1. acceptor present, donor present
    2. acceptor present, donor empty
    3. acceptor empty, donor present

    Atom occupancies are modified in place.
    """

    for atom_name, altloc_set in residue['acceptor']['atom_altlocs'].items():

        if residue['complete']: # case 1
            donor_altloc_set = residue['donor']['atom_altlocs'][atom_name]
            original_acceptor_set = altloc_set - donor_altloc_set
            cardinality_original_acceptor_set = len(original_acceptor_set)
            cardinality_donor_set = len(donor_altloc_set)
            acceptor_occ = bdc*(1/cardinality_original_acceptor_set)
            donor_occ = (1-bdc)*(1/cardinality_donor_set)
            if altloc_set == {'\x00'}:
                update_occ = 1
                update_atom_occupancy(residue['acceptor'], atom_name, altloc_set.pop(), update_occ)
            else:
                for al in altloc_set:
                    if al in altloc_set - donor_altloc_set:
                        update_occ = acceptor_occ
                    elif al in donor_altloc_set:
                        update_occ = donor_occ
                    else:
                        raise ValueError(f'Error updating occ, could not assign altloc {al} in {residue}')
                    update_atom_occupancy(residue['acceptor'], atom_name, al, update_occ)

        else: # case 2, 3
            if residue['donor'] is None and residue['acceptor'] is not None:
                update_occ = bdc*(1/len(altloc_set))
            elif residue['donor'] is not None:
                update_occ = (1-bdc)*(1/len(altloc_set))
            else:
                raise ValueError(f'unknown condition encountered when assigning occupancy')
            for al in altloc_set:
                update_atom_occupancy(residue['acceptor'], atom_name, al, update_occ)

# occupancy restraint generation

def cluster_altloc_atoms(st: gemmi.Structure, **kwargs):
    """
    Cluster all altloc atoms in a structure using DBSCAN, with
    symmetry-aware periodic distance. Validates that atoms from a single
    residue are never split across multiple clusters.

    Returns
    -------
    atoms : list[tuple]
        (chain, seqid, resname, atom_name, altloc, gemmi.Atom)
    labels : np.ndarray
        Cluster label for each atom
    atom_clusters : dict[int, set]
        cluster_label -> set of atom tuples
    residue_clusters : dict[int, set]
        cluster_label -> set of (chain, seqid, resname)
    """

    # -------------------------------------------------------------------
    # 1. Define symmetry-aware distance function
    # -------------------------------------------------------------------
    def periodic_dist(atom1: gemmi.Atom, atom2: gemmi.Atom) -> float:
        nearest = st.cell.find_nearest_pbc_image(atom1.pos, atom2.pos)
        return nearest.dist()

    # -------------------------------------------------------------------
    # 2. Collect all atoms that have altloc identifiers
    # -------------------------------------------------------------------
    atoms = [
        (chain.name, residue.seqid.num, residue.name,
         atom.name, atom.altloc, atom)
        for model in st
        for chain in model
        for residue in chain
        for atom in residue
        if atom.has_altloc()
    ]

    # Pre-compute N×N distance matrix
    atom_objs = [record[-1] for record in atoms]
    D = pairwise_distances(atom_objs, metric=periodic_dist)

    # Cluster using DBSCAN
    labels = DBSCAN(metric="precomputed", **kwargs).fit_predict(D)
    unique_labels = set(labels)

    # -------------------------------------------------------------------
    # 3. Build cluster maps
    # -------------------------------------------------------------------
    residue_clusters = {label: set() for label in unique_labels}
    atom_clusters = {label: set() for label in unique_labels}

    for idx, label in enumerate(labels):
        chain, seqid, resname, _, _, atom = atoms[idx]
        residue_id = (chain, seqid, resname)

        residue_clusters[label].add(residue_id)
        atom_clusters[label].add(atoms[idx])

    # -------------------------------------------------------------------
    # 4. Validate: NO residue may appear in >1 cluster
    # -------------------------------------------------------------------
    all_residues = set().union(*residue_clusters.values())
    total_counts = sum(len(res) for res in residue_clusters.values())

    if len(all_residues) != total_counts:
        # Find exactly which residues are split
        split = set()
        for c1, c2 in combinations(residue_clusters.values(), 2):
            overlap = c1.intersection(c2)
            split.update(overlap)
        raise ValueError(
            f"Residues split across clusters (increase eps?): {split}"
        )

    # -------------------------------------------------------------------
    # 5. Return results
    # -------------------------------------------------------------------
    return atoms, labels, atom_clusters, residue_clusters


def generate_occupancy_restraints(merged_residues: list[dict], **kwargs) -> str:
    """
    Generate REFMAC-style occupancy restraints for a merged ensemble.

    Parameters
    ----------
    merged_residues : list[dict]
        Output from `map_residues_to_merge`, each element containing
        `acceptor_structure_ref`, `acceptor`, `donor`, and `complete`.
        
    **kwargs :
        Additional keyword arguments forwarded to
        ``cluster_altloc_atoms`` (e.g. DBSCAN parameters such as
        ``eps`` or ``min_samples``).

    Returns
    -------
    str
        Multi-line string of occupancy restraints.
    """

    def residue_key(res: dict) -> tuple[str, int, str]:
        """(chain, seqid, resname) key for completeness lookup."""
        return (res['chain'], res['seqid'], res['resname'])

    # ------------------------------------------------------------------#
    # 1. Build completeness map: for each residue, say whether the merge
    #    was complete (acceptor AND donor present) or not.
    # ------------------------------------------------------------------#
    completeness_map: dict[tuple[str, int, str], bool] = {}

    for mr in merged_residues:
        complete = mr['complete']
        donor = mr['donor']
        acceptor = mr['acceptor']

        if complete or (not complete and donor is None):
            ref_residue = acceptor
        elif not complete and donor is not None:
            ref_residue = donor
        else:
            raise ValueError(f"Unexpected merge state encountered for {mr}")

        completeness_map[residue_key(ref_residue)] = complete

    # ------------------------------------------------------------------#
    # 2. Get the (single) acceptor structure and cluster its altloc atoms.
    # ------------------------------------------------------------------#
    
    acceptor_structures = {mr['acceptor_structure_ref'] for mr in merged_residues}
    if len(acceptor_structures) != 1:
        raise ValueError(
            "Expected exactly one acceptor_structure_ref, "
            f"found {len(acceptor_structures)}"
        )
    acceptor_structure = acceptor_structures.pop()

    atoms, labels, atom_clusters, residue_clusters = cluster_altloc_atoms(
        acceptor_structure, **kwargs
    )

    # atoms: (chain, seqid, resname, atom_name, altloc, gemmi.Atom)

    # ------------------------------------------------------------------#
    # 3. Assign a numeric group ID to each altloc within each cluster.
    #    These will be used for refmac unique group IDs.
    # ------------------------------------------------------------------#
    residue_altloc_keys = sorted(
        {(label, atom[4]) for label, atom in zip(labels, atoms)}  # cluster, altloc
    )
    group_id_for_residue_altloc = {
        key: i for i, key in enumerate(residue_altloc_keys, start=1)
    }

    # ------------------------------------------------------------------#
    # 4. Build the restraint string cluster-by-cluster.
    # ------------------------------------------------------------------#
    lines: list[str] = []

    for cluster_label, atom_set in atom_clusters.items():
        # Add all atom lines for this cluster
        group_ids = set()
        for chain, seqid, resname, atom_name, altloc, _ in atom_set:
            residue_altloc_key = (cluster_label, altloc)
            group_id = group_id_for_residue_altloc[residue_altloc_key]
            group_ids.add(group_id)
            lines.append(
                f"occupancy group id {group_id} "
                f"chain {chain} resi {seqid:>4} "
                f"alte {altloc} atom {atom_name}"
            )

        # Decide if this cluster is a complete or incomplete group
        # (based on whether any residue in it was merged completely).
        residues_in_cluster = {(a[0], a[1], a[2]) for a in atom_set}  # chain, seqid, resname
        is_complete_group = any(
            completeness_map[res_key] for res_key in residues_in_cluster
        )

        cluster_completeness = "complete" if is_complete_group else "incomplete"
        lines.append(f'occupancy group alts {cluster_completeness} {" ".join([str(s) for s in sorted(group_ids)])}')

    restraint_string = "\n".join(lines) + "\n"
    return restraint_string


def generate_refmac_restraints(merged_residues: list[dict], **kwargs) -> str:
    """
    Generate a complete REFMAC restraint block consisting of both distance
    restraints and occupancy restraints. Distance restraints are generated
    using the altloc pair correspondence information contained in the residue
    dictionaries. Occupancy restraints use spatial information directly
    from the gemmi structure object combined with information from the merged
    residue dictionary about the completeness of partial residues.

    Parameters
    ----------
    merged_residues : list of dict
        A list of merged residue dictionaries, typically produced by
        ``map_residues_to_merge``. Each dictionary must contain an
        ``"acceptor"`` entry with an ``"atom_altloc_refmac_restraints"``
        field.

    Returns
    -------
    str
        A single multi-line string containing:
        
        * All per-atom distance restraints collected from the acceptor
          residues
        * All occupancy-group restraints generated by
          ``generate_occupancy_restraints``
        * A final ``"occupancy refine"`` directive for REFMAC

        The output is formatted as a ready-to-write REFMAC external
        restraints file.

    Notes
    -----
    This function collects restraint lines from each residue’s
    ``atom_altloc_refmac_restraints`` mapping and appends the occupancy
    restraints generated across all merged residues. No validation is
    performed here; it assumes all previous steps (merging, altloc
    expansion, pair checking, etc.) have already ensured consistency.
    """

    dist_lines = []

    for residue in merged_residues:
        atom_restraints = residue["acceptor"].get("atom_altloc_refmac_restraints", {})
        for restraints in atom_restraints.values():
            dist_lines.extend(restraints)

    distance_restraints = "\n".join(dist_lines) + ("\n" if dist_lines else "")
    occupancy_restraints = generate_occupancy_restraints(merged_residues, **kwargs)

    return f"{distance_restraints}{occupancy_restraints}occupancy refine"


class EnsembleMerger:
    """
    Build a merged structural ensemble from acceptor and donor
    `gemmi.Structure` objects corresponding to ground and bound/changed
    states respectively. The acceptor object is iteratively modified in place
    into a single multi-conformer representation of the ensemble.
    
    Performs residue mapping, altloc merging, conformer equivalence
    detection, altloc telescoping, occupancy assignment, and generation of
    REFMAC-compatible occupancy and distance restraints.

    The workflow implemented in :meth:`run` corresponds to the steps:

        1. Map all altlocs in acceptor and donor structures.
        2. Determine residue-to-residue correspondence.
        3. Merge donor conformers into acceptor residues.
        4. Compute altloc-equivalence relationships.
        5. Collapse redundant altloc IDs ("telescope" them).
        6. Apply ensemble occupancy rules.
        7. Generate per-atom REFMAC distance restraints.
        8. Generate cross-conformer REFMAC occupancy restraints.

    Parameters
    ----------
    acceptor : gemmi.Structure
        Base structure that will receive merged donor conformers.
    donor : gemmi.Structure
        Complementary structure providing additional altloc states.
    bdc : float, optional
        Background density correction (BDC) value used in some occupancy
        policies. Default is ``0.8``.
    occupancy_kwargs : dict, optional
        Extra keyword arguments passed to the occupancy-restraint generator.
        These parameters typically control clustering behavior for
        `cluster_altloc_atoms`, such as ``eps`` or ``min_samples``.
        If ``None``, defaults to ``{'eps': 2.5, 'min_samples': 1}``.

    Attributes
    ----------
    acceptor : gemmi.Structure
        The acceptor model, modified in-place during merging.
    donor : gemmi.Structure
        The donor model used to supplement missing altlocs.
    _acceptor_residues : dict
        Altloc-mapped residue dictionary for the acceptor.
    _donor_residues : dict
        Altloc-mapped residue dictionary for the donor.
    _residues_to_merge : list[dict]
        List of residue merge descriptors produced during mapping.
    _max_acceptor_altlocs : set[str]
        Set of all altloc identifiers observed in the merged acceptor.
    _refmac_restraints : str or None
        Final REFMAC occupancy restraint block generated after running
        the full pipeline.

    Notes
    -----
    - All modifications are applied in-place to the ``acceptor`` structure.
    - The class does not return intermediate structures; the primary outputs
      are the updated acceptor model and the generated REFMAC restraint text.
    """
    
    def __init__(self, acceptor: gemmi.Structure,
                 donor: gemmi.Structure,
                 xtal_id: str=None,
                 bdc: np.float64=0.8,
                 occupancy_kwargs: dict={'eps': 2.5, 'min_samples':1},
                 sync_solvent=True):

        #donor model is merged into acceptor model which is modified in place
        self.xtal_id = xtal_id
        self.acceptor = acceptor
        self.donor = donor
        self.bdc = bdc
        self.occupancy_kwargs = occupancy_kwargs or {}

        #check model, raise exception if possible mismatch detected
        self._model_precheck()

        if sync_solvent:
            self._sync_solvent_labels()
        
        #populated during merging
        self._acceptor_residues = None
        self._donor_residues = None
        self._residues_to_merge = None
        self._max_acceptor_altlocs = None
        self._merged_residues = None
        self._refmac_restraints: str | None = None

    def _model_precheck(self):
        if not self.acceptor.cell.approx(self.donor.cell, 0.001):
            raise ValueError(f"cell mismatch, acceptor: {self.acceptor.cell} vs. donor: {self.donor.cell}")

        if not self.acceptor.spacegroup_hm == self.donor.spacegroup_hm:
            raise ValueError(f"mismatching spacegroups, acceptor: {self.acceptor.spacegroup_hm} vs. donor: {self.donor.spacegroup_hm}")

        if not self.acceptor.resolution == self.donor.resolution:
            raise ValueError(f"mismatching resolution, acceptor: {self.acceptor.resolution} vs. donor: {self.donor.resolution}")

    def _sync_solvent_labels(self):
        sync_solvent_labels(self.acceptor, self.donor)

    def _map_altlocs(self):
        
        self._acceptor_residues = map_altlocs(self.acceptor)
        self._donor_residues = map_altlocs(self.donor)

    def _map_residues_to_merge(self):
        self._residues_to_merge = map_residues_to_merge(self._donor_residues, self._acceptor_residues)

    def _merge_residues(self):
        for residue in self._residues_to_merge:
            merge_residue(residue)

    def _check_pairs(self):
        for residue in self._residues_to_merge:
            check_pairs(residue['acceptor'])
            
    def _telescope_altlocs(self):
        max_acceptor_altlocs = {
            atom.altloc
            for model in self.acceptor
            for chain in model
            for residue in chain
            for atom in residue
        }
        for residue in self._residues_to_merge:
            telescope_altlocs(residue['acceptor'], max_acceptor_altlocs)

    def _apply_occupancy(self):
        for residue in self._residues_to_merge:
            apply_occupancy_policy(residue)

    def _generate_refmac_dist_restraints(self):
        for residue in self._residues_to_merge:
            generate_refmac_dist_restraints(residue['acceptor'])

    def _generate_refmac_restraints(self):
        self._refmac_restraints = generate_refmac_restraints(self._residues_to_merge, **self.occupancy_kwargs)

    def refmac_restraints_to_string(self):
        if self._refmac_restraints:
            return self._refmac_restraints
        else:
            raise ValueError("missing refmac restraints for presumably merged structure")

    def run(self):
        try:
            self._map_altlocs()
            self._map_residues_to_merge()
            self._merge_residues()
            self._check_pairs()
            self._telescope_altlocs()
            self._apply_occupancy()
            self._generate_refmac_dist_restraints()
            self._generate_refmac_restraints()
        except Exception as e:
            print(f"caught {e} while merging {self.acceptor.name}")
            
#----
# need to validate structure, ensure consistent altloc labels per residue 
#create residue dicts
#map residues to merge
#merge residues
#evaluate updated acceptor structure/residue dict
#consolidate altlocs to ground state
#apply occupancy policy
#generate restraints

#restraints
# "exte dist first chain resi alte atom second chain resi atom value sigma type"
# occupancy group id 5 chain A resi  111 alte A
# occupancy group id 5 chain A resi  111 alte A atom CA