import gemmi
import numpy as np
import string


def periodic_dist(atom1: gemmi.Atom, atom2: gemmi.Atom, cell) -> float:
    nearest = cell.find_nearest_pbc_image(atom1.pos, atom2.pos)
    return nearest.dist()


def map_solvent_residues(st: gemmi.Structure):
    solvent_residues = {}
    for model in st:
        for chain in model:
            for res in chain:
                if res.name == "HOH":
                    solvent_residues[(chain.name, res.seqid.num, "HOH")] = list(res)
    return solvent_residues


def generate_acceptor_solvent_alias(
    st: gemmi.Structure, solvent_residues: dict, chain_name="S"
):
    for model in st:
        for chain in model:
            if chain.name == chain_name:
                raise ValueError(f"existing chain {chain_name} found in {model}")
    new_seqid = 1
    solvent_alias = {}
    for res in solvent_residues.keys():
        _, old_seqid, res_name = res
        solvent_alias[res] = (chain_name, new_seqid, res_name)
        new_seqid += 1
    return solvent_alias


def generate_donor_solvent_alias(
    st: gemmi.Structure,
    acceptor_residues: dict,
    acceptor_solvent_alias: dict,
    chain_name="S",
    threshold=2,
):
    for model in st:
        for chain in model:
            if chain.name == chain_name:
                raise ValueError(f"existing chain {chain_name} found in {model}")

    donor_residues = map_solvent_residues(st)
    donor_solvent_alias = {}
    acceptor_residues_set = set(acceptor_residues.keys())

    for donor_residue in donor_residues.keys():
        donor_atoms = donor_residues[donor_residue]

        # check for existing label in acceptor residues
        if donor_residue in acceptor_residues_set:
            # need to check position to ensure proximity
            acceptor_atoms = acceptor_residues[donor_residue]
            distance_checks = [
                periodic_dist(a_i, a_j, st.cell) < threshold
                for a_i in acceptor_atoms
                for a_j in donor_atoms
            ]
            if all(distance_checks):
                # can safely use label
                donor_solvent_alias[donor_residue] = acceptor_solvent_alias[
                    donor_residue
                ]

            elif not any(distance_checks):
                # all atoms are far away, we need to use fresh label
                new_seqid = (
                    max(
                        (
                            seqid
                            for _, seqid, _ in set(donor_solvent_alias.values())
                            | set(acceptor_solvent_alias.values())
                        ),
                        default=0,
                    )
                    + 1
                )
                donor_solvent_alias[donor_residue] = (chain_name, new_seqid, "HOH")

            else:
                raise ValueError(
                    f"Found some possible mixed overlap solvent altlocs in {st} donor_res: {donor_residue}"
                )
        else:
            # the solvent was placed in the donor model, but not acceptor
            # need to check positions of all acceptor atoms
            # if proximal atom is found use that resid label
            # else use fresh label
            matched = False
            for acceptor_residue, acceptor_atoms in acceptor_residues.items():
                distance_checks = [
                    periodic_dist(a_i, a_j, st.cell) < threshold
                    for a_i in acceptor_atoms
                    for a_j in donor_atoms
                ]
                if all(distance_checks):
                    # can safely use label
                    donor_solvent_alias[donor_residue] = acceptor_solvent_alias[
                        acceptor_residue
                    ]
                    # also need to update label in donor solvent dict

                    matched = True
                    break

                elif any(distance_checks) and not all(distance_checks):
                    raise ValueError(
                        f"Found some possible mixed overlap solvent altlocs in {st} donor_res: {donor_residue}"
                    )

            if not matched:
                # all atoms are far away, we need to use fresh label
                new_seqid = (
                    max(
                        (
                            seqid
                            for _, seqid, _ in set(donor_solvent_alias.values())
                            | set(acceptor_solvent_alias.values())
                        ),
                        default=0,
                    )
                    + 1
                )
                donor_solvent_alias[donor_residue] = (chain_name, new_seqid, "HOH")

    return donor_solvent_alias


def insert_solvent_chain(model: gemmi.Model, solvent_residues: dict, solvent_alias: dict, chain_name='S'):
    new_solvent_chain = gemmi.Chain(chain_name)
    inverse_solvent_alias = {v: [] for v in set(solvent_alias.values())}
    for k, v in solvent_alias.items():
        inverse_solvent_alias[v].append(k)
    
    for residue in sorted(inverse_solvent_alias.keys(), key=lambda r: r[1]):
        _, new_seqid, _ = residue
        solvent_atoms = []
        for old_residue in inverse_solvent_alias[residue]:
            solvent_atoms.extend([a.clone() for a in solvent_residues[old_residue]])
        add_solvent_residue_to_chain(solvent_atoms, new_seqid, new_solvent_chain)

    # normalize altlocs, necessary for cases where more than one solvent are assigned
    # to the same solvent residue label
    for residue in new_solvent_chain:
        if len(residue) > 1:
            for atom, altloc_label in zip(residue, list(string.ascii_uppercase)):
                atom.altloc = altloc_label
                atom.occ = 1 / len(residue)
    model.add_chain(new_solvent_chain)


def add_solvent_residue_to_chain(
    atoms: list[gemmi.Atom], seqid: int, chain: gemmi.Chain
):
    res = gemmi.Residue()
    res.seqid.num = seqid
    res.name = "HOH"
    for atom in atoms:
        res.add_atom(atom)
    chain.add_residue(res)
    return


def prune_solvents(st: gemmi.Structure, chains_to_keep=["S"]):
    for model in st:
        for chain in model:
            for res in chain:
                if res.name == "HOH" and chain.name not in chains_to_keep:
                    for i in range(len(res) - 1, -1, -1):
                        del res[i]
    st[0]  # refresh references?


def check_for_one_atom_res_clash(st: gemmi.Structure, resname: str, cutoff=0.01):
    # Map each solvent atom -> (chain_name, residue_number)
    atom_to_res = {
        a: (c.name, r.seqid.num)
        for m in st
        for c in m
        for r in c
        if r.name == resname
        for a in r
    }

    if not atom_to_res:
        return

    atoms = np.array(list(atom_to_res), dtype=object)
    n = len(atoms)

    # Compute only the upper triangle (including diagonal)
    overlap = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i, n):
            overlap[i, j] = periodic_dist(atoms[i], atoms[j], st.cell) < cutoff

    col_counts = overlap.sum(axis=0)
    clashing = col_counts > 1

    if np.any(clashing):
        clashing_residues = [atom_to_res[a] for a in atoms[clashing]]
        raise ValueError(
            f"overlapping residues detected structure: {st}, residues: {clashing_residues}:{resname}"
        )

def check_for_solvent_clash(st: gemmi.Structure):
    check_for_one_atom_res_clash(st, 'HOH')

def check_for_nonpolymer_clashes(st: gemmi.Structure):
    for resname in ['Cl','K','Na','Zn','Br','Ca','Mg','DMS','GOL','SO4','PO4','ATP','UNL','LIG']:
        check_for_one_atom_res_clash(st, resname)


def sync_solvent_labels(acceptor: gemmi.Structure, donor: gemmi.Structure):
    check_for_solvent_clash(acceptor)
    check_for_solvent_clash(donor)
    acceptor_solvent = map_solvent_residues(acceptor)
    donor_solvent = map_solvent_residues(donor)
    acceptor_alias = generate_acceptor_solvent_alias(acceptor, acceptor_solvent)
    # information about acceptor structure is needed to relabel donor solvents
    donor_alias = generate_donor_solvent_alias(donor, acceptor_solvent, acceptor_alias)
    insert_solvent_chain(acceptor[0], acceptor_solvent, acceptor_alias)
    insert_solvent_chain(donor[0], donor_solvent, donor_alias)
    prune_solvents(acceptor)
    prune_solvents(donor)
    check_for_solvent_clash(acceptor)
    check_for_solvent_clash(donor)
