import numpy as np
import gemmi


def event_map_stats(
    st: gemmi.Structure,
    density_map: gemmi.Ccp4Map,
    res: gemmi.Residue,
    tol: np.float64=0.01,
):

    grid = density_map.grid
    a,b,c = grid.unit_cell.a, grid.unit_cell.b, grid.unit_cell.c
    a_, b_, c_ = st.cell.a, st.cell.b, st.cell.c
    if not (np.isclose(a,a_,atol=tol) and np.isclose(b,b_,atol=tol) and np.isclose(c,c_,atol=tol)):
        raise Exception('mismatching cell constants for density and structure!')
    density_map.setup(np.nan)
    grid.spacegroup = gemmi.SpaceGroup("P 1")
    density_map.update_ccp4_header()
    dencalc = gemmi.DensityCalculatorX()
    dencalc.d_min = st.resolution
    dencalc.grid.setup_from(st)
    dencalc.put_model_density_on_grid(st[0])
    dencalc.grid.spacegroup = gemmi.SpaceGroup("P 1")
    rho = np.zeros((len(res), 2))
    for idx, atom in enumerate(res):
        rho[idx, 0] = grid.interpolate_value(atom.pos)
        rho[idx, 1] = dencalc.grid.interpolate_value(atom.pos)
    return np.mean(rho[:, 0]), np.var(rho[:, 0]), np.corrcoef(rho.T)[0, 1]

def real_space_map_stats(
    st: gemmi.Structure,
    rblock: gemmi.ReflnBlock,
    res: gemmi.Residue,
    amplitude_label: str = "pdbx_FWT",
    phase_label: str = "pdbx_PHWT",
):
    """
    Computes real-space electron density statistics for a given residue.

    Parameters
    ----------
    st : gemmi.Structure
        The refined ensemble model.
    rblock : gemmi.ReflnBlock
        Reflection data block containing refined structure factors.
    res : gemmi.Residue
        The residue for which statistics will be computed.
    amplitude_label : str, optional
        The label for structure factor amplitudes in the reflection block.
    phase_label : str, optional
        The label for phases in the reflection block.

    Returns
    -------
    tuple of float
        A tuple containing:
        - Mean density from the experimental map over residue atoms.
        - Variance of the experimental map values over residue atoms.
        - Pearson correlation coefficient between experimental and model density.

    Notes
    -----
    Default labels assume phenix.refine generated structure/reflections.
    """

    density_map = rblock.transform_f_phi_to_map(amplitude_label, phase_label)
    dencalc = gemmi.DensityCalculatorX()
    dencalc.d_min = np.min(rblock.make_d_array())
    dencalc.grid.setup_from(st)
    dencalc.put_model_density_on_grid(st[0])
    rho = np.zeros((len(res), 2))
    for idx, atom in enumerate(res):
        rho[idx, 0] = density_map.interpolate_value(atom.pos)
        rho[idx, 1] = dencalc.grid.interpolate_value(atom.pos)
    return np.mean(rho[:, 0]), np.var(rho[:, 0]), np.corrcoef(rho.T)[0, 1]


def lig_prot_b_iso_ratio(
    st: gemmi.Structure,
    res: gemmi.Residue,
):
    """
    Occ-weighted ratio of avg isotropic B-factors between ligand/protein.

    Parameters
    ----------
    st : gemmi.Structure
        Structure containing both the ligand and protein chains.
    res : gemmi.Residue
        Ligand residue whose B-factors will be averaged.

    Returns
    -------
    float
        The ratio of average ligand B_iso to average protein B_iso (weighted by occupancy).
    """

    lig_b_iso = []
    for atom in res:
        lig_b_iso.append(atom.b_iso)
    prot_b_iso = []
    sel = gemmi.Selection(";polymer")
    for model in st:
        for chain in model:
            for residue in chain:
                for atom in sel.atoms(residue):
                    prot_b_iso.append(atom.b_iso * atom.occ)

    return np.mean(np.array(lig_b_iso)) / np.mean(np.array(prot_b_iso))


def real_space_diff_zscore(
    st: gemmi.Structure,
    rblock: gemmi.ReflnBlock,
    res: gemmi.Residue,
    amplitude_label: str = "pdbx_DELFWT",
    phase_label: str = "pdbx_DELPHWT",
):
    """
    Computes the average Z-score of difference electron density over a given residue.

    Parameters
    ----------
    st : gemmi.Structure
        The structure containing the atomic model.
    rblock : gemmi.ReflnBlock
        Reflection block with difference map coefficients (e.g., Fo-Fc).
    res : gemmi.Residue
        The residue to evaluate.
    amplitude_label : str, optional
        Label for the amplitude column in the reflection block.
    phase_label : str, optional
        Label for the phase column in the reflection block.

    Returns
    -------
    float
        Mean Z-score of interpolated difference density values over atoms in the residue.
    """
    density_map = rblock.transform_f_phi_to_map(amplitude_label, phase_label)
    sigma = np.std(density_map.array)
    rho = []
    for atom in res:
        rho.append(density_map.interpolate_value(atom.pos))
    return np.mean(np.array(rho) / sigma)
