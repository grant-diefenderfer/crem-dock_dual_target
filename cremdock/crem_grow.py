import logging

import numpy as np
from crem.crem import grow_mol
from rdkit import Chem, rdBase
from rdkit.Chem.Crippen import MolLogP
from rdkit.Chem.rdMolDescriptors import CalcTPSA
from scipy.spatial import distance_matrix

from cremdock.auxiliary import calc_rtb
from cremdock.molecules import neutralize_atoms
from cremdock.user_protected_atoms import get_atom_idxs_for_canon


def get_protected_ids(mol, protein_xyz, dist_threshold, mol2=None, protein_xyz2=None):
    """
    Returns list of ids of heavy atoms ids which have ALL hydrogen atoms close to the protein
    :param mol: molecule
    :param protein_xyz: coordinates of heavy atoms of a protein
    :param dist_threshold: minimum distance to hydrogen atoms
    :return:
    """
    hids = np.array([a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1])
    xyz = mol.GetConformer().GetPositions()[hids]
    min_xyz = xyz.min(axis=0) - dist_threshold
    max_xyz = xyz.max(axis=0) + dist_threshold
    # select protein atoms which are within a box of min-max coordinates of ligand hydrogen atoms
    pids = (protein_xyz >= min_xyz).any(axis=1) & (protein_xyz <= max_xyz).any(axis=1)
    pxyz = protein_xyz[pids]
    m = distance_matrix(xyz, pxyz)  # get matrix (ligandH x protein)
    ids = set(hids[(m <= dist_threshold).any(axis=1)].tolist())  # ids of H atoms close to a protein

    mapping = None
    if mol2 is not None and protein_xyz2 is not None:
        hids2 = np.array([a.GetIdx() for a in mol2.GetAtoms() if a.GetAtomicNum() == 1])
        xyz2 = mol2.GetConformer().GetPositions()[hids2]
        min_xyz2 = xyz2.min(axis=0) - dist_threshold
        max_xyz2 = xyz2.max(axis=0) + dist_threshold
        # select protein atoms which are within a box of min-max coordinates of ligand hydrogen atoms
        pids2 = (protein_xyz2 >= min_xyz2).any(axis=1) & (protein_xyz2 <= max_xyz2).any(axis=1)
        pxyz2 = protein_xyz2[pids2]
        m2 = distance_matrix(xyz2, pxyz2)
        ids2= set(hids2[(m2 <= dist_threshold).any(axis=1)].tolist())

        mol1_no_h = Chem.RemoveHs(mol)
        mol2_no_h = Chem.RemoveHs(mol2)
        match = mol2_no_h.GetSubstructMatch(mol1_no_h)
        if match and len(match) == mol1_no_h.GetNumAtoms():
            m1_map = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
            m2_map = [a.GetIdx() for a in mol2.GetAtoms() if a.GetAtomicNum() > 1]
            mapping = {m1_map[i]: m2_map[j] for i, j in enumerate(match)}
    

    output_ids = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() > 1:
            mol_h_neighbors = set(n.GetIdx() for n in a.GetNeighbors() if n.GetAtomicNum() == 1)
            # all hydrogens of a heavy atom are close to protein
            mol2_h_neighbors = None
            if mol2 is not None and mapping is not None:
                try:
                    a2_idx = mapping.get(a.GetIdx())
                    if a2_idx is not None:
                        a2 = mol2.GetAtomWithIdx(a2_idx)
                        mol2_h_neighbors = set(n.GetIdx() for n in a2.GetNeighbors() if n.GetAtomicNum() == 1)
                except Exception:
                    pass
            if mol2 is not None and protein_xyz2 is not None and mol2_h_neighbors is not None:
                # only consider protection if there are hydrogens to protect in the corresponding model
                prot1 = not (mol_h_neighbors - ids) if mol_h_neighbors else False
                prot2 = not (mol2_h_neighbors - ids2) if mol2_h_neighbors else False
                if prot1 or prot2:
                    output_ids.append(a.GetIdx())
            else:
                if mol_h_neighbors and not (mol_h_neighbors - ids):
                    output_ids.append(a.GetIdx())

    return output_ids


def __get_protein_heavy_atom_xyz(protein):
    """
    Returns coordinates of heavy atoms
    :param protein: protein file (pdb or pdbqt), explicit hydrogens are not necessary
    :return: 2d_array (n_atoms x 3)
    """
    with open(protein) as f:
        pdb_block = f.read()
    return get_protein_heavy_atoms_xyz_from_string(pdb_block)


def get_protein_heavy_atoms_xyz_from_string(pdb_block):
    protein = Chem.MolFromPDBBlock('\n'.join([line[:66] for line in pdb_block.split('\n')]), sanitize=False)
    if protein is None:
        raise ValueError("Protein structure is incorrect. Please check protein pdbqt file.")
    xyz = protein.GetConformer().GetPositions()
    xyz = xyz[[a.GetAtomicNum() > 1 for a in protein.GetAtoms()], ]
    return xyz


def grow_mol_crem(mol, protein_xyz, max_mw, max_rtb, max_logp, max_tpsa, protein_xyz2=None, h_dist_threshold=2, ncpu=1, **kwargs):
    mol_0 = neutralize_atoms(mol)  # add neutralize_atoms to calc correct logp and tpsa
    mw = max_mw - Chem.Descriptors.MolWt(mol_0)
    if mw <= 0:
        return []
    rtb = max_rtb - calc_rtb(mol_0) - 1  # it is necessary to take into account the formation of bonds during the growth of the molecule
    if rtb == -1:
        rtb = 0
    logp = max_logp - MolLogP(mol_0) + 0.5
    tpsa = max_tpsa - CalcTPSA(mol_0)

    mol = Chem.AddHs(mol, addCoords=True)

    mol2=None
    if protein_xyz2 is not None and mol.HasProp('mol_block_2'):
        mol2 = Chem.MolFromMolBlock(mol.GetProp('mol_block_2'), removeHs=False)
        if mol2 is not None:
            mol2 = Chem.AddHs(mol2, addCoords=True)

    _protected_user_ids = set()
    if mol.HasProp('protected_user_canon_ids'):
        _protected_user_ids = set(
            get_atom_idxs_for_canon(mol, list(map(int, mol.GetProp('protected_user_canon_ids').split(',')))))
    _protected_alg_ids = set(get_protected_ids(mol, protein_xyz, h_dist_threshold, mol2=mol2, protein_xyz2=protein_xyz2))
    protected_ids = _protected_alg_ids | _protected_user_ids

    # remove explicit hydrogen and charges and redefine protected atom ids
    for i in protected_ids:
        mol.GetAtomWithIdx(i).SetIntProp('__tmp', 1)
    mol = neutralize_atoms(mol)
    protected_ids = []
    for a in mol.GetAtoms():
        if a.HasProp('__tmp') and a.GetIntProp('__tmp'):
            protected_ids.append(a.GetIdx())
            a.ClearProp('__tmp')

    blocker = rdBase.BlockLogs()  # suppress CReM warnings, https://github.com/rdkit/rdkit/issues/2683
    logging.debug(f'conformer positions: {mol.GetConformer().GetPositions()[:3]}')
    logging.debug(f'protected_ids: {protected_ids}')
    logging.debug(f'num free attachment points: {mol.GetNumAtoms() - len(protected_ids)}')
    try:
        res = list(grow_mol(mol, protected_ids=protected_ids, return_rxn=False, return_mol=True, ncores=ncpu,
                            symmetry_fixes=True, mw=(1, mw), rtb=(0, rtb), logp=(-100, logp), tpsa=(0, tpsa), **kwargs))

    except Exception as e:
        logging.error(f'grow error, {mol.GetProp("_Name")} {Chem.MolToSmiles(mol)}, {e}',
                      stack_info=True, exc_info=True)
        res = []

    res = tuple(m for smi, m in res)

    return res


def grow_mols_crem(mols, protein_xyz, max_mw, max_rtb, max_logp, max_tpsa,protein_xyz2=None, h_dist_threshold=2, ncpu=1, **kwargs):
    """

    :param mols: list of molecules
    :param protein_xyz: 2D array of heavy atoms coordinates
    :param max_mw:
    :param max_rtb:
    :param max_logp:
    :param max_tpsa:
    :param h_dist_threshold: maximum distance from H atoms to the protein to mark them as protected from further grow
    :param ncpu: number of cpu
    :param kwargs: arguments passed to crem function grow_mol
    :return: dict of parent mols and lists of corresponding generated mols
    """
    res = dict()
    for mol in mols:
        tmp = grow_mol_crem(mol, protein_xyz, max_mw=max_mw, max_rtb=max_rtb, max_logp=max_logp, max_tpsa=max_tpsa,
                            protein_xyz2=protein_xyz2, h_dist_threshold=h_dist_threshold, ncpu=ncpu, **kwargs)
        if tmp:
            res[mol] = tmp
    return res