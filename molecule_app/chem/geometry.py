from rdkit.Chem import AllChem, rdMolTransforms


def generate_3d_conformer(mol):
    """Generate and optimize a 3D conformer."""
    params = AllChem.ETKDGv3()
    if AllChem.EmbedMolecule(mol, params) != 0:
        raise ValueError("3D embedding failed.")
    AllChem.UFFOptimizeMolecule(mol)
    return mol


def compute_bond_angles(mol, max_angles=30):
    """Compute sample bond angles (i-j-k)."""
    conf = mol.GetConformer()
    angles = []
    for j in range(mol.GetNumAtoms()):
        nbrs = [a.GetIdx() for a in mol.GetAtomWithIdx(j).GetNeighbors()]
        for a in range(len(nbrs)):
            for b in range(a + 1, len(nbrs)):
                i, k = nbrs[a], nbrs[b]
                angle = rdMolTransforms.GetAngleDeg(conf, i, j, k)
                angles.append((i, j, k, round(angle, 2)))
                if len(angles) >= max_angles:
                    return angles
    return angles
