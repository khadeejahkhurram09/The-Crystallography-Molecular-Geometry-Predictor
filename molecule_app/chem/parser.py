from rdkit import Chem


def parse_molecule(user_input: str):
    """Parse SMILES/InChI into an RDKit molecule."""
    mol = Chem.MolFromSmiles(user_input)
    if mol is None and user_input.startswith("InChI="):
        mol = Chem.MolFromInchi(user_input)
    if mol is None:
        raise ValueError("Invalid molecule input. Use SMILES or InChI.")
    mol = Chem.AddHs(mol)
    return mol
