from rdkit import Chem
from rdkit.Chem import AllChem

def get_proton_acceptor_atom(molecule):
    """
    Get the most likely atom that will accept a proton in the molecule.

    Parameters:
    - molecule (RDKit Mol object): The input molecule.

    Returns:
    - proton_acceptor_atom_idx (int): Index of the proton acceptor atom.
    """
    # Calculate Gasteiger-Marsili charges
    AllChem.ComputeGasteigerCharges(molecule)
    
    # Get the atom with the lowest charge as a potential proton acceptor
    proton_acceptor_atom_idx = min(range(molecule.GetNumAtoms()), key=lambda x: molecule.GetAtomWithIdx(x).GetDoubleProp("_GasteigerCharge"))
    
    return proton_acceptor_atom_idx

# Example usage with water molecule
water_smiles = 'O'
water_molecule = Chem.MolFromSmiles(water_smiles)

# Find the most likely proton acceptor atom in water molecule
proton_acceptor_atom_idx = get_proton_acceptor_atom(water_molecule)

print(f"Most likely proton acceptor atom index: {proton_acceptor_atom_idx}")
