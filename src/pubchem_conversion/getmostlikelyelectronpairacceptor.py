from rdkit import Chem
from rdkit.Chem import AllChem

def find_acceptor_atom_after_reaction(acid_mol, base_mol):
    # Combine acid and base molecules to simulate the reaction
    reaction_result = Chem.CombineMols(acid_mol, base_mol)

    # Generate 3D coordinates for the reaction product
    AllChem.EmbedMolecule(reaction_result)

    # Calculate partial charges for the reaction product
    AllChem.ComputeGasteigerCharges(reaction_result)

    # Find the acceptor atom in the reaction product
    acceptor_candidates = []

    for atom in reaction_result.GetAtoms():
        # Check if the atom has a positive charge
        if atom.GetFormalCharge() > 0:
            acceptor_candidates.append((atom, atom.GetTotalNumHs()))

    if acceptor_candidates:
        # Select the atom with the highest positive charge and available electron pairs
        acceptor_atom, num_hydrogens = max(acceptor_candidates, key=lambda x: x[0].GetDoubleProp("_GasteigerCharge"))
        return acceptor_atom, num_hydrogens

    return None, None

# Define your acid and base molecules as SMILES strings
acid_smiles = "CC(=O)O"
base_smiles = "CCN(CC)CC"

# Create RDKit molecule objects
acid_mol = Chem.MolFromSmiles(acid_smiles)
base_mol = Chem.MolFromSmiles(base_smiles)

# Generate 3D coordinates for the molecules
AllChem.EmbedMolecule(acid_mol)
AllChem.EmbedMolecule(base_mol)

# Calculate partial charges
AllChem.ComputeGasteigerCharges(acid_mol)
AllChem.ComputeGasteigerCharges(base_mol)

# Find the acceptor atom after the reaction
acceptor_atom, num_hydrogens = find_acceptor_atom_after_reaction(acid_mol, base_mol)

# Print results
if acceptor_atom:
    print(f"The most likely atom to accept an electron pair after the reaction is atom {acceptor_atom.GetIdx()} "
          f"with charge {acceptor_atom.GetDoubleProp('_GasteigerCharge')} and {num_hydrogens} available electron pairs.")
else:
    print("No suitable acceptor atom found after the reaction.")
