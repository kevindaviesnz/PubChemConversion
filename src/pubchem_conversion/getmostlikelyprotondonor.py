from rdkit import Chem
from rdkit.Chem import AllChem

def get_proton_donor_atom(molecule):
    """
    Get the most likely atom that will donate a proton in the molecule.

    Parameters:
    - molecule (RDKit Mol object): The input molecule.

    Returns:
    - proton_donor_atom_idx (int): Index of the proton donor atom.
    """
    # Calculate Gasteiger-Marsili charges
    AllChem.ComputeGasteigerCharges(molecule)
    
    # Get the atom with the highest charge as a potential proton donor
    proton_donor_atom_idx = max(range(molecule.GetNumAtoms()), key=lambda x: molecule.GetAtomWithIdx(x).GetDoubleProp("_GasteigerCharge"))
    
    return proton_donor_atom_idx

# Example usage with HCl molecule
hcl_smiles = 'Cl[H]'
hcl_molecule = Chem.MolFromSmiles(hcl_smiles)

# Find the most likely proton donor atom in HCl molecule
proton_donor_atom_idx = get_proton_donor_atom(hcl_molecule)

print(f"Most likely proton donor atom index: {proton_donor_atom_idx}")
from rdkit import Chem
from rdkit.Chem import AllChem

def find_donor_atom_after_reaction(acid_mol, base_mol):
    # Combine acid and base molecules to simulate the reaction
    reaction_result = Chem.CombineMols(acid_mol, base_mol)

    # Generate 3D coordinates for the reaction product
    AllChem.EmbedMolecule(reaction_result)

    # Calculate partial charges for the reaction product
    AllChem.ComputeGasteigerCharges(reaction_result)

    # Find the electron pair donor atom in the reaction product
    donor_candidates = []

    for atom in reaction_result.GetAtoms():
        # Check if the atom has available electron pairs
        num_hydrogens = atom.GetTotalNumHs()
        if num_hydrogens > 0:
            donor_candidates.append((atom, num_hydrogens))

    if donor_candidates:
        # Select the atom with the most available electron pairs
        donor_atom, num_hydrogens = max(donor_candidates, key=lambda x: x[1])
        return donor_atom, num_hydrogens

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

# Find the electron pair donor atom after the reaction
donor_atom, num_hydrogens = find_donor_atom_after_reaction(acid_mol, base_mol)

# Print results
if donor_atom:
    print(f"The most likely electron pair donor atom after the reaction is atom {donor_atom.GetIdx()} "
          f"with {num_hydrogens} available electron pairs.")
else:
    print("No suitable electron pair donor atom found after the reaction.")


