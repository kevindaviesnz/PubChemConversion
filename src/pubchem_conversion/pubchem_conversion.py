from models import Chemical, ChemicalAtoms

# $ pip install rdkit
# https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem, AllChem, Draw

def determine_proton_donor_atom_index():
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
    

def determine_proton_acceptor_atom_index():
    
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


def determine_electron_donor_acceptor_atom_index():
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



def determine_electron_pair_acceptor_atom_index():
    
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


def to_object(chemical_from_pubchem:dict):

    # Create mol from SMILES
    m = Chem.MolFromSmiles(chemical_from_pubchem["PropertyTable"]["Properties"][0]["CanonicalSMILES"], )
    # img = Draw.MolToImage(m)

    electronegativity = 0
    
    chemical_atoms = ChemicalAtoms(
        proton_donor_atom_index = 0,
        proton_acceptor_atom_index = 0,
        electon_pair_donor_atom_index = 0,
        electron_pair_acceptor_atom_index = 0
    )
    
    chemical = Chemical(
        **chemical_from_pubchem["PropertyTable"]["Properties"][0], 
        electronegativity=electronegativity,
        chemical_atoms=chemical_atoms
    )
    return chemical
