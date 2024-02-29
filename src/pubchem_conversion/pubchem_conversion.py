from models import Chemical, ChemicalAtoms

# $ pip install rdkit
# https://www.rdkit.org/docs/GettingStartedInPython.html
from rdkit import Chem, AllChem, Draw




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
