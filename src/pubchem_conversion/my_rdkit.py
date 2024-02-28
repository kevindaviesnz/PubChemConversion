from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw

# https://www.rdkit.org/docs/GettingStartedInPython.html

# Create mol from SMILES
m = Chem.MolFromSmiles('Cc1ccccc1')
print(m)

# Draw molecule
#img = Draw.MolToImage(m)
#print(img)

# Mol to smiles
print(Chem.MolToSmiles(m))

# Mol block
print(Chem.MolToMolBlock(m))

# Looping over atoms
for atom in m.GetAtoms():
    print(atom.GetAtomicNum())
    
print(m.GetBonds()[0].GetBondType())

# Atom
print(m.GetAtomWithIdx(0).GetSymbol())
print(m.GetAtomWithIdx(0).GetExplicitValence())

# Bonds
print(m.GetBondWithIdx(0).GetBeginAtomIdx())
print(m.GetBondWithIdx(0).GetEndAtomIdx())
print(m.GetBondBetweenAtoms(0,1).GetBondType())

# Remove substructure
m = Chem.MolFromSmiles('CC(=O)O')
patt = Chem.MolFromSmarts('C(=O)[OH]')
rm = AllChem.DeleteSubstructs(m,patt)
print(Chem.MolToSmiles(rm))

# Chemical reactions
rxn = AllChem.ReactionFromSmarts('[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]')
rxn.GetNumProductTemplates()
ps = rxn.RunReactants((Chem.MolFromSmiles('CC(=O)O'),Chem.MolFromSmiles('NC')))
print(Chem.MolToSmiles(ps[0][0])) # CNC(C)=O

rxn = AllChem.ReactionFromSmarts('[C:1]=[C:2].[C:3]=[*:4][*:5]=[C:6]>>[C:1]1[C:2][C:3][*:4]=[*:5][C:6]1')
ps = rxn.RunReactants((Chem.MolFromSmiles('OC=C'), Chem.MolFromSmiles('C=CC(N)=C')))
print(Chem.MolToSmiles(ps[0][0]))

"""
from rdkit import Chem
from rdkit.Chem import AllChem

def get_acidic_and_basic_sites(molecule):
    acidic_site = AllChem.GetBestRDKFunc(molecule, usePh=True, isAcid=True)
    basic_site = AllChem.GetBestRDKFunc(molecule, usePh=True, isAcid=False)
    return acidic_site, basic_site

# Example: Ethanol (CCO)
molecule_smiles = 'CCO'
mol = Chem.MolFromSmiles(molecule_smiles)

acidic_site, basic_site = get_acidic_and_basic_sites(mol)
print(f'Most acidic site: {acidic_site}')
print(f'Most basic site: {basic_site}')

"""

"""
HCl + water
"""

# Create HCl molecule using SMILES
hcl_smiles = 'Cl[H]'
hcl_molecule = Chem.MolFromSmiles(hcl_smiles)

# Create water molecule using SMILES
water_smiles = 'O'
water_molecule = Chem.MolFromSmiles(water_smiles)

# Find in HCl molecule the atom that will donate a proton (acid atom)
# proton_donor_atom_idx = [atom.GetIdx() for atom in hcl_molecule.GetAtoms() if atom.GetAtomicNum() == 1][0]
proton_donor_atom_idx = AllChem.GetBestRDKFunc(hcl_molecule, usePh=True, isAcid=True)

# Find in water molecule atom that will accept the proton (base atom)
proton_acceptor_atom_idx = [atom.GetIdx() for atom in water_molecule.GetAtoms() if atom.GetAtomicNum() == 8][0]

# Remove proton from HCl molecule
hcl_molecule = Chem.DeleteSubstructs(hcl_molecule, Chem.MolFromSmiles('[H]'))

# Add proton to water molecule
water_molecule = Chem.CombineMols(water_molecule, Chem.MolFromSmiles('[H]'))

# Render SMILES of water molecule after protonation
protonated_water_smiles = Chem.MolToSmiles(water_molecule)
print(f"Protonated water SMILES: {protonated_water_smiles}")

# Render SMILES of HCl molecule after deprotonation
deprotonated_hcl_smiles = Chem.MolToSmiles(hcl_molecule)
print(f"Deprotonated HCl SMILES: {deprotonated_hcl_smiles}")

# Optional: Draw the structures using RDKit's drawing functions
Draw.MolToImage(hcl_molecule).show()
Draw.MolToImage(water_molecule).show()
