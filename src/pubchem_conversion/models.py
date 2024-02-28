from pydantic import BaseModel
from typing import List
from pydantic.dataclasses import dataclass

class ChemicalAtoms(BaseModel):
    proton_donor_atom_index:int
    proton_acceptor_atom_index:int
    electon_pair_donor_atom_index:int
    electron_pair_acceptor_atom_index:int
    
class Chemical(BaseModel):
    CID:str
    MolecularFormula:str
    MolecularWeight:int
    CanonicalSMILES:str
    IsomericSMILES:str
    InChI:str
    InChIKey:str
    IUPACName:str
    XLogP:int
    ExactMass:int
    MonoisotopicMass:int
    TPSA:int
    Complexity:int
    Charge:int
    HBondDonorCount:int
    HBondAcceptorCount:int
    RotatableBondCount:int
    HeavyAtomCount:int
    IsotopeAtomCount:int
    AtomStereoCount:int
    DefinedAtomStereoCount:int
    UndefinedAtomStereoCount:int
    BondStereoCount:int
    DefinedBondStereoCount:int
    UndefinedBondStereoCount:int
    CovalentUnitCount:int
    Volume3D:int
    XStericQuadrupole3D:int
    YStericQuadrupole3D:int
    ZStericQuadrupole3D:int
    FeatureCount3D:int
    FeatureAcceptorCount3D:int
    FeatureDonorCount3D:int
    FeatureAnionCount3D:int
    FeatureCationCount3D:int
    FeatureRingCount3D:int
    FeatureHydrophobeCount3D:int
    ConformerModelRMSD3D:int
    EffectiveRotorCount3D:int
    ConformerCount3D:int
    Fingerprint2D:str
    electronegativity:int
    atoms: ChemicalAtoms

