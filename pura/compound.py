from pura.units import Mass, Amount, Volume
from rdkit import Chem
from pydantic import BaseModel
from typing import Optional, List, Any
from enum import Enum


class CompoundIdentifierType(Enum):
    UNSPECIFIED = 0
    CUSTOM = 1
    # Simplified molecular-input line-entry system
    SMILES = 2
    # IUPAC International Chemical Identifier
    INCHI = 3
    # Molblock from a MDL Molfile V3000
    MOLBLOCK = 4
    # Chemical name following IUPAC nomenclature recommendations
    IUPAC_NAME = 5
    # Any accepted common name, trade name, etc.
    NAME = 6
    # Chemical Abstracts Service Registry Number (with hyphens)
    CAS_NUMBER = 7
    # PubChem Compound ID number
    PUBCHEM_CID = 8
    # ChemSpider ID number
    CHEMSPIDER_ID = 9
    # ChemAxon extended SMILES
    CXSMILES = 10
    # IUPAC International Chemical Identifier key
    INCHI_KEY = 11
    # XYZ molecule file
    XYZ = 12
    # UniProt ID (for enzymes)
    UNIPROT_ID = 13
    # Protein data bank ID (for enzymes)
    PDB_ID = 14
    # Amino acid sequence (for enzymes).
    AMINO_ACID_SEQUENCE = 15
    # HELM; https:#www.pistoiaalliance.org/helm-notation/.
    HELM = 16
    # SMILES arbitrary target specification
    SMARTS = 17


class Data:
    pass


class Analysis:
    pass


class Source:
    pass


class CompoundIdentifier(BaseModel):
    identifier_type: CompoundIdentifierType
    value: str
    details: Optional[str] = None

    def __eq__(self, other: Any) -> bool:
        if not other.identifier_type == self.identifier_type:
            raise TypeError(
                f"Not the same identifier type ({other.identifier_type} != {self.identifier_type} is not"
            )
        return self.value == other.value


class Compound(BaseModel):
    identifiers: List[CompoundIdentifier]
    # amount: Union[Amount, Mass, Volume]
    amount: Amount = None
    mass: Mass = None
    volume: Volume = None
    # source: Source = None
    # data: Dict[str, Data] = None
    # analyses: Dict[str, Analysis] = None


def standardize_identifier(identifier: CompoundIdentifier):
    if identifier.identifier_type == CompoundIdentifierType.SMILES:
        mol = Chem.MolFromSmiles(identifier.value)
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache(strict=False)
        if mol is not None:
            identifier.value = Chem.MolToSmiles(mol)
