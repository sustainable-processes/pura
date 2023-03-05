from pura.units import *
from rdkit import Chem
from pydantic import BaseModel
from typing import Optional, List, Any, Union
from enum import Enum
import re
import warnings


class TextureType(str, Enum):
    UNSPECIFIED = "UNSPECIFIED"
    CUSTOM = "CUSTOM"
    POWDER = "POWDER"
    CRYSTAL = "CRYSTAL"
    OIL = "OIL"
    AMORPHOUS_SOLID = "AMORPHOUS_SOLID"
    FOAM = "FOAM"
    WAX = "WAX"
    SEMI_SOLID = "SEMI_SOLID"
    SOLID = "SOLID"
    LIQUID = "LIQUID"


class Texture(BaseModel):
    texture_type: TextureType
    details: Optional[str] = None


class CompoundIdentifierType(str, Enum):
    UNSPECIFIED = "UNSPECIFIED"
    CUSTOM = "CUSTOM"
    # Simplified molecular-input line-entry system
    SMILES = "SMILES"
    # IUPAC International Chemical Identifier
    INCHI = "INCHI"
    # Molblock from a MDL Molfile V3000
    MOLBLOCK = "MOLBLOCK"
    # Chemical name following IUPAC nomenclature recommendations
    IUPAC_NAME = "IUPAC_NAME"
    # Any accepted common name, trade name, etc.
    NAME = "NAME"
    # Chemical Abstracts Service Registry Number (with hyphens)
    CAS_NUMBER = "CAS_NUMBER"
    # PubChem Compound ID number
    PUBCHEM_CID = "PUBCHEM_CID"
    # ChemSpider ID number
    CHEMSPIDER_ID = "CHEMSPIDER_ID"
    # ChemAxon extended SMILES
    CXSMILES = "CXSMILES"
    # IUPAC International Chemical Identifier key
    INCHI_KEY = "INCHI_KEY"
    # XYZ molecule file
    XYZ = "XYZ"
    # UniProt ID (for enzymes)
    UNIPROT_ID = "UNIPROT_ID"
    # Protein data bank ID (for enzymes)
    PDB_ID = "PDB_ID"
    # Amino acid sequence (for enzymes).
    AMINO_ACID_SEQUENCE = "AMINO_ACID_SEQUENCE"
    # HELM https:#www.pistoiaalliance.org/helm-notation/.
    HELM = "HELM"
    # SMILES arbitrary target specification
    SMARTS = "SMARTS"


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
        return (
            other.identifier_type == self.identifier_type and self.value == other.value
        )


pura_to_rdkit_converters = {
    CompoundIdentifierType.SMILES: Chem.MolFromSmiles,
    CompoundIdentifierType.INCHI: Chem.MolFromInchi,
}

rdkit_to_pura_converters = {
    CompoundIdentifierType.SMILES: Chem.MolToSmiles,
    CompoundIdentifierType.INCHI: Chem.MolToInchi,
}


class Compound(PintModel):
    identifiers: List[CompoundIdentifier]
    quantity: Optional[Union[Amount, Mass, Volume]] = None

    @classmethod
    def from_rdkit_mol(
        cls,
        mol: Chem.Mol,
        identifier_types: List[CompoundIdentifierType] = [
            CompoundIdentifierType.SMILES,
            CompoundIdentifierType.INCHI,
        ],
    ):
        """Construct a Compound from an RDKit Mol object."""
        identifiers = []
        for identifier_type in identifier_types:
            if identifier_type in pura_to_rdkit_converters:
                converter = rdkit_to_pura_converters[identifier_type]
                value = converter(mol)
                identifiers.append(
                    CompoundIdentifier(identifier_type=identifier_type, value=value)
                )
        return cls(identifiers=identifiers)

    def to_rdkit_mol(self) -> Chem.Mol:
        """Convert a Compound to an RDKit Mol object."""
        for identifier_type, converter in pura_to_rdkit_converters.items():
            for identifier in self.identifiers:
                if identifier.identifier_type == identifier_type:
                    return converter(identifier.value)
        raise ValueError("No valid identifier found.")

    @classmethod
    def from_smiles(cls, smiles: str):
        mol = Chem.MolFromSmiles(smiles)
        return cls.from_rdkit_mol(mol)

    def to_smiles(self) -> str:
        for identifier in self.identifiers:
            if identifier.identifier_type == CompoundIdentifierType.SMILES:
                return identifier.value
        mol = self.to_rdkit_mol()
        return Chem.MolToSmiles(mol)


def unique_identifiers(identifiers: List[CompoundIdentifier]):
    unique_identifiers = []
    for identifier in identifiers:
        if identifier not in unique_identifiers:
            unique_identifiers.append(identifier)
    return unique_identifiers


def standardize_identifier(identifier: CompoundIdentifier):
    if identifier.identifier_type == CompoundIdentifierType.SMILES:
        smi = identifier.value
        # check smi and raise warnings
        if "*" in smi:
            # SMILES example: '*', '*CCC'
            warnings.warn("Warning: * in SMILES string.")
        if "." in smi:
            smi_segs = smi.split(".")
            smi_segs = [re.sub("-\w", "", _) for _ in smi_segs]
            if all([(("+" not in _) & ("-" not in _)) for _ in smi_segs]):
                # SMILES example: '[HH].[HH].[HH].[HH].[HH].[Ir].[MgH2].[MgH2]'
                warnings.warn(
                    "Warning: SMILES of a mixture, rather than a pure compound, was found."
                )
        if ("+" in smi) | ("-" in smi):
            # calculate charge of the compound based on its SMILES
            p = re.findall("\+(.*?)\]", smi)
            n = re.findall("\-(.*?)\]", smi)
            p_total = []
            n_total = []
            for _ in p:
                if str.isdigit(_):
                    p_total.append(int(_))
                elif _ == "":
                    p_total.append(1)
                else:
                    continue
            for _ in n:
                if str.isdigit(_):
                    n_total.append(int(_))
                elif _ == "":
                    n_total.append(1)
                else:
                    continue
            p_total = sum(p_total)
            n_total = -sum(n_total)
            comp_charge = p_total + n_total
            if comp_charge != 0:
                # SMILES example: 'O=C([O-])[O-].[K+]', '[NH4+]', '[H]O[H].[O-][I+2]([O-])[O-]'
                warnings.warn(
                    "Warning: Compound is not electrically neutral based on its SMILES string."
                )
                # NOTE: did not add the function to neutralize charged SMILES, the balance function could do sth like this:
                # ['[Cs+].[H-].[PdH4-2]', 'O=[V+]([O-])O.[NH4+]', '[H]C(=O)[O-].[H]O[H].[Rh+2]', '[O-][I+2]([O-])[O-]'] ->
                # ['[Cs+].[Cs+].[Cs+].[H-].[PdH4-2]', 'O=[V+]([O-])O.[NH4+].[OH-]', '[H]C(=O)[O-].[H]C(=O)[O-].[H]O[H].[Rh+2]', '[O-][I+2]([O-])[O-].[H+]']
                # if this feature is necessary we can add later
        mol = Chem.MolFromSmiles(smi)
        if mol is not None:
            for a in mol.GetAtoms():
                if a.GetIsotope() != 0:
                    # SMILES example: '[2H]P([2H])[2H]'
                    warnings.warn("Warning: SMILES string contains isotopes.")
            Chem.SanitizeMol(mol)
            mol.UpdatePropertyCache(strict=False)
            identifier.value = Chem.MolToSmiles(mol)
