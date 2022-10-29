from pura.units import Mass, Amount, Volume
from rdkit import Chem
from pydantic import BaseModel
from typing import Optional, List, Any
from enum import Enum
import re
import warnings


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
        return (
            other.identifier_type == self.identifier_type and self.value == other.value
        )


class Compound(BaseModel):
    identifiers: List[CompoundIdentifier]
    # amount: Union[Amount, Mass, Volume]
    amount: Amount = None
    mass: Mass = None
    volume: Volume = None
    # source: Source = None
    # data: Dict[str, Data] = None
    # analyses: Dict[str, Analysis] = None


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
        for a in mol.GetAtoms():
            if a.GetIsotope() != 0:
                # SMILES example: '[2H]P([2H])[2H]'
                warnings.warn("Warning: SMILES string contains isotopes.")
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache(strict=False)
        if mol is not None:
            identifier.value = Chem.MolToSmiles(mol)
