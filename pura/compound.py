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
        if mol is not None:
            identifier.value = Chem.MolToSmiles(mol)
            
            
### tidy up SMILES ###
class Erratic_smiles_catcher:
    """
    Background: in Reaxys, some SMILES are problematic or not useful, this code help to identify these SMILES
    dependencies: re, rdkit
    """
    def __init__(self, smis):
        self.smiles = smis
    def get_erratic_smiles_str(self, err_str):
        # identify smiles containing string strs
        # err_str: str, example of erratic smiles strings: '*'
        res = []
        idxs = []
        for idx, smi in enumerate(self.smiles):
            if err_str in smi:
                res.append(smi)
                idxs.append(idx)
        return res, idxs
    def get_mix_smiles(self):
        # some SMILES of mixtures, e.g., '[HH].[HH].[HH].[HH].[HH].[Ir].[MgH2].[MgH2]', 'O.O.[Pt+2]'
        # this function identify these mixtures but avoid to include salts, e.g., not salt 'O=C([O-])[O-].[K+]'
        res = []
        idxs = []
        for idx, smi in enumerate(self.smiles):
            if '.' in smi:
                smi_segs = smi.split('.')
                if any([(('+' not in _)&('-' not in _)) for _ in smi_segs]):
                    res.append(smi)
                    idxs.append(idx)
        return {'smis': res, 'idxs': idxs}
    def get_isotope_smiles(self):
        # identify SMILES containing isotopes, e.g., '[2H]P([2H])[2H]'
        res = []
        idxs = []
        for idx, smi in enumerate(self.smiles):
            if self.is_isotope_exist(smi):
                res.append(smi)
                idxs.append(idx)
        return {'smis': res, 'idxs': idxs}
    def get_smiles_by_heavyatomcount(self, heavy_atom_range):
        # small or big molecules could be problematic, e.g., '17718916': '[H][H-]', 4122981:'[CH-].[H+]'
        # heavy_atom_range: [int_lowerbound, int_upperbond], <=upper and >=lower bounds of heavyatomcount
        # dependency: self.get_heavyatomcount 
        idxs = []
        smis = []
        heavy_atom_counts = []
        for idx, smi in enumerate(self.smiles):
            heavy_atom_count = self.get_heavyatomcount(smi)
            if (heavy_atom_count>=heavy_atom_range[0])&(heavy_atom_count<=heavy_atom_range[1]):
                smis.append(smi)
                idxs.append(idx)
                heavy_atom_counts.append(heavy_atom_count)
        return {'idxs':idxs, 'smis':smis, 'n_heavy_atoms':heavy_atom_counts}
    def get_charged_smiles(self):
        # many unneutral salts, even ions, e.g.,: 'O=C([O-])[O-].[K+]', '[NH4+]'
        res = []
        idxs = []
        for idx, smi in enumerate(self.smiles):
            smi_charge = self.get_smi_charge(smi)
            if smi_charge != 0:
                res.append((smi, smi_charge))
                idxs.append(idx)
        return {'smis': res, 'idxs': idxs}
    def get_duplicate_smiles(self):
        # some compound share the same SMILES but different IDs, though they could be different compounds!
        smis_idxs = {} # {smi: [idx, ..], ..}
        smis_dup = []
        for i,smi in enumerate(self.smiles):
            if smi not in smis_idxs:
                smis_idxs[smi] = [i]
            else:
                smis_idxs[smi].append(i)
                smis_dup.append(smi)
        res = {smi_dup:smis_idxs[smi_dup] for smi_dup in smis_dup}
        return res
    def get_heavyatomcount(self, smi):
        mol = Chem.MolFromSmiles(smi)
        n_heavy = mol.GetNumHeavyAtoms()
        return n_heavy
    def get_smi_charge(self, smi):
        # smi: list of SMILES
        # return: e.g., ['CCC', 'O=[V+]([O-])O', '[NH4+]', '[H]O[H].[O-][I+2]([O-])[O-]'] -> [0, 0, 1, -1]
        if ('+' not in smi)&('-' not in smi):
            res = 0
        else:
            p = re.findall('\+(.*?)\]', smi)
            n = re.findall('\-(.*?)\]', smi)
            p_total = sum([int(_) if _ != '' else 1 for _ in p])
            n_total = -sum([int(_) if _ != '' else 1 for _ in n])
            res = p_total + n_total
        return res
    def is_isotope_exist(self, smi):
        res = False
        mol = Chem.MolFromSmiles(smi)
        for a in mol.GetAtoms():
            if a.GetIsotope()!=0:
                res = True
                break
        return res

class Smiles_healer:
    """
    Background: after identify erratic SMILES, probaly we want to fix them in order to perserve as much info as possible
    dependencies: rdkit, re, from pulp import *
    """
    def __init__(self, smis):
        self.smiles = smis
    def heal_err_smiles(self):
        """
        this is just a proposed subproject, here only nuetralize SMILES for salts and ions was done, see nuetralize_salts
        """
        pass
    def nuetralize_salts(self):
        # return: e.g.,
        #['[Cs+].[H-].[PdH4-2]', 'O=[V+]([O-])O.[NH4+]', '[H]C(=O)[O-].[H]O[H].[Rh+2]', '[O-][I+2]([O-])[O-]'] -> 
        # ['[Cs+].[Cs+].[Cs+].[H-].[PdH4-2]', 'O=[V+]([O-])O.[NH4+].[OH-]', '[H]C(=O)[O-].[H]C(=O)[O-].[H]O[H].[Rh+2]', '[O-][I+2]([O-])[O-].[H+]']
        # dependencies: self.get_smi_charge, self.get_ortho_vector
        len_smis = len(self.smiles)
        smis_nuetral = ['']*len_smis
        for i in range(len_smis):
            smi = self.smiles[i]
            smi_charge = self.get_smi_charge(smi)
            if smi_charge==0:
                smis_nuetral[i] = smi
            else:
                frags = smi.split('.')
                frags_charges = [self.get_smi_charge(_) for _ in frags ]
                if all([_>=0 for _ in frags_charges]):
                    frags.append('[OH-]')
                    frags_charges.append(-1)
                elif all([_<=0 for _ in frags_charges]):
                    frags.append('[H+]')
                    frags_charges.append(1)
                frags_charged = [(frags[idx], _ ) for idx,_ in enumerate(frags_charges) if _!=0]
                coeffs = self.get_ortho_vector([_[1] for _ in frags_charged]) 
                frags_charged_coeffs = dict(zip([_[0] for _ in frags_charged], coeffs))
                smi_new = ''
                for frag in frags:
                    if frag in frags_charged_coeffs:
                        for j in range(frags_charged_coeffs[frag]):
                            smi_new = smi_new + frag + '.'
                    else:
                        smi_new = smi_new + frag + '.'
                smis_nuetral[i] = smi_new[:-1]
        return smis_nuetral
    def get_smi_charge(self, smi):
        if ('+' not in smi)&('-' not in smi):
            res = 0
        else:
            # example resultes from re.findall: ['', '2', '']
            # '+]' and '-]' -> '' -> 1+, and 1-
            p = re.findall('\+(.*?)\]', smi)
            n = re.findall('\-(.*?)\]', smi)
            p_total = sum([int(_) if (_ != '')&str.isdigit(_) else 1 for _ in p])
            n_total = -sum([int(_) if (_ != '')&str.isdigit(_) else 1 for _ in n])
            res = p_total + n_total
        return res
    def get_ortho_vector(self, v):
        # v: vetcor of ints, should not contain 0
        # return: e.g., [1,2,-2] -> [2,1,2]
        # dependencies: pulp
        variables = []
        for i in range(len(v)):
            variables.append(LpVariable('x'+str(i), cat='Integer', lowBound=1))
        prob = LpProblem("Neutralization", LpMinimize)
        prob += 0,"Objective Function"
        constraint = None
        for i in range(len(v)):
            constraint += variables[i]*v[i]
        prob += constraint == 0,f'Constraint_for_a' # the name canbe any words: Constraint_for_b or _c ...
        # to mute log: prob.solve(pulp.PULP_CBC_CMD(msg=False))
        prob.solve()
        coeffs = [int(_.value()) for _ in variables]
        return coeffs
