from typing import List,Tuple,Dict
from .reaction import ReactionIdentifier, ReactionIdentifierType
from rdkit import Chem
from rdkit.Chem import PeriodicTable, GetPeriodicTable
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import ray

class CustomError(Exception):
    '''
    For error handling *Change/adapt as needed*
    '''
    pass


def initray(restart=True,num_cpus=16,log_to_driver=False):
    """
    Initializes a Ray CPU cluster for parallel execution of functions

    Args:
        restart (bool, optional): Restart existing cluster. Defaults to True.
        num_cpus (int, optional): Number of CPUs for parallel execution. Defaults to 16.
        log_to_driver (bool, optional): Controls whether results are logged to driver. Defaults to False.
    """
    if restart:
        ray.shutdown()
    ray.init(num_cpus=num_cpus,log_to_driver=log_to_driver)
    
def molfromsmiles(SMILES):
    """
    Converts a smiles string into a mol object for RDKit use. Returns a mol object. *Unsure how to integrate in compound.py*

    Args:
        SMILES (str): SMILES string to convert to mol object

    Returns:
        mol: RDKit mol object
    """
    mol=Chem.MolFromSmiles(SMILES)
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)
    return mol
    
def atomtypes(mol)->Tuple[Dict,int]:
    """
    Generates an atom type dictionary with counts for a given mol file. Also
    returns overall change of mol. 

    Args:
        mol (RDKit mol): RDKit mol object to generate atom type dictionary for

    Returns:
        Tuple[Dict,int]: Dictionary with keys as atom element and value as count, overall charge
    """

    typedict={}
    mol2=Chem.AddHs(mol) #Assumes hydrogens are added
    charge=0
    for atom in mol2.GetAtoms():
        elem=PeriodicTable.GetElementSymbol(GetPeriodicTable(),atom.GetAtomicNum())
        charge+=atom.GetFormalCharge()
        if elem not in typedict.keys():
            typedict[elem]=1
        else:
            typedict[elem]+=1
    return typedict, charge



def getcompdict(ID=1,mol=None,smiles=None,formula=None)->Dict:
    """
    Wrapper for atomtypes. ID must be input whether Reaxys ID or random number. Default is 1. 
    Smiles and mol can also be specified. Can be converted to an object instance if intended.

    Args:
        ID (int, optional): ID of compound/mol. Defaults to 1.
        mol (RDKit mol, optional): RDKit mol. Defaults to None.
        smiles (str, optional): SMILES string. Defaults to None.
        formula (str, optional): Molecular formula. Defaults to None.

    Raises:
        CustomError: Raises error if smiles specified
        is not valid/cannot be processed

    Returns:
        Dict: Compound dictionary with ID as key and another dictionary with the
        following keys: 'atomdict' (Output of atomtypes, dictionary with atom
        elements and count), 'charge' (Overall charge of mol),'smiles'
        (smiles of mol),'formula' (chemical formula),'count' (instance number, always 1)
    """

    if smiles is not None:
        try:
            mol=molfromsmiles(smiles)
        except Exception as e:
            raise CustomError("Compound "+str(ID)+" smiles is invalid")
    if formula is None:
        formula=CalcMolFormula(mol)
    atomdata=atomtypes(mol)
    compddict={ID:{'atomdict':atomdata[0],'charge':atomdata[1],'smiles':smiles,'formula':formula,'count':1}}
    return compddict



def balance_reactions(
    reactions: List[ReactionIdentifier],
    ignore_hydrogens: bool = True,
    use_ray: bool = False,
    
) -> List[Tuple[ReactionIdentifier, ReactionIdentifier]]:
    if use_ray:

        @ray.remote
        def balance_reaction_ray(**kwargs):
            balance_reaction(**kwargs)

        for reaction in reactions:
            ray.remote(balance_reaction_ray)
    else:
        [balance_reaction(reaction) for reaction in reactions]


def balance_reaction(
    input_reaction_identifier: ReactionIdentifier,
) -> Tuple[ReactionIdentifier, ReactionIdentifier]:
    pass

    # do balancing
    output_reaction_identifier = None
    return input_reaction_identifier, output_reaction_identifier
