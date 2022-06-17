from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import AllChem
from functools import partial
import types
import psutil
from psutil._common import bytes2human
from typing import List, Tuple, Union, Optional, Dict
import pandas as pd
from tqdm import tqdm
import numpy as np
from pathlib import Path
from sklearn.model_selection import GroupShuffleSplit
import logging
from collections import namedtuple
import json
from IPython.display import clear_output, Image
from pint import UnitRegistry



__all__ = [
    "get_mem",
    "remove_frame",
    "RecursiveNamespace",
    "canonicalize_smi",
    "clean_rxn_smiles",
    "SmilesCleaner",
    "get_pistachio_dataset",
    "ureg"
]

# Unit registry
ureg = UnitRegistry()
ureg.define('micro- = 1e-6 = μ-')

def get_mem():
    """Get memory available in a human readable format"""
    mem_usage = psutil.virtual_memory().available
    return bytes2human(mem_usage)


def remove_frame(ax, sides=["top", "left", "right"]):
    """Remove the frame of a matplotlib plot"""
    for side in sides:
        ax_side = ax.spines[side]
        ax_side.set_visible(False)

class RecursiveNamespace(types.SimpleNamespace):
    # def __init__(self, /, **kwargs):  # better, but Python 3.8+
    def __init__(self, **kwargs):
        """Create a SimpleNamespace recursively"""
        self.__dict__.update({k: self.__elt(v) for k, v in kwargs.items()})

    def __elt(self, elt):
        """Recurse into elt to create leaf namepace objects"""
        if type(elt) is dict:
            return type(self)(**elt)
        if type(elt) in (list, tuple):
            return [self.__elt(i) for i in elt]
        return elt

# ------ SMILES utilities ------------------

class SmilesCleaner:
    def __init__(self, cache: Optional[Dict] = None):
        self._cache = cache if cache is not None else {}
    
    def __call__(self, unclean_smiles: Optional[str] = None) -> str:
        if unclean_smiles is None:
            return 
        smiles = self._cache.get(unclean_smiles)
        if smiles is None:
            try:
                smiles = canonicalize_smi(unclean_smiles)
            except: 
                return
            self._cache[unclean_smiles] = smiles
        return smiles

### Code below from rxnfp by Philippe Schwaller
def canonicalize_smi(smi, remove_atom_mapping=False):
    """
    Canonicalize SMILES
    """
    mol = Chem.MolFromSmiles(smi)
    if not mol:
        raise ValueError("Molecule not canonicalizable")
    if remove_atom_mapping:
        for atom in mol.GetAtoms():
            if atom.HasProp("molAtomMapNumber"):
                atom.ClearProp("molAtomMapNumber")
    return Chem.MolToSmiles(mol)


def clean_rxn_smiles(rxn_smiles):
    rxn_components = rxn_smiles.split(">")
    reactants = [rxn_component.split(".") for rxn_component in rxn_components[:2]]
    try:
        reactants = [canonicalize_smi(smi, remove_atom_mapping=True) + "."
                     for rxn_component in reactants
                     for smi in rxn_component]
    except ValueError:
        return ""
    reactants = "".join(reactants).rstrip(".")
    products = rxn_components[2].split(" ")[0].split(".")
    try:
        products = [canonicalize_smi(smi, remove_atom_mapping=True) + "." for smi in products]
    except ValueError:
        return ""
    products = "".join(products).rstrip(".")
    return f"{reactants}>>{products}"

# Data preprocessing ------------------------------------

class Component:
    """ Representation of a reaction component"""
    def __init__(self, 
        role: str,
        smiles: str,
        amount: Union[str, ureg.Quantity] = None,
        mass: Union[str, ureg.Quantity] = None,
        molecular_weight: Union[str, ureg.Quantity] = None,
        volume: Union[str, ureg.Quantity] = None,
        concentration: Union[str, ureg.Quantity] = None 
    ):
        self._role = role
        self._smiles = smiles
        self._mass = mass
        self._molecular_weight = molecular_weight
        self._amount  = amount
        self._volume = volume
        self._concentration = concentration

    @property
    def amount(self) -> ureg.Quantity:
        # If nothing, return nothing
        if self._amount is None and self._mass is None and  self._molecular_weight is None:
            return
        # Return directly if possible
        elif type(self._amount) is ureg.Quantity:
            return self._amount
        elif type(self._amount) is str:
            return ureg(self._amount)
        # If can, calculate, do so
        elif self._amount is None and self._mass is not None and self._molecular_weight is not None:
            return self._mass / self._molecular_weight


    @property
    def mass(self) -> ureg.Quantity:
        # If nothing, return nothing
        if self._mass is None and self._amount is None and  self._molecular_weight is None:
            return
        # Return directly if possible
        elif type(self._mass) is ureg.Quantity:
            return self._mass
        elif type(self._mass) is str:
            return ureg(self._mass)
        # If can, calculate, do so
        elif self._mass is None and self._amount is not None and  self._molecular_weight is not None:
            return self._amount * self._molecular_weight

    @property
    def volume(self) -> ureg.Quantity:
        if self._volume is None:
            return
        elif type(self._volume) == ureg.Quantity:
            return self._volume
        else:
            return ureg(self._volume)

    @property
    def concentration(self):
        if self._concentration is None:
            return
        elif type(self._concentration) == ureg.Quantity:
            return self._concentration
        else:
            return ureg(self._concentration)

    @property
    def molecular_weight(self):
        if self._molecular_weight is None:
            return
        elif type(self._molecular_weight) == ureg.Quantity:
            return self._molecular_weight
        else:
            return ureg(self._molecular_weight)

    @property
    def smiles(self):
        return self._smiles

    @property
    def role(self):
        return self._role
    
      
class ReactionDatapoint:
  def __init__(
  	self,
    reaction_id: int,
    components: List[Component],
    reaction_yield: float,
    description: Optional[str]=None,
  ):
     self.reaction_id = reaction_id
     self.components = components
     self.reaction_yield = reaction_yield
     self.description = description


class SuzukiReactionDatapoint(ReactionDatapoint):
    def __init__(
        self,
        reaction_id: int,
        components: List[Component],
        reaction_yield: float,
        description: str,
        catalyst_loading: Union[str, ureg.Quantity] = None,
        base_loading: Union[str, ureg.Quantity] = None,
        reaction_scale: Union[str, ureg.Quantity] = None,
        total_volume: Union[str, ureg.Quantity] = None
    ):
        super().__init__(reaction_id, components, reaction_yield, description)
        self._catalyst_loading = catalyst_loading
        self._base_loading = base_loading
        self._reaction_scale = reaction_scale
        self._total_volume = total_volume
    
    @property
    def catalyst_loading(self):
        if self._catalyst_loading is None:
            return
        elif type(self._catalyst_loading) is ureg.Quantity:
            return self._catalyst_loading
        elif type(self._catalyst_loading) is str:
            return ureg(self._catalyst_loading)

    @property
    def base_loading(self):
        if self._base_loading is None:
            return
        elif type(self._base_loading) is ureg.Quantity:
            return self._base_loading
        elif type(self._base_loading) is str:
            return ureg(self._base_loading)
    
    @property
    def reaction_scale(self):
        if self._reaction_scale is None:
            return
        elif type(self._reaction_scale) is ureg.Quantity:
            return self._reaction_scale
        elif type(self._reaction_scale) is str:
            return ureg(self._reaction_scale)

    @property
    def total_volume(self):
        if self._total_volume is None:
            return
        elif type(self._total_volume) is ureg.Quantity:
            return self._total_volume
        elif type(self._total_volume) is str:
            return ureg(self._total_volume)

    @property
    def reaction_smiles(self):
        reactants = "".join(
            [component.smiles + "." for component in self.components if component.role != SuzukiRole.product]
        ).rstrip(".")
        products = "".join(
            [component.smiles for component in self.components if component.role == SuzukiRole.product]
        ).rstrip(".")
        return f"{reactants}>>{products}"

    def _repr_png_(self):
        reaction = AllChem.ReactionFromSmarts(self.reaction_smiles, useSmiles=True)
        return reaction._repr_png_()


class ReactionDataset:
    def __init__(self, data: List[ReactionDatapoint]):
        self._data = data


    @property
    def data(self):
        return self._data

    def to_frame(self, component_representations: Optional[dict]=None):
        """To dataframe
        Can specify a dictionary with the representations that should be used for the components
        (i.e., SMILES, one-hot encoding, fingerprint)
        """
        raise NotImplementedError()

class SuzukiReactionDataset(ReactionDataset):
    def to_frame(self):
        reaction_dicts = []
        for reaction in self._data:
            reaction_dict = {}
            n_solvents = 1
            for component in reaction.components:
                if component.role == SuzukiRole.boronic_acid_reactant:
                    reaction_dict["reactant_1_smiles"] = component.smiles
                elif component.role == SuzukiRole.organohalide_reactant:
                    reaction_dict["reactant_2_smiles"] = component.smiles
                elif component.role == SuzukiRole.pre_catalyst:
                    pass
                elif component.role == SuzukiRole.ligand:
                    reaction_dict["ligand_smiles"] = component.smiles
                elif component.role == SuzukiRole.base:
                    reaction_dict["base_smiles"] = component.smiles
                elif component.role == SuzukiRole.solvent:
                    reaction_dict[f"solvent_{n_solvents}_smiles"] = component.smiles
                    n_solvents += 1
                elif component.role == SuzukiRole.product:
                    reaction_dict[f"product_smiles"] = component.smiles
            reaction_dict["catalyst_loading"] = reaction.catalyst_loading.magnitude if reaction.catalyst_loading else None
            reaction_dict["base_loading"] = reaction.base_loading.magnitude if reaction.base_loading else None
            reaction_dict["reaction_scale (g)"] = reaction.reaction_scale.to(ureg.gram).magnitude if reaction.reaction_scale else None
            reaction_dict["total_volume (mL)"] = reaction.total_volume.to(ureg.ml).magnitude if reaction.total_volume else None
            reaction_dict["yield"] = reaction.reaction_yield
            reaction_dict["reaction_id"] = reaction.reaction_id
            reaction_dicts.append(reaction_dict)
        return pd.DataFrame(reaction_dicts)
        
    
def read_batch_details(path: Path, batch_num: int):
    with open(path / f"suzuki_details_batch_{batch_num}.json", "r") as f:
        details_data = json.load(f)
    return details_data


def get_pistachio_dataset(
    reactions_df: pd.DataFrame, 
    reaction_details_path: Union[str, Path],
    base_smiles_list: List[str],
    ligand_smiles_list: Optional[List[str]]= None,
    solvent_smiles_list: Optional[List[str]] = None,
    reaction_id_col: str = "reaction_id",
    paragraph_col: str = "paragraphText",
    yield_col: str = "yield",
    smiles_cleaner: Optional[SmilesCleaner] = None,
    n_batches: int = 298,
    start_batch: int = 0
):
    """
    reactions_df: pd.DataFrame
        Dataframe withr eactions
    reaction_details_path 
    """
    # Drop rows in reaction df without / duplicate reaction_id, paragraph or numerical yield
    reactions_df = reactions_df.dropna()
    reactions_df = reactions_df.drop_duplicates(subset=reaction_id_col)
    reactions_df[yield_col] = reactions_df[yield_col].astype(float)
    reactions_df = reactions_df.set_index(reaction_id_col)

    # Create smiles cleaner
    if smiles_cleaner is None:
        smiles_cleaner = SmilesCleaner()

    # Ligand molecules list
    if ligand_smiles_list is not None:
        ligand_mols = {smiles: Chem.MolFromSmiles(smiles) for smiles in ligand_smiles_list}
    else: 
        ligand_mols = None
    
    # Canonicalize smiles lists
    if base_smiles_list is not None:
        base_smiles_list = [smiles_cleaner(smiles) for smiles in base_smiles_list]
    if solvent_smiles_list is not None:
        solvent_smiles_list = [smiles_cleaner(smiles) for smiles in solvent_smiles_list]

    data = []
    reaction_details_path = Path(reaction_details_path)
    for batch in range(start_batch, n_batches):
        # Read in batch of reaction details
        print(f"Batch {batch+1} of {n_batches}")
        new_data = read_batch_details(reaction_details_path, batch)
        # Loop through reactions
        for data_dict in new_data:
            for reaction_id, d in data_dict.items():
                try:
                    reaction_metadata = reactions_df.loc[int(reaction_id)]
                except KeyError:
                    # Skip reactions that were already dropped
                    continue
                components = parse_suzuki_components(
                    d["components"], smiles_cleaner, ligand_mols, base_smiles_list, solvent_smiles_list
                )
                good = check_components(components)
                if not good: 
                    continue
                catalyst_loading, base_loading, reaction_scale, total_volume = calculate_conditions(components)
                data.append(
                    SuzukiReactionDatapoint(
                        reaction_id=reaction_id,
                        components=components,
                        description=d["data"].get("paragraphText"),
                        catalyst_loading=catalyst_loading,
                        base_loading=base_loading,
                        reaction_scale=reaction_scale,
                        total_volume=total_volume,
                        reaction_yield=reaction_metadata[yield_col],
                    )
                )
        clear_output(wait=True)
    return SuzukiReactionDataset(data)

class PistachioRole:
    reactant = "Reactant"
    product = "Product"
    catalyst = "Catalyst"
    solvent = "Solvent"
    agent = "Agent"

class SuzukiRole:
    boronic_acid_reactant = "Boronic Acid Reactant"
    organohalide_reactant = "Organohalide Reactant"
    product = "Product"
    solvent = "Solvent"
    base = "Base"
    ligand = "Ligand"
    pre_catalyst = "Pre-catalyst"
    agent = "Agent"

def check_components(components: List[Component])-> bool:
    """Check for presence of correct components"""
    boronic_acid_present = False
    organohalide_present = False
    base_present = False
    # TODO: Check there is a Pd pre-catalyst
    for component in components:
        if component.role == SuzukiRole.boronic_acid_reactant:
            boronic_acid_present = True
        elif component.role == SuzukiRole.organohalide_reactant:
            organohalide_present = True
        elif component.role == SuzukiRole.base:
            base_present = True
    return boronic_acid_present and organohalide_present and base_present

def calculate_conditions(components: List[Component]) -> (ureg.Quantity, ureg.Quantity, ureg.Quantity):
    """Calculate continuous conditions"""
    # Catalyst loading, base loading and reaction scale
    ligand_amount = product_mass = base_amount = None
    limiting_reactant_amount = 1e6 * ureg.mol
    reaction_scale = 0.0 * ureg.gram
    for component in components:
        if component.role == SuzukiRole.boronic_acid_reactant or component.role == SuzukiRole.organohalide_reactant:
            amount = component.amount
            if amount is None:
                continue
            elif amount < limiting_reactant_amount:
                limiting_reactant_amount = amount
            if component.mass is None:
                continue
            reaction_scale += component.mass
        elif component.role == SuzukiRole.base:
            amount = component.amount
            if amount is None:
                continue
            base_amount = amount
        elif component.role == SuzukiRole.ligand:
            ligand_amount = component.amount

    if ligand_amount is not None and limiting_reactant_amount is not None:
        catalyst_loading = ligand_amount / limiting_reactant_amount
    else:
        catalyst_loading = None
    if base_amount is not None and limiting_reactant_amount is not None:
        base_loading = base_amount / limiting_reactant_amount
    else:
        base_loading = None

    # Total volume (reaction scale)
    total_volume = sum(
        [component.volume for component in components if component.volume is not None]
    )

    return catalyst_loading, base_amount, reaction_scale, total_volume


def get_quantities(quantity_dicts: List[dict])-> dict:
    new_quantities = {}
    for quantity_dict in quantity_dicts:
        if not quantity_dict.get("value"):
            continue
        if quantity_dict["type"] == "Amount":
            new_quantities["amount"] = quantity_dict["value"] * ureg.mole
        elif  quantity_dict["type"] == "Volume":
            new_quantities["volume"] = quantity_dict["value"] * ureg.liter
        elif quantity_dict["type"] == "Molarity":
            new_quantities["concentration"] = quantity_dict["value"] * ureg.mol / ureg.liter
        elif quantity_dict["type"] == "Mass":
            new_quantities["mass"] = quantity_dict["value"] * ureg.gram
    return new_quantities  


def parse_suzuki_components(
    components_list: List[Dict],
    smiles_cleaner: SmilesCleaner,
    base_smiles: List[str],
    ligand_mols: Optional[Dict[str, Chem.Mol]] = None,
    solvent_smiles: Optional[List[str]] = None,
) -> List[Component]:
    """Parse components for Suzuki reaction"""
    components = []
    should_continue = True
    logger = logging.getLogger(__name__)
    for component_dict in components_list:
        role = component_dict["role"] 
        if role == PistachioRole.reactant:
            component = handle_reactant(component_dict, smiles_cleaner)
        elif role == PistachioRole.product:
            component = handle_product(component_dict, smiles_cleaner)
        elif role == PistachioRole.agent or role == PistachioRole.catalyst:
            component = handle_agent_catalyst(
                component_dict, smiles_cleaner, ligand_mols, base_smiles
            )
        elif role == PistachioRole.solvent:
            component = handle_solvent(component_dict, smiles_cleaner, solvent_smiles)
            
        if component is None:
            continue
        elif type(component) is Component:
            components.append(component)
        elif type(component) is list:
            if len(component) > 0:
                components.extend(component)
    return components


def handle_reactant(component: Dict, smiles_cleaner: SmilesCleaner) -> Union[Component, None]:
    smiles = smiles_cleaner(component.get("smiles"))
    if smiles is None:
        return
    mol = Chem.MolFromSmiles(smiles)

    # Quantities
    quantity_dicts = component.get("quantities")
    if quantity_dicts is not None:
        quantities = get_quantities(quantity_dicts)
    else:
        quantities = {}

    # Check if reactant is a Boronic acid
    boronic_acid = Chem.MolFromSmiles("B(O)O")
    boronic_acid_matches = mol.GetSubstructMatches(boronic_acid)
    if len(boronic_acid_matches) > 0:
        return Component(
            role=SuzukiRole.boronic_acid_reactant,
            smiles=smiles,
            **quantities
        )
    
    # Check if reactant is a halide
    halides = [Chem.MolFromSmiles(e) for e in ["Br", "Cl", "I", "F"]]
    for halide in halides:
        halide_mataches = mol.GetSubstructMatches(halide)
        if len(halide_mataches) > 0:
            return Component(
                role=SuzukiRole.organohalide_reactant,
                smiles=smiles,
                **quantities
            )

def handle_product(component: Dict, smiles_cleaner: SmilesCleaner) -> Component:
    smiles = smiles_cleaner(component.get("smiles"))
    # Quantities
    quantity_dicts = component.get("quantities")
    if quantity_dicts is not None:
        quantities = get_quantities(quantity_dicts)
    else:
        quantities = {}
    
    if smiles is not None:
        return Component(
            role=SuzukiRole.product,
            smiles=smiles,
            **quantities
        )

def handle_agent_catalyst(
    component: Dict, 
    smiles_cleaner: SmilesCleaner, 
    base_smiles: List[str],
    solvent_smiles: Optional[List[str]] = None,
    ligand_mols: Optional[Dict[str, Chem.Mol]] = None,
) -> List[Component]:
    # Split string on period
    smiles = component.get("smiles")
    if smiles is None:
        return []
    smiles_list = smiles.split(".")
    components = []

    # Quantities
    quantity_dicts = component.get("quantities")
    if quantity_dicts is not None:
        quantities = get_quantities(quantity_dicts)
    else:
        quantities = {}
    
    # Handle ligands and bases
    phosphorus = Chem.MolFromSmiles("P")
    palladium = Chem.MolFromSmiles("[Pd]")
    for smiles in smiles_list:
        smiles = smiles_cleaner(smiles)
        if smiles is None:
            continue
        mol = Chem.MolFromSmiles(smiles)
        # Bases
        if smiles in base_smiles:
            components.append(
                Component(
                    role=SuzukiRole.base,
                    smiles=smiles,
                    **quantities
                )
            )
        # P-Ligands
        elif len(mol.GetSubstructMatches(phosphorus)) > 0:
            if ligand_mols is not None:
                smiles, counts = get_ligand_smiles_counts(smiles, ligand_mols)
            components.append(
                Component(
                    role=SuzukiRole.ligand,
                    smiles=smiles,
                    **quantities
                )
            )
        elif len(mol.GetSubstructMatches(palladium)) > 0:
            components.append(
                Component(
                    role=SuzukiRole.pre_catalyst,
                    smiles=smiles,
                    **quantities
                )
            )
        else:
            c = handle_solvent(component, smiles_cleaner, solvent_smiles)
            if c is None:
                components.append(c)
            else:
                components.append(
                    Component(
                        role=SuzukiRole.agent,
                        smiles=smiles,
                        **quantities
                    )
                )

    return components


def get_ligand_smiles_counts(smiles: str, ligand_mols: Dict[str, Chem.Mol])-> dict:
    # Get rdkit mol for clean smiles
    mol = Chem.MolFromSmiles(smiles)

    # Loop through each perera ligand and check for matches
    ligand_matches = {}
    for ligand_smiles, ligand_mol in ligand_mols.items():
        ligand_matches[ligand_smiles] = len(mol.GetSubstructMatches(ligand_mol))

    # If there are multiple matches, take the substructure with the highest molecular weight
    n = 0
    wt = 0
    final_smiles = ""
    for ligand_smiles, n_matches in ligand_matches.items():
        if n_matches > 0 and n==0:
            final_smiles = ligand_smiles
            n = n_matches
            wt = Chem.rdMolDescriptors.CalcExactMolWt(ligand_mols[ligand_smiles])
        elif n_matches > 0 and n > 0:
            wt_temp = Chem.rdMolDescriptors.CalcExactMolWt(ligand_mols[ligand_smiles])
            if wt_temp > wt:
                final_smiles = ligand_smiles
                n = n_matches
                wt = wt_temp
    return final_smiles, n

def handle_solvent(
    component: Dict, smiles_cleaner: SmilesCleaner, solvent_smiles: Optional[List[str]] = None
) -> Union[Component, None]:
    smiles = smiles_cleaner(component.get("smiles"))
    # Quantities
    quantity_dicts = component.get("quantities")
    volume = None
    if quantity_dicts is not None:
        quantities = get_quantities(quantity_dicts)
    else: 
        quantities = {}
    if solvent_smiles is not None:
        if smiles in solvent_smiles:
            return Component(role=SuzukiRole.solvent, smiles=smiles, **quantities)
    elif smiles is not None:
        return Component(role=SuzukiRole.solvent, smiles=smiles, **quantities)


# Split  -------------------------------------------------------------

def make_mol(smiles_input):
    mol_output = Chem.MolFromSmiles(smiles_input)

    return mol_output


def scaffold_split(df, column_to_split:str, test_size=0.1):
  """" Split a dataframe by scaffolds
  
  Parameters
  ----------
  df : pd.DataFrame
    Dataframe with dataset
  column_to_split : str
    Column with smiles strings to split by scaffold
  test_size : float
    Size of the test set
    
  Returns
  -------
  train_df : pd.DataFrame
    The dataframe with the training set
  test_df : pd.DataFrame
    The dataframe with the training set
  train_img: 
    An image of the scaffolds in the training set
  test_img:
    An image of the scaffolds in the test set
    
  Examples
  --------
  >>> train_df, test_df, train_img, test_img = scaffold_split(df, "reactant_1_smiles")
  """
  logger = logging.getLogger(__name__)

  # Get scaffolds
  unique_mols = pd.unique(df[column_to_split])
  scaffolds, scaffold_indices = scaffold_groups(unique_mols)
  logger.info(f"Number of scaffolds: {len(scaffolds)}")

  # Merge scaffolds into dataframe
  scaffold_df = pd.DataFrame(
      {"scaffold_smiles": unique_mols, "scaffold_index": scaffold_indices}
  )
  df = df.merge(scaffold_df, left_on=column_to_split, right_on="scaffold_smiles", how="left")

  # Train-test split
  splitter = GroupShuffleSplit(n_splits=2, test_size=test_size, random_state=1995)
  train_indices, test_indices = next(
      splitter.split(X=df, groups=df["scaffold_index"])
  )
  train_df = df.iloc[train_indices]
  test_df = df.iloc[test_indices]
  
  # Log statistics
  logger.info(f"Training set size: ~{train_df.shape[0]/1e3:.0f}k")
  logger.info(f"Test set size: ~{test_df.shape[0]/1e3:.0f}k")
  #null_scaffold = scaffolds[""]
  scaffold_counts = df["scaffold_index"].value_counts()
  #logger.info(
   #   f"Number of records that match the null scaffold: ~{scaffold_counts.loc[null_scaffold]/1e3:.0f}k"
  #)

  # Visualize scaffolds
  train_img = visualize_scaffolds(train_df, scaffolds)
  test_img = visualize_scaffolds(test_df, scaffolds)

  return train_df, test_df, train_img, test_img

def generate_scaffold(mol: Union[str, Chem.Mol, Tuple[Chem.Mol, Chem.Mol]], include_chirality: bool = False) -> str:
    """
    Computes the Bemis-Murcko scaffold for a SMILES string.
    :param mol: A SMILES or an RDKit molecule.
    :param include_chirality: Whether to include chirality in the computed scaffold..
    :return: The Bemis-Murcko scaffold for the molecule.
    """
    if isinstance(mol, str):
        #mol = make_mol(mol, keep_h = False)
        # Alexander deleted the keep_h argument because it gave errors
        mol = make_mol(mol)
    if isinstance(mol, tuple):
        mol = mol[0]
    scaffold = MurckoScaffold.MurckoScaffoldSmiles(mol = mol, includeChirality = include_chirality)

    return scaffold


def scaffold_groups(mols: List[str]):
    """Find all the scaffolds and reference each one by index
    Parameters
    ---------
    mols: list of str
        The list of smiles strings
    Returns
    -------
    scaffolds, scaffold_indices
        scaffolds is a dictionary mapping scaffold to index.
        scaffold_indices gives the index of the scaffold for each molecule
    """
    logger = logging.getLogger(__name__)
    scaffolds = dict()
    scaffold_num = 0
    scaffold_list = [""] * len(mols)
    n_mols = len(mols)
    logger.info("Finding scaffolds")
    for i, mol in tqdm(enumerate(mols), total=len(mols)):
        scaffold = generate_scaffold(mol)
        if scaffold not in scaffolds:
            scaffolds[scaffold] = scaffold_num
            scaffold_num += 1
        scaffold_list[i] = scaffold
    scaffold_indices = [scaffolds[s] for s in scaffold_list]
    return scaffolds, scaffold_indices
  
def visualize_scaffolds(df: pd.DataFrame, scaffolds: dict):
    """Visualize scaffolds with counts of scaffolds below each molecule"""
    scaffold_counts = df["scaffold_index"].value_counts()
    #null_scaffold = scaffolds[""]
    scaffold_idx_to_smiles = {idx: scaffold for scaffold, idx in scaffolds.items()}
    img = Chem.Draw.MolsToGridImage(
        [
            Chem.MolFromSmiles(scaffold_idx_to_smiles[idx])
            for idx in scaffold_counts.index
        ],
        molsPerRow=6,
        subImgSize=(200, 200),
        legends=[
            str(scaffold_counts.iloc[i])
            for i in range(len(scaffold_counts))
        ],
        returnPNG=False,
    )
    return img
