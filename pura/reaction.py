from .compound import Compound

# from .units import Time, Temperature
from typing import List, Tuple, Optional, Dict, Union
import pandas as pd
from enum import Enum
from dataclasses import dataclass, asdict


# class ReactionIdentifierType(Enum):
#     """Generic identifiers for a reaction"""

#     UNSPECIFIED = 0

#     CUSTOM = 1
#     # Reaction SMILES
#     REACTION_SMILES = 2
#     # Reaction_CXSMILES
#     REACTION_CXSMILES = 6
#     # Extended SMILES.
#     RDFILE = 3
#     # Reaction data file.
#     RINCHI = 4
#     # Reaction InChI.
#     NAME = 5
#     # Named reaction or reaction category.


class ReactionRole(Enum):
    UNSPECIFIED = 0
    # A reactant is any compound that contributes atoms to a desired or
    # observed product.
    REACTANT = 1
    # An agent is any compound/species that enables the reaction but does not contribute atoms to the product.
    # Reagent, solvents, and catalysts are all examples of agents.
    AGENT = 2
    REAGENT = 3
    SOLVENT = 4
    CATALYST = 5
    # The workup role is used when defining quenches, buffer additives for
    # liquid-liquid separations, etc.
    WORKUP = 5
    # Internal standards can be included as part of a reaction input (when
    # added prior to the start of the reaction) or as part of a workup
    # step of addition.
    INTERNAL_STANDARD = 6
    AUTHENTIC_STANDARD = 7
    # A product can be any species produced by the reaction, whether desired
    # or undesired.
    PRODUCT = 8


@dataclass
class ReactionInput:
    compound: Compound
    role: ReactionRole
    addition_order: Optional[int] = None
    # addition_time: Optional[Time] = None
    # addition_duration: Optional[Time] = None
    # flow_rate: Optional[float] = None
    # addition_temperature: Optional[Temperature] = None


@dataclass
class ReactionOutcome:
    products: Compound
    # Specify the role the product played in the reaction. This
    # aids in constructing proper reaction SMILES and product lists.
    # Species produced by the reaction should be identified as PRODUCT,
    # whether desired or undesired. Recovered starting materials can be
    # specified as REACTANT, REAGENT, CATALYST, etc. as appropriate
    role: ReactionRole = ReactionRole.PRODUCT


class ReactionConditions:
    temperature: Optional[float] = None
    pressure: Optional[float] = None
    reflux: Optional[bool] = None
    ph: Optional[float] = None


@dataclass
class Reaction:
    inputs: List[ReactionInput]
    conditions: ReactionConditions
    outcomes: List[ReactionOutcome]
    rxn_yield: Optional[float] = None

    def rxn_smiles(self):
        """Construct a reaction SMILES"""
        pass

    def to_dict(self):
        """Convert to a dictionary"""
        return asdict(self)

    def from_dict(self, d: Dict):
        """Convert from a dictionary"""
        raise NotImplementedError()

    @property
    def reactants(self):
        return [i for i in self.inputs if i.role == ReactionRole.REACTANT]

    @property
    def products(self):
        return [o for o in self.outcomes if o.role == ReactionRole.PRODUCT]

    @property
    def agents(self):
        return [i for i in self.inputs if i.role == ReactionRole.AGENT]


def to_frame(reactions: List[Reaction], columns: List[str]) -> pd.DataFrame:
    reaction_dicts = [reaction.to_row() for reaction in reactions]
    return pd.DataFrame(reaction_dicts)


def reaction_from_smiles(
    smiles: str,
    role_lookup: Optional[Dict[ReactionRole, List[Union[str, Compound]]]] = None,
) -> Reaction:
    """Construct a reaction from a reaction SMILES

    Parameters
    ----------
    smiles : str
        Reaction SMILES
    role_lookup : Optional[Dict[ReactionRole, List[Union[str, Compound]]]], optional
        A dictionary where keys are ReactionRole and values are lists of compounds or SMILES strings, by default None.

    Returns
    -------
    Reaction

    Examples
    --------
    >>> role_lookup = {ReactionRole.SOLVENT: ['CCO', 'CCO']}
    >>> reaction_from_smiles('CC(=O)O.[Na+].[O-]S(=O)(=O)C1=CC=CC=C1>>CC(=O)O.[Na+].[O-]S(=O)(=O)C1=CC=CC=C1', role_lookup=role_lookup)

    """
    raise NotImplementedError()
