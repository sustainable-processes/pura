from .compound import Compound

# from .units import Time, Temperature
from typing import List, Tuple, Optional, Dict
import pandas as pd
from enum import Enum
from pydantic import BaseModel


class ReactionIdentifierType(BaseModel):
    """Generic identifiers for a reaction"""

    UNSPECIFIED = 0
    CUSTOM = 1
    REACTION_SMILES = 2
    REACTION_CXSMILES = 6
    # Extended SMILES.
    RDFILE = 3
    # Reaction data file.
    RINCHI = 4
    # Reaction InChI.
    NAME = 5
    # Named reaction or reaction category.


class ReactionRole(BaseModel):
    UNSPECIFIED = 0
    # A reactant is any compound that contributes atoms to a desired or
    # observed product.
    REACTANT = 1
    REAGENT = 2
    SOLVENT = 3
    CATALYST = 4
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


class ReactionIdentifier(BaseModel):
    identifier_type: ReactionIdentifierType
    value: str
    details: Optional[str] = None


class ReactionInput(BaseModel):
    # This is different than ORD. We don't include ReactionRole inside
    # the Compound object because compounds can exist independent of reactions
    components: Dict[ReactionRole, Compound]
    addition_order: Optional[int] = None
    # addition_time: Optional[Time] = None
    # addition_duration: Optional[Time] = None
    flow_rate: Optional[float] = None
    # addition_temperature: Optional[Temperature] = None


class ReactionConditions(BaseModel):
    pass


class Reaction(BaseModel):
    """

    Reaction identifiers are the whole reaction and should be convertible
    to inputs + outcomes.
    """

    identifiers: List[ReactionIdentifier]

    # List of pure substances or mixtures that were added to the reaction vessel.
    # This is a map instead of a repeated field to simplify reaction templating
    # through the use of keys. String keys are simple descriptions and are
    # present only for convenience.
    inputs: List[ReactionInput]
    conditions: ReactionConditions

    def to_row(self):
        return


def to_frame(reactions: List[Reaction], columns: List[str]) -> pd.DataFrame:
    reaction_dicts = [reaction.to_row() for reaction in reactions]
    return pd.DataFrame(reaction_dicts)
