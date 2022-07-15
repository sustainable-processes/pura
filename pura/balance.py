from typing import List, Tuple, Optional, Dict
import pandas as pd
from enum import Enum

from pydantic import BaseModel


class ReactionIdentifierType:
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


class ReactionIdentifier(BaseModel):
    identifier_type: ReactionIdentifierType
    value: str
    details: Optional[str] = None


class Reaction(BaseModel):
    def __init__(
        self, reaction_identifier: List[ReactionIdentifier], conditions: Dict
    ) -> None:
        pass

    def to_row(self):
        return


def to_frame(
    reactions: List[Reaction], columns: List[str], max_rows: int = 2
) -> pd.DataFrame:
    reaction_dicts = [reaction.to_row() for reaction in reactions]
    return pd.DataFrame(reaction_dicts)


# def balance_reactions(
#     reactions: List[ReactionIdentifier],
#     ignore_hydrogens: bool = True,
#     use_ray: bool = False,
# ) -> List[Tuple[ReactionIdentifier, ReactionIdentifier]]:
#     if use_ray:
#         for reaction in reactions:
#             ray.remote(balance_reaction_ray)
#     else:
#         [balance_reaction(reaction) for reaction in reactions]


# @ray.remote
# def balance_reaction_ray(**kwargs):
#     balance_reaction(**kwargs)


def balance_reaction(
    input_reaction_identifier: ReactionIdentifier,
) -> Tuple[ReactionIdentifier, ReactionIdentifier]:
    pass

    # do balancing
    output_reaction_identifier = None
    return input_reaction_identifier, output_reaction_identifier
