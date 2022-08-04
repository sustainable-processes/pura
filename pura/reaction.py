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


def to_frame(reactions: List[Reaction], columns: List[str]) -> pd.DataFrame:
    reaction_dicts = [reaction.to_row() for reaction in reactions]
    return pd.DataFrame(reaction_dicts)
