from typing import List, Tuple
from .reaction import ReactionIdentifier, ReactionIdentifierType


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
