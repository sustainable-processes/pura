# from pura.reaction import reaction_from_smiles, ReactionRole, Reaction
import json
from pura.reaction import Reaction, reaction_from_smiles, ReactionRole
from pura.units import *
import pytest


def test_reaction_from_smiles():
    rxn_smiles = "Brc1cncnc1.OB(O)c1ccsc1.COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O.O=P([O-])([O-])[O-].[K+].[K+].[K+].C[Mg+].[Cl-]>>c1ncc(-c2ccsc2)cn1"
    role_lookup = {ReactionRole.REAGENT: ["[K+]", "[Cl-]"]}
    rxn = reaction_from_smiles(
        rxn_smiles,
        reaction_time=60 * ureg.minute,
        reaction_yield=50.0,
        desired_product_check=lambda c: True,
        role_lookup=role_lookup,
    )
    for c in rxn.reagent_compounds:
        assert c.to_smiles() in role_lookup[ReactionRole.REAGENT]
    assert rxn.reaction_yield == 50.0
    assert rxn.reaction_time == 1 * ureg.hour

    # Should raise value error for negative yield
    with pytest.raises(ValueError):
        rxn = reaction_from_smiles(
            rxn_smiles,
            reaction_time=60 * ureg.minute,
            reaction_yield=-1.0,
            desired_product_check=lambda c: True,
            role_lookup=role_lookup,
        )

    # Should raise value error for negative time
    with pytest.raises(ValueError):
        rxn = reaction_from_smiles(
            rxn_smiles,
            reaction_time=-60 * ureg.minute,
            reaction_yield=50.0,
            desired_product_check=lambda c: True,
            role_lookup=role_lookup,
        )

    # Should raise value error if reaction yield passed but not desired product check
    with pytest.raises(ValueError):
        rxn = reaction_from_smiles(
            rxn_smiles,
            reaction_time=60 * ureg.minute,
            reaction_yield=50.0,
            desired_product_check=None,
            role_lookup=role_lookup,
        )

    # Should warn for yield greater than 100
    with pytest.warns():
        rxn = reaction_from_smiles(
            rxn_smiles,
            reaction_time=60 * ureg.minute,
            reaction_yield=101.0,
            desired_product_check=lambda c: True,
            role_lookup=role_lookup,
        )

    # test serialization
    output = rxn.json()
    rxn == Reaction.parse_obj(json.loads(output))


if __name__ == "__main__":
    print(Reaction.schema_json(indent=2))
