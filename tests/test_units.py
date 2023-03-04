from pura.units import ureg
from pura.reaction import ReactionConditions
from pura.compound import Compound, CompoundIdentifier, CompoundIdentifierType
import json

Q_ = ureg.Quantity

print("Schema:")
print(Compound.schema_json(indent=2))

c = Compound(
    identifiers=[
        CompoundIdentifier(identifier_type=CompoundIdentifierType.SMILES, value="CO")
    ],
    quantity=25 * ureg.mol,
)
print("Example:")
js = c.json(indent=2)
print(js)


print("After serialization and deserialization:")
print(Compound.parse_obj(json.loads(js)).json(indent=2))
# rc = ReactionConditions(temperature=Q_(25, ureg.degC))
# print(rc)
