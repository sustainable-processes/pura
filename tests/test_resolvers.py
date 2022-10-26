import pytest
from pura.resolvers import resolve_identifiers, CompoundResolver
from pura.compound import Compound, CompoundIdentifier, CompoundIdentifierType
from pura.services import CIR
from pura.services.pubchem import PubChem, OUTPUT_IDENTIFIER_MAP
from rdkit import Chem


@pytest.mark.parametrize("identifier_type", OUTPUT_IDENTIFIER_MAP)
def test_resolve_identifiers_no_agreement(identifier_type):
    resolved = resolve_identifiers(
        ["Josiphos SL-J001-1"],
        input_identifer_type=CompoundIdentifierType.NAME,
        output_identifier_type=identifier_type,
        services=[
            PubChem(),
            CIR(),
        ],
        agreement=0,
    )
    print(resolved)
    # correct_smiles = (
    #     "CC(C1CCCC1P(C2=CC=CC=C2)C3=CC=CC=C3)P(C4CCCCC4)C5CCCCC5.C1CCCC1.[Fe]"
    # )
    # correct_inchi = "ULNUEACLWYUUMO-UHFFFAOYSA-N"
    # # correct_smiles = Chem.CanonSmiles(correct_smiles)
    # input_compound, resolved_identifiers = resolved[0]
    # assert resolved_identifiers[0].value == correct_inchi
