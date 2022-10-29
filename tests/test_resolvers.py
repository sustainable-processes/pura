import asyncio
import pytest
from pura.resolvers import resolve_identifiers, CompoundResolver
from pura.compound import Compound, CompoundIdentifier, CompoundIdentifierType
from pura.services import CIR
from pura.services.pubchem import PubChem, OUTPUT_IDENTIFIER_MAP
from rdkit import Chem
from aiohttp import *


@pytest.mark.parametrize("identifier_type", OUTPUT_IDENTIFIER_MAP)
def test_resolve_identifiers_no_agreement(identifier_type):
    resolved = resolve_identifiers(
        ["Josiphos SL-J001-1"],
        input_identifer_type=CompoundIdentifierType.NAME,
        output_identifier_type=identifier_type,
        services=[
            PubChem(),
            # CIR(),
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


def test_resolve_backup_identifiers():
    resolved = resolve_identifiers(
        ["Josiphos SL-J001-1"],
        input_identifer_type=CompoundIdentifierType.NAME,
        output_identifier_type=CompoundIdentifierType.SMILES,
        backup_identifier_types=[CompoundIdentifierType.INCHI_KEY],
        services=[
            PubChem(),
            # CIR(),
        ],
        agreement=0,
    )
    print(resolved)


def test_pubchem():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(async_test_pubchem())


async def async_test_pubchem():
    async with ClientSession() as session:
        service = PubChem()
        resolved = await service.resolve_compound(
            session=session,
            input_identifier=CompoundIdentifier(
                identifier_type=CompoundIdentifierType.NAME, value="Josiphos SL-J001-1"
            ),
            output_identifier_types=[
                CompoundIdentifierType.SMILES,
                CompoundIdentifierType.INCHI_KEY,
            ],
        )
        print(resolved)


if __name__ == "__main__":
    # test_resolve_backup_identifiers()
    test_pubchem()
