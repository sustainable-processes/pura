import pytest
from unittest.mock import AsyncMock
from aiohttp import ClientSession
from pura.services.pubchem import PubChem
from pura.compound import CompoundIdentifier, CompoundIdentifierType

@pytest.mark.asyncio
async def test_pubchem_resolve_compound():
    # Create a mock ClientSession
    session = AsyncMock(spec=ClientSession)

    # Create an instance of PubChem
    service = PubChem()

    # Define the test data
    input_identifier = CompoundIdentifier(
        identifier_type=CompoundIdentifierType.NAME, value="Josiphos SL-J001-1"
    )
    output_identifier_types = [
        CompoundIdentifierType.SMILES,
        CompoundIdentifierType.INCHI_KEY,
    ]

    # Call the resolve_compound method with the test data
    resolved = await service.resolve_compound(
        session=session,
        input_identifier=input_identifier,
        output_identifier_types=output_identifier_types,
    )

    # Define the expected output
    expected_output = [
        CompoundIdentifier(
            identifier_type=CompoundIdentifierType.SMILES,
            value="CC(C1CCCC1P(C2=CC=CC=C2)C3=CC=CC=C3)P(C4CCCCC4)C5CCCCC5.C1CCCC1.[Fe]"
        ),
        CompoundIdentifier(
            identifier_type=CompoundIdentifierType.INCHI_KEY,
            value="ULNUEACLWYUUMO-UHFFFAOYSA-N"
        ),
    ]

    # Assert that the resolved output matches the expected output
    assert resolved == expected_output
```

