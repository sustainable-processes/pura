import asyncio
import pytest
from pura.resolvers import resolve_identifiers, CompoundResolver
from pura.compound import Compound, CompoundIdentifier, CompoundIdentifierType
from pura.services import CIR, Opsin, ChemSpider, CAS
from pura.services.pubchem import PubChem, OUTPUT_IDENTIFIER_MAP, autocomplete
from rdkit import Chem
from aiohttp import *
from dotenv import load_dotenv

load_dotenv()


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
        [
            "Pd(OAc)2",
            "Josiphos SL-J001-1",
            "Rh(NBD)2BF4",
            "Dichloro(p-cymene)ruthenium(II) dimer",
            "DuPhos",
        ],
        input_identifer_type=CompoundIdentifierType.NAME,
        output_identifier_type=CompoundIdentifierType.SMILES,
        backup_identifier_types=[
            CompoundIdentifierType.INCHI_KEY,
            CompoundIdentifierType.CAS_NUMBER,
        ],
        services=[PubChem(autocomplete=True), CIR(), CAS(), ChemSpider()],
        agreement=1,
        silent=True,
    )
    print("\nResults\n")
    for input_compound, resolved_identifiers in resolved:
        print(input_compound.identifiers[0].value, resolved_identifiers)
        print()


# Josiphos SL-J001-1 [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='C1CCCC1.CC(C1CCCC1P(c1ccccc1)c1ccccc1)P(C1CCCCC1)C1CCCCC1.[Fe]', details=None)]

# Rh(NBD)2BF4 [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='C1=CC2C=CC1C2.C1=CC2C=CC1C2.F[B-](F)(F)F.[Rh]', details=None)]

# Dichloro(p-cymene)ruthenium(II) dimer [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='Cc1ccc(C(C)C)cc1.Cc1ccc(C(C)C)cc1.Cl[Ru]Cl.Cl[Ru]Cl', details=None)]

# DuPhos [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='CC(C)C1CCC(C(C)C)P1c1ccccc1P1C(C(C)C)CCC1C(C)C', details=None)]


def test_opsin():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(async_test_opsin())


async def async_test_opsin():
    async with ClientSession() as session:
        service = Opsin()
        resolved = await service.resolve_compound(
            session=session,
            input_identifier=CompoundIdentifier(
                identifier_type=CompoundIdentifierType.NAME, value="methylbenzene"
            ),
            output_identifier_types=[
                CompoundIdentifierType.SMILES,
                CompoundIdentifierType.INCHI_KEY,
            ],
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


def test_chempsider():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(async_test_chempsider())


async def async_test_chempsider():
    async with ClientSession() as session:
        service = ChemSpider()
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


def test_cas():
    loop = asyncio.get_event_loop()
    loop.run_until_complete(async_test_cas())


async def async_test_cas():
    async with ClientSession() as session:
        service = CAS()
        resolved = await service.resolve_compound(
            session=session,
            input_identifier=CompoundIdentifier(
                identifier_type=CompoundIdentifierType.CAS_NUMBER, value="108-88-3"
            ),
            output_identifier_types=[
                CompoundIdentifierType.SMILES,
                CompoundIdentifierType.INCHI_KEY,
            ],
        )
        print(resolved)


if __name__ == "__main__":
    test_resolve_backup_identifiers()
    # test_pubchem()
    # test_opsin()
    # test_chempsider()
    # test_cas()
