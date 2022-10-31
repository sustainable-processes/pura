from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from aiohttp.web_exceptions import HTTPNotFound
from typing import List, Union
from urllib.parse import quote
import logging


BASE_API = "https://opsin.ch.cam.ac.uk/opsin/"

IDENTIFIER_MAP = {
    CompoundIdentifierType.SMILES: "smiles",
    CompoundIdentifierType.INCHI: "stdinchi",
    CompoundIdentifierType.INCHI_KEY: "stdinchikey",
}

logger = logging.getLogger(__name__)


class Opsin(Service):
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[CompoundIdentifier]:
        representations = [
            IDENTIFIER_MAP.get(output_identifier_type)
            for output_identifier_type in output_identifier_types
        ]
        if not any(representations):
            raise ValueError(
                f"{output_identifier_types} is not one of the valid identifier types for Opsin."
            )

        output_compound_identifiers = []
        for representation, output_identifier_type in zip(
            representations, output_identifier_types
        ):
            result = await opsin_request(
                session, input_identifier.value, representation
            )

            if result is not None:
                output_compound_identifiers += [
                    CompoundIdentifier(
                        identifier_type=output_identifier_type,
                        value=result,
                    )
                ]
        return output_compound_identifiers


async def opsin_request(
    session: ClientSession, name: str, output_representation: str
) -> str:
    api_url = BASE_API + quote(name)

    logger.debug(f"Request: {api_url}")
    async with session.get(api_url) as resp:
        logging.debug(f"Response status: {resp.status}")
        if resp.status == 404:
            raise HTTPNotFound(text=f"Failed request for {name}")
        response = await resp.json()
    return response.get(output_representation)
