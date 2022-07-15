from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from typing import List, Union
from aiohttp.web_exceptions import HTTPNotFound
import logging

BASE_API = "https://rboq1qukh0.execute-api.us-east-2.amazonaws.com"

logger = logging.getLogger(__name__)

INPUT_IDENTIFIER_MAP = {
    CompoundIdentifierType.CAS_NUMBER: "cas_rn",
    # CompoundIdentifierType.SMILES: "smiles",
    # CompoundIdentifierType.INCHI: "stdinchi",
    # CompoundIdentifierType.INCHI_KEY: "stdinchikey",
}

OUTPUT_IDENTIFIER_MAP = {
    CompoundIdentifierType.CAS_NUMBER: "rn",
    CompoundIdentifierType.SMILES: "canonicalSmile",
    CompoundIdentifierType.INCHI: "stdinchi",
    CompoundIdentifierType.INCHI_KEY: "inchiKey",
}


class CAS(Service):
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_type: CompoundIdentifierType,
    ) -> List[Union[CompoundIdentifier, None]]:
        input_identifier_cas = INPUT_IDENTIFIER_MAP.get(
            input_identifier.identifier_type
        )
        if input_identifier_cas is None:
            raise ValueError(
                f"{input_identifier.identifier_type} is not one of the valid identifier types for the Pubchem."
            )

        output_identifier_cas = OUTPUT_IDENTIFIER_MAP.get(output_identifier_type)
        if output_identifier_cas is None:
            raise ValueError(
                f"{output_identifier_type} is not one of the valid identifier types for the chemical identifier resolver."
            )

        result = await common_chemistry_request(
            session, input_identifier.value, input_identifier_cas, output_identifier_cas
        )

        if result is not None:
            return [
                CompoundIdentifier(
                    identifier_type=output_identifier_type,
                    value=result,
                )
            ]
        else:
            return []


async def common_chemistry_request(
    session: ClientSession,
    identifier: str,
    input_representation: str,
    output_representation: str,
) -> str:
    api_url = BASE_API + "/default/detail"
    params = {input_representation: identifier}

    logger.debug(f"Request: {api_url}")
    async with session.get(api_url, params=params) as resp:
        logger.debug(f"Response status: {resp.status}")
        if resp.status == 404:
            raise HTTPNotFound(text=f"Failed request for {identifier}")
        response = await resp.json()
    return response.get(output_representation)
