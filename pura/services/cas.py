from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from pura.utils import inverse_map
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
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[Union[CompoundIdentifier, None]]:
        input_identifier_cas = INPUT_IDENTIFIER_MAP.get(
            input_identifier.identifier_type
        )
        if input_identifier_cas is None:
            raise ValueError(
                f"{input_identifier.identifier_type} is not one of the valid identifier types for the CAS."
            )

        output_representations = [
            OUTPUT_IDENTIFIER_MAP.get(output_identifier_type)
            for output_identifier_type in output_identifier_types
        ]
        if not any(output_representations):
            raise ValueError(
                f"{output_identifier_types} are not one of the valid identifier types for CAS."
            )

        results = await common_chemistry_request(
            session,
            input_identifier.value,
            input_identifier_cas,
            output_representations,
        )

        inv_identifier_map = inverse_map(OUTPUT_IDENTIFIER_MAP)
        return [
            CompoundIdentifier(
                identifier_type=inv_identifier_map[output_representation],
                value=result,
            )
            for output_representation, result in results.items()
        ]


async def common_chemistry_request(
    session: ClientSession,
    identifier: str,
    input_representation: str,
    output_representations: List[str],
) -> str:
    api_url = BASE_API + "/default/detail"
    params = {input_representation: identifier}

    logger.debug(f"Request: {api_url}")
    async with session.get(api_url, params=params) as resp:
        logger.debug(f"Response status: {resp.status}")
        if resp.status == 404:
            raise HTTPNotFound(text=f"Failed request for {identifier}")
        response = await resp.json()
    return {
        output_representation: response.get(output_representation)
        for output_representation in output_representations
    }
    # return response.get(output_representation)
