from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from typing import Optional, List
import os
import logging


class ChemSpider(Service):
    """

    Arguments
    ---------
    token : str, optional
        Chemspider API token. If not passed, will try to read from the
        environmental variable `CHEMSPIDER_TOKEN`.
    """

    def __init__(self, token: Optional[str] = None) -> None:
        if token is None:
            env_token = os.environ.get("CHEMSPIDER_TOKEN")
            if env_token is None:
                raise ValueError(
                    "No Chemspider API token passed or found in the environment."
                )
        self.token = token

    async def post(
        self,
        session: ClientSession,
        api: str,
        namespace: str,
        endpoint: str,
        params: dict = None,
        json: dict = None,
    ):
        logger = logging.getLogger(__name__)
        # Construct request URL
        url = "{}/{}/{}/{}/{}".format(
            self.api_url, api, self.api_version, namespace, endpoint
        )

        # Set apikey header
        headers = {"apikey": self.api_key}

        logger.debug("{} : {} : {} : {}".format(url, headers, params, json))

        # Make request
        async with session.post(url, data=json, params=params, headers=headers) as resp:
            response = await resp.json()

        return response

    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_type: CompoundIdentifierType,
    ) -> List[CompoundIdentifierType]:

        json = {
            "name": name,
            "orderBy": ORDERS.get(order),
            "orderDirection": DIRECTIONS.get(direction),
        }
        response = self.post(
            api="compounds", namespace="filter", endpoint="name", json=json
        )
        return response["queryId"]
