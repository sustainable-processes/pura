"""
Based heavily on Chemspipy by Matt Swain
https://github.com/mcs07/ChemSpiPy
"""
from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
import asyncio
from typing import Optional, List
import os
import logging
from dotenv import load_dotenv

load_dotenv()
logger = logging.getLogger(__name__)

#: 2D coordinate dimensions
MOL2D = "2d"
#: 3D coordinate dimensions
MOL3D = "3d"
#: Both coordinate dimensions
BOTH = "both"

#: Ascending sort direction
ASCENDING = "ascending"
#: Descending sort direction
DESCENDING = "descending"

#: Record ID sort order
RECORD_ID = "record_id"
#: CSID sort order (same as RECORD_ID, kept for backwards compatibility)
CSID = "csid"
#: Mass defect sort order
MASS_DEFECT = "mass_defect"
#: Molecular weight sort order
MOLECULAR_WEIGHT = "molecular_weight"
#: Reference count sort order
REFERENCE_COUNT = "reference_count"
#: Datasource count sort order
DATASOURCE_COUNT = "datasource_count"
#: Pubmed count sort order
PUBMED_COUNT = "pubmed_count"
#: RSC count sort order
RSC_COUNT = "rsc_count"

#: Map sort directions to strings required by REST API.
DIRECTIONS = {ASCENDING: "ascending", DESCENDING: "descending"}

#: Map sort orders to strings required by REST API.
ORDERS = {
    RECORD_ID: "recordId",
    CSID: "recordId",
    MASS_DEFECT: "massDefect",
    MOLECULAR_WEIGHT: "molecularWeight",
    REFERENCE_COUNT: "referenceCount",
    DATASOURCE_COUNT: "dataSourceCount",
    PUBMED_COUNT: "pubMedCount",
    RSC_COUNT: "rscCount",
}

#: All available compound details fields.
FIELDS = [
    "SMILES",
    "Formula",
    "AverageMass",
    "MolecularWeight",
    "MonoisotopicMass",
    "NominalMass",
    "CommonName",
    "ReferenceCount",
    "DataSourceCount",
    "PubMedCount",
    "RSCCount",
    "Mol2D",
    "Mol3D",
]
IDENTIFIER_TYPES = [
    CompoundIdentifierType.CHEMSPIDER_ID,
    CompoundIdentifierType.SMILES,
    CompoundIdentifierType.INCHI,
    CompoundIdentifierType.INCHI_KEY,
]


class ChemSpider(Service):
    api_url = "https://api.rsc.org"
    api_version = "v1"

    """

    Arguments
    ---------
    token : str, optional
        Chemspider API token. If not passed, will try to read from the
        environmental variable `CHEMSPIDER_TOKEN`.
    """

    def __init__(self, token: Optional[str] = None) -> None:
        if token is None:
            token = os.environ.get("CHEMSPIDER_TOKEN")
            if token is None:
                raise ValueError(
                    "No Chemspider API token passed or found in the environment."
                )
        self.token = token

    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[CompoundIdentifierType]:
        if not any(
            [
                output_identifier_type in IDENTIFIER_TYPES
                for output_identifier_type in output_identifier_types
            ]
        ):
            raise ValueError(
                f"{output_identifier_types} are invalid output identifier type(s) for ChemSpider"
            )

        query_id = await self.filter_name(
            session, input_identifier.value, order=None, direction=ASCENDING
        )
        record_ids = await self.filter_results(session, query_id=query_id)

        resolved_identifiers = []
        if CompoundIdentifierType.CHEMSPIDER_ID in output_identifier_types:
            resolved_identifiers += [
                CompoundIdentifier(
                    identifier_type=CompoundIdentifierType.CHEMSPIDER_ID,
                    value=record_id,
                )
                for record_id in record_ids
            ]
        for record_id in record_ids:
            details = await self.get_details(session, record_id)
            if CompoundIdentifierType.SMILES in output_identifier_types:
                resolved_identifiers += [
                    CompoundIdentifier(
                        identifier_type=CompoundIdentifierType.SMILES,
                        value=details.get("smiles"),
                    )
                ]
            if CompoundIdentifierType.INCHI in output_identifier_types:
                mol_2d = details["mol2D"]
                value = await self.convert(session, mol_2d, "Mol", "InChI")
                resolved_identifiers += [
                    CompoundIdentifier(
                        value=value, identifier_type=CompoundIdentifierType.INCHI
                    )
                ]
            if CompoundIdentifierType.INCHI_KEY in output_identifier_types:
                mol_2d = details["mol2D"]
                value = await self.convert(session, mol_2d, "Mol", "InChIKey")
                resolved_identifiers += [
                    CompoundIdentifier(
                        identifier_type=CompoundIdentifierType.INCHI_KEY, value=value
                    )
                ]

        return resolved_identifiers

    async def request(
        self,
        session: ClientSession,
        method: str,
        api: str,
        namespace: str,
        endpoint: str,
        params: dict = None,
        json: dict = None,
    ):

        # if params is None:
        #     params = {}
        # if json is None:
        #     json = {}

        # Construct request URL
        url = "{}/{}/{}/{}/{}".format(
            self.api_url, api, self.api_version, namespace, endpoint
        )

        # Set apikey header
        headers = {
            "apikey": self.token,
        }

        logger.debug("{} : {} : {} : {}".format(url, headers, params, json))

        # Make request
        async with session.request(
            method=method, url=url, json=json, params=params, headers=headers
        ) as resp:
            logger.debug(f"Response status: {resp.status}")
            response = await resp.json()

        return response

    async def filter_name(
        self, session: ClientSession, input_identifier, order=None, direction=None
    ):
        json = {
            "name": input_identifier,
            "orderBy": ORDERS.get(order),
            "orderDirection": DIRECTIONS.get(direction),
        }
        response = await self.request(
            session,
            method="post",
            api="compounds",
            namespace="filter",
            endpoint="name",
            json=json,
        )
        return response["queryId"]

    async def filter_results(self, session: ClientSession, query_id):
        """Get filter results using a query ID that was returned by a previous filter request.

        :param string query_id: Query ID from a previous filter request.
        :param int start: Zero-based results offset.
        :param int count: Number of results to return.
        :return: List of results.
        :rtype: list[int]
        """
        endpoint = "{}/results".format(query_id)
        response = await self.request(
            session,
            method="get",
            api="compounds",
            namespace="filter",
            endpoint=endpoint,
        )
        return response["results"]

    async def convert(self, session: ClientSession, input, input_format, output_format):
        """Convert a chemical from one format to another.
        Format: ``SMILES``, ``InChI``, ``InChIKey`` or ``Mol``.
        Allowed conversions: from InChI to InChIKey, from InChI to Mol file, from InChI to SMILES, from InChIKey to
        InChI, from InChIKey to Mol file, from Mol file to InChI, from Mol file to InChIKey, from SMILES to InChI.
        :param string input: Input chemical.
        :param string input_format: Input format.
        :param string output_format: Output format.
        :return: Input chemical in output format.
        :rtype: string
        """
        json = {
            "input": input,
            "inputFormat": input_format,
            "outputFormat": output_format,
        }
        response = await self.request(
            session,
            method="post",
            api="compounds",
            namespace="tools",
            endpoint="convert",
            json=json,
        )
        return response["output"]

    async def get_details(self, session: ClientSession, record_id, fields=FIELDS):
        """Get details for a compound record.
        The available fields are listed in :data:`~chemspipy.api.FIELDS`.
        :param int record_id: Record ID.
        :param list[string] fields: (Optional) List of fields to include in the result.
        :return: Record details.
        :rtype: dict
        """
        params = {"fields": ",".join(fields)}
        endpoint = "{}/details".format(record_id)
        response = await self.request(
            session,
            method="get",
            api="compounds",
            namespace="records",
            endpoint=endpoint,
            params=params,
        )
        return response
