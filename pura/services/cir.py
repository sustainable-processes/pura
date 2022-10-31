"""
Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.
Based heavily on CIRpy by Matt Swain
https://github.com/mcs07/CIRpy

"""
from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from pura.utils import inverse_map
from aiohttp import ClientSession
from typing import List, Union
import logging

try:
    from urllib.error import HTTPError
    from urllib.parse import quote, urlencode
    from urllib.request import urlopen
except ImportError:
    from urllib import urlencode
    from urllib2 import quote, urlopen, HTTPError

try:
    from lxml import etree
except ImportError:
    try:
        import xml.etree.cElementTree as etree
    except ImportError:
        import xml.etree.ElementTree as etree


logger = logging.getLogger(__name__)

API_BASE = "https://cactus.nci.nih.gov/chemical/structure"


FILE_FORMATS = {
    "alc",
    "cdxml",
    "cerius",
    "charmm",
    "cif",
    "cml",
    "ctx",
    "gjf",
    "gromacs",
    "hyperchem",
    "jme",
    "maestro",
    "mol",
    "mol2",
    "mrv",
    "pdb",
    "sdf3000",
    "sln",
    "xyz",
}
IDENTIFIER_MAP = {
    CompoundIdentifierType.SMILES: "smiles",
    CompoundIdentifierType.MOLBLOCK: "mol",
    CompoundIdentifierType.INCHI: "stdinchi",
    CompoundIdentifierType.IUPAC_NAME: "iupac_name",
    CompoundIdentifierType.CAS_NUMBER: "cas",
    CompoundIdentifierType.INCHI_KEY: "stdinchikey",
    CompoundIdentifierType.XYZ: "xyz",
}


class CIR(Service):
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[Union[CompoundIdentifierType, None]]:
        representations = [
            IDENTIFIER_MAP.get(output_identifier_type)
            for output_identifier_type in output_identifier_types
        ]
        if not any(representations):
            raise ValueError(
                f"{output_identifier_types} is not one of the valid identifier types for the chemical identifier resolver."
            )

        output_compound_identifiers = []
        for representation, output_identifier_type in zip(
            representations, output_identifier_types
        ):
            if representation is None:
                continue
            values = await resolve(
                session, input=input_identifier.value, representation=representation
            )
            if isinstance(values, str):
                values = [values]

            if len(values) > 0:
                output_compound_identifiers += [
                    CompoundIdentifier(
                        identifier_type=output_identifier_type, value=value
                    )
                    for value in values
                ]

        return output_compound_identifiers


async def resolve(
    session: ClientSession,
    input: str,
    representation: str,
    resolvers: List[str] = None,
    get3d: bool = False,
    **kwargs,
):
    """Resolve input to the specified output representation.

    :param aiohttp.ClientSession() session: Aiohttp session
    :param string input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :returns: Output representation or None
    :rtype: string or None
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    # Take first result from XML query
    results = await query(
        session, input, representation, resolvers, False, get3d, **kwargs
    )
    result = results[0].value if results else None
    if result is None:
        return []
    elif not isinstance(result, list):
        return [result]
    else:
        return result


class Result(object):
    """A single result returned by CIR."""

    def __init__(self, input, notation, input_format, resolver, representation, value):
        """

        :param string input: Originally supplied input identifier that produced this result
        :param string notation: Identifier matched by the resolver or tautomer ID
        :param string input_format: Format of the input as interpreted by the resolver
        :param string resolver: Resolver used to produce this result
        :param string representation: Requested output representation
        :param value: Actual result value
        :type value: string or list(string)
        """
        self.input = input
        self.representation = representation
        self.resolver = resolver
        self.input_format = input_format
        self.notation = notation
        self.value = value

    def __repr__(self):
        return (
            "Result(input=%r, representation=%r, resolver=%r, input_format=%r, notation=%r, value=%r)"
            % (
                self.input,
                self.representation,
                self.resolver,
                self.input_format,
                self.notation,
                self.value,
            )
        )

    def __str__(self):
        return self.value

    def __eq__(self, other):
        return isinstance(other, type(self)) and self.__dict__ == other.__dict__

    def __getitem__(self, prop):
        """Allow dict-style access to attributes to ease transition from when results were dicts."""
        if prop in self.__dict__:
            return getattr(self, prop)
        raise KeyError(prop)

    def __setitem__(self, prop, val):
        """Allow dict-style setting of attributes to ease transition from when results were dicts."""
        setattr(self, prop, val)

    def __contains__(self, prop):
        """Allow dict-style checking of attributes to ease transition from when results were dicts."""
        return prop in self.__dict__

    def to_dict(self):
        """Return a dictionary containing Result data."""
        return self.__dict__


async def query(
    session,
    input,
    representation,
    resolvers=None,
    get3d=False,
    tautomers=False,
    **kwargs,
):
    """Get all results for resolving input to the specified output representation.

    :param aiohttp.ClientSession() session: Aiohttp session
    :param string input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :returns: List of resolved results
    :rtype: list(Result)
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    tree = await request(
        session, input, representation, resolvers, get3d, tautomers, **kwargs
    )

    results = []
    if tree is None:
        return results
    for data in tree.findall(".//data"):
        value = [item.text for item in data.findall("item")]
        result = Result(
            input=tree.attrib["string"],
            representation=tree.attrib["representation"],
            resolver=data.attrib["resolver"],
            input_format=data.attrib["string_class"],
            notation=data.attrib["notation"],
            value=value[0] if len(value) == 1 else value,
        )
        results.append(result)
    logger.debug("Received %s query results", len(results))
    return results


async def request(
    session: ClientSession,
    input: str,
    representation: str,
    resolvers: List[str] = None,
    get3d: bool = False,
    tautomers: bool = False,
    **kwargs,
):
    """Make a request to CIR and return the XML response.

    :param aiohttp.ClientSession() session: Aiohttp session
    :param string input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :returns: XML response from CIR
    :rtype: Element
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    """
    url = construct_api_url(
        input, representation, resolvers, get3d, tautomers, **kwargs
    )
    logger.debug("Making request: %s", url)
    async with session.get(url) as resp:
        response = await resp.text()
    try:
        feed = etree.fromstring(response.encode("ascii"))
        tree = etree.ElementTree(feed)
        return tree.getroot()
    except (UnicodeDecodeError, UnicodeEncodeError) as e:
        logger.error(e)
        return


def construct_api_url(
    input,
    representation,
    resolvers=None,
    get3d=False,
    tautomers=False,
    xml=True,
    **kwargs,
):
    """Return the URL for the desired API endpoint.

    :param string input: Chemical identifier to resolve
    :param string representation: Desired output representation
    :param list(str) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :param bool tautomers: (Optional) Whether to return all tautomers
    :param bool xml: (Optional) Whether to return full XML response
    :returns: CIR API URL
    :rtype: str
    """
    # File formats require representation=file and the format in the querystring
    if representation in FILE_FORMATS:
        kwargs["format"] = representation
        representation = "file"
    # Prepend input with 'tautomers:' to return all tautomers
    if tautomers:
        input = "tautomers:%s" % input
    url = "%s/%s/%s" % (API_BASE, quote(input), representation)
    if xml:
        url += "/xml"
    if resolvers:
        kwargs["resolver"] = ",".join(resolvers)
    if get3d:
        kwargs["get3d"] = True
    if kwargs:
        url += "?%s" % urlencode(kwargs)
    return url
