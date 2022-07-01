"""
Python interface for the Chemical Identifier Resolver (CIR) by the CADD Group at the NCI/NIH.
Based heavily on CIRpy by Matt Swain
https://github.com/mcs07/CIRpy

"""
from .service import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from typing import List
import functools
import inspect
import logging
import os

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


log = logging.getLogger("cirpy")
log.addHandler(logging.NullHandler())

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
        output_identifier_type: CompoundIdentifierType,
    ) -> List[CompoundIdentifier]:
        representation = IDENTIFIER_MAP.get(output_identifier_type)
        if representation is None:
            raise ValueError(
                f"{output_identifier_type} is not one of the valid identifier types for the chemical identifier resolver."
            )

        values = await resolve(
            session, input=input_identifier.value, representation=representation
        )
        if isinstance(values, str):
            values = [values]
        return [
            CompoundIdentifier(identifier_type=output_identifier_type, value=value)
            for value in values
        ]


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
    log.debug("Received %s query results", len(results))
    return results


async def request(
    session,
    input,
    representation,
    resolvers=None,
    get3d=False,
    tautomers=False,
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
    log.debug("Making request: %s", url)
    async with session.get(url) as resp:
        response = await resp.text()
    feed = etree.fromstring(response.encode("ascii"))
    tree = etree.ElementTree(feed)
    return tree.getroot()


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


def resolve_image(
    input,
    resolvers=None,
    fmt="png",
    width=300,
    height=300,
    frame=False,
    crop=None,
    bgcolor=None,
    atomcolor=None,
    hcolor=None,
    bondcolor=None,
    framecolor=None,
    symbolfontsize=11,
    linewidth=2,
    hsymbol="special",
    csymbol="special",
    stereolabels=False,
    stereowedges=True,
    header=None,
    footer=None,
    **kwargs,
):
    """Resolve input to a 2D image depiction.

    :param string input: Chemical identifier to resolve
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param string fmt: (Optional) gif or png image format (default png)

    :param int width: (Optional) Image width in pixels (default 300)
    :param int height: (Optional) Image height in pixels (default 300)
    :param bool frame: (Optional) Whether to show border frame (default False)
    :param int crop: (Optional) Crop image with specified padding

    :param int symbolfontsize: (Optional) Atom label font size (default 11)
    :param int linewidth: (Optional) Bond line width (default 2)

    :param string bgcolor: (Optional) Background color
    :param string atomcolor: (Optional) Atom label color
    :param string hcolor: (Optional) Hydrogen atom label color
    :param string bondcolor: (Optional) Bond color
    :param string framecolor: (Optional) Border frame color

    :param bool hsymbol: (Optional) Hydrogens: all, special or none (default special)
    :param bool csymbol: (Optional) Carbons: all, special or none (default special)
    :param bool stereolabels: (Optional) Whether to show stereochemistry labels (default False)
    :param bool stereowedges: (Optional) Whether to show wedge/dash bonds (default True)
    :param string header: (Optional) Header text above structure
    :param string footer: (Optional) Footer text below structure

    """
    # Aggregate all arguments into kwargs
    args, _, _, values = inspect.getargvalues(inspect.currentframe())
    for arg in args:
        if values[arg] is not None:
            kwargs[arg] = values[arg]
    # Turn off anti-aliasing for transparent background
    if kwargs.get("bgcolor") == "transparent":
        kwargs["antialiasing"] = False
    # Renamed parameters
    if "stereolabels" in kwargs:
        kwargs["showstereo"] = kwargs.pop("stereolabels")
    if "fmt" in kwargs:
        kwargs["format"] = kwargs.pop("fmt")
    # Toggle stereo wedges
    if "stereowedges" in kwargs:
        status = kwargs.pop("stereowedges")
        kwargs.update({"wedges": status, "dashes": status})
    # Constant values
    kwargs.update({"representation": "image", "xml": False})
    url = construct_api_url(**kwargs)
    log.debug("Making image request: %s", url)
    response = urlopen(url)
    return response.read()


def download(
    input,
    filename,
    representation,
    overwrite=False,
    resolvers=None,
    get3d=False,
    **kwargs,
):
    """Convenience function to save a CIR response as a file.

    This is just a simple wrapper around the resolve function.

    :param string input: Chemical identifier to resolve
    :param string filename: File path to save to
    :param string representation: Desired output representation
    :param bool overwrite: (Optional) Whether to allow overwriting of an existing file
    :param list(string) resolvers: (Optional) Ordered list of resolvers to use
    :param bool get3d: (Optional) Whether to return 3D coordinates (where applicable)
    :raises HTTPError: if CIR returns an error code
    :raises ParseError: if CIR response is uninterpretable
    :raises IOError: if overwrite is False and file already exists
    """
    result = resolve(input, representation, resolvers, get3d, **kwargs)
    # Just log and return if nothing resolved
    if not result:
        log.debug("No file to download.")
        return
    # Only overwrite an existing file if explicitly instructed to.
    if not overwrite and os.path.isfile(filename):
        raise IOError(
            "%s already exists. Use 'overwrite=True' to overwrite it." % filename
        )
    # Ensure file ends with a newline
    if not result.endswith("\n"):
        result += "\n"
    with open(filename, "w") as f:
        f.write(result)


def memoized_property(fget):
    """Decorator to create memoized properties."""
    attr_name = "_{0}".format(fget.__name__)

    @functools.wraps(fget)
    def fget_memoized(self):
        if not hasattr(self, attr_name):
            setattr(self, attr_name, fget(self))
        return getattr(self, attr_name)

    return property(fget_memoized)
