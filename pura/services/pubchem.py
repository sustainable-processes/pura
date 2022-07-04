# -*- coding: utf-8 -*-
"""
PubChemPy

Python interface for the PubChem PUG REST service.
https://github.com/mcs07/PubChemPy
"""


from pura.services import Service
from pura.compound import Compound, CompoundIdentifier, CompoundIdentifierType
import json
import logging
import sys
from aiohttp import ClientSession
from typing import List


log = logging.getLogger("pubchempy")
log.addHandler(logging.NullHandler())


text_types = str, bytes


API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"

IDENTIFIER_MAP = {
    CompoundIdentifierType.SMILES: "CanonicalSMILES",
    CompoundIdentifierType.MOLBLOCK: "mol",
    CompoundIdentifierType.INCHI: "stdinchi",
    CompoundIdentifierType.IUPAC_NAME: "iupac_name",
    CompoundIdentifierType.CAS_NUMBER: "cas",
    CompoundIdentifierType.INCHI_KEY: "stdinchikey",
    CompoundIdentifierType.XYZ: "xyz",
    CompoundIdentifierType.NAME: "name",
}


class Pubchem(Service):
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_type: CompoundIdentifierType,
    ) -> List[CompoundIdentifierType]:

        namespace = IDENTIFIER_MAP.get(input_identifier.identifier_type)
        if namespace is None:
            raise ValueError(
                f"{input_identifier.identifier_type} is not one of the valid identifier types for the Pubchem."
            )

        representation = IDENTIFIER_MAP.get(output_identifier_type)
        if representation is None:
            raise ValueError(
                f"{output_identifier_type} is not one of the valid identifier types for the chemical identifier resolver."
            )

        results = await get_properties(
            session, representation, input_identifier.value, namespace, searchtype=None
        )
        return [
            CompoundIdentifier(
                identifier_type=output_identifier_type,
                value=result[representation],
            )
            for result in results
        ]


class PubChemPyDeprecationWarning(Warning):
    """Warning category for deprecated features."""

    pass


class PubChemPyError(Exception):
    """Base class for all PubChemPy exceptions."""

    pass


class ResponseParseError(PubChemPyError):
    """PubChem response is uninterpretable."""

    pass


class PubChemHTTPError(PubChemPyError):
    """Generic error class to handle all HTTP error codes."""

    def __init__(self, e):
        self.code = e.code
        self.msg = e.reason
        try:
            self.msg += ": %s" % json.loads(e.read().decode())["Fault"]["Details"][0]
        except (ValueError, IndexError, KeyError):
            pass
        if self.code == 400:
            raise BadRequestError(self.msg)
        elif self.code == 404:
            raise NotFoundError(self.msg)
        elif self.code == 405:
            raise MethodNotAllowedError(self.msg)
        elif self.code == 504:
            raise TimeoutError(self.msg)
        elif self.code == 501:
            raise UnimplementedError(self.msg)
        elif self.code == 500:
            raise ServerError(self.msg)

    def __str__(self):
        return repr(self.msg)


class BadRequestError(PubChemHTTPError):
    """Request is improperly formed (syntax error in the URL, POST body, etc.)."""

    def __init__(self, msg="Request is improperly formed"):
        self.msg = msg


class NotFoundError(PubChemHTTPError):
    """The input record was not found (e.g. invalid CID)."""

    def __init__(self, msg="The input record was not found"):
        self.msg = msg


class MethodNotAllowedError(PubChemHTTPError):
    """Request not allowed (such as invalid MIME type in the HTTP Accept header)."""

    def __init__(self, msg="Request not allowed"):
        self.msg = msg


class TimeoutError(PubChemHTTPError):
    """The request timed out, from server overload or too broad a request.

    See :ref:`Avoiding TimeoutError <avoiding_timeouterror>` for more information.
    """

    def __init__(self, msg="The request timed out"):
        self.msg = msg


class UnimplementedError(PubChemHTTPError):
    """The requested operation has not (yet) been implemented by the server."""

    def __init__(self, msg="The requested operation has not been implemented"):
        self.msg = msg


class ServerError(PubChemHTTPError):
    """Some problem on the server side (such as a database server down, etc.)."""

    def __init__(self, msg="Some problem on the server side"):
        self.msg = msg


async def request(
    session: ClientSession,
    identifier,
    namespace="cid",
    domain="compound",
    operation=None,
    output="JSON",
    searchtype=None,
    **kwargs,
):
    """
    Construct API request from parameters and return the response.

    Full specification at http://pubchem.ncbi.nlm.nih.gov/pug_rest/PUG_REST.html
    """
    if not identifier:
        raise ValueError("identifier/cid cannot be None")
    # If identifier is a list, join with commas into string
    if isinstance(identifier, int):
        identifier = str(identifier)
    if not isinstance(identifier, text_types):
        identifier = ",".join(str(x) for x in identifier)
    # Filter None values from kwargs
    kwargs = dict((k, v) for k, v in kwargs.items() if v is not None)
    # Build API URL
    urlid, postdata = None, {}
    if namespace == "sourceid":
        identifier = identifier.replace("/", ".")
    if (
        namespace in ["listkey", "formula", "sourceid"]
        or searchtype == "xref"
        or (searchtype and namespace == "cid")
        or domain == "sources"
    ):
        urlid = quote(identifier.encode("utf8"))
    else:
        # postdata = urlencode([(namespace, identifier)]).encode("utf8")
        postdata = {namespace: identifier}
    comps = filter(
        None, [API_BASE, domain, searchtype, namespace, urlid, operation, output]
    )
    apiurl = "/".join(comps)
    if kwargs:
        apiurl += "?%s" % urlencode(kwargs)

    # Make request
    log.debug("Request URL: %s", apiurl)
    log.debug("Request data: %s", postdata)
    async with session.post(apiurl, data=postdata) as resp:
        response = await resp.json()
    return response


async def get_json(
    session: ClientSession,
    identifier,
    namespace="cid",
    domain="compound",
    operation=None,
    searchtype=None,
    **kwargs,
):
    """Request wrapper that automatically parses JSON response and supresses NotFoundError."""
    try:
        return await request(
            session=session,
            identifier=identifier,
            namespace=namespace,
            domain=domain,
            operation=operation,
            output="JSON",
            searchtype=searchtype,
            **kwargs,
        )

    except NotFoundError as e:
        log.info(e)
        return None


# Allows properties to optionally be specified as underscore_separated, consistent with Compound attributes
PROPERTY_MAP = {
    "molecular_formula": "MolecularFormula",
    "molecular_weight": "MolecularWeight",
    "canonical_smiles": "CanonicalSMILES",
    "isomeric_smiles": "IsomericSMILES",
    "inchi": "InChI",
    "inchikey": "InChIKey",
    "iupac_name": "IUPACName",
    "xlogp": "XLogP",
    "exact_mass": "ExactMass",
    "monoisotopic_mass": "MonoisotopicMass",
    "tpsa": "TPSA",
    "complexity": "Complexity",
    "charge": "Charge",
    "h_bond_donor_count": "HBondDonorCount",
    "h_bond_acceptor_count": "HBondAcceptorCount",
    "rotatable_bond_count": "RotatableBondCount",
    "heavy_atom_count": "HeavyAtomCount",
    "isotope_atom_count": "IsotopeAtomCount",
    "atom_stereo_count": "AtomStereoCount",
    "defined_atom_stereo_count": "DefinedAtomStereoCount",
    "undefined_atom_stereo_count": "UndefinedAtomStereoCount",
    "bond_stereo_count": "BondStereoCount",
    "defined_bond_stereo_count": "DefinedBondStereoCount",
    "undefined_bond_stereo_count": "UndefinedBondStereoCount",
    "covalent_unit_count": "CovalentUnitCount",
    "volume_3d": "Volume3D",
    "conformer_rmsd_3d": "ConformerModelRMSD3D",
    "conformer_model_rmsd_3d": "ConformerModelRMSD3D",
    "x_steric_quadrupole_3d": "XStericQuadrupole3D",
    "y_steric_quadrupole_3d": "YStericQuadrupole3D",
    "z_steric_quadrupole_3d": "ZStericQuadrupole3D",
    "feature_count_3d": "FeatureCount3D",
    "feature_acceptor_count_3d": "FeatureAcceptorCount3D",
    "feature_donor_count_3d": "FeatureDonorCount3D",
    "feature_anion_count_3d": "FeatureAnionCount3D",
    "feature_cation_count_3d": "FeatureCationCount3D",
    "feature_ring_count_3d": "FeatureRingCount3D",
    "feature_hydrophobe_count_3d": "FeatureHydrophobeCount3D",
    "effective_rotor_count_3d": "EffectiveRotorCount3D",
    "conformer_count_3d": "ConformerCount3D",
}


async def get_properties(
    session: ClientSession,
    properties,
    identifier,
    namespace="cid",
    searchtype=None,
    **kwargs,
):
    """Retrieve the specified properties from PubChem.

    :param identifier: The compound, substance or assay identifier to use as a search query.
    :param namespace: (optional) The identifier type.
    :param searchtype: (optional) The advanced search type, one of substructure, superstructure or similarity.
    :param as_dataframe: (optional) Automatically extract the properties into a pandas :class:`~pandas.DataFrame`.
    """
    if isinstance(properties, text_types):
        properties = properties.split(",")
    properties = ",".join([PROPERTY_MAP.get(p, p) for p in properties])
    properties = "property/%s" % properties

    results = await get_json(
        session=session,
        identifier=identifier,
        namespace=namespace,
        domain="compound",
        operation=properties,
        searchtype=searchtype,
        **kwargs,
    )

    results = results["PropertyTable"]["Properties"] if results else []

    return results
