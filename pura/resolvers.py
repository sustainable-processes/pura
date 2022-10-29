from pura.compound import (
    CompoundIdentifier,
    CompoundIdentifierType,
    Compound,
    standardize_identifier,
    unique_identifiers,
)
from pura.services import *
from tqdm import tqdm
from aiohttp import *
from aiohttp.web_exceptions import (
    HTTPClientError,
    HTTPServerError,
    HTTPServiceUnavailable,
)
import asyncio
from typing import Optional, List, Union, Tuple, Dict, Callable
from itertools import combinations
from functools import reduce
import logging
import queue

logger = logging.getLogger(__name__)

aiohttp_errors = (
    HTTPServiceUnavailable,
    ClientConnectionError,
    TimeoutError,
    ClientConnectorCertificateError,
    ClientConnectorError,
    ClientConnectorSSLError,
    ClientError,
    ClientHttpProxyError,
    ClientOSError,
    ClientPayloadError,
    ClientProxyConnectionError,
    ClientResponseError,
    ClientSSLError,
    ContentTypeError,
    InvalidURL,
    ServerConnectionError,
    ServerDisconnectedError,
    ServerFingerprintMismatch,
    ServerTimeoutError,
    WSServerHandshakeError,
    asyncio.TimeoutError,
)


def base_check_agreement(
    identifiers_list: List[List[CompoundIdentifier]],
    agreement: int,
    service: Service,
) -> Tuple[List[List[CompoundIdentifier]], bool]:
    """
    Check if identifier lists from multiple services agree.

    Parameters
    ---------

    identifiers_list :  list of list of :class:`~pura.compound.CompoundIdentifier`
        List of lists of CompoudIdentifier objects. Each item in the first list is
        considered results from a  different service.
    agreement : int
        The number of services that must agree.
    service : Service
        The service used at the time of calling this method. Can be used
        by methods that override for further filtering.

    Returns
    -------

    A tuple of the following two objects:
    1. list of lists of CompoundIdentifier objects representing the unique identifiers
    2. boolean that is true if there is sufficient agreement between services

    Notes
    ------

    Algorithm:
    1. Find all combinatons of services that can satisfy agreement (combinations)
    2. Find the intersection of each combination
    3. If the intersection is greater than zero, then you have sufficient agreement.

    """
    if agreement <= 0:
        raise ValueError("Agreement must be greater than 0.")
    identifiers_list_new = []
    for identifiers in identifiers_list:
        if len(identifiers) > 0:
            identifiers_list_new.append(
                [identifier.value for identifier in identifiers]
            )
            identifier_type = identifiers[0].identifier_type
    identifiers_sets = [set(ident) for ident in identifiers_list_new]
    options = list(range(len(identifiers_list_new)))
    intersection = []
    for combo in combinations(options, agreement):
        intersection = reduce(
            set.intersection, [identifiers_sets[combo_i] for combo_i in combo]
        )
        if len(intersection) > 0:
            break
    if len(intersection) == 0:
        return identifiers_list, False
    else:
        return [
            [
                CompoundIdentifier(identifier_type=identifier_type, value=identifier)
                for identifier in intersection
            ]
        ], True


class CompoundResolver:
    """Resolve compound identifier types using external services such as PubChem.

    Parameters
    ----------

    services : list of :class:`~pura.service.Service`
        The services used for resolution.
    silent : bool, optional
        If True, logs errors but does not raise them. Default is False
    agreement_check : callable, optional
        A function that checks for agreement. See :class:`~pura.resolvers.base_check_agreement`
        for the API and the default agreement_check function.

    Examples
    --------

    The example below demonstrates how multiple identifiers of the same compound
    can be used to help increase the chances of resolution.  Each pairing of input identifier
    and service will be tried until sufficient agreement is reached or the pairings are exhausted.

    >>> services = [CAS(), CIR(), Opsin(), PubChem()]
    >>> compounds = [
    ...    Compound(
    ...            identifiers=[
    ...                CompoundIdentifier(
    ...                    identifier_type=CompoundIdentifierType.CAS_NUMBER,
    ...                    value="6737-11-7",
    ...                ),
    ...                CompoundIdentifier(
    ...                    identifier_type=CompoundIdentifierType.NAME,
    ...                    value="i-Butyric acid 1,2,3-propanetriyl ester",
    ...                ),
    ...            ]
    ...        )
    ... ]
    >>> resolver = CompoundResolver(services=services, silent=True)
    >>> compound_identifiers_list = resolver.resolve(
    ...        input_compounds=compounds,
    ...        output_identifier_type=CompoundIdentifierType.SMILES,
    ...        agreement=2,
    ...        batch_size=1,
    ... )
    >>> print(compound_identifiers_list)

    """

    def __init__(
        self,
        services: List[Service],
        silent: Optional[bool] = False,
        agreement_check: Optional[Callable] = None,
    ):
        self._services = services
        self.silent = silent
        self.agreement_check = (
            agreement_check if agreement_check is not None else base_check_agreement
        )

    def resolve(
        self,
        input_compounds: List[Compound],
        output_identifier_type: CompoundIdentifierType,
        backup_identifier_types: Optional[List[CompoundIdentifierType]] = None,
        agreement: Optional[int] = 1,
        batch_size: Optional[int] = None,
        n_retries: Optional[int] = 3,
    ) -> List[Tuple[Compound, Union[List[CompoundIdentifier], None]]]:
        """Resolve a list of compound identifiers to another identifier type(s).

        Parameters
        ---------

        input_compounds : list of :class:`~pura.compund.Compound`
            The list of Compounds that should be resolved.
        output_identifiers_types : list of :class:`~pura.compund.CompoundIdentifier`
            The compound identifier type that should be outputted.
        backup_identifier_types : list of :class:`~pura.compound.CompoundIdentifierType`
            A list of identifier types that can be looked up and then used to resolve
            to the desired output identifier type. Default is None.
        agreement : int, optional
            The number of services that must give the same resolved
            CompoundIdentifier for the resolution to be considered correct.
            If set to zero, then all services will be tried and the unique
            results returned. Default is 1.
        batch_size : int, optional
            The batch size sets the number of requests to send simultaneously.
            Defaults to 10 or the length input_idententifier, whichever is smaller.
        n_retries : int, optional
            The number of times a request should be retried if there is a problem
            in contacting a service (e.g., an internet outage). Defaults to 3.


        Notes
        -----

        The retries option uses an exponential backoff, so a sleep of 2^n seconds will occur before
        the next retry, where n is the number of retries thus far.

        Returns
        -------

        A list of tuples where the first element of each tuple is the input compound and the second
        element is a list of compound identifier(s).

        """
        # Make sure output identifier type is different than backup identifier types
        if (
            backup_identifier_types is not None
            and output_identifier_type in backup_identifier_types
        ):
            raise ValueError(
                "Output identifier type cannot be in backup identifier types."
            )
        try:
            loop = asyncio.get_event_loop()
        except RuntimeError as e:
            if str(e).startswith("There is no current event loop in thread"):
                loop = asyncio.new_event_loop()
                asyncio.set_event_loop(loop)
            else:
                raise
        return loop.run_until_complete(
            self._resolve(
                input_compounds=input_compounds,
                output_identifier_type=output_identifier_type,
                backup_identifier_types=backup_identifier_types,
                agreement=agreement,
                batch_size=batch_size,
                n_retries=n_retries,
            )
        )

    async def _resolve(
        self,
        input_compounds: List[Compound],
        output_identifier_type: CompoundIdentifierType,
        backup_identifier_types: Optional[List[CompoundIdentifierType]] = None,
        agreement: Optional[int] = 1,
        batch_size: Optional[int] = None,
        n_retries: Optional[int] = 3,
    ) -> List[Tuple[Compound, Union[List[CompoundIdentifier], None]]]:
        """This is the async function with the same API as resolve"""

        n_identifiers = len(input_compounds)
        if batch_size is None:
            batch_size = 10 if n_identifiers >= 10 else n_identifiers
        n_batches = n_identifiers // batch_size
        n_batches += 0 if n_identifiers % batch_size == 0 else 1
        resolved_identifiers = []
        backup_identifier_types = (
            backup_identifier_types if backup_identifier_types is not None else []
        )
        # Iterate through batches
        for batch in tqdm(range(n_batches), position=0, desc="Batch"):
            # Get subset of data
            start = batch * batch_size
            batch_identifiers = input_compounds[start : start + batch_size]

            # Start aiohttp session
            async with ClientSession() as session:
                # Create series of tasks to run in parallel
                tasks = [
                    self._resolve_one_compound(
                        session,
                        compound_identifier,
                        output_identifier_type,
                        backup_identifier_types,
                        agreement,
                        n_retries=n_retries,
                    )
                    for compound_identifier in batch_identifiers
                ]
                batch_bar = tqdm(
                    asyncio.as_completed(tasks),
                    total=len(tasks),
                    desc=f"Batch {batch} Progress",
                    position=1,
                    leave=True,
                )
                resolved_identifiers.extend([await f for f in batch_bar])
                batch_bar.clear()

        return resolved_identifiers

    async def _resolve_one_compound(
        self,
        session: ClientSession,
        input_compound: Compound,
        output_identifier_type: CompoundIdentifierType,
        backup_identifier_types: List[CompoundIdentifierType],
        agreement: int,
        n_retries: Optional[int],
    ) -> Tuple[Compound, Union[List[CompoundIdentifier], None]]:
        """Resolve one compound"""
        resolved_identifiers_list = []
        agreement_satisfied = False
        # Create input identifier queue
        input_identifiers_queue = queue.Queue()
        for identifier in input_compound.identifiers:
            input_identifiers_queue.put(identifier)

        # Main loop
        while not input_identifiers_queue.empty() and not agreement_satisfied:
            input_identifier = input_identifiers_queue.get()
            for service in self._services:
                for j in range(n_retries):
                    try:
                        output_identifier_types = [
                            output_identifier_type
                        ] + backup_identifier_types
                        if input_identifier.identifier_type in output_identifier_types:
                            output_identifier_types.remove(
                                input_identifier.identifier_type
                            )
                        resolved_identifiers = await service.resolve_compound(
                            session,
                            input_identifier=input_identifier,
                            output_identifier_types=output_identifier_types,
                        )
                        logger.debug(
                            f"{service} | {input_identifier.value}->{resolved_identifiers}"
                        )
                        # Standardize identifiers (e.g., SMILES canonicalization)
                        for identifier in resolved_identifiers:
                            if identifier is not None:
                                standardize_identifier(identifier)
                        # Add backup identifiers to queue
                        for identifier in resolved_identifiers:
                            if (
                                identifier.identifier_type in backup_identifier_types
                                and input_identifier.identifier_type
                                not in backup_identifier_types
                            ):
                                input_identifiers_queue.put(identifier)
                                resolved_identifiers.remove(identifier)
                        resolved_identifiers_list.append(resolved_identifiers)
                        break
                    except aiohttp_errors as e:  # type: ignore
                        # If server is busy, use exponential backoff
                        logger.debug(f"Sleeping for {2**j} ({e})")
                        await asyncio.sleep(2**j)
                    except (HTTPClientError, HTTPServerError) as e:
                        # Log/raise on all other HTTP errors
                        if self.silent:
                            logger.error(msg=f"{service}, {input_identifier}: {e}")
                            break
                        else:
                            raise e
                    except ValueError as e:
                        if self.silent:
                            logger.error(e)
                            break
                        else:
                            raise e

                # Check agreement between services
                if len(resolved_identifiers_list) > 0 and agreement > 1:
                    (
                        resolved_identifiers_list,
                        agreement_satisfied,
                    ) = self.agreement_check(
                        resolved_identifiers_list,
                        agreement,
                        service,
                    )
                elif (
                    agreement == 1
                    and len(resolved_identifiers_list) > 0
                    and len(resolved_identifiers_list[0]) > 0
                ):
                    agreement_satisfied = True

                if agreement_satisfied:
                    break

        # If agreement is 0, then we just want to dedup and return
        if agreement == 0:
            resolved_identifiers_list = unique_identifiers(resolved_identifiers_list)
            agreement_satisfied = True

        if not agreement_satisfied:
            error_txt = f"Not sufficient agreement for {input_identifier} (outputs: {resolved_identifiers_list})"
            if self.silent:
                logger.error(error_txt)
                resolved_identifiers_list = []
            else:
                raise ResolverError(error_txt)
        return input_compound, flatten_list(resolved_identifiers_list)


def flatten_list(l: List):
    """Flatten a multidimensional list to one dimension"""
    lnew = []
    for li in l:
        if type(li) != list:
            lnew.append(li)
        else:
            lnew.extend(flatten_list(li))
    return lnew


class ResolverError(Exception):
    """General error problems with resolving names."""

    pass


def resolve_identifiers(
    names: List[str],
    output_identifier_type: CompoundIdentifierType,
    input_identifer_type: CompoundIdentifierType = CompoundIdentifierType.NAME,
    backup_identifier_types: Optional[List[CompoundIdentifierType]] = None,
    agreement: int = 1,
    batch_size: int = 100,
    services: Optional[List[Service]] = None,
    silent: Optional[bool] = False,
) -> List[Tuple[Compound, Union[List[CompoundIdentifier], None]]]:
    """Resolve a list of names (or any other identifier) to an identifier type.

    Parameters
    -----------

    names : list of str
        The list of compound names that should be resolved
    output_identifier_type : :class:`~pura.compund.CompoundIdentifierType`
        The list of compound identifier types to resolve to.
    input_identifier_type : :class:`~pura.compund.CompoundIdentifierType`, optional
        The input identifier type, Defaults to name
    backup_identifier_types : list of :class:`~pura.compound.CompoundIdentifierType`
        A list of identifier types that can be looked up and then used to resolve
        to the desired output identifier type. Default is None.
    agreement : int, optional
        The number of services that must give the same resolved
        `CompoundIdentifier` for the resolution to be considered correct.
        If set to zero, then all services will be tried and the unique
        results returned. Default is 1.
    batch_size : int, optional
        The batch size sets the number of requests to send simultaneously.
        Defaults to 100 or the length input_idententifier, whichever is smaller.
    services : list of :class:`~pura.service.Service`, optional
        Services used to do resolution. Defaults to PubChem and CIR.
    silent : bool, optional
        If True, logs errors but does not raise them. Default is False

    Example
    -------

    >>> from pura.services import  Pubchem, CIR
    >>> smiles = resolve_identifiers(
    ...     ["aspirin", "ibuprofen", "toluene"],
    ...     output_identifier_type=CompoundIdentifierType.SMILES,
    ...     services=[Pubchem(), CIR()],
    ...     agreement=2,
    ... )

    Notes
    -----

    This is a convenience function for quickly resolving a list of strings without having
    to create Compound objects.

    """
    if services is None:
        services = [PubChem(), CIR()]
    compounds = [
        Compound(
            identifiers=[
                CompoundIdentifier(identifier_type=input_identifer_type, value=name)
            ]
        )
        for name in names
    ]
    resolver = CompoundResolver(services=services, silent=silent)
    return resolver.resolve(
        input_compounds=compounds,
        output_identifier_type=output_identifier_type,
        backup_identifier_types=backup_identifier_types,
        agreement=agreement,
        batch_size=batch_size,
    )
