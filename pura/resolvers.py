from multiprocessing.sharedctypes import Value
from pura.compound import (
    CompoundIdentifier,
    CompoundIdentifierType,
    Compound,
    standardize_identifier,
)
from pura.services import Service, CIR, PubChem, ChemSpider, Opsin
from tqdm import tqdm
from aiohttp import *
from aiohttp.web_exceptions import (
    HTTPClientError,
    HTTPServerError,
    HTTPServiceUnavailable,
)
import asyncio
from typing import Optional, List, Union, Tuple, Dict
from itertools import combinations
from functools import reduce
import logging

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


class CompoundResolver:
    """Resolve compound identifier types using external services such as PubChem.

    Parameters
    ----------
    services : list of Service
        The services used for resolution.
    silent : bool, optional
        If True, logs errors but does not raise them. Default is False

    Examples
    --------


    """

    def __init__(self, services: List[Service], silent: Optional[bool] = False):
        self._services = services
        self.silent = silent

    def resolve(
        self,
        input_identifiers: List[CompoundIdentifier],
        output_identifier_type: List[CompoundIdentifierType],
        agreement: Optional[int] = 1,
        batch_size: Optional[int] = None,
    ) -> List[List[CompoundIdentifier]]:
        """Resolve a list of compound identifiers to another identifier type(s).

        Arguments
        ---------
        input_identifiers : List[CompoundIdentifier]
            The list of compound identifiers that should be resolved
        output_identifiers_types: List[CompoundIdentifierType]
            The list of compound identifier types to resolve to
        agreement : int, optional
            The number of services that must give the same resolved
            compoundidentifier for the resolution to be considered correct.
            Default is 1.
        batch_size : int, optional
            The batch size sets the number of requests to send simultaneously.
            Defaults to 100 or the length input_idententifier, whichever is smaller.

        Returns
        -------
        A list of tuples where the first element of each tuple is the input identifier and the second
        element is the output identifier.

        """
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
                input_identifiers=input_identifiers,
                output_identifier_type=output_identifier_type,
                agreement=agreement,
                batch_size=batch_size,
            )
        )

    async def _resolve(
        self,
        input_identifiers: List[CompoundIdentifier],
        output_identifier_type: CompoundIdentifierType,
        agreement: Optional[int] = 1,
        batch_size: Optional[int] = None,
    ) -> List[Tuple[CompoundIdentifier, Union[List[CompoundIdentifier], None]]]:
        """Resolve a list of compound identifiers to another identifier type(s).

        Arguments
        ---------
        input_identifiers : List[CompoundIdentifier]
            The list of compound identifiers that should be resolved
        output_identifiers_types: List[CompoundIdentifierType]
            The list of compound identifier types to resolve to
        agreement : int, optional
            The number of services that must give the same resolved
            compoundidentifier for the resolution to be considered correct.
            Default is 1.
        batch_size : int, optional
            The batch size sets the number of requests to send simultaneously.
            Defaults to 100 or the length input_idententifier, whichever is smaller.

        Returns
        -------
        A list of tuples where the first element of each tuple is the input identifier and the second
        element is the output identifier.

        """

        n_identifiers = len(input_identifiers)
        if batch_size is None:
            batch_size = 100 if n_identifiers >= 100 else n_identifiers
        n_batches = n_identifiers // batch_size
        n_batches += 0 if n_identifiers % batch_size == 0 else 1
        resolved_identifiers = []
        # Iterate through batches
        for batch in tqdm(range(n_batches), position=0, desc="Batch"):
            # Get subset of data
            start = batch * batch_size
            batch_identifiers = input_identifiers[start : start + batch_size]

            # Start aiohttp session
            async with ClientSession() as session:
                # Create series of tasks to run in parallel
                tasks = [
                    self.resolve_one_identifier(
                        session,
                        compound_identifier,
                        output_identifier_type,
                        agreement,
                        n_retries=7,
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

    async def resolve_one_identifier(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_type: CompoundIdentifierType,
        agreement: int,
        n_retries: Optional[int] = 7,
    ) -> Tuple[CompoundIdentifier, Union[List[CompoundIdentifier], None]]:

        resolved_identifiers_list = []
        agreement_satisfied = False
        for i, service in enumerate(self._services):
            for j in range(n_retries):
                try:
                    resolved_identifiers = await service.resolve_compound(
                        session,
                        input_identifier=input_identifier,
                        output_identifier_type=output_identifier_type,
                    )
                    logger.debug(
                        f"{service} | {input_identifier.value}-> {resolved_identifiers}"
                    )
                    # Standardize identifiers (e.g., SMILES canonicalization)
                    for identifier in resolved_identifiers:
                        if identifier is not None:
                            standardize_identifier(identifier)
                    resolved_identifiers_list.append(resolved_identifiers)

                    break
                except aiohttp_errors as e:
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

            # Chceck agreement between services
            if i > 0:
                resolved_identifiers_list, agreement_satisfied = self.check_agreement(
                    resolved_identifiers_list, agreement
                )

            if agreement_satisfied:
                break

        if not agreement_satisfied:
            error_txt = f"Not sufficient agreement for {input_identifier} (outputs: {resolved_identifiers_list})"
            if self.silent:
                logger.error(error_txt)
                resolved_identifiers_list = []
            else:
                raise ResolverError(error_txt)

        return input_identifier, resolved_identifiers_list

    def check_agreement(
        self,
        identifiers_list: List[List[Union[CompoundIdentifier, None]]],
        agreement: int,
    ) -> Tuple[List[CompoundIdentifier], bool]:
        """
        Check if identifier lists from multiple services agree.

        Arguments
        ---------
        identifiers_list
            List of lists of CompoudIdentifier objects. Each item in the first list is
            considered results from a  different service.
        agreement : int
            The number of services that must agreee.
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
                CompoundIdentifier(identifier_type=identifier_type, value=identifier)
                for identifier in intersection
            ], True


class ResolverError(Exception):
    pass


def resolve_names(
    names: List[str],
    output_identifier_type: CompoundIdentifierType,
    agreement: int = 1,
    batch_size: int = 100,
    services: Optional[List[Service]] = None,
    silent: Optional[bool] = False,
) -> List[CompoundIdentifier]:
    """Resolve a list of names to an identifier type.

    Arguments
    ---------
    names : list of str
        The list of compound names that should be resolved
    output_identifiers_type : CompoundIdentifierType
        The list of compound identifier types to resolve to
    agreement : int, optional
        The number of services that must give the same resolved
        compoundidentifier for the resolution to be considered correct.
        Default is 1.
    batch_size : int, optional
        The batch size sets the number of requests to send simultaneously.
        Defaults to 100 or the length input_idententifier, whichever is smaller.
    services : list of `Service`, optional
        Services used to do resolution. Defaults to PubChem and CIR.
    silent : bool, optional
        If True, logs errors but does not raise them. Default is False

    Example
    -------
    >>> from pura.services import  Pubchem, CIR
    >>> smiles = resolve_names(
    ...     ["aspirin", "ibuprofen", "toluene"],
    ...     output_identifier_type=CompoundIdentifierType.SMILES,
    ...     services=[Pubchem(), CIR()],
    ...     agreement=2,
    ... )

    """
    if services is None:
        services = [PubChem(), CIR()]
    name_identifiers = [
        CompoundIdentifier(identifier_type=CompoundIdentifierType.NAME, value=name)
        for name in names
    ]
    resolver = CompoundResolver(services=services, silent=silent)
    return resolver.resolve(
        input_identifiers=name_identifiers,
        output_identifier_type=output_identifier_type,
        agreement=agreement,
        batch_size=batch_size,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    compound_identifiers_list = resolve_names(
        ["aspirin", "ibuprofen", "[Ru(p-cymene)I2]2", "hay", "2,4,6-trinitrotoluene"],
        output_identifier_type=CompoundIdentifierType.SMILES,
        services=[PubChem(), CIR(), Opsin()],
        agreement=1,
        batch_size=10,
        silent=True,
    )
    print(compound_identifiers_list)
    resolved = [
        {"name": input_identifier.value, "smiles": output_identifiers[0].value}
        if len(output_identifiers) > 0
        else {"name": input_identifier.value, "smiles": None}
        for input_identifier, output_identifiers in compound_identifiers_list
    ]
    print(resolved)


# (
#   CompoundIdentifier(
#       identifier_type=<CompoundIdentifierType.NAME: 6>,
#        value='[Ru(p-cymene)I2]2', details=None
#   ),
#   [[], []]
#  )
