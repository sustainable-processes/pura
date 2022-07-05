from pura.compound import (
    CompoundIdentifier,
    CompoundIdentifierType,
    Compound,
    standardize_identifier,
)
from pura.services import Service, CIR, Pubchem, ChemSpider
from tqdm import tqdm
from aiohttp import ClientSession
from aiohttp.web_exceptions import (
    HTTPClientError,
    HTTPServerError,
    HTTPServiceUnavailable,
)
import asyncio
from typing import Optional, List, Union
from itertools import combinations
from functools import reduce
import logging

logger = logging.getLogger(__name__)


class CompoundResolver:
    def __init__(self, services: List[Service]):
        self._services = services

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
        """
        loop = asyncio.get_event_loop()
        logging.info("Running download")
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
    ) -> List[List[Union[CompoundIdentifier, None]]]:
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
        silent: bool = False,
    ) -> Union[List[CompoundIdentifier], None]:

        agreement_count = 0
        resolved_identifiers_list = []
        for i, service in enumerate(self._services):
            for j in range(n_retries):
                try:
                    resolved_identifiers = await service.resolve_compound(
                        session,
                        input_identifier=input_identifier,
                        output_identifier_type=output_identifier_type,
                    )
                    # Standardize identifiers (e.g., SMILES canonicalization)
                    for identifier in resolved_identifiers:
                        if identifier is not None:
                            standardize_identifier(identifier)
                    resolved_identifiers_list.append(resolved_identifiers)
                    break
                except HTTPServiceUnavailable:
                    # If server is busy, use exponential backoff
                    asyncio.sleep(2**j)
                except (HTTPClientError, HTTPServerError) as e:
                    # Log/raise on all other HTTP errors
                    if silent:
                        logger.error(e)
                        return
                    else:
                        raise TypeError(e)

            # Chceck agreement between services
            if i > 0 and len(resolved_identifiers_list) > 0:
                resolved_identifiers = self.reduce_options(
                    resolved_identifiers_list, agreement
                )
                agreement_count += 1
            elif len(resolved_identifiers_list) > 0:
                agreement_count += 1
            if agreement_count >= agreement:
                break

        if agreement_count < agreement:
            error_txt = f"Not sufficient agreeement ({agreement_count}) for {resolved_identifiers_list}"
            if silent:
                logger.error(error_txt)
                return
            else:
                raise TypeError(error_txt)

        return resolved_identifiers

    def reduce_options(
        self, identifiers_list: List[List[CompoundIdentifier]], agreement: int
    ) -> List[CompoundIdentifier]:
        """
        Reduce and deduplcate options

        Notes
        ------
        Algorithm:
        1. Find all combinatons of services that can satisfy agreement (combinations)
        2. Find the intersection of each combination
        3. If the intersection is greater than zero, then you have sufficient agreement.
        """
        options = list(range(len(identifiers_list)))
        intersection = []
        identifier_type = identifiers_list[0][0].identifier_type
        identifiers_list = [
            [identifier.value for identifier in identifiers]
            for identifiers in identifiers_list
        ]
        identifiers_sets = [set(ident) for ident in identifiers_list]
        for combo in combinations(options, agreement):
            intersection = reduce(
                set.intersection, [identifiers_sets[combo_i] for combo_i in combo]
            )
            if len(intersection) > 0:
                break
        if len(intersection) == 0:
            return []
        else:
            return [
                CompoundIdentifier(identifier_type=identifier_type, value=identifier)
                for identifier in intersection
            ]


""""
[A, B, C]
[B, C]
[]
|
|
V
[B]

"""


def resolve_names(
    names: List[str],
    output_identifier_type: CompoundIdentifierType,
    agreement: int = 1,
    batch_size: int = 100,
    services: Optional[List[Service]] = None,
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
    services : list of `Service`
        Services used to do resolution

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
        services = [Pubchem(), CIR()]
    resolver = CompoundResolver(services=services)
    name_identifiers = [
        CompoundIdentifier(identifier_type=CompoundIdentifierType.NAME, value=name)
        for name in names
    ]
    return resolver.resolve(
        input_identifiers=name_identifiers,
        output_identifier_type=output_identifier_type,
        agreement=agreement,
        batch_size=batch_size,
    )


if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)

    smiles = resolve_names(
        [
            "Pentylbenzene",
            "2-Phenylpentane",
            "3-Phenylpentane",
            "Oxalic acid",
            "3-Methyl-2-phenylbutane",
        ],
        output_identifier_type=CompoundIdentifierType.SMILES,
        services=[Pubchem(), CIR(), ChemSpider()],
        agreement=2,
        batch_size=5,
    )
    print(smiles)
