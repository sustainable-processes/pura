from pura.compound import (
    CompoundIdentifier,
    CompoundIdentifierType,
    Compound,
    standardize_identifiers,
)
from pura.services import Service, CIR, Pubchem
from tqdm import tqdm
from aiohttp import *
import asyncio
from typing import Optional, List, Union
import logging

aiohttp_errors = (
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
    ) -> List[Union[CompoundIdentifier, None]]:
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
                        n_retries=1,
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
        logger = logging.getLogger(__name__)

        agreement_count = 0
        resolved_identifiers_list = []
        for i, service in enumerate(self._services):

            for j in range(n_retries):
                # try:
                resolved_identifiers = await service.resolve_compound(
                    session,
                    input_identifier=input_identifier,
                    output_identifier_type=output_identifier_type,
                )
                # Standardize identifiers (e.g., SMILES canonicalization)
                standardize_identifiers(resolved_identifiers)
                resolved_identifiers_list.append(resolved_identifiers)
                break
                # except aiohttp_errors:
                #     # Increasing back off by 2^n with each retry
                #     # This should deal with the internet going out temporarily, etc.
                #     asyncio.sleep(2**j)
                # except TypeError:
                #     error_txt = f"Could not construct request for {input_identifier}"
                #     if silent:
                #         logger.error(error_txt)
                #         return
                #     else:
                #         raise TypeError(error_txt)

            if i > 0:
                if self._check_agreement(resolved_identifiers_list):
                    agreement_count += 1
            else:
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

    def _check_agreement(
        self, identifiers_list: List[List[CompoundIdentifier]]
    ) -> bool:
        # TODO: make this more robust
        return all([ident[0] for ident in identifiers_list])


def resolve_names(
    names: List[str],
    output_identifier_type: CompoundIdentifierType,
    agreement: int = 1,
    batch_size: int = 100,
    services: Optional[List[Service]] = None,
) -> List[CompoundIdentifier]:
    """Resolve a list of names to an identifier type(s).

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
    """
    if services is None:
        servicess = [Pubchem(), CIR()]
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
    # import pint

    # ureg = pint.UnitRegistry()
    # aspirin = Compound(
    #     identifiers=[
    #         CompoundIdentifier(
    #             identifier_type=CompoundIdentifierType.SMILES,
    #             value="O=C(C)Oc1ccccc1C(=O)O",
    #         )
    #     ],
    # )

    smiles = resolve_names(
        ["aspirin", "ibuprofen", "toluene"],
        output_identifier_type=CompoundIdentifierType.SMILES,
        services=[Pubchem(), CIR()],
        agreement=2,
    )
    print(smiles)
