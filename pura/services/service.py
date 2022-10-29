from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from abc import ABC, abstractmethod
from typing import List, Union


class Service(ABC):
    def __init__(self) -> None:
        pass

    @abstractmethod
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[Union[CompoundIdentifier, None]]:
        pass
