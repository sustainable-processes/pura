from abc import ABC, abstractmethod
from typing import List
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession


class Service(ABC):
    def __init__(self) -> None:
        pass

    @abstractmethod
    def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_type: CompoundIdentifierType,
    ) -> List[CompoundIdentifierType]:
        pass
