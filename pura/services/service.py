from abc import ABC, abstractmethod


class Service(ABC):
    def __init__(self) -> None:
        pass

    # @abstractmethod
    # def resolve_compound(
    #     session: ClientSession,
    #     input_identifier: cp.CompoundIdentifier,
    #     output_identifier_type: cp.CompoundIdentifierType,
    # ):
    #     pass
