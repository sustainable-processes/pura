from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from abc import ABC, abstractmethod
from typing import List, Union
from functools import wraps
import time


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        # first item in the args, ie `args[0]` is `self`
        print(f"Function {func.__name__}{args} {kwargs} Took {total_time:.4f} seconds")
        return result

    return timeit_wrapper


class Service(ABC):
    def __init__(self) -> None:
        pass

    async def setup(self):
        pass

    async def teardown(self):
        pass

    @abstractmethod
    @timeit
    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[Union[CompoundIdentifier, None]]:
        pass
