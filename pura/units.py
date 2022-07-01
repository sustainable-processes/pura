from pydantic import BaseModel
import pint
from typing import List, TypeVar, Generic


T = TypeVar('T')
# from https://github.com/hgrecco/pint/issues/1166#issuecomment-1116309404
class PintType(Generic[T]):

    Q = pint.Quantity

    def __init__(self, q_check: str):
        self.q_check = q_check

    def __get_validators__(self):
        yield self.validate

    def validate(self, v):
        q = self.Q(v)
        assert q.check(self.q_check), f"Dimensionality must be {self.q_check}"
        return q

Mass = PintType("[mass]")
Amount = PintType("[substance]")
Volume = PintType("[volume]")