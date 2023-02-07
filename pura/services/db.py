from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from typing import List, Optional, Union
from databases import Database
import sqlalchemy


metadata = sqlalchemy.MetaData()

# Schema

# compound
# --------------
# id: INT, KEY
# inchi_key: STR
# canonical_smiles: STR
# other_identifiers: Foreign Key to identifiers
#
# identifiers
# --------------
# id: INT, KEY
# identifier_type: ENUM
# identifier_value: STR

compound = sqlalchemy.Table(
    "compound",
    metadata,
    sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True),
    sqlalchemy.Column("identifiers", sqlalchemy.Integer),
)

identifiers = sqlalchemy.Table(
    "identifiers",
    metadata,
    sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True),
    sqlalchemy.Column("identifier_type", sqlalchemy.Enum(CompoundIdentifierType)),
    sqlalchemy.Column("value", sqlalchemy.String),
    sqlalchemy.Column("compound_id", sqlalchemy.Integer),
)


class Database(Service):
    def __init__(self, db_path: Optional[str] = None) -> None:
        db_path = db_path or "/path/to/db.json"
        self.db = Database(db_path)

    async def setup(self):
        await self.db.connect()

    async def teardown(self):
        await self.db.disconnect()

    async def resolve_compound(
        self,
        session: ClientSession,
        input_identifier: CompoundIdentifier,
        output_identifier_types: List[CompoundIdentifierType],
    ) -> List[Union[CompoundIdentifier, None]]:

        query = identifiers.select()
        values = {
            "value": input_identifier.value,
            "identifier_value": input_identifier.identifier_type,
        }
        rows = await self.db.fetch_all(query, values)
        for row in rows:
            query = compound.select()
            values = {"id": row.compound_id}
            compound = await self.db.fetch_one(query, values)
            if compound:
                query = identifiers.where(identifiers.c.compound_id == compound.id)
                await self.db.fetch_all(
        return [
            CompoundIdentifier(identifier_type=row.identifer_type, value=row.value)
            for row in rows
            if row.identifier_type in output_identifier_types
        ]
