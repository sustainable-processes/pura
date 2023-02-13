from pura.services import Service
from pura.compound import CompoundIdentifier, CompoundIdentifierType
from aiohttp import ClientSession
from typing import List, Optional, Union
from databases import Database as AsyncDatabase
import sqlalchemy
import pandas as pd
import asyncio
import nest_asyncio
from rdkit import Chem
from sqlalchemy.dialects.sqlite import insert

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

metadata = sqlalchemy.MetaData()

compound_table = sqlalchemy.Table(
    "compound",
    metadata,
    sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True),
    sqlalchemy.Column("inchi", sqlalchemy.String, nullable=False),
    sqlalchemy.UniqueConstraint("inchi", sqlite_on_conflict="IGNORE"),
)

identifiers_table = sqlalchemy.Table(
    "identifiers",
    metadata,
    sqlalchemy.Column("id", sqlalchemy.Integer, primary_key=True),
    sqlalchemy.Column("identifier_type", sqlalchemy.Enum(CompoundIdentifierType)),
    sqlalchemy.Column("value", sqlalchemy.String),
    sqlalchemy.Column(
        "compound_id",
        sqlalchemy.Integer,
        sqlalchemy.ForeignKey(compound_table.c.id),
        nullable=False,
    ),
    sqlalchemy.UniqueConstraint(
        "identifier_type", "value", "compound_id", sqlite_on_conflict="IGNORE"
    ),
)


class Database(Service):
    def __init__(self, db_path: Optional[str] = None) -> None:
        db_path = db_path or "/path/to/db.json"
        self.db = AsyncDatabase(db_path)

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
        # Find the compound id
        query = compound_table.select()
        query = query.where(
            (identifiers_table.c.value == input_identifier.value)
            & (identifiers_table.c.identifier_type == input_identifier.identifier_type)
        )
        id = await self.db.fetch_one(query=query)
        if id is None:
            return []
        id = id[0]

        # Find all matching ouptut identifiers
        query = sqlalchemy.select(
            identifiers_table.c.identifier_type, identifiers_table.c.value
        )
        query = query.where(
            (identifiers_table.c.compound_id == id)
            & (identifiers_table.c.identifier_type.in_(output_identifier_types))
        )
        rows = await self.db.fetch_all(query)

        return [
            CompoundIdentifier(
                identifier_type=row[0],
                value=row[1],
            )
            for row in rows
        ]


async def create_tables(db_path: str):
    """Create tables in the database"""
    db = AsyncDatabase(db_path)
    await db.connect()
    for table in metadata.tables.values():
        # Set `if_not_exists=False` if you want the query to throw an
        # exception when the table already exists
        schema = sqlalchemy.schema.CreateTable(table, if_not_exists=True)
        query = str(schema.compile())
        await db.execute(query=query)
    await db.disconnect()


async def load_into_database(
    data: pd.DataFrame,
    db_path: str,
    input_identifier_type: CompoundIdentifierType,
    input_identifier_column: str,
    output_identifier_type: CompoundIdentifierType,
    output_identifier_column: str,
    smiles_column: Optional[str] = None,
    inchi_column: Optional[str] = None,
    update_on_conflict: bool = False,
):
    """Load data into the database

    Parameters
    ----------
    data : pd.DataFrame
        Data to load into the database
    db_path : str
        Path to the database
    identifier_type : CompoundIdentifierType
        Type of identifier
    identifier_column : str
        Column containing the identifier
    inchi_column : Optional[str], optional
        Column containing InChI, by default None
    smiles_column : Optional[str], optional
        Column containing SMILES, used for converting to InChI
    update_on_conflict : bool, optional
        Whether to update the database if the compound and identifier already exists, by default False

    """
    # connect to the database
    db = AsyncDatabase(db_path)
    await db.connect()

    # Get inchi
    if inchi_column is not None:
        inchi = data[inchi_column]
    elif smiles_column is not None:
        inchi = [
            Chem.MolToInchi(Chem.MolFromSmiles(smiles))
            for smiles in data[smiles_column]
        ]
    else:
        raise ValueError("Must provide inchi or SMILES column.")

    # Upsert compounds
    query = insert(compound_table)
    if update_on_conflict:
        query = query.on_conflict_do_update(
            index_elements=[compound_table.c.inchi],
        )
    else:
        query = query.on_conflict_do_nothing(
            index_elements=[compound_table.c.inchi],
        )
    values = [{"inchi": inchi} for inchi in inchi]
    await db.execute_many(query=query, values=values)

    # Get compound ids
    query = "SELECT id FROM compound WHERE inchi = :inchi"
    ids = []
    for value in values:
        id = await db.fetch_one(query=query, values=value)
        ids.append(id[0])

    # Upsert identifiers
    query = insert(identifiers_table)
    index_elements = [
        identifiers_table.c.compound_id,
        identifiers_table.c.identifier_type,
        identifiers_table.c.value,
    ]
    if update_on_conflict:
        query = query.on_conflict_do_update(
            index_elements=index_elements,
        )
    else:
        query = query.on_conflict_do_nothing(
            index_elements=index_elements,
        )
    values = []
    for col, identifier_type in zip(
        [input_identifier_column, output_identifier_column],
        [input_identifier_type, output_identifier_type],
    ):
        values += [
            {
                "value": row[col],
                "identifier_type": identifier_type,
                "compound_id": ids[i],
            }
            for i, row in data.reset_index(drop=True).iterrows()
        ]
    await db.execute_many(query=query, values=values)

    # Disconnect from the database
    await db.disconnect()


async def test_query():
    db_path = "sqlite+aiosqlite:///pura.db"
    db = AsyncDatabase(db_path)
    await db.connect()

    query = compound_table.select()
    query = query.where(
        (identifiers_table.c.value == "fluticasone furoate")
        & (identifiers_table.c.identifier_type == CompoundIdentifierType.NAME)
    )
    id = await db.fetch_one(query=query)
    if id is not None:
        id = id[0]
        print(f"ID: {id}")

        query = sqlalchemy.select(
            identifiers_table.c.value, identifiers_table.c.identifier_type
        )
        query = query.where(
            (identifiers_table.c.compound_id == id)
            & (
                identifiers_table.c.identifier_type.in_(
                    [CompoundIdentifierType.SMILES, CompoundIdentifierType.NAME],
                )
            )
        )
        rows = await db.fetch_all(query)
        print([row for row in rows])

    await db.disconnect()


async def _main():
    db_path = "sqlite+aiosqlite:///pura.db"
    print("Creating tables...")
    await create_tables(db_path=db_path)
    data = pd.read_csv("smiles.csv")
    print("Loading data...")
    await load_into_database(
        data=data,
        db_path=db_path,
        input_identifier_column="SMILES",
        input_identifier_type=CompoundIdentifierType.SMILES,
        output_identifier_column="Name",
        output_identifier_type=CompoundIdentifierType.NAME,
        smiles_column="SMILES",
    )
    await test_query()


def main():
    loop = asyncio.new_event_loop()
    asyncio.set_event_loop(loop)

    loop.run_until_complete(_main())


if __name__ == "__main__":
    main()
