# Dictionary for manual canonicalization

# initialize a dict that maps catalysts to the humanly cleaned smiles
catalyst_replacements = {}

# Add a catalyst to the catalyst_replacements df
# catalyst_replacements[
#     "C1CCC(P(C2CCCCC2)C2CCCCC2)CC1.O=P([O-])([O-])[O-].[K+].[K+].[K+]"
# ] = "C1CCC(P(C2CCCCC2)C2CCCCC2)CC1.O=P([O-])([O-])O.[K+].[OH-].[Pd+2]"
# catalyst_replacements[
#     "C1CCC(P(C2CCCCC2)C2CCCCC2)CC1.CC(=O)[O-].CC(=O)[O-].O=P([O-])([O-])[O-].[K+].[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "C1CCC(P(C2CCCCC2)C2CCCCC2)CC1.[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "CC(=O)O.CC(=O)O.COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O=P([O-])([O-])[O-].[K+].[Pd+2]"
# ] = "CC(=O)O.CC(=O)O.COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O.O=P([O-])([O-])[O-].[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)O.CC(=O)O.O=P([O-])([O-])[O-].[K+].[Pd+2]"
# ] = "CC(=O)O.CC(=O)O.O.O=P([O-])([O-])[O-].[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O.O=P([O-])(O)O.[K+].[Pd+2]"
# ] = "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O=P([O-])(O)O.[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.[Pd+2]"
# ] = "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O=P([O-])(O)O.[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O=P([O-])([O-])[O-].[K+].[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "CC(=O)[O-].CC(=O)[O-].COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].O.O=P([O-])([O-])[O-].[K+].[K+].[K+].[Pd+2].c1ccc(P(C2CCCCC2)C2CCCCC2)c(-n2c3ccccc3c3ccccc32)c1"
# ] = " CC(=O)[O-].CC(=O)[O-].O.O=P([O-])(O)O.[K+].[Pd+2].c1ccc(P(C2CCCCC2)C2CCCCC2)c(-n2c3ccccc3c3ccccc32)c1"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].O.O.O.O.O.O.O.O=P([O-])(O)O.[K+].[Pd+2]"
# ] = "CC(=O)[O-].CC(=O)[O-].O=P([O-])(O)O.[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].O.O=P([O-])(O)O.[K+].[Pd+2]"
# ] = "CC(=O)[O-].CC(=O)[O-].O=P([O-])(O)O.[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].O.O=P([O-])([O-])[O-].[K+].[K+].[K+].[Pd+2]"
# ] = "CC(=O)[O-].CC(=O)[O-].O=P([O-])(O)O.[K+].[Pd+2]"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "CC(=O)[O-].CC(=O)[O-].[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "CC(C)(C)P(C(C)(C)C)C(C)(C)C"
# ] = "CC(C)(C)P(C(C)(C)C)C(C)(C)C.Cl[Pd]Cl"
# catalyst_replacements[
#     "CC(C)(C)P(C(C)(C)C)C(C)(C)C.CC(C)(C)P(C(C)(C)C)C(C)(C)C.[Pd]"
# ] = "CC(C)(C)P(C(C)(C)C)C(C)(C)C.Cl[Pd]Cl"
# catalyst_replacements[
#     "CC(C)(C)P(CCCS(=O)(=O)O)C(C)(C)C"
# ] = "CC(C)(C)P(CCCS(=O)(=O)[O-])C(C)(C)C"
# catalyst_replacements[
#     "CC1(C)c2cccc(P(c3ccccc3)c3ccccc3)c2Oc2c(P(c3ccccc3)c3ccccc3)cccc21"
# ] = "CC1(C)c2cccc(P(c3ccccc3)c3ccccc3)c2Oc2c(P(c3ccccc3)c3ccccc3)cccc21.Cl[Pd]Cl"
# catalyst_replacements[
#     "CN(C)c1ccc(P(C(C)(C)C)C(C)(C)C)cc1.CN(C)c1ccc(P(C(C)(C)C)C(C)(C)C)cc1.[Cl-].[Cl-].[Pd+2]"
# ] = "CN(C)c1ccc(P(C(C)(C)C)C(C)(C)C)cc1.Cl[Pd]Cl"
# catalyst_replacements[
#     "CN(C)c1ccc([P+](C(C)(C)C)(C(C)(C)C)[Pd](Cl)(Cl)[P+](c2ccc(N(C)C)cc2)(C(C)(C)C)C(C)(C)C)cc1.O.O.O.O.O.O.O.O.O.O.O.O.O=P([O-])([O-])[O-].[Na+]"
# ] = "CN(C)c1ccc([P+](C(C)(C)C)(C(C)(C)C)[Pd](Cl)(Cl)[P+](c2ccc(N(C)C)cc2)(C(C)(C)C)C(C)(C)C)cc1"
# catalyst_replacements[
#     "CN(C)c1ccc([PH](C(C)(C)C)(C(C)(C)C)[Pd-2](Cl)(Cl)[PH](c2ccc(N(C)C)cc2)(C(C)(C)C)C(C)(C)C)cc1"
# ] = "CN(C)c1ccc([P+](C(C)(C)C)(C(C)(C)C)[Pd](Cl)(Cl)[P+](c2ccc(N(C)C)cc2)(C(C)(C)C)C(C)(C)C)cc1"
# catalyst_replacements[
#     "Cc1ccccc1P(c1ccccc1C)c1ccccc1C"
# ] = "Cc1ccccc1P(c1ccccc1C)c1ccccc1C.Cc1ccccc1P(c1ccccc1C)c1ccccc1C.Cl[Pd]Cl"
# catalyst_replacements[
#     "Cc1ccccc1P(c1ccccc1C)c1ccccc1C.Cl[Pd]Cl"
# ] = "Cc1ccccc1P(c1ccccc1C)c1ccccc1C.Cc1ccccc1P(c1ccccc1C)c1ccccc1C.Cl[Pd]Cl"
# catalyst_replacements[
#     "O.O.O.O.O.O.O.O=P([O-])([O-])[O-].[K+].[K+].[K+]"
# ] = "O.O=P([O-])([O-])[O-].[K+].[K+].[K+]"
# catalyst_replacements[
#     "O.O.O.O=P([O-])([O-])[O-].[K+].[K+].[K+]"
# ] = "O.O=P([O-])([O-])[O-].[K+].[K+].[K+]"
# catalyst_replacements[
#     "O=C(C=Cc1ccccc1)C=Cc1ccccc1"
# ] = "O=C(C=Cc1ccccc1)C=Cc1ccccc1.O=C(C=Cc1ccccc1)C=Cc1ccccc1.[Pd]"
# catalyst_replacements[
#     "O=C(C=Cc1ccccc1)C=Cc1ccccc1.O=C(C=Cc1ccccc1)C=Cc1ccccc1.O=C(C=Cc1ccccc1)C=Cc1ccccc1.[Pd].[Pd]"
# ] = "O=C(C=Cc1ccccc1)C=Cc1ccccc1.O=C(C=Cc1ccccc1)C=Cc1ccccc1.[Pd]"
# catalyst_replacements[
#     "O=P([O-])([O-])[O-].[Ag+].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "O=P([O-])([O-])[O-].[Ag+].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "O=P([O-])([O-])[O-].[K+].[Pd+2].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "O=P([O-])([O-])[O-].[K+].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "O=P([O-])([O-])[O-].[K+].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "O=P([O-])([O-])[O-].[K+].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[Cl-].[Cl-].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "[Cl-].[Cl-].[Pd+2].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[Pd].[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# ] = "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Rh+3].[Rh+3]"
# ] = "CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Rh+2].[Rh+2]"
# catalyst_replacements[
#     "[CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Rh+3]]"
# ] = "CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].CC(=O)[O-].[Rh+2].[Rh+2]"
# catalyst_replacements[
#     "[CC(C)(C)[P]([Pd][P](C(C)(C)C)(C(C)(C)C)C(C)(C)C)(C(C)(C)C)C(C)(C)C]"
# ] = "CC(C)(C)[PH]([Pd][PH](C(C)(C)C)(C(C)(C)C)C(C)(C)C)(C(C)(C)C)C(C)(C)C"
# catalyst_replacements[
#     "CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.CCCC[N+](CCCC)(CCCC)CCCC.[Br-].[Br-].[Br-]"
# ] = "CCCC[N+](CCCC)(CCCC)CCCC.[Br-]"
# catalyst_replacements["[CCO.CCO.CCO.CCO.[Ti]]"] = "CCO[Ti](OCC)(OCC)OCC"
# catalyst_replacements["[CC[O-].CC[O-].CC[O-].CC[O-].[Ti+4]]"] = "CCO[Ti](OCC)(OCC)OCC"
# catalyst_replacements[
#     "[Cl[Ni]Cl.c1ccc(P(CCCP(c2ccccc2)c2ccccc2)c2ccccc2)cc1]"
# ] = "Cl[Ni]1(Cl)[P](c2ccccc2)(c2ccccc2)CCC[P]1(c1ccccc1)c1ccccc1"
# catalyst_replacements[
#     "[Cl[Pd](Cl)([P](c1ccccc1)(c1ccccc1)c1ccccc1)[P](c1ccccc1)(c1ccccc1)c1ccccc1]"
# ] = "Cl[Pd](Cl)([PH](c1ccccc1)(c1ccccc1)c1ccccc1)[PH](c1ccccc1)(c1ccccc1)c1ccccc1"
# catalyst_replacements["[Cl[Pd+2](Cl)(Cl)Cl.[Na+].[Na+]]"] = "Cl[Pd]Cl"

# catalyst_replacements["[O=C([O-])[O-].[Ag+2]]"] = "O=C([O-])[O-].[Ag+].[Ag+]"
# catalyst_replacements["[O=S(=O)([O-])[O-].[Ag+2]]"] = "O=S(=O)([O-])[O-].[Ag+].[Ag+]"
# catalyst_replacements["[O=[Ag-]]"] = "O=[Ag]"
# catalyst_replacements["[O=[Cu-]]"] = "O=[Cu]"

# catalyst_replacements[
#     "[[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1]"
# ] = "[Pd].c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1.c1ccc(P(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[c1ccc([PH](c2ccccc2)(c2ccccc2)[Pd-4]([PH](c2ccccc2)(c2ccccc2)c2ccccc2)([PH](c2ccccc2)(c2ccccc2)c2ccccc2)[PH](c2ccccc2)(c2ccccc2)c2ccccc2)cc1]"
# ] = "c1ccc([PH](c2ccccc2)(c2ccccc2)[Pd]([PH](c2ccccc2)(c2ccccc2)c2ccccc2)([PH](c2ccccc2)(c2ccccc2)c2ccccc2)[PH](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[c1ccc([P]([Pd][P](c2ccccc2)(c2ccccc2)c2ccccc2)(c2ccccc2)c2ccccc2)cc1]"
# ] = "c1ccc([PH](c2ccccc2)(c2ccccc2)[Pd]([PH](c2ccccc2)(c2ccccc2)c2ccccc2)([PH](c2ccccc2)(c2ccccc2)c2ccccc2)[PH](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"
# catalyst_replacements[
#     "[c1ccc([P](c2ccccc2)(c2ccccc2)[Pd]([P](c2ccccc2)(c2ccccc2)c2ccccc2)([P](c2ccccc2)(c2ccccc2)c2ccccc2)[P](c2ccccc2)(c2ccccc2)c2ccccc2)cc1]"
# ] = "c1ccc([PH](c2ccccc2)(c2ccccc2)[Pd]([PH](c2ccccc2)(c2ccccc2)c2ccccc2)([PH](c2ccccc2)(c2ccccc2)c2ccccc2)[PH](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"
catalyst_replacements["Karstedt catalyst"] = "C[Si](C)(C=C)O[Si](C)(C)C=C.[Pt]"
# catalyst_replacements["Karstedt's catalyst"] = "C[Si](C)(C=C)O[Si](C)(C)C=C.[Pt]"
catalyst_replacements["Pd on-carbon"] = "[C].[Pd]"
catalyst_replacements["TEA"] = "OCCN(CCO)CCO"
# catalyst_replacements["Ti-superoxide"] = "O=[O-].[Ti]"
catalyst_replacements["sulfated tin oxide"] = "O=S(O[Sn])(O[Sn])O[Sn]"
catalyst_replacements[
    "tereakis(triphenylphosphine)palladium(0)"
] = "c1ccc([PH](c2ccccc2)(c2ccccc2)[Pd]([PH](c2ccccc2)(c2ccccc2)c2ccccc2)([PH](c2ccccc2)(c2ccccc2)c2ccccc2)[PH](c2ccccc2)(c2ccccc2)c2ccccc2)cc1"
catalyst_replacements["zeolite"] = "O=[Al]O[Al]=O.O=[Si]=O"

import pandas as pd
from pura.services import load_into_database, create_tables, LocalDatabase
from pura.compound import CompoundIdentifierType
from pura.resolvers import resolve_identifiers
import asyncio

# db_path = "pura.db"


async def _main():
    data = [{"name": key, "smiles": val} for key, val in catalyst_replacements.items()]
    # data.append(
    #     {
    #         "name": "zeolite",
    #         "smiles_ident": "O=[Al]O[Al]=O.O=[Si]=O",
    #         "smiles": "O=[Al]O[Al]=O.O=[Si]",
    #     }
    # )

    df = pd.DataFrame(data)

    await create_tables(error_if_exists=False)
    await load_into_database(
        df,
        identifier_columns=[
            ("name", CompoundIdentifierType.NAME, False),
            ("smiles", CompoundIdentifierType.SMILES, True),
        ],
        smiles_column="smiles",
    )


def main():
    # loop = asyncio.new_event_loop()
    # asyncio.set_event_loop(loop)

    # loop.run_until_complete(_main())

    results = resolve_identifiers(
        list(catalyst_replacements.keys())[:2],
        output_identifier_type=CompoundIdentifierType.SMILES,
        input_identifer_type=CompoundIdentifierType.NAME,
        services=[LocalDatabase(return_canonical_only=True)],
    )
    for res in results:
        print(res)


if __name__ == "__main__":
    # import logging

    # logging.basicConfig(level=logging.DEBUG)
    main()
    # from rdkit import Chem

    # for name, smi in catalyst_replacements.items():
    #     print(name, Chem.CanonSmiles(smi))
