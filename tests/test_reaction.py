from pura.reaction import reaction_from_smiles
import pytest


def test_reaction_from_smiles():
    reaction_from_smiles(
        "Brc1cncnc1.OB(O)c1ccsc1.COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O.O=P([O-])([O-])[O-].[K+].[K+].[K+].C[Mg+].[Cl-]>>c1ncc(-c2ccsc2)cn1"
    )
