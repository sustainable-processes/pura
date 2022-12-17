import streamlit as st
from typing import List, Optional, Tuple


st.title("Molecule name resolver")

with st.spinner("Loading..."):

    from pura.resolvers import CompoundResolver
    from pura.compound import CompoundIdentifierType
    from pura.compound import (
        CompoundIdentifier,
        CompoundIdentifierType,
        Compound,
        standardize_identifier,
        unique_identifiers,
    )
    from pura.services import *
    from rdkit import Chem
    from rdkit.Chem import Draw
    import numpy as np
    import pandas as pd


st.markdown(
    "This app resolves names of molecules to their corresponding SMILES strings using [pura](https://github.com/sustainable-processes/pura)."
)


@st.cache(suppress_st_warning=True)
def get_predictions(
    names: List[str],
    agreement: Optional[int] = 1,
) -> List[Tuple[str, str]]:
    output_identifier_type = CompoundIdentifierType.SMILES
    input_identifer_type: CompoundIdentifierType = CompoundIdentifierType.NAME
    backup_identifier_types = [
        CompoundIdentifierType.INCHI_KEY,
        CompoundIdentifierType.CAS_NUMBER,
    ]
    batch_size: int = 100
    services = [PubChem(autocomplete=True), CIR(), CAS()]
    silent = True

    compounds = [
        Compound(
            identifiers=[
                CompoundIdentifier(identifier_type=input_identifer_type, value=name)
            ]
        )
        for name in names
    ]
    resolver = CompoundResolver(services=services, silent=silent)
    results = resolver.resolve(
        input_compounds=compounds,
        output_identifier_type=output_identifier_type,
        backup_identifier_types=backup_identifier_types,
        agreement=agreement,
        batch_size=batch_size,
        progress_bar_type="streamlit",
    )
    return [
        (
            input_compound.identifiers[0].value,
            [identifier.value for identifier in resolved_identifiers],
        )
        for input_compound, resolved_identifiers in results
        if len(resolved_identifiers) > 0
    ]


upload_csv = st.checkbox("Upload a CSV")
# Upload CSV
names = None
if upload_csv:
    name_column = st.text_input("Column with molecule names", value="Name")
    csv = st.file_uploader("Upload CSV", type="csv")
    if csv is not None:
        df = pd.read_csv(csv)
        if name_column not in df.columns:
            st.error(
                f"""Column "{name_column}" not found in CSV. Change column name above."""
            )
        else:
            names = ",".join(df[name_column].tolist())
else:
    sample_names = "aspirin, ibuprofen, acetaminophen"
    names = st.text_input("Enter names", value=sample_names)

# Get and display predictions
if names:
    names = names.split(",")
    with st.spinner("Resolving names..."):
        results = get_predictions(names)
    smiles = [smi[0] for _, smi in results]
    names = [name for name, _ in results]

    # CSV
    df = pd.DataFrame({"Name": names, "SMILES": smiles})
    st.download_button(
        "Download SMILES",
        data=df.to_csv(index=False).encode("utf-8"),
        file_name="smiles.csv",
        mime="text/csv",
    )

    # Dipslay
    if len(smiles) > 10:
        st.write("Showing first 10 results")
    img = Draw.MolsToGridImage(
        [Chem.MolFromSmiles(smi) for smi in smiles[:10]],
        legends=names[:10],
        molsPerRow=2,
        subImgSize=(600, 400),
    )
    st.image(img)
