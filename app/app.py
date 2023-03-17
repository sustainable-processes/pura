import streamlit as st
from typing import List, Optional, Tuple

from PIL import Image


st.set_page_config(
    page_title="Molecule name resolver",
    page_icon=Image.open("app/favicon-32x32.png"),
    # layout="wide",
    initial_sidebar_state="collapsed",
)

st.title("Molecule resolver")

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

with st.sidebar:
    st.markdown("## Options")
    agreement_ui = st.number_input(
        "Number of services that must agree on SMILES.",
        value=1,
        min_value=0,
        max_value=3,
        step=1,
    )
    n_retries_ui = st.number_input(
        "Number of retries. Decrease if taking a long time",
        value=3,
        min_value=0,
        max_value=10,
        step=1,
    )
    batch_size_ui = st.number_input(
        "Batch size. Increase if taking a long time",
        value=50,
        min_value=1,
        max_value=100,
        step=10,
    )
    st.markdown("---")
    st.markdown(
        "Thanks to [PubChem](https://pubchem.ncbi.nlm.nih.gov/), [CIR](https://cactus.nci.nih.gov/chemical/structure), [Opsin](https://opsin.ch.cam.ac.uk/), and [CAS](https://commonchemistry.cas.org/) for providing the data."
    )


st.markdown(
    "This app resolves names of molecules to their corresponding SMILES strings using [pura](https://github.com/sustainable-processes/pura). Click the arrow in the top left corner for more options."
)


# @st.cache(suppress_st_warning=True, show_spinner=False)
def get_predictions(
    names: List[str],
    agreement: int = 1,
    n_retries: int = 3,
    batch_size: int = 50,
) -> List[Tuple[str, str]]:
    output_identifier_type = CompoundIdentifierType.SMILES
    input_identifer_type: CompoundIdentifierType = CompoundIdentifierType.NAME
    backup_identifier_types = [
        CompoundIdentifierType.INCHI_KEY,
        CompoundIdentifierType.CAS_NUMBER,
    ]
    services = [
        PubChem(autocomplete=True),
        CIR(),
        Opsin(),
    ]
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
        progress_bar_type="streamlit",
        n_retries=n_retries,
        batch_size=batch_size,
    )
    return [
        (
            input_compound.identifiers[0].value,
            [identifier.value for identifier in resolved_identifiers],
        )
        for input_compound, resolved_identifiers in results
        if len(resolved_identifiers) > 0
    ]


names = None
container = st.container()
columns = st.columns(5)
df = None
name_column = "Name"
with columns[0]:
    do_resolve = st.button("Get SMILES", type="primary")

with container:
    sample_names = "aspirin, ibuprofen, acetaminophen"
    names = st.text_input("Enter names", value=sample_names)
    names = names.split(",")

    st.markdown("**OR**...")

    # Upload CSV
    csv = st.file_uploader("Upload CSV", type="csv")
    if csv is not None:
        df = pd.read_csv(csv)
        name_column = st.selectbox("Select column with molecule names", df.columns)
        names = df[name_column].astype(str).tolist()
        if df.shape[0] > 5:
            st.write(f"Showing first 5 of {df.shape[0]} rows")
        st.table(df.head(5))

if len(names) > 500:
    st.warning(
        f"Too many names {len(names)}. Please enter in a smaller batch or use pura from python."
    )

# Get and display predictions
if len(names) > 0 and do_resolve:
    names = [name.lstrip(" ").rstrip(" ") for name in names]
    with st.spinner("Resolving names..."):
        results = get_predictions(
            names,
            agreement=agreement_ui,
            n_retries=n_retries_ui,
            batch_size=batch_size_ui,
        )
    smiles = [smi[0] for _, smi in results]
    names = [name for name, _ in results]
    labels = [f"{name}\n({smi[0]})" for name, smi in results]

    # CSV
    final_df = pd.DataFrame({name_column: names, "SMILES": smiles})
    if df is not None:
        final_df = df.merge(final_df, on=name_column, how="left")
    st.download_button(
        "Download CSV",
        data=final_df.to_csv(index=False).encode("utf-8"),
        file_name="smiles.csv",
        mime="text/csv",
    )

    # Dipslay
    if len(smiles) > 10:
        st.write("Showing first 10 results")
    img = Draw.MolsToGridImage(
        [Chem.MolFromSmiles(smi) for smi in smiles[:10]],
        legends=labels[:10],
        molsPerRow=2,
        subImgSize=(600, 400),
    )
    st.image(img)
