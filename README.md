[![Open app](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://moleculeresolver.streamlit.app/)

# pura
Pura helps you clean chemical and reaction data. Right now, you can use it to resolve common chemical names (e.g., aspirin) to standard cheminformatics identifiers like SMILES.

You can now access pura using our [web app](https://moleculeresolver.streamlit.app/)!


## Installation

```bash
pip install pura
```

## Resolve compound identifiers

Compounds are often recorded as common names instead of a machine readable identifier like SMILES.

There are several services that can do name resolution (PubChem, Chemical Identity Resolver, ChemSpider), and they sometimes disagree. Pura enables you to check several services asynchronously and ensure that a certain number agree on the resolved identifier. You can then discard or manually check the names that could not be resolved.

You can find a full list of services [here](https://github.com/sustainable-processes/pura/tree/main/pura/services).

```python
# Import pura
from pura.resolvers import resolve_identifiers
from pura.compound import CompoundIdentifierType
from pura.services import PubChem, CIR, CAS

# Resolve names to SMILES
resolved = resolve_identifiers(
    ["Josiphos SL-J001-1", "Rh(NBD)2BF4", "DuPhos"],
    input_identifer_type=CompoundIdentifierType.NAME,
    output_identifier_type=CompoundIdentifierType.SMILES,
    backup_identifier_types=[
        CompoundIdentifierType.INCHI_KEY,
        CompoundIdentifierType.CAS_NUMBER,
    ],
    services=[PubChem(autocomplete=True), CIR(), CAS()],
    agreement=1,
    silent=True,
)
print("\nResults\n")
for input_compound, resolved_identifiers in resolved:
    print(input_compound, resolved_identifiers, "\n")
#Josiphos SL-J001-1 [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, #value='C1CCCC1.CC(C1CCCC1P(c1ccccc1)c1ccccc1)P(C1CCCCC1)C1CCCCC1.[Fe]', details=None)]

# Rh(NBD)2BF4 [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='C1=CC2C=CC1C2.C1=CC2C=CC1C2.F[B-](F)(F)F.[Rh]', details=None)]

# Dichloro(p-cymene)ruthenium(II) dimer [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='Cc1ccc(C(C)C)cc1.Cc1ccc(C(C)C)cc1.Cl[Ru]Cl.Cl[Ru]Cl', details=None)]

# DuPhos [CompoundIdentifier(identifier_type=<CompoundIdentifierType.SMILES: 2>, value='CC(C)C1CCC(C(C)C)P1c1ccccc1P1C(C(C)C)CCC1C(C)C', details=None)]
```

## Development

<!-- ### Roadmap -->

<!-- - [x] Name resolution (July 2022)
- [x] Reaction representations (July 2022)
- [ ] Reaction balancing (July - August 2022)
- [ ] Reaction mapping (reaction mapper initially) (July - August 2022)
- [ ] Reports on quality (August 2022)
- [ ] Comparison quality of balancing and mapping on reaxys, USPTO and pistachio (September 2022)
- [ ] Write and submit paper to Neurips science workshop (September - October 2022)
- [ ] Publish package on pypi (September 2022)
- [ ] Documentation and website (November 2022)
- [ ] Template extraction (December 2022)
- [ ] Agreement/consensus algorithms for multiple representations of compounds -->

<!-- ### Getting set up -->

1. Install poetry using the following or via the instructions [here](https://python-poetry.org/docs/#installation):

    ```bash
    curl -sSL https://raw.githubusercontent.com/python-poetry/poetry/master/get-poetry.py | python -
    ```

2. Clone the repository:

    ```bash
    git clone https://github.com/sustainable-processes/pura.git
    cd pura
    ```

3. Install the dependencies from `poetry.lock`:

    ```bash
    poetry install
    ```

    Once you make some changes, commit and push:

    ```bash
    git commit -am <YOUR COMMIT MESSAGE>
    git push
    ```

# Subtleties of name resolution

Molecules will often be referred to with an English name, however, the same molecule can have many different names, and different molecules can have very similar (and sometimes even the same?) name! As an example, consider these two very similar names that refer to two different molecules:
- Phenyl acetate is the ester of phenol and acetic acid (CC(=O)Oc1ccccc1)
- Phenylacetate is an organic compound containing a phenyl functional group and a carboxylic acid functional group (O=C(O)Cc1ccccc1)

Furthermore, the name resolution can sometimes be further complicated by formal charges. Phenylacetate (a.k.a phenylacetic acid) is a carboxylic acid, so in water it will both be found as O=C(O)Cc1ccccc1 and O=C([O-])Cc1ccccc1, and indeed when querying services, both the charged and uncharged molecule was returned, which led to lack of agreement between services, despite the services having the same idea about what the molecule was.

Finally the presence/absense of stereochemical information can again cause disagreement between different services(Discussed in [Issue #45](https://github.com/sustainable-processes/pura/issues/45)). An example would be:
- Given the molecule: (e)-2-butenenitrile
- PubChem will resolve to: ['C/C=C/C#N']
- CIR will resolve to: ['CC=CC#N']

Using agreement=2 will require (at least) 2 data providers to be in agreement with each other, which would flag cases with ambiguity (since Pura would return None, so you avoid getting the wrong result).

The way these disagreements should be resolved will depend on the context, so it's probably not possible to apply a standardised way of resolving conflict - rather, researchers should be aware of these subtleties, and make informed decisions that fit with the goals of their own projects.

<!-- 
## Resources

- [Reaction Data Curation I: Chemical Structures and Transformations Standardization](https://doi.org/10.1002/minf.202100119)
- [RDchiral](https://github.com/connorcoley/rdchiral)
- [Selfies](https://github.com/aspuru-guzik-group/selfies)
- [CGRTools](https://doi.org/10.1021/acs.jcim.9b00102)
- [ChemDataExtractor](https://github.com/mcs07/ChemDataExtractor) -->
