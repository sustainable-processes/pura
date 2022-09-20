# pura
Pura helps you clean chemical and reaction data. It fills the gap of making chemical data useable for machine learning algorithms. You can use pura to:

- Resolve common chemical names (e.g., aspirin) to standard cheminformatics identifiers like SMILES
- Balance and atom-map reactions (future)
- Extract reaction templates (future)

## Installation

```bash
pip install pura
```

## What you can do with Pura

Pura can help with both compounds and reactions. Below are examples of its key features.

### Resolve common names to SMILES

Compounds are often recorded as common names instead of a machine readable identifier like SMILES.

There are several services that can do name resolution (PubChem, Chemical Identity Resolver, ChemSpider), and they sometimes disagree. Pura enables you to check several services asynchronously and ensure that a certain number agree on the resolved identifier. You can then discard or manually check the names that could not be resolved.

You can find a full list of services [here](https://github.com/sustainable-processes/pura/tree/main/pura/services).

```python
# Import pura
from pura.resolvers import resolve_names
from pura.compound import CompoundIdentifierType
from pura.services import Pubchem, CIR, Opsin

# Resolve names to SMILES
# Agreement=2 ensures that at least two services agree on the SMILES
smiles = resolve_names(
    names=["aspirin", "ibuprofen", "[Ru(p-cymene)I2]2"],
    output_identifier_type=CompoundIdentifierType.SMILES,
    services=[Pubchem(), CIR(), Opsin()],
    agreement=2,
)
#  Output (resolves the first two names, but not the third)
# [
#   [
#       CompoundIdentifier(
#           identifier_type=<CompoundIdentifierType.SMILES: 2>, 
#           value='CC(C)Cc1ccc(C(C)C(=O)O)cc1', 
#           details=None
#       )
#   ],
#   [
#       CompoundIdentifier(
#           identifier_type=<CompoundIdentifierType.SMILES: 2>,
#           value='CC(=O)Oc1ccccc1C(=O)O', 
#           details=None
#       )
#   ],
#   [],
# ]
```
## Concepts behind Pura

Pura is a package for the transform part of extract-transform-load (ETL) workflows in cheminformatics.

```mermaid
flowchart LR
G(Extract) ---> T(Transform) ---> L(Load)
```

## Development

### Roadmap

- [x] Name resolution (July 2022)
- [x] Reaction representations (July 2022)
- [ ] Reaction balancing (July - August 2022)
- [ ] Reaction mapping (reaction mapper initially) (July - August 2022)
- [ ] Reports on quality (August 2022)
- [ ] Comparison quality of balancing and mapping on reaxys, USPTO and pistachio (September 2022)
- [ ] Write and submit paper to Neurips science workshop (September - October 2022)
- [ ] Publish package on pypi (September 2022)
- [ ] Documentation and website (November 2022)
- [ ] Template extraction (December 2022)
- [ ] Agreement/consensus algorithms for multiple representations of compounds

### Getting set up

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


## Resources

- [Reaction Data Curation I: Chemical Structures and Transformations Standardization](https://doi.org/10.1002/minf.202100119)
- [RDchiral](https://github.com/connorcoley/rdchiral)
- [Selfies](https://github.com/aspuru-guzik-group/selfies)
- [CGRTools](https://doi.org/10.1021/acs.jcim.9b00102)
- [ChemDataExtractor](https://github.com/mcs07/ChemDataExtractor)
