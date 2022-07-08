# pura
Pura helps you clean your messy chemical data. It fills the gap of making chemical data useable for machine learning algorithms.

## Installation

```bash
pip install pura
```

## What you can do with Pura

Pura can help with both compounds and reactions. Below are examples of its key features.

### Resolve common names to SMILES

Compounds are often recorded as common names instead of a machine readable identifier like SMILES.

There are several services that can do name resolution (PubChem, Chemical Identity Resolver, ChemSpider), and they sometimes disagree. Pura enables you to check several services asynchronously and ensure that they all agree on the resolved identifier. You can then discard or  manually check the names that could not be resolved.

```python
# Import pura
from pura.resolvers import resolve_names
from pura.compound import CompoundIdentifierType
from pura.services import Pubchem, CIR

# Resolve names to SMILES
# Agreement=2 ensures that at least two services agree on the SMILES
smiles = resolve_names(
    names=["aspirin", "ibuprofen", "[Ru(p-cymene)I2]2"],
    output_identifier_type=CompoundIdentifierType.SMILES,
    services=[Pubchem(), CIR()],
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


## Roadmap

- [x] Name resolution (July 2022)
- [ ] Reaction representations (July 2022)
- [ ] Reaction balancing (rxnmapper to start) (July - August 2022)
- [ ] Reports on quality (August 2022)
- [ ] Comparison quality of balancing and mapping on reaxys, USPTO and pistachio (September 2022)
- [ ] Submit paper to Neurips science workshop (September/October 2022_
- [ ] Documentation and website (November 2022)
- [ ] Template extraction (December 2022)
- [ ] Agreement/consensus algorithms for multiple representations of compounds



## Resources

- [Reaction Data Curation I: Chemical Structures and Transformations Standardization](https://doi.org/10.1002/minf.202100119)
- [RDchiral](https://github.com/connorcoley/rdchiral)
- [Selfies](https://github.com/aspuru-guzik-group/selfies)