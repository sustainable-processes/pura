# pura
Clean reaction data quickly

## Installation

```bash
pip install pura
```

## Resolve molecules

```python
# Import pura
from pura.resolvers import resolve_names
from pura.compound import CompoundIdentifierType
from pura.services import Pubchem, CIR

smiles = resolve_names(
    ["aspirin", "ibuprofen", "toluene"],
    output_identifier_type=CompoundIdentifierType.SMILES,
    services=[Pubchem(), CIR()],
    agreement=2,
)
```

