# pura
Clean reaction data quickly

## Installation

```bash
pip install pura
```

## Resolve molecules

```python
# Import pura
from pura.compound import Resolver

names = ["aspirin", "water", "methane", "1-methyl-propane"]
res = Resolver()
smiles = res.resolve(names, output=Resolver.SMILES)
```

