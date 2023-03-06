import warnings
from .compound import Compound, Texture

from .units import *
from typing import Any, Callable, List, Optional, Dict, Union
import pandas as pd
from enum import Enum
from pydantic import BaseModel, validator
from rdkit.Chem import rdChemReactions
from rdkit import Chem


class ReactionRole(Enum):
    UNSPECIFIED = "UNSPECIFIED"
    REACTANT = "REACTANT"
    """
    A reactant is any compound that contributes atoms to a desired or
    observed product.
    """

    AGENT = "AGENT"
    """
    An agent is any compound/species that enables the reaction but does not contribute atoms to the product.
    Reagent, solvents, and catalysts are all examples of agents.
    """
    REAGENT = "REAGENT"
    SOLVENT = "SOLVENT"
    CATALYST = "CATALYST"

    WORKUP = "WORKUP"
    """
    The workup role is used when defining quenches, buffer additives for
    liquid-liquid separations, etc.
    """

    INTERNAL_STANDARD = "INTERNAL_STANDARD"
    """
    Internal standards can be included as part of a reaction input (when
    added prior to the start of the reaction) or as part of a workup
    step of addition.
    """
    AUTHENTIC_STANDARD = "AUTHENTIC_STANDARD"

    PRODUCT = "PRODUCT"
    """
    A product can be any species produced by the reaction, whether desired
    or undesired.
    """


class ReactionInput(PintModel):
    """Inputs to a reaction"""

    compound: Compound
    """Compound information for the input"""

    role: ReactionRole
    """The role of the compound in the reaction. Cannot be ReactionRole.PRODUCT"""

    addition_order: Optional[int] = None
    """
    The order in which the compound was added to the reaction starting at 1
    If order is specified for one input, it must be specified for all inputs.
    """

    addition_time: Optional[Time] = None
    """The time at which the compound was added to the reaction"""

    addition_duration: Optional[Time] = None
    """The duration of the addition"""

    flow_rate: Optional[Union[MassFlow, MolarFlow, VolumeFlow]] = None
    """The flow rate of the compound. Useful in continuous flow reactions"""

    addition_temperature: Optional[Temperature] = None
    """The temperature at which the compound was added
    
    Note that when passing pint temperatures, you need to use the quantity object:
    `addition_temperature=Quantity(25, "degC")`

    See: https://pint.readthedocs.io/en/stable/user/nonmult.html
    
    """

    @validator("role")
    def check_role(cls, v):
        if v == ReactionRole.PRODUCT:
            raise ValueError("Reaction input role cannot be {ReactionRole.PRODUCT}")
        return v


class ReactionProduct(PintModel):
    compound: Compound
    """Compound information for the product"""

    role: ReactionRole = ReactionRole.PRODUCT
    """Reaction time (for flow, equivalent to residence time or spacetime)."""

    product_yield: Optional[float] = None
    """The yield of the product as a percentage between 0 and 100"""

    selectivity: Optional[float] = None
    """The selectivity of the product as a percentage between 0 and 100"""

    peak_area: Optional[float] = None
    """The peak area of the product in a chromatogram"""

    counts: Optional[float] = None
    """The raw counts of the product in a mass spectrum"""

    is_desired_product: Optional[bool] = None
    """If this is the desired reaction product"""

    isolated_color: Optional[str] = None
    """A description of the color of the isolated product"""

    texture: Optional[Texture] = None
    """Specify the texture type and details such as 'fine crystals'"""

    @validator("product_yield")
    def yield_limits(cls, v):
        if v and v < 0:
            raise ValueError(f"Yield must be greater than 0% (currently {v}%)")
        if v and v > 100:
            warnings.warn(f"Yield is greater than 100% (currently {v}%)")
        return v

    @validator("selectivity")
    def selectivity_limits(cls, v):
        if v and v < 0:
            raise ValueError(f"Selectivity must be greater than 0% (currently {v}%)")
        if v and v > 100:
            warnings.warn(f"Selectivity is greatter than 100% (currently {v}%)")
        return v


class ReactionOutcome(PintModel):
    """Reaction outcomes are measurements of reaction results at a specific time"""

    products: List[ReactionProduct]
    """The product(s) of the reaction measured at the reaction time (if known)
    If the quantities of leftover starting materials are measured, these 
    starting materials may also be defined as product.
    """

    reaction_time: Optional[Time] = None
    """The time at which the reaction outcome was measured"""

    @validator("products")
    def one_desired_product(cls, v):
        if sum([1 if p.is_desired_product else 0 for p in v]) > 1:
            raise ValueError(
                "Only one desired product can be specified per reaction outcome"
            )
        return v

    @validator("reaction_time")
    def reaction_time_positive(cls, v):
        if v is not None and v < 0:
            raise ValueError("Reaction time must be positive")
        return v


class ReactionConditions(PintModel):
    temperature: Optional[Temperature] = None
    """Reaction temperature as a pint Quantity"""

    pressure: Optional[Pressure] = None
    """Reaction pressure as a pint Quantity"""

    reflux: bool = False
    """Whether the reaction is a reflux reaction. Default is False"""

    ph: Optional[float] = None
    """The pH of the reaction"""


class Reaction(PintModel):
    """A reaction is specified by inputs, conditions and outcomes."""

    inputs: List[ReactionInput]
    """The inputs to the reaction. These include reactants, reagents, solvents, and catalysts"""

    conditions: ReactionConditions
    """The reaction conditions. If no reaction conditions are known, pass an empty ReactionConditions object."""

    outcomes: List[ReactionOutcome]
    """
    The reaction outcomes are measurements of reaction results at a specific time.
    Each reaction can have multiple outcomes, each measured at a different time.
    """

    @validator("inputs")
    def addition_order_linear(cls, v):
        if not all([i.addition_order is None for i in v]):
            if not all([i.addition_order is not None for i in v]):
                raise ValueError(
                    "If any addition order is specified, all addition orders must be specified"
                )
            if not sorted(list(set([i.addition_order for i in v]))) == list(
                range(1, len(v) + 1)
            ):
                raise ValueError(
                    "Addition order must be a linear sequence starting at 1"
                )
        return v

    def reaction_smiles(self, split_agents: bool = False):
        """Construct a reaction SMILES"""
        reactants = ".".join([i.compound.to_smiles() for i in self.inputs])
        agents = ".".join([i.compound.to_smiles() for i in self.inputs])
        reactants = (
            f"{reactants}.{agents}" if not split_agents else f"{reactants}>{agents}"
        )
        products = ".".join(
            [p.compound.to_smiles() for p in self.outcomes[-1].products]
        )
        if split_agents:
            return f"{reactants}>{products}"
        else:
            return f"{reactants}>>{products}"

    @property
    def reactant_compounds(self):
        """Return a list of reactant compounds"""
        return [i.compound for i in self.inputs if i.role == ReactionRole.REACTANT]

    @property
    def product_compounds(self):
        """Return a list of product compounds from the last reaction outcome"""
        return [
            p.compound
            for p in self.outcomes[-1].products
            if p.role == ReactionRole.PRODUCT
        ]

    @property
    def agent_compounds(self):
        """Return a list of agents. This includes reagents, solvents, and catalysts"""
        return [
            i.compound
            for i in self.inputs
            if (i.role == ReactionRole.AGENT)
            or (i.role == ReactionRole.SOLVENT)
            or (i.role == ReactionRole.CATALYST)
            or (i.role == ReactionRole.REAGENT)
        ]

    @property
    def reagent_compounds(self):
        """Return a list of reagents"""
        return [i.compound for i in self.inputs if i.role == ReactionRole.REAGENT]

    @property
    def solvent_compounds(self):
        """Return a list of solvents"""
        return [i.compounds for i in self.inputs if i.role == ReactionRole.SOLVENT]
    
    @property
    def catalyst_compounds(self):
        """Return a list of catalysts"""
        return [i.compounds for i in self.inputs if i.role == ReactionRole.CATALYST]

    @property
    def reaction_time(self):
        """Return the reaction time for the last reaction outcome"""
        return self.outcomes[-1].reaction_time

    @property
    def reaction_yield(self):
        """Return the reaction yield for the desired product at the last reaction outcome"""
        for p in self.outcomes[-1].products:
            if p.is_desired_product:
                return p.product_yield
        warnings.warn("No desired product found.")


def to_frame(reactions: List[Reaction], columns: List[str]) -> pd.DataFrame:
    reaction_dicts = [reaction.to_dict() for reaction in reactions]
    return pd.DataFrame(reaction_dicts)


def reaction_from_smiles(
    reaction_smiles: str,
    conditions: Optional[ReactionConditions] = None,
    reaction_time: Optional[Time] = None,
    reaction_yield: Optional[float] = None,
    role_lookup: Optional[Dict[ReactionRole, List[Union[str, Compound]]]] = None,
    desired_product_check: Optional[Callable[[Compound], bool]] = None,
) -> Reaction:
    """Construct a reaction from a reaction SMILES

    Parameters
    ----------
    reaction_smiles : str
        Reaction SMILES
    conditions : ReactionConditions, optional
        Reaction conditions, by default None
    role_lookup : Optional[Dict[ReactionRole, List[Union[str, Compound]]]], optional
        A dictionary of compounds that should be treated a specific role.
        This only applies to reactants and agents.s
    reaction_time : Time, optional
        Reaction time, by default None
    reaction_yield : float, optional
        Reaction yield, by default None
    desired_product_check : Callable[[Compound], bool], optional
        A function that checks if a compound is the desired product, by default None

    Returns
    -------
    Reaction

    Notes
    -----
    The role_lookup dictionary is used to assign roles to reactants, agents, and products.
    The keys of the dictionary are ReactionRole and the values are lists of compounds or SMILES strings.

    The desired_product_check function is used to check if a compound is the desired product. It must
    be passed if reaction_yield is passed. If only one product is always passed the following can b
    used to always set it as the desired product:
    ```python
    reaction_from_smiles(..., desired_product_check=lambda c: True)
    ```

    Examples
    --------
    >>> role_lookup = {ReactionRole.SOLVENT: ['CCO', 'CCO']}
    >>> rxn_smiles = "Brc1cncnc1.OB(O)c1ccsc1.COc1cccc(OC)c1-c1ccccc1P(C1CCCCC1)C1CCCCC1.O.O=P([O-])([O-])[O-].[K+].[K+].[K+].C[Mg+].[Cl-]>>c1ncc(-c2ccsc2)cn1"
    >>> reaction_from_smiles(rxn_smiles, role_lookup=role_lookup)

    """
    # Read reaction SMILES using RDKit
    rxn = rdChemReactions.ReactionFromSmarts(reaction_smiles, useSmiles=True)

    # Get reactant compounds
    inputs = [
        ReactionInput(compound=Compound.from_rdkit_mol(mol), role=ReactionRole.REACTANT)
        for mol in rxn.GetReactants()
    ]

    # Get agent compounds
    inputs.extend(
        [
            ReactionInput(
                compound=Compound.from_rdkit_mol(mol), role=ReactionRole.AGENT
            )
            for mol in rxn.GetAgents()
        ]
    )
    if role_lookup:
        for role, compounds in role_lookup.items():
            for compound in compounds:
                if isinstance(compound, str):
                    compound = Compound.from_smiles(compound)
                for input in inputs:
                    if input.compound == compound:
                        input.role = role

    # Get reaction outcome
    if reaction_yield is not None and desired_product_check is None:
        raise ValueError(
            "If reaction yield is specified, desired product check must also be specified."
        )
    products = [
        ReactionProduct(compound=Compound.from_rdkit_mol(mol))
        for mol in rxn.GetProducts()
    ]
    if desired_product_check is not None:
        for product in products:
            if desired_product_check(product.compound):
                product.is_desired_product = True
                product.product_yield = reaction_yield
    outcome = ReactionOutcome(products=products, reaction_time=reaction_time)

    # Construct reaction
    conditions = conditions or ReactionConditions()
    return Reaction(inputs=inputs, conditions=conditions, outcomes=[outcome])
