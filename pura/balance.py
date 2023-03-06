from collections import Counter
from typing import List, Optional,Tuple,Dict, Union

import pandas as pd
import modin.pandas as mpd

from .compound import Compound, CompoundIdentifier, CompoundIdentifierType, standardize_identifier
from .reaction import Reaction, ReactionIdentifier, ReactionIdentifierType, ReactionInput, ReactionRole
from rdkit import Chem
from rdkit.Chem import PeriodicTable, GetPeriodicTable
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import ray

from .helpCompound import hc_Dict # Help compounds


#%% Maps
reaction_identifier_map = {
    "REACTION_SMILES": [ReactionIdentifier, {"identifier_type": ReactionIdentifierType.REACTION_SMILES}],
    "REACTION_CXSMILES": [ReactionIdentifier, {"identifier_type": ReactionIdentifierType.REACTION_CXSMILES}],
    "RDFILE": [ReactionIdentifier, {"identifier_type": ReactionIdentifierType.RDFILE}],
    "RINCHI": [ReactionIdentifier, {"identifier_type": ReactionIdentifierType.RINCHI}],
    "NAME": [ReactionIdentifier, {"identifier_type": ReactionIdentifierType.NAME}],
    "REACTANT": [ReactionInput, {"role": ReactionRole.REACTANT}]}
reaction_role_map={
    "AGENT": [ReactionRole.AGENT],
    "PRODUCT": [ReactionRole.PRODUCT],
    "REAGENT": [ReactionRole.REAGENT],
    "SOLVENT": [ReactionRole.SOLVENT],
    "CATALYST": [ReactionRole.CATALYST]}
    
compound_identifier_map={"SMILES":[CompoundIdentifier,{"identifier_type":CompoundIdentifierType.SMILES}],
                         "INCHI":[CompoundIdentifier,{"identifier_type":CompoundIdentifierType.INCHI}],
                         "NAME":[CompoundIdentifier,{"identifier_type":CompoundIdentifierType.NAME}],
                         "CAS":[CompoundIdentifier,{"identifier_type":CompoundIdentifierType.CAS_NUMBER}]}

#%% Functions
                         
def parse_smiles_list(reaction_smiles_list:List[str])->List[Reaction]:
    """ Parse a list of reaction smiles into a list of Reaction objects.

    Args:
        reaction_smiles_list (List[str]): List of reaction smiles

    Returns:
        List[Reaction]: List of Reaction objects
    """
    reactions=[]
    for reaction_smiles in reaction_smiles_list:
        reactions.append(parse_input({"REACTION_SMILES":reaction_smiles}))
    return reactions

    
# reaction_smiles=[r1,r2],reactants=[...],agents=[[a1,a2,a3],[a3,a4,a5]],products=[...]

# General list
# a.b>c.d>e.f, solvents, catalysts, reagents
# a.b>>c.d, solvents, catalysts, agents, reagents (either compound or smiles)

# No general list
# List of dictionaries
# List of nested dictionaries
# Dataframe



# Subject to change

def parse_input(input_identifiers:Union(Dict[str,str],List[Dict[str,str]]),canonicalize=True)->Reaction:
    """
    Parse an input identifier dictionary into a Reaction object. Note the format is very specific. Need to generalize,
    make simpler
    
    Args:
        
        reaction_identifiers (Dict[str,str]): Reaction identifiers. Defaults to None.
                                            eg. {'REACTION_SMILES': 'CCO.CC(=O)O.[O-][N+](=O)c1ccc(Cl)cc1>>CCO.CC(=O)O.[N+](=O)[O-]c1ccc(Cl)cc1'}
                                            
                                            eg. {'REACTION_SMILES': 'CCO.CC(=O)O.[O-][N+](=O)c1ccc(Cl)cc1>>CCO.CC(=O)O.[N+](=O)[O-]c1ccc(Cl)cc1', 
                                                'REACTANT': [{"SMILES":'CCO'},{"SMILES":...}], 'AGENT': [{"SMILES":'CC(=O)O'}], 'PRODUCT': {"SMILES":'[O-][N+](=O)c1ccc(Cl)cc1'}}
                                                
                                            eg. {'REACTANT': [{"SMILES":'CCO'}], 'AGENT': [{"SMILES":'CC(=O)O'}], 'PRODUCT': [{"SMILES":'[O-][N+](=O)c1ccc(Cl)cc1'}]}
                                            
                                            eg. {'REACTION_SMILES':'[O-][N+](=O)c1ccc(Cl)cc1>CCO.CC(=O)O>CCO.CC(=O)O.[N+](=O)[O-]c1ccc(Cl)cc1'}
                                            
                                            eg. {'REACTANT': [{"SMILES":'CCO'}], 'REAGENT': [{"SMILES:'CC(=O)O'}], 'PRODUCT': [{"SMILES:'[O-][N+](=O)c1ccc(Cl)cc1'}]}
    Returns:
        Reaction: Reaction object
    """
    reaction_identifiers=[]
    reaction_inputs=[]
    for elem,val in input_identifiers.items():
        if elem in reaction_identifier_map:
            args=reaction_identifier_map[elem][1]
            args.update({"value":val})
            reaction_identifiers.append(reaction_identifier_map[elem][0](**args))
        elif elem in reaction_role_map:
            for compddict in val:
                compound_identifiers=[]
                for compound_identifier_type,compound_val in compddict.items():
                    if compound_identifier_type in compound_identifier_map:
                        args=compound_identifier_map[compound_identifier_type][1]
                        args.update({"value":compound_val})
                        compound_identifier=compound_identifier_map[compound_identifier_type][0](**args)
                        if canonicalize:
                            compound_identifier=standardize_identifier(compound_identifier)
                        compound_identifiers.append(compound_identifier)
                compound=Compound(identifiers=compound_identifiers)
                reaction_input=ReactionInput(components={reaction_role_map[elem][0]:compound})
                reaction_inputs.append(reaction_input)
    reaction=Reaction(identifiers=reaction_identifiers,inputs=reaction_inputs)
    return reaction
        
            
            
#%% Balancing functions (not yet integrated)
class CustomError(Exception):
    '''
    For error handling *Change/adapt as needed*
    '''
    pass

def initray(restart=True,num_cpus=16,log_to_driver=False):
    """
    Initializes a Ray CPU cluster for parallel execution of functions

    Args:
        restart (bool, optional): Restart existing cluster. Defaults to True.
        num_cpus (int, optional): Number of CPUs for parallel execution. Defaults to 16.
        log_to_driver (bool, optional): Controls whether results are logged to driver. Defaults to False.
    """
    if restart:
        ray.shutdown()
    ray.init(num_cpus=num_cpus,log_to_driver=log_to_driver)
    
def molfromsmiles(SMILES):
    """
    Converts a smiles string into a mol object for RDKit use. Returns a mol object. *Unsure how to integrate in compound.py*

    Args:
        SMILES (str): SMILES string to convert to mol object

    Returns:
        mol: RDKit mol object
    """
    mol=Chem.MolFromSmiles(SMILES)
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)
    return mol
    
def atomtypes(mol)->Tuple[Dict,int]:
    """
    Generates an atom type dictionary with counts for a given mol file. Also
    returns overall change of mol. 

    Args:
        mol (RDKit mol): RDKit mol object to generate atom type dictionary for

    Returns:
        Tuple[Dict,int]: Dictionary with keys as atom element and value as count, overall charge
    """

    typedict={}
    mol2=Chem.AddHs(mol) #Assumes hydrogens are added
    charge=0
    for atom in mol2.GetAtoms():
        elem=PeriodicTable.GetElementSymbol(GetPeriodicTable(),atom.GetAtomicNum())
        charge+=atom.GetFormalCharge()
        if elem not in typedict.keys():
            typedict[elem]=1
        else:
            typedict[elem]+=1
    return typedict, charge



def getcompdict(ID=1,mol=None,smiles=None,formula=None)->Dict:
    """
    Wrapper for atomtypes. ID must be input whether Reaxys ID or random number. Default is 1. 
    Smiles and mol can also be specified. Can be converted to an object instance if intended.

    Args:
        ID (int, optional): ID of compound/mol. Defaults to 1.
        mol (RDKit mol, optional): RDKit mol. Defaults to None.
        smiles (str, optional): SMILES string. Defaults to None.
        formula (str, optional): Molecular formula. Defaults to None.

    Raises:
        CustomError: Raises error if smiles specified
        is not valid/cannot be processed

    Returns:
        Dict: Compound dictionary with ID as key and another dictionary with the
        following keys: 'atomdict' (Output of atomtypes, dictionary with atom
        elements and count), 'charge' (Overall charge of mol),'smiles'
        (smiles of mol),'formula' (chemical formula),'count' (instance number, always 1)
    """

    if smiles is not None:
        try:
            mol=molfromsmiles(smiles)
        except Exception as e:
            raise CustomError("Compound "+str(ID)+" smiles is invalid")
    if formula is None:
        formula=CalcMolFormula(mol)
    atomdata=atomtypes(mol)
    compddict={ID:{'atomdict':atomdata[0],'charge':atomdata[1],'smiles':smiles,'formula':formula,'count':1}}
    return compddict






# def balance_reactions(
#     reactions: List[ReactionIdentifier],
#     ignore_hydrogens: bool = True,
#     use_ray: bool = False,
    
# ) -> List[Tuple[ReactionIdentifier, ReactionIdentifier]]:
#     if use_ray:

#         @ray.remote
#         def balance_reaction_ray(**kwargs):
#             balance_reaction(**kwargs)

#         for reaction in reactions:
#             ray.remote(balance_reaction_ray)
#     else:
#         [balance_reaction(reaction) for reaction in reactions]


# def balance_reaction(
#     input_reaction_identifier: ReactionIdentifier,
# ) -> Tuple[ReactionIdentifier, ReactionIdentifier]:
#     pass

#     # do balancing
#     output_reaction_identifier = None
#     return input_reaction_identifier, output_reaction_identifier




IP = {
    "reagents": [],  # Only if one reaction is inputted
    "solvents": [],  # Only if one reaction is inputted
    "coefflim": 6,  # Maximum tolerable stoichiometric coefficient
    "usemapper": True,  # If mapper needs to be used
    "addrctonly": False,  # If only reactants should be included for balancing
    "ignoreH": False,  # If hydrogens are to be ignored when balancing
    "hc_prod": hc_Dict,  # Help compound dictionary,
    "coefflim": 6,
    "hc_react": None,
    "first": True,
    "ncpus": 1,
    "restart": True,
    "shutdown_after": False,
}  # Help reactant dictionary

def balance_rxn_uspto_df(df: pd.DataFrame, IP=IP, **kwargs):
    if kwargs:
        IP = {**IP, **kwargs}
    if IP["ncpus"] > 1:
        initray(num_cpus=IP["ncpus"], restart=IP["restart"])
        dfdis = mpd.DataFrame(df)
    else:
        dfdis = df
    dfbal = dfdis.apply(
        balance_rxn_uspto_row,
        IP=IP,
        axis=1,
        result_type="reduce",
    )
    dfser = pd.Series(
        data=dfbal.values, index=dfbal.index
    )  # Optional convert modin back to pandas
    dfbal0 = pd.DataFrame(
        data=dfser.tolist(),
        index=dfser.index,
        columns=[
            "rxnsmiles0",
            "Rdata",
            "Pdata",
            "balrxnsmiles",
            "msg",
            "LHSids",
            "RHSids",
            "hcrct",
            "hcprod",
            "LHSdata",
            "RHSdata",
            "Rgtdata",
            "Solvdata",
            "mappedrxn",
            "conf",
            "msg1",
        ],
    )
    if IP["shutdown_after"] and IP["ncpus"] > 1:
        ray.shutdown()
    return dfbal0


def balance_rxn_uspto_row(row: pd.Series, IP=IP, **kwargs):
    if kwargs:
        IP = {**IP, **kwargs}
    (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    ) = balance_rxn_dp(rxn=row["uspto_reaction_smiles"], IP=IP)
    return (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    )


def parse_uspto(
    rxnsmiles: str,
):  # NOTE: MIXTURES NOT WELL REPRESENTED HERE AS REACTION SMILES IS INSUFFICIENT
    reagents = []
    if ">>" not in rxnsmiles:  # Reagents are present
        rxnsmilesgroup = rxnsmiles.split(">")
        reagents = rxnsmilesgroup[1].split(".")
        rxnsmiles = ">>".join([rxnsmilesgroup[0], rxnsmilesgroup[2]])
    return rxnsmiles, reagents

def getfragments(chemlist, smiles=False, ref=None):
    """
    Concatenates smiles strings in a list separated by '.' If ID's are given instead,
    looks up IDs in reference dictionary or database.

    Parameters
    ----------
    chemlist : list
        List of smiles strings OR substance IDs
    ref : dict/database (Optional if list of smiles, mandatory if list of IDs)
        If database provided should be multiindexed with smiles string attached to the ID as key.
        At the moment only fragment database allowed.
    smiles: bool
        Indicate if chemlist is smiles or not

    Returns
    -------
    frag : str
        String joining all smiles given in/calculated from chemlist with '.' operator

    """
    if chemlist == []:
        raise CustomError("ERROR: NO CHEMICALS IN LIST PROVIDED")
    if smiles:  # List of smiles provided
        return ".".join(chemlist)
    else:  # List of IDs given
        if ref is None:
            raise CustomError(
                "Please supply reference (fragment) database multiindexed with ID as second level and containing a Smiles column"
            )
        try:
            frag = ".".join([locrecord(ID, ref, smiles=True) for ID in chemlist])
        except Exception:
            raise CustomError("One or more compounds missing from database")
        else:
            return frag

def locrecord(ID, DB, smiles=False, fragsmarts=False, fragsmiles=False, mixture=False):

    """
    Returns desired information by searching fragment database by ID

    Parameters
    ----------
    ID : int
        Substance ID to search for
    DB : pandas dataframe
        Fragment database with primary index as fragment smiles and secondary index as substance ID
    smiles : bool, optional
        Specify True if only smiles needs to be returned. The default is False.
    fragsmarts : bool, optional
        Specify True if only fragment smarts need to be returned. The default is False.
    fragsmiles : bool, optional
        Specify True if only fragment smiles needs to be returned . The default is False.
    mixture : bool, optional
        Specify True if only mixture presence needs to be returned. The default is False.

    Returns
    -------
    Output : pandas row or str depending on input preference. If everything is False, row will be returned

    """

    record = DB.xs(ID, level=1)
    if smiles:
        return record.Smiles.unique()[0]
    if fragsmarts:
        return record.FragmentSmarts.values
    if fragsmiles:
        return record.index.values
    if mixture:
        return record[">1 Compound"].unique()[0]
    return record


def getspecdat_rxn(
    rxnsmiles: str, reagents: List = [], solvents: List = [], dicts=['Rdata','Pdata','Rgtdata','Solvdata']
) -> Tuple[Dict, Dict, Dict, Dict]:
    """
    A more general function of getspecdat that retrieves species data from a given reaction SMILES.
    Optionally, reagents and/or solvents can be given in a list [SMILES,SMILES etc.]

    Args:
        rxnsmiles (str): Reaction SMILES
        reagents (List, optional): List of reagent SMILES. Defaults to [].
        solvents (List, optional): List of solvent SMILES. Defaults to [].
        dicts (list, optional): List of dictionaries to be returned. Defaults to ['Rdata','Pdata','Rgtdata','Solvdata'].

    Raises:
        Exception: If invalid reaction SMILES is given

    Returns:
        Tuple[Dict, Dict,Dict,Dict]: Tuple of dictionaries of reactant, product, reagent and solvent data
    """
    Rdata={}
    Pdata={}
    Rgtdata={}
    Solvdata={}
    if "Rdata" in dicts or "Pdata" in dicts:
        splitrxn = rxnsmiles.split(">>")
        if len(splitrxn) == 1:  # Only reactants specified
            raise Exception
        if "Rdata" in dicts:
            rcts = splitrxn[0].split(".")
            rcts = [Chem.MolToSmiles(molfromsmiles(rct)) for rct in rcts]
            rcts = Counter(rcts)
            for i, rct in enumerate(rcts):
                Rdata.update(getcompdict(ID=i, smiles=rct))
                Rdata[i]["count"] = rcts[rct]
        if "Pdata" in dicts: 
            prods = splitrxn[1].split(".")
            prods = [Chem.MolToSmiles(molfromsmiles(prod)) for prod in prods]
            prods = Counter(prods)
            for j, prod in enumerate(prods):
                Pdata.update(getcompdict(ID=j, smiles=prod))
                Pdata[j]["count"] = prods[prod]
            
    if reagents and 'Rgtdata' in dicts:
        reagents=[Chem.MolToSmiles(molfromsmiles(rgt)) for rgt in reagents]
        for k, rgt in enumerate(reagents):
            if Rdata:
                k = max(list(Rdata.keys())) + k + 1
            Rgtdata.update(getcompdict(ID=k, smiles=rgt))
        
    if solvents and 'Solvdata' in dicts:
        solvents=[Chem.MolToSmiles(molfromsmiles(solv)) for solv in solvents]
        for f, solv in enumerate(solvents):
            if Rdata and Rgtdata:
                f = max(list(Rdata.keys()) + list(Rgtdata.keys())) + f + 1
            elif Rdata:
                f=max(list(Rdata.keys())) + f + 1
            elif Rgtdata:
                f=max(list(Rgtdata.keys())) + f + 1
            Solvdata.update(getcompdict(ID=f, smiles=solv))
    return Rdata, Pdata, Rgtdata, Solvdata

def balance_rxn_dp(rxn: str = "", IP=IP, rctstorgts=True, **kwargs):
    """_summary_

    Args:
        rxn (str): _description_

    Returns:
        _type_: _description_
    """
    if kwargs:
        IP = {**IP, **kwargs}
    if IP["hc_prod"] is None:
        IP["hc_prod"] = {}
    Rdata = {}
    Rgtdata = {}
    Solvdata = {}
    Pdata = {}
    dicts = [
        "Rdata",
        "Pdata",
        "Rgtdata",
        "Solvdata",
    ]  # Allowable inputs, either {} or {..} output of getcompdict
    try:  # First try with user inputs
        if "Rgtdata" in IP:
            Rgtdata = IP["Rgtdata"]
            dicts.remove("Rgtdata")
        if "Solvdata" in IP:
            Solvdata = IP["Solvdata"]
            dicts.remove("Solvdata")
        if "Rdata" in IP:
            Rdata = IP["Rdata"]
            dicts.remove("Rdata")
        if "Pdata" in IP:
            Pdata = IP["Pdata"]
            dicts.remove("Pdata")
        if dicts:
            if rxn:  # Reaction string is specified
                if ">>" not in rxn:  # with reagents eg. USPTO)
                    rxn, reagents = parse_uspto(rxn)
                    IP["reagents"] = reagents

                output = getspecdat_rxn(
                    rxn, reagents=IP["reagents"], solvents=IP["solvents"], dicts=dicts
                )
                if "Rdata" in dicts:
                    Rdata = output[0]
                if "Pdata" in dicts:
                    Pdata = output[1]
                if "Rgtdata" in dicts:
                    Rgtdata = output[2]
                if "Solvdata" in dicts:
                    Solvdata = output[3]

    except Exception as e:
        msg = "Invalid. Species missing. " + str(e)
        return (
            "Error",
            "Error",
            "Error",
            "Error",
            msg,
            [],
            [],
            [],
            [],
            "Error",
            "Error",
            {},
            {},
            "Error",
            "Error",
            msg,
        )
    # IP["addedspecies"] = [i for i in Rdata]
    rxnsmiles0 = ">>".join(
        [
            getfragments(
                [Rdata[r]["smiles"] for r in Rdata for _ in range(Rdata[r]["count"])],
                smiles=True,
            ),
            getfragments(
                [Pdata[p]["smiles"] for p in Pdata for _ in range(Pdata[p]["count"])],
                smiles=True,
            ),
        ]
    )
    # New
    if "msg" in IP:
        if IP["msg"] and "With hydrogen carriers" in IP["msg"]:
            IP["msg"] = "With hydrogen carriers: " + ",".join(
                [
                    hcarrier
                    for hcarrier in IP["msg"]
                    .split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                    if hcarrier.isdigit()
                ]
            )
        else:
            IP["msg"] = ""
    IP["rxnsmiles0"] = rxnsmiles0
    input = {
        key: IP[key]
        for key in [
            "rxnsmiles0",
            "first",
            "usemapper",
            "addedspecies",
            "hc_prod",
            "hc_react",
            "coefflim",
            "addrctonly",
            "ignoreH",
            "mandrcts",
            "msg",
        ]
        if key in IP
    }

    (
        rxnsmiles0,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
    ) = balancerxn(
        Rdata,
        Pdata,
        Rgtdata=Rgtdata,
        Solvdata=Solvdata,
        **input,
    )
    # print(msg)
    # print(balrxnsmiles)
    # print(Rdata)
    # print(LHSdata)
    mappedrxn_ = maprxn([balrxnsmiles])[0]
    if mappedrxn_ != "Error":
        mappedrxn = mappedrxn_.get("mapped_rxn")
        conf = mappedrxn_.get("confidence")
    else:
        mappedrxn = "Error"
        conf = 0
    if mappedrxn != "Error":
        rseries = pd.DataFrame(
            [
                {
                    "mapped_rxn": mappedrxn,
                    "LHSdata": LHSdata,
                    "RHSdata": RHSdata,
                    "msg": msg,
                }
            ]
        ).iloc[0]
        LHSdata, RHSdata, msg1 = checkrxnrow(rseries)
        if ("Unmapped" in msg1 or "unmapped" in msg1) and (
            "Smiles discrepancy" not in msg1
        ):  # Need to skip if mandatory reactants carries over
            if (
                rctstorgts and "unmapped" in msg1
            ):  # Adding unmapped mandatory reactants as reagents
                Rgtdata.update(
                    {
                        int(ID): copy.deepcopy(Rdata[int(ID)])
                        for ID in msg1.split("unmapped from LHS: ")[1]
                        .split(", ")[0]
                        .split(",")
                        if ID.isdigit()
                    }
                )

            inputdict = {
                "Rdata": Rdata,
                "LHSdata": LHSdata,
                "RHSdata": RHSdata,
                "msg": msg,
                "msg1": msg1,
                "hcprod": hcprod,
                "hcrct": hcrct,
                "Rgtdata": Rgtdata,
                "Solvdata": Solvdata,
            }
            if (
                "mandrcts" in input and input["mandrcts"]
            ):  # External mandatory reactants supplied
                inputdict["Rdata"] = input["mandrcts"]
            # print(inputdict)

            (
                mappedrxn,
                conf,
                balrxnsmiles,
                msg,
                LHSids,
                RHSids,
                hcrct,
                hcprod,
                LHSdata,
                RHSdata,
                msg1,
            ) = updaterxns_(pd.DataFrame([inputdict]).iloc[0], hc_prod=IP["hc_prod"])
    else:
        msg1 = "Mapping error"
    LHSids = [ID for ID in LHSdata for _ in range(LHSdata[ID]["count"])]
    RHSids = [ID for ID in RHSdata for _ in range(RHSdata[ID]["count"])]
    return (
        rxnsmiles0,
        Rdata,
        Pdata,
        balrxnsmiles,
        msg,
        LHSids,
        RHSids,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        Rgtdata,
        Solvdata,
        mappedrxn,
        conf,
        msg1,
    )


# balance_rxn_dp('CCCCCCc1ccc(C(=O)Cl)cc1.CCCCc1ccc(C#Cc2ccc(CNc3ccc4c(c3)C(=O)OC(C)(C)O4)cc2)cc1.Cl>>CCCCCCc1ccc(C(=O)N(Cc2ccc(C#Cc3ccc(CCCC)cc3)cc2)c2ccc3c(c2)C(=O)OC(C)(C)O3)cc1')


def rxn_center_dp(mappedrxn: str, LHSdata: Dict, RHSdata: Dict, expand: int = 1):
    # print(LHSdata)
    # print(RHSdata)
    for rctid in LHSdata:
        userinput = LHSdata[rctid]["smiles"]
        fraglist = getCarrierFrags0(userinput, resFormat="smiles", expand=expand)
        fragloc = {}
        nofg = set()
        for idx, cleanmol in enumerate(LHSdata[rctid]["cleanmol"]):
            # cleanmol = Chem.AddHs(cleanmol)
            for frag in fraglist:
                fragloc, nofg = update_matches(
                    Chem.AddHs(Chem.RemoveAllHs(cleanmol)),
                    frag,
                    fragloc=fragloc,
                    nofg=nofg,
                    idx=idx,
                    rctid=rctid,
                )
        LHSdata[rctid]["fragloc"] = fragloc
    specmap, rnbmap, rxncentermapnum, msg = getrxncenterrow(
        pd.DataFrame(
            [{"mapped_rxn": mappedrxn, "LHSdata": LHSdata, "RHSdata": RHSdata}]
        ).iloc[0],
    )
    if msg:
        LHSdata, msg, outfrag, outfg, outneighbor, unusedanalogue = validrxncenterrow(
            pd.DataFrame(
                [
                    {
                        "specmap": specmap,
                        "rxncentermapnum": rxncentermapnum,
                        "LHSdata": LHSdata,
                        "rnbmap": rnbmap,
                    }
                ]
            ).iloc[0]
        )
    else:
        outfrag = "Error"
        outfg = "Error"
        outneighbor = "Error"
        unusedanalogue = "Error"
    return (
        specmap,
        rnbmap,
        rxncentermapnum,
        LHSdata,
        msg,
        outfrag,
        outfg,
        outneighbor,
        unusedanalogue,
    )


def gen_template_ip(
    LHSdata: Dict,
    RHSdata: Dict,
    specmap: Dict,
    outfrag: Dict = {},
    rnbmap: Dict = {},
    unusedanalogue: List = [],
    specificity="loose",
    processall=True,
):
    template, LHSdata, RHSdata, msg4, farfg, unusedprod = gen_template_row(
        pd.DataFrame(
            [
                {
                    "LHSdata": LHSdata,
                    "RHSdata": RHSdata,
                    "specmap": specmap,
                    "outfrag": outfrag,
                    "rnbmap": rnbmap,
                    "unusedanalogue": unusedanalogue,
                }
            ]
        ).iloc[0],
        specificity=specificity,
        processall=processall,
    )
    return template, LHSdata, RHSdata, msg4, farfg, unusedprod


def updaterxns_(row, hc_prod={}):
    """
    Updates reactions if there are unmapped species and balances if there are changes (optional). Assumes both
    balancerxn, maprxn and checkrxn have all been called already

    rctstorgts: Move unmapped mandatory reactants to reagents

    """
    #     breakpoint()
    msg1 = copy.deepcopy(row["msg1"])
    msg = copy.deepcopy(row["msg"])  # Balanced output message
    LHSdata = copy.deepcopy(row["LHSdata"])
    RHSdata = copy.deepcopy(row["RHSdata"])
    hcprod = copy.deepcopy(row["hcprod"])
    hcrct = copy.deepcopy(row["hcrct"])
    if "Rgtdata" in row.keys():
        Rgtdata = row["Rgtdata"]
    else:
        Rgtdata = {}
    if "Solvdata" in row.keys():
        Solvdata = row["Solvdata"]
    else:
        Solvdata = {}
    if "Rdata" in row.keys():
        mandrcts = row["Rdata"]
    else:
        mandrcts = LHSdata
    addedspecies = list(set(LHSdata.keys()) - set(mandrcts.keys()))
    storemsg = ""
    i = 0
    if "Smiles discrepancy" in msg1:
        raise Exception("Smiles discrepancy")
    while (
        "Unmapped" in msg1 or "unmapped" in msg1
    ) or i == 0:  # Unmapped species exist not reflected
        if "from RHS" or "from LHS" in msg1:  # Mandatory products/reactants unmapped
            if "from LHS" in msg1 and i != 0:  # To avoid infinite loop
                break
            storemsg = msg1
        if hcprod is not None:
            hcprod = [hcprod_ for hcprod_ in hcprod if hcprod_ in RHSdata]
        hcrct = [hcrct_ for hcrct_ in hcrct if hcrct_ in LHSdata]
        balrxnsmiles, msg, LHS, RHS, hcrct, hcprod, LHSdata, RHSdata = balancerxn(
            LHSdata,
            RHSdata,
            first=False,
            Rgtdata=Rgtdata,
            Solvdata=Solvdata,
            addedspecies=addedspecies,
            hc_prod=hc_prod,
            coefflim=6,
            mandrcts=mandrcts,
            usemapper=False,
            ignoreH=False,
        )
        #         if 'Hydrogen carriers' in msg or not hc_prod: #No point balancing again as hydrogen deficit always present
        #             balrxnsmiles,_,LHSids,RHSids,_,_,_,_=update_rxn(LHSdata,RHSdata,hc_prod=hc_prod,hcprod=hcprod,hcrct=hcrct,msg=msg)
        mappedrxn = maprxn([balrxnsmiles])[0]
        if mappedrxn == "Error":
            mapped_rxn = "Error"
            conf = "Error"
            msg1 = "Mapping error"
            break
        else:
            mapped_rxn = mappedrxn.get("mapped_rxn")
            conf = mappedrxn.get("confidence")
            if "With hydrogen carriers" in msg:
                hcarriers = [
                    int(hcarrier)
                    for hcarrier in msg.split("With hydrogen carriers: ")[1]
                    .split(", ")[0]
                    .split(",")
                ]
            else:
                hcarriers = []
            LHSdata, RHSdata, msg1 = checkrxn(
                mapped_rxn,
                Rdata=LHSdata,
                Pdata=RHSdata,
                updateall=True,
                removeunmapped=True,
                mandrcts=mandrcts,
                hcarriers=hcarriers,
            )
        if storemsg:
            if msg1 == "Valid":
                msg1 = storemsg
            elif msg1 != storemsg:
                msg1 = storemsg + ", " + msg1
            break
        i += 1
    return (
        mapped_rxn,
        conf,
        balrxnsmiles,
        msg,
        LHS,
        RHS,
        hcrct,
        hcprod,
        LHSdata,
        RHSdata,
        msg1,
    )