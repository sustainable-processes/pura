from collections import Counter, OrderedDict
import copy
import json
from math import ceil
from typing import List, Optional,Tuple,Dict, Union

import pandas as pd
import modin.pandas as mpd
from rxnmapper import RXNMapper

from .compound import Compound, CompoundIdentifier, CompoundIdentifierType, standardize_identifier
from .reaction import Reaction, ReactionInput, ReactionRole, reaction_from_smiles
from rdkit import Chem
from rdkit.Chem import PeriodicTable, GetPeriodicTable, GetFormalCharge, rdChemReactions
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import ray
from .helpCompound import hc_Dict # Help compounds
from .units import *
from pydantic.tools import parse_obj_as
from chempy import balance_stoichiometry

# This codebase aims to balance a reaction given a reaction object 
# (it is highly recommended not to pass in just the reactions smiles unless mixtures are guaranteed to not be present)




#%% There are many raw data formats of reaction data. This section aims to parse these formats into a standard format given a certain mapping dictionary.

# To be done: Add more parsers for different raw data formats

#%% For parallel execution (TBD)

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
    
    


#%%Base functions for balancing reactions
def gencompdicts(compound_list:List[Compound],start_idx:int=0)->Dict:
    """
    Takes a list of Compound objects and generates a dictionary of compound dictionaries with keys as compound index and values as compound dictionaries

    Args:
        compound_list (List[Compound]): List of Compound objects to generate dictionaries for
        start_idx (int, optional): Starting index for compound dictionaries. Defaults to 0.

    Returns:
        compddicts: Dictionary of compound dictionaries with keys as compound index and values as compound dictionaries
    """
    compddicts={}
    idx=start_idx
    errordict={}
    smiles_list=[compound.to_smiles() for compound in compound_list]
    smiles_counter=Counter(smiles_list)
    smiles_done=set()
    for i,compound in enumerate(compound_list):
        compdsmiles=smiles_list[i]
        if compdsmiles not in smiles_done:
            compddict=getcompdict(compound,smiles=compdsmiles,count=smiles_counter[compdsmiles])
            smiles_done.add(compdsmiles)
            if isinstance(compddict,dict):
                compddicts.update({idx:compddict})
                idx+=1
            else:
                errordict.update({i:compddict})
    return compddicts,errordict
        

def getcompdict(compound:Compound,**kwargs)->Dict: #Can be a method in the compound function; keys in dict can be attributes or properties (TBD). Need to deal with mixtures.
    """
    Generates a compound dictionary from a Compound object by converting to RDKit mol and applying atomtypes function.

    Args:
        compound (Compound): Compound object to generate dictionary for

    Returns:
        Dict: Compound dictionary with atomdict (Output of atomtypes, dictionary with atom elements and count), charge (Overall charge of mol), 
              formula (chemical formula), count (instance number, always 1)
    """
    try:
        compddict={}
        if "mol" not in kwargs:
            mol=compound.to_rdkit_mol()
            compddict.update({'mol':mol})
        if "formula" not in kwargs:
            formula=CalcMolFormula(mol) #Function in RDKit. Can be a method in compound class
            compddict.update({'formula':formula})
        if "atomdata" not in kwargs:
            atomdata=atomtypes(mol)
            compddict.update({'atomdict':atomdata})
        if "charge" not in kwargs:
            charge=GetFormalCharge(mol)
            compddict.update({'charge':charge})
        if "smiles" not in kwargs:
            smiles=compound.to_smiles()
            compddict.update({'smiles':smiles})
        if "count" not in kwargs:
            compddict.update({'count':1})
        if kwargs:
            compddict={**compddict,**kwargs}
        return compddict
    except Exception as e:
        return e
    

def atomtypes(mol:Chem.Mol)->Dict[str,int]:
    """
    Generates an atom type dictionary with counts for a given mol file. Also
    returns overall change of mol. 

    Args:
        mol (RDKit mol): RDKit mol object to generate atom type dictionary for

    Returns:
        Tuple[Dict,int]: Dictionary with keys as atom element and value as count
    """

    typedict={}
    mol2=Chem.AddHs(mol) #Assumes hydrogens are added
    for atom in mol2.GetAtoms():
        elem=PeriodicTable.GetElementSymbol(GetPeriodicTable(),atom.GetAtomicNum())
        if elem not in typedict.keys():
            typedict[elem]=1
        else:
            typedict[elem]+=1
    return typedict

#%% Base functions for reaction and compound handling
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
        raise ValueError("ERROR: NO CHEMICALS IN LIST PROVIDED")
    if smiles:  # List of smiles provided
        return ".".join(chemlist)
    else:  # List of IDs given
        if ref is None:
            raise ValueError(
                "Please supply reference (fragment) database multiindexed with ID as second level and containing a Smiles column"
            )
        try:
            frag = ".".join([locrecord(ID, ref, smiles=True) for ID in chemlist])
        except Exception:
            raise ValueError("One or more compounds missing from database")
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

def buildrxn(Rdata, Pdata):
    """
    Takes reaction and product data, and concatenates smiles strings, forming
    reaction smiles/smarts

    """
    LHS = [Rdata[ID]["smiles"] for ID in Rdata for _ in range(Rdata[ID]["count"])]
    RHS = [Pdata[ID]["smiles"] for ID in Pdata for _ in range(Pdata[ID]["count"])]
    return ">>".join([getfragments(LHS, smiles=True), getfragments(RHS, smiles=True)])

def molfromsmiles(SMILES: str) -> Chem.rdchem.Mol:
    """
    Converts a smiles string into a mol object for RDKit use.

    Args:
        SMILES (str): Smiles string of molecule

    Returns:
        Chem.rdchem.Mol: RDKit mol object
    """
    mol = Chem.MolFromSmiles(SMILES)
    Chem.SanitizeMol(mol)
    mol.UpdatePropertyCache(strict=False)
    return mol

def gensmilesfreq(specdict, validate=True):
    """
    Generates smile frequency dictionary (Sometimes species have the same SMILES with different IDs eg. mixture vs pure)

    """
    smilesfreq = {}
    for ID0 in specdict:
        specsmiles = specdict[ID0]["smiles"]
        specsmiles = specsmiles.split(".")
        for specsmile in specsmiles:
            if validate:
                specsmile = Chem.MolToSmiles(molfromsmiles(specsmile))
            if specsmile in smilesfreq:
                smilesfreq[specsmile].extend([ID0])
            else:
                smilesfreq.update({specsmile: [ID0]})
    return smilesfreq

#%% Main functions for mapping and checking reactions

def maprxn(rxns: List[str]):
    """
    For a given list of reactions, rxns, returns mapped reactions with confidence scores.
    Uses IBM transformer model.

    Parameters
    ----------
    rxns : List[str]
        List of reaction SMILES (no reactant/reagent split)

    Returns
    -------
    Output : list
        Mapped reactions with confidence scores

           mapped_rxn: str
               Mapped reaction SMARTS

           confidence: str
               Model confidence in the mapping rxn

    ['Error']: list
        If code doesn't run or mapper doesn't work
    """

    rxn_mapper = RXNMapper()
    try:
        return rxn_mapper.get_attention_guided_atom_maps(rxns)
    except Exception:
        return ["Error"]

def checkrxnrow(row, updateall=True, removeunmapped=True):
    #     breakpoint()
    mappedrxn = row["mapped_rxn"]
    Rdata = row["LHSdata"]
    Pdata = row["RHSdata"]
    msg = row["msg"]
    if "with species" in msg:
        mandrcts = set(Rdata.keys()) - set(
            [
                int(addedspec)
                for addedspec in msg.rsplit("with species: ", 1)[1]
                .split(" with help product(s): ")[0]
                .split(",")
            ]
        )
    else:
        mandrcts = set(Rdata.keys())
    if "With hydrogen carriers" in msg:
        hcarriers = [
            int(hcarrier)
            for hcarrier in msg.split("With hydrogen carriers: ")[1]
            .split(", ")[0]
            .split(",")
        ]
    else:
        hcarriers = []
    if "Mandatory" in msg:
        mandrcts = mandrcts.union(
            {
                int(mandrct)
                for mandrct in msg.split("Mandatory species unmapped from LHS: ")[1]
                .split(", ")[0]
                .split(",")
            }
        )
    res = checkrxn(
        mappedrxn,
        Rdata=Rdata,
        Pdata=Pdata,
        updateall=updateall,
        removeunmapped=removeunmapped,
        mandrcts=mandrcts,
        hcarriers=hcarriers,
    )

    return res

def updatespecdict(
    refdict,
    smilesfreq,
    cleanmol,
    mappedmol,
    updateddict=OrderedDict({}),
    mixtures={},
    hcarriers=[],
    updateall=True,
):
    """
    Updates species dictionary based on given reactant and cleaned molecule from a reaction
    hcarriers (list of hydrogen containing species involved in reaction but not mapped)
    """
    foundmatch = False
    specsmiles = Chem.MolToSmiles(
        molfromsmiles(Chem.MolToSmiles(cleanmol))
    )  # Ensuring RDKit smiles
    #     breakpoint()
    if specsmiles not in smilesfreq:
        ID0 = ""
        msg = "Smiles discrepancy for species"
    else:
        idx = ""
        msg = "Valid"
        IDlist = smilesfreq[specsmiles]
        mixtures_ = ["." in refdict[ID0]["smiles"] for ID0 in IDlist]
        if len(IDlist) > 1:  # Try mixtures first
            pure = [
                ID0 for i, ID0 in enumerate(IDlist) if not mixtures_[i]
            ]  # Pure matches
            if any(mixtures):
                for i, ID0 in enumerate(IDlist):
                    if mixtures_[i] and ID0 in mixtures:
                        if specsmiles in mixtures[ID0]:
                            loccount = len(mixtures[ID0][specsmiles])
                        else:
                            loccount = 0
                        if any(
                            [
                                len(mixtures[ID0][specsmiles_]) > loccount
                                for specsmiles_ in mixtures[ID0]
                            ]
                        ):
                            idx = i
                            break
                if not idx and not pure:
                    for i, ID0 in enumerate(IDlist):
                        if mixtures_[i] and ID0 not in mixtures:
                            idx = i
                            break
            if not idx and pure:
                for i, ID0 in enumerate(IDlist):
                    if not mixtures_[i] and ID0 not in updateddict:
                        idx = i
                        break
            if not idx:
                idx = 0
        else:
            idx = 0
        ID0 = IDlist[idx]
        # New
        if (
            hcarriers
            and ID0 not in hcarriers
            and not any(
                [atom.HasProp("molAtomMapNumber") for atom in mappedmol.GetAtoms()]
            )
        ):
            return updateddict, mixtures, msg, ID0
        # New
        mixture = mixtures_[idx]
        if updateall:
            mappedsmiles = Chem.MolToSmiles(mappedmol)
            if mixture:
                if ID0 not in mixtures:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    mixtures.update({ID0: {specsmiles: [(mappedsmiles, cleanmol)]}})
                elif specsmiles not in mixtures[ID0]:
                    mixtures[ID0].update({specsmiles: [(mappedsmiles, cleanmol)]})
                else:
                    mixtures[ID0][specsmiles].extend([(mappedsmiles, cleanmol)])
            else:
                if ID0 not in updateddict:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    updateddict[ID0]["count"] = 1
                    updateddict[ID0].update(
                        {"mappedsmiles": [mappedsmiles], "cleanmol": [cleanmol]}
                    )
                else:
                    updateddict[ID0]["count"] += 1
                    updateddict[ID0]["mappedsmiles"].extend([mappedsmiles])
                    updateddict[ID0]["cleanmol"].extend([cleanmol])
        else:
            if mixture:
                if ID0 not in mixtures:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    mixtures.update({ID0: {specsmiles: [()]}})
                elif specsmiles not in mixtures[ID0]:
                    mixtures[ID0].update({specsmiles: [()]})
                else:
                    mixtures[ID0][specsmiles].extend([()])
            else:
                if ID0 not in updateddict:
                    updateddict.update({ID0: copy.deepcopy(refdict[ID0])})
                    updateddict[ID0]["mixture"] = mixture
                    updateddict[ID0]["count"] = 1
                else:
                    updateddict[ID0]["count"] += 1
    return updateddict, mixtures, msg, ID0

def checkrxn(
    mappedrxn,
    Rdata={},
    Pdata={},
    ordered=True,
    updateall=True,
    removeunmapped=True,
    mandrcts=[],
    mandprods=[],
    hcarriers=[],
):  # Assume same rxn smiles stored next to each other
    """
    Checks reaction, updating mapped and clean molecules, removing unmapped species

    """

    #     breakpoint()
    rdrxn = rdChemReactions.ReactionFromSmarts(mappedrxn, useSmiles=True)
    cleanrxn = copy.copy(rdrxn)
    rdChemReactions.RemoveMappingNumbersFromReactions(cleanrxn)
    if ordered:
        LHSdata = OrderedDict({})
        RHSdata = OrderedDict({})
    else:
        LHSdata = {}
        RHSdata = {}
    msgr = []
    msgp = []
    # Updating LHSdata
    if Rdata:
        rsmilesfreq = gensmilesfreq(Rdata)
        smilesmismatch = []
        rmixtures = {}
        mismatch = False
        for ID, rct in enumerate(cleanrxn.GetReactants()):
            mappedmol = rdrxn.GetReactants()[ID]
            formula = Chem.rdMolDescriptors.CalcMolFormula(rct)
            if removeunmapped:
                if (
                    any(
                        [
                            atom.HasProp("molAtomMapNumber")
                            for atom in mappedmol.GetAtoms()
                        ]
                    )
                    or formula == "H2"
                    or hcarriers
                ):  # Confirmed, mapped reactant
                    LHSdata, rmixtures, msg_, ID0 = updatespecdict(
                        Rdata,
                        rsmilesfreq,
                        rct,
                        mappedmol,
                        updateddict=LHSdata,
                        mixtures=rmixtures,
                        updateall=updateall,
                        hcarriers=hcarriers,
                    )
                    if msg_ != "Valid":
                        #                         if ID0 not in smilesmismatch:
                        #                             smilesmismatch+=[ID0]
                        mismatch = True
            else:
                LHSdata, rmixtures, msg_, ID0 = updatespecdict(
                    Rdata,
                    rsmilesfreq,
                    rct,
                    mappedmol,
                    updateddict=LHSdata,
                    mixtures=rmixtures,
                    updateall=updateall,
                    hcarriers=hcarriers,
                )
                if msg_ != "Valid":
                    #                     if ID0 not in smilesmismatch:
                    #                         smilesmismatch+=[ID0]
                    mismatch = True
        if mismatch:
            smilesmismatch = [
                ID0 for ID0 in Rdata if ID0 not in LHSdata and ID0 not in rmixtures
            ]
            if smilesmismatch:
                msgr += [
                    "Smiles discrepancy for LHS species: "
                    + ",".join([str(ID) for ID in smilesmismatch])
                ]
            else:
                msgr += ["Smiles discrepancy for LHS species"]
        #         breakpoint()
        if rmixtures:
            msgr += [
                "Mixture detected for LHS species: "
                + ",".join([str(ID) for ID in rmixtures])
            ]
            for ID0 in rmixtures:
                localcounts = [
                    int(Counter(rsmilesfreq[mixsmiles])[ID0])
                    for mixsmiles in rmixtures[ID0]
                ]
                numinsts = [
                    len(rmixtures[ID0][mixsmiles]) for mixsmiles in rmixtures[ID0]
                ]
                lb = [0 for j in range(len(numinsts))]
                ub = [0 for j in range(len(numinsts))]
                div = [
                    int(ceil(numinst / localcount))
                    for numinst, localcount in zip(numinsts, localcounts)
                ]
                count = max(div)
                if updateall:
                    LHSdata[ID0].update(
                        {"mappedsmiles": [], "cleanmol": [], "unmappedmix": []}
                    )
                    for i in range(count):
                        mappedsmiles = []
                        unmappedsmiles = []
                        cleanmol = []
                        for j, mixsmiles in enumerate(rmixtures[ID0]):
                            ub_ = min(numinsts[j], localcounts[j])
                            ub[j] = ub_ + lb[j]
                            mappedlist = rmixtures[ID0][mixsmiles][lb[j] : ub[j]]
                            mappedsmiles += [comb[0] for comb in mappedlist]
                            cleanmol += [comb[1] for comb in mappedlist]
                        mappedsmiles = tuple(mappedsmiles)
                        cleanmol = tuple(cleanmol)
                        unmappedsmiles = tuple(
                            [
                                mixsmiles
                                for mixsmiles in rsmilesfreq
                                for k in range(
                                    int(Counter(rsmilesfreq[mixsmiles])[ID0])
                                )
                                if mixsmiles in Rdata[ID0]["smiles"]
                                if mixsmiles not in rmixtures[ID0]
                            ]
                        )
                        LHSdata[ID0]["mappedsmiles"].extend([mappedsmiles])
                        LHSdata[ID0]["cleanmol"].extend([cleanmol])
                        LHSdata[ID0]["unmappedmix"].extend([unmappedsmiles])
                        lb = copy.deepcopy(ub)
                LHSdata[ID0]["count"] = count
        if removeunmapped:
            rem = Counter()
            Rspecies = Counter(
                [ID0 for ID0 in Rdata for _ in range(Rdata[ID0]["count"])]
            )
            LHSspecies = Counter(
                [ID0 for ID0 in LHSdata for _ in range(LHSdata[ID0]["count"])]
            )
            rem.update(Rspecies)
            rem.subtract(LHSspecies)
            #             removedmandrcts=[ID0 for ID0 in rem if ID0 in mandrcts if ID0 not in LHSspecies]
            #             removedrcts=[ID0 for ID0 in rem if ID0 not in removedmandrcts if ID0 not in smilesmismatch if ID0 in Rspecies if ID0 not in LHSspecies]
            removedrcts = [
                ID0
                for ID0 in rem
                if ID0 in Rspecies
                if ID0 not in LHSspecies
                if ID0 not in smilesmismatch
            ]
            removedmandrcts = [ID0 for ID0 in removedrcts if ID0 in mandrcts]
            # New
            removedmandrcts = list(
                set(removedmandrcts).union(
                    {mandrct for mandrct in mandrcts if mandrct not in LHSspecies}
                )
            )
            # New
            runmapped = [
                ID0
                for ID0 in rem
                if rem[ID0] > 0
                if ID0 not in removedrcts
                if ID0 not in smilesmismatch
            ]  #!=0
            if removedmandrcts:
                msgr += [
                    "Mandatory species unmapped from LHS: "
                    + ",".join([str(ID) for ID in removedmandrcts])
                ]
                removedrcts = [ID0 for ID0 in removedrcts if ID0 not in removedmandrcts]
            if removedrcts:
                msgr += [
                    "Unmapped species from LHS: "
                    + ",".join([str(ID) for ID in removedrcts])
                ]
            if runmapped:
                msgr += [
                    "Unmapped species instances from LHS: "
                    + ",".join([str(ID) for ID in runmapped])
                ]
    if Pdata:
        #         breakpoint()
        psmilesfreq = gensmilesfreq(Pdata)
        smilesmismatch = []
        pmixtures = {}
        mismatch = False
        for ID, prod in enumerate(cleanrxn.GetProducts()):
            mappedmol = rdrxn.GetProducts()[ID]
            #             formula=rdkit.Chem.rdMolDescriptors.CalcMolFormula(prod)
            if removeunmapped:
                if any(
                    [atom.HasProp("molAtomMapNumber") for atom in mappedmol.GetAtoms()]
                ):  # Confirmed, mapped reactant
                    RHSdata, pmixtures, msg_, ID0 = updatespecdict(
                        Pdata,
                        psmilesfreq,
                        prod,
                        mappedmol,
                        updateddict=RHSdata,
                        mixtures=pmixtures,
                        updateall=updateall,
                    )
                    if msg_ != "Valid":
                        #                         if ID0 not in smilesmismatch:
                        #                             smilesmismatch+=[ID0]
                        mismatch = True
            else:
                RHSdata, pmixtures, msg_, ID0 = updatespecdict(
                    Pdata,
                    psmilesfreq,
                    prod,
                    mappedmol,
                    updateddict=RHSdata,
                    mixtures=pmixtures,
                    updateall=updateall,
                )
                if msg_ != "Valid":
                    #                     if ID0 not in smilesmismatch:
                    #                         smilesmismatch+=[ID0]
                    mismatch = True
        if mismatch:
            smilesmismatch = [
                ID0 for ID0 in Pdata if ID0 not in RHSdata and ID0 not in pmixtures
            ]
            if smilesmismatch:
                msgp += [
                    "Smiles discrepancy for RHS species: "
                    + ",".join([str(ID) for ID in smilesmismatch])
                ]
            else:
                msgp += ["Smiles discrepancy for RHS species"]
        if pmixtures:
            msgp += [
                "Mixture detected for RHS species: "
                + ",".join([str(ID) for ID in pmixtures])
            ]
            #             breakpoint()
            for ID0 in pmixtures:
                localcounts = [
                    int(Counter(psmilesfreq[mixsmiles])[ID0])
                    for mixsmiles in pmixtures[ID0]
                ]
                numinsts = [
                    len(pmixtures[ID0][mixsmiles]) for mixsmiles in pmixtures[ID0]
                ]
                lb = [0 for j in range(len(numinsts))]
                ub = [0 for j in range(len(numinsts))]
                div = [
                    int(ceil(numinst / localcount))
                    for numinst, localcount in zip(numinsts, localcounts)
                ]
                count = max(div)
                if updateall:
                    RHSdata[ID0].update(
                        {"mappedsmiles": [], "cleanmol": [], "unmappedmix": []}
                    )
                    for i in range(count):
                        mappedsmiles = []
                        unmappedsmiles = []
                        cleanmol = []
                        for j, mixsmiles in enumerate(pmixtures[ID0]):
                            ub_ = min(numinsts[j], localcounts[j])
                            ub[j] = ub_ + lb[j]
                            mappedlist = pmixtures[ID0][mixsmiles][lb[j] : ub[j]]
                            mappedsmiles += [comb[0] for comb in mappedlist]
                            cleanmol += [comb[1] for comb in mappedlist]
                        mappedsmiles = tuple(mappedsmiles)
                        cleanmol = tuple(cleanmol)
                        unmappedsmiles = tuple(
                            [
                                mixsmiles
                                for mixsmiles in psmilesfreq
                                for k in range(
                                    int(Counter(psmilesfreq[mixsmiles])[ID0])
                                )
                                if mixsmiles in Pdata[ID0]["smiles"]
                                if mixsmiles not in pmixtures[ID0]
                            ]
                        )
                        RHSdata[ID0]["mappedsmiles"].extend([mappedsmiles])
                        RHSdata[ID0]["cleanmol"].extend([cleanmol])
                        RHSdata[ID0]["unmappedmix"].extend([unmappedsmiles])
                        lb = copy.deepcopy(ub)
                RHSdata[ID0]["count"] = count
        if removeunmapped:
            rem = Counter()
            Pspecies = Counter(
                [ID0 for ID0 in Pdata for _ in range(Pdata[ID0]["count"])]
            )
            RHSspecies = Counter(
                [ID0 for ID0 in RHSdata for _ in range(RHSdata[ID0]["count"])]
            )
            rem.update(Pspecies)
            rem.subtract(RHSspecies)
            #             removedmandprod=[ID0 for ID0 in rem if ID0 in mandprods if ID0 not in RHSspecies]
            #             removedprods=[ID0 for ID0 in rem if ID0 not in removedmandprod if ID0 not in smilesmismatch if ID0 in Pspecies if ID0 not in RHSspecies]
            removedprods = [
                ID0
                for ID0 in rem
                if ID0 in Pspecies
                if ID0 not in RHSspecies
                if ID0 not in smilesmismatch
            ]
            removedmandprods = [ID0 for ID0 in removedprods if ID0 in mandprods]
            punmapped = [
                ID0
                for ID0 in rem
                if rem[ID0] > 0
                if ID0 not in removedprods
                if ID0 not in smilesmismatch
            ]  #!=0
            if removedmandprods:
                msgp += [
                    "Mandatory species unmapped from RHS: "
                    + ",".join([str(ID) for ID in removedmandprods])
                ]
                removedprods = [
                    ID0 for ID0 in removedprods if ID0 not in removedmandprods
                ]
            if removedprods:
                msgp += [
                    "Unmapped species from RHS: "
                    + ",".join([str(ID) for ID in removedprods])
                ]
            if punmapped:
                msgp += [
                    "Unmapped species instances from RHS: "
                    + ",".join([str(ID) for ID in punmapped])
                ]
    msg = msgr + msgp
    if not msg:
        msg = "Valid"
    else:
        msg = ", ".join(msg)
    return LHSdata, RHSdata, msg


#%% Main functions for balancing

IP = {
    "reagents": [],  # Only if one reaction is inputted
    "solvents": [],  # Only if one reaction is inputted
    "coefflim": 6,  # Maximum tolerable stoichiometric coefficient
    "usemapper": True,  # If mapper needs to be used
    "addrctonly": False,  # If only reactants should be included for balancing
    "ignoreH": False,  # If hydrogens are to be ignored when balancing
    "hc_prod": hc_Dict,  # Help compound dictionary
    "coefflim": 6,      # Maximum tolerable stoichiometric coefficient
    "hc_react": None,   # Help compound dictionary for reactants
    "first": True,      # If first iteration
    "ncpus": 1,         # Number of CPUs for parallel execution
    "restart": True,    # Restart existing Ray cluster
    "shutdown_after": False,  # Shutdown Ray cluster after execution
} 

def balance_reaction(reaction:Reaction, IP:Optional[Dict]=IP, rctstorgts=True, **kwargs):
    """_summary_

    Args:
        reaction (Reaction): Reaction object to balance
        IP (Optional[Dict], optional): Input parameters. Defaults to IP.
        **kwargs: Additional input parameters. IP will be updated with these.

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
    Agtdata={}
    rgterrordict={}
    solerrordict={}
    agterrordict={}
    
    # Generate reactant and product dictionaries (raise error if SMILES issue or missing reactants/products)
    if getattr(reaction,"reactant_compounds"):
        Rdata,rerrordict=gencompdicts(reaction.reactant_compounds)
        if rerrordict:
            errormsg="Error in reactant compounds: "+", ".join([f"Compound {key}: {val}" for key,val in rerrordict.items()])
            raise ValueError(errormsg)
    else:
        raise ValueError("No reactants found in reaction")
    if getattr(reaction,"product_compounds"):
        Pdata,perrordict=gencompdicts(reaction.product_compounds,start_idx=len(Rdata))
        if perrordict:
            errormsg="Error in product compounds: "+", ".join([f"Compound {key}: {val}" for key,val in perrordict.items()])
            raise ValueError(errormsg)
    else:
        raise ValueError("No products found in reaction")
    # Generate agent dictionaries
    if getattr(reaction,"reagent_compounds"):
        Rgtdata,rgterrordict=gencompdicts(reaction.reagent_compounds,start_idx=len(Rdata)+len(Pdata))
    if getattr(reaction,"solvent_compounds"):
        Solvdata,solerrordict=gencompdicts(reaction.solvent_compounds,start_idx=len(Rdata)+len(Pdata)+len(Rgtdata))
    Agtdata={**Rgtdata,**Solvdata}
    if getattr(reaction,"agent_compounds"):
        Agtdata2,agterrordict=gencompdicts(reaction.agent_compounds,start_idx=len(Rdata)+len(Pdata)+len(Rgtdata)+len(Solvdata))
        agtsmiles_list=[agtdict['smiles'] for agtdict in Agtdata.values()]
        Agtdata.update({id:agtdict for id,agtdict in Agtdata2.items() if agtdict['smiles'] not in agtsmiles_list})
    
    # return Rdata,Pdata,Agtdata,rgterrordict,solerrordict,agterrordict
                                                  
    rxnsmiles0=reaction.reaction_smiles()
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

def update_stoich(stoich, compdict, hcID=None, hc_Dict=None):
    """
    Based on balanced stoichiometry output of balance_stoichiometry function from chempy, and given a
    dictionary of help compounds and relevant help IDs, updates species dictionary
    """
    usedid = []
    formdict = {}
    msg = ""
    for ID in compdict:
        form = compdict[ID]["formula"]
        if form not in formdict:
            formdict.update({form: [ID]})
        else:
            formdict[form].extend([ID])
    #     breakpoint()
    if hcID:
        if hc_Dict is None:
            raise ValueError(
                "Please supply help compound reference dictionary/dataframe"
            )
        for hcid in hcID:
            form = hc_Dict[hcid]["formula"]
            if form not in formdict:
                formdict.update({form: [hcid]})
            else:
                formdict[form].extend([hcid])
    for formula, coeff in stoich.items():
        if formula in formdict:
            for ID in formdict[formula]:
                if ID not in compdict:
                    compdict.update({ID: hc_Dict[ID]})
                compdict[ID]["count"] = coeff
                usedid += [ID]
        else:
            msg = "Invalid balancing. Formula indicated in stoich outside compound dictionary"
            break
    if msg:
        return "Error", msg, formdict
    else:
        valid = True
        unusedid = set(compdict.keys()) - set(usedid)
        if unusedid:
            for ID in unusedid:
                if "rxs_ids" not in compdict[ID]:  # Identifying help compounds
                    valid = False
                    break
            if valid:
                for ID in unusedid:
                    del compdict[
                        ID
                    ]  # 'Reactants' that aren't reactants..these are removed
        if valid and compdict:
            return compdict, msg, formdict
        else:
            #             breakpoint()
            msg = "Invalid balancing. Species missing: " + ",".join(
                [str(unused) for unused in unusedid]
            )
            return compdict, msg, formdict

def update_rxn(
    Rdata,
    Pdata,
    reac=None,
    prod=None,
    hc_prod=None,
    hcprod=[],
    hcrct=[],
    rxnsmiles0=None,
    msg=None,
):
    """
    Wrapper for calling update_stoich function

    """
    stoichupdated = False
    addmsgr = ""
    addmsgp = ""
    if reac is not None:  # Switched order
        stoichupdated = True
        Rdata1, addmsgr, formdictr = update_stoich(reac, Rdata)
    if prod is not None:
        stoichupdated = True
        Pdata1, addmsgp, formdictp = update_stoich(
            prod, Pdata, hcID=hcprod, hc_Dict=hc_prod
        )
    if addmsgr and addmsgr != "Valid":
        if msg is not None:
            msg = addmsgr + " from LHS" + ", " + msg
        else:
            msg = addmsgr + " from LHS"
    if addmsgp and addmsgp != "Valid":
        if msg is not None:
            msg = addmsgp + " from RHS" + ", " + msg
        else:
            msg = addmsgp + " from RHS"
    if stoichupdated:
        if Rdata1 != "Error" and Pdata1 != "Error":
            Rdata = Rdata1
            Pdata = Pdata1
            try:
                balrxnsmiles = buildrxn(Rdata, Pdata)
            except Exception:
                balrxnsmiles = "Error"
                msg = (
                    "Invalid balancing. Species missing in database" + ", " + msg
                )  # Just to make sure, although this error should never happen
        else:
            balrxnsmiles = "Error"
    else:
        if (
            ("LHS species insufficient" in msg)
            | ("Invalid" in msg)
            | ("Mapping error" in msg)
            | ("discrepancy" in msg)
        ):
            balrxnsmiles = "Error"
        else:
            balrxnsmiles = buildrxn(Rdata, Pdata)

    if hcrct:
        LHSids = [
            ID for ID in Rdata if ID not in hcrct for _ in range(Rdata[ID]["count"])
        ]
    else:
        LHSids = [ID for ID in Rdata for _ in range(Rdata[ID]["count"])]
    if hcprod:
        RHSids = [
            ID for ID in Pdata if ID not in hcprod for _ in range(Pdata[ID]["count"])
        ]
    else:
        RHSids = [ID for ID in Pdata for _ in range(int(Pdata[ID]["count"]))]
    if rxnsmiles0 is not None:
        return (
            rxnsmiles0,
            balrxnsmiles,
            msg,
            LHSids,
            RHSids,
            hcrct,
            hcprod,
            Rdata,
            Pdata,
        )  # Same number and type of atoms reactant and product side and same charge ie. perfectly balanced reaction. Pretty much impossible.
    else:
        return balrxnsmiles, msg, LHSids, RHSids, hcrct, hcprod, Rdata, Pdata

def findmatch(atomdeficit, atomdict, strict=True, returnmultdict=True):
    """
    Calculates proportion of atoms mapped, based on given atom deficit dictionary and
    atom dictionary of a candidate
    """
    #     breakpoint()
    if not set(atomdict.keys()).intersection(
        set(atomdeficit.keys())
    ):  # No match at all
        return False, False
    rem2 = Counter()
    rem2.update(atomdict)
    rem2.subtract(Counter(atomdeficit))
    if any([val < 0 for val in rem2.values()]):
        multdict = {
            k: abs(atomdeficit[k] / atomdict[k]) for k in atomdeficit if k in atomdict
        }
        if (atomdeficit.keys() == atomdict.keys()) & (
            len(set(Counter(multdict).values())) == 1
        ):
            return 1.0, 1  # Exact multiple
        elif returnmultdict:
            return False, multdict
        else:
            if strict:
                mult = int(ceil(max(multdict.values())))
            else:
                mult = int(ceil(min(multdict.values())))
            return False, mult
    mapprop = (sum(atomdeficit.values())) / (sum(atomdict.values()))
    #     mapprop=1-(sum(rem2.values())/sum(atomdict.values()))
    return round(mapprop, 1), 1

def resolvecandidates(
    postype,
    Rdata,
    specdict,
    candidates,
    Pdata,
    update=True,
    validate=True,
    rctonly=False,
    coefflim=6,
    ignoreH=False,
):
    """
    Resolves candidates based on atom deficit (postype) and supplied candidate matches


    """

    #     breakpoint()
    msg = ""
    if validate:
        combinedkeys = {
            key for candi in candidates for key in specdict[candi]["atomdict"].keys()
        }
        if not set(postype.keys()).issubset(
            combinedkeys
        ):  # Reactants/reagents/solvents cannot account for atom imbalance
            msg = "LHS species insufficient"
            return Rdata, candidates, msg
    if len(candidates) > 1:
        matches = []
        mult = []
        for candi in candidates:
            match, mult_ = findmatch(postype, specdict[candi]["atomdict"])
            matches += [match]
            mult += [mult_]
        if 1 in matches:
            index_max = [i for i, match in enumerate(matches) if match == 1]
            candidates = [candidates[idx] for idx in index_max]
            mult = [mult[idx] for idx in index_max]
        else:
            if len(candidates) > 1:
                if len(postype) == 1 and "H" in postype and ignoreH:
                    msg = "Hydrogen carriers: " + ",".join(
                        [str(candi) for candi in candidates]
                    )
                    return Rdata, candidates, msg
                if all(match is not False for match in matches):
                    if (
                        len(postype) == 1 and "H" in postype
                    ):  # Deficit is only hydrogen, so mapper will not help
                        counter = Counter(matches)
                        maxmatch = max(matches)
                        if counter[maxmatch] == 1:
                            index_max = matches.index(maxmatch)
                            candidates = [candidates[index_max]]
                            mult = [mult[index_max]]
                    elif len(set(matches)) == len(
                        matches
                    ):  # Possibility of mapping wrong species still there
                        index_max = max(
                            range(len(matches)), key=matches.__getitem__
                        )  # Species with maximum atoms mappable
                        candidates = [candidates[index_max]]
                        mult = [mult[index_max]]
                else:  # Higher stoichiometric coefficients needed or more than one species needed
                    #                     breakpoint()
                    atompop = {k: [] for k in postype}
                    for candi, mult_ in zip(candidates, mult):
                        if type(mult_) == dict:
                            for k in mult_:
                                if k in atompop:
                                    atompop[k].extend([candi])
                        else:
                            for k in specdict[candi]["atomdict"]:
                                if k in atompop:
                                    atompop[k].extend([candi])
                    #                     breakpoint()
                    for k in sorted(
                        atompop, key=lambda k: (len(atompop[k]), postype[k])
                    ):
                        if rctonly:
                            extrarct = [
                                rctid
                                for rctid in Rdata
                                if k in Rdata[rctid]["atomdict"]
                                if rctid not in candidates
                            ]
                        else:
                            extrarct = [
                                rctid
                                for rctid in specdict
                                if k in specdict[rctid]["atomdict"]
                                if rctid not in candidates
                            ]
                        if extrarct:
                            for rctid in extrarct:
                                mult_ = int(
                                    ceil(postype[k] / specdict[rctid]["atomdict"][k])
                                )
                                if mult_ <= 1:
                                    mult_ = 1
                                    matches += [True]
                                else:
                                    matches += [False]
                                mult += [mult_]
                            candidates += extrarct
                            atompop[k].extend(extrarct)
                        if len(atompop[k]) == 1:  # Only one candidate for atom
                            candi = atompop[k][0]
                            mult_ = mult[candidates.index(candi)]
                            if type(mult_) == dict:
                                mult = [int(ceil(mult_[k]))]
                            else:
                                mult = [mult_]
                            candidates = [candi]
                            break
                        else:
                            if len(postype) == 1 and "H" in postype:
                                msg = "Hydrogen carriers: " + ",".join(
                                    [str(candi) for candi in candidates]
                                )
                                return Rdata, candidates, msg
                            matches2 = []
                            candidates2 = []
                            for candi in atompop[k]:
                                idx = candidates.index(candi)
                                mult_ = mult[idx]
                                if type(mult_) == dict:
                                    mult_ = int(ceil(mult_[k]))
                                mult[idx] = mult_
                                match = matches[idx]
                                if match == False:
                                    atomdict = copy.deepcopy(
                                        specdict[candi]["atomdict"]
                                    )
                                    totatomdict = {
                                        j: atomdict[j] * mult_ for j in atomdict
                                    }
                                    totmapped = {
                                        j: min(totatomdict[j], postype[j])
                                        for j in postype
                                        if j in atomdict
                                    }
                                    mapprop = round(
                                        (sum(totmapped.values()))
                                        / (sum(totatomdict.values())),
                                        1,
                                    )
                                    matches[idx] = mapprop
                                    matches2 += [mapprop]
                                else:
                                    matches2 += [match]
                            counter = Counter(matches2)
                            maxmatch = max(matches2)
                            if rctonly:
                                candidates2 = [
                                    candi for candi in atompop[k] if candi in Rdata
                                ]
                            if not candidates2:
                                candidates2 = [
                                    candi
                                    for candi in atompop[k]
                                    if matches[candidates.index(candi)] == maxmatch
                                    or candi in Rdata
                                ]
                            mult2 = [
                                mult[candidates.index(candi)] for candi in candidates2
                            ]
                            candidates = candidates2
                            mult = mult2
                            break
    else:
        match, _ = findmatch(postype, specdict[candidates[0]]["atomdict"])
        matches=[match]
        mult = [1]
    print(postype)
    print(candidates)
    if (len(candidates) > 1 or (ignoreH and len(candidates)==1 and 1 not in matches)) and "Hydrogen carriers" not in msg:  # or ignoreH
        if len(postype) == 1 and "H" in postype:  # Still multiple options
            msg = "Hydrogen carriers: " + ",".join([str(candi) for candi in candidates])
            return Rdata, candidates, msg
    #     breakpoint()
    msg = "Valid"
    if update:
        for candi, mult_ in zip(candidates, mult):
            if candi in Rdata.keys():
                Rdata[candi]["count"] += (mult_*specdict[candi]['count']) #Added * specdict[count]
            else:
                Rdata.update({candi: specdict[candi]})
                Rdata[candi]["count"] = mult_*specdict[candi]['count'] #Added * specdict[count]
    return Rdata, candidates, msg

def balance(
    Rdata,
    Pdata,
    hc_prod={},
    balbefore=True,
    coefflim=6,
    addedspecies=[],
    addedhc=[],
    hc_react={},
    mandrcts={},
):
    """
    Balances reaction given LHS and RHS species by invoking the balance_stoichiometry function from ChemPy

    """
    #     breakpoint()
    chempyr = [Rdata[ID]["formula"] for ID in Rdata]
    chempyp = [Pdata[ID]["formula"] for ID in Pdata]
    #     chempyr=[Rdata[ID]['formula'] for ID in Rdata for _ in range(Rdata[ID]['count'])]
    #     chempyp=[Pdata[ID]['formula'] for ID in Pdata for _ in range(Pdata[ID]['count'])]
    if len(set(chempyr + chempyp)) < len(
        chempyr + chempyp
    ):  # Isomers/same formulae present
        raise ValueError("Isomers detected. Balancing will not work")
    chempyr = set(chempyr)
    chempyp = set(chempyp)
    #     highcoeff=False
    reac0 = {}
    prod0 = {}
    msg = "Balanced"
    addedstr = ""
    if addedspecies:
        addedstr = ",".join(
            [str(species) for species in set(addedspecies) if species not in mandrcts]
        )
        if addedstr:
            addedstr = " with species: " + addedstr
    if addedhc:
        addedstr2 = " with help reactant(s): " + (
            ",".join([hc_react[species]["formula"] for species in addedhc])
        )
        if addedstr:
            addedstr=addedstr + ", " + addedstr2
        else:
            addedstr = addedstr2
    msg += addedstr

    if (
        balbefore or not hc_prod
    ):  # Either user has indicated that first pass at balancing needs to be done or fails to specify help compounds
        try:
            reac, prod = balance_stoichiometry(
                chempyr, chempyp, underdetermined=None, allow_duplicates=True
            )  # Try balancing once without adding compounds
            if any(idx < 0 for idx in reac.values()) or any(
                idx < 0 for idx in prod.values()
            ):  # Don't want negative stoich coefficients
                raise Exception
            elif any(
                [
                    idx > coefflim
                    for tup in zip(reac.values(), prod.values())
                    for idx in tup
                ]
            ):
                raise Exception

        except Exception as e:
            if not hc_prod:  # User has not supplied help compounds
                raise ValueError(
                    "Reaction does not balance. Help compounds not provided"
                )
            pass
        else:  # Reaction successfully balanced
            return reac, prod, None, msg
    #     breakpoint()
    if hc_prod:  # Can try help compounds
        try:
            reac, prod, hcid = tryhelp(hc_prod, chempyr, chempyp, coefflim=coefflim)
        except Exception as e:
            raise ValueError("Reaction does not balance even with help compounds")
        else:
            hclist = ",".join([hc_prod[hc]["formula"] for hc in hcid])
            return reac, prod, hcid, msg + " with help product(s): " + hclist

def tryhelp(hc_atomtype, chempyr, chempyp, coefflim=6):
    """
    Attempts to balance reaction with addition of help compounds in helpCompounds.py.

    hc_atomtype is a dictionary or list of help compounds
    chempyr is a set of LHS species formulae
    chemyp is a set of RHS species formulae
    coefflim is the maximum acceptable stoichiometric coefficient value after balancing

    Returns reac,prod (outputs of balance_stoichiometry function in chempy) and a list of help compounds added

    """
    #     breakpoint()
    reac = {}
    prod = {}
    hcid = None
    lim = len(hc_atomtype)
    keylist = [hcid for hcid in hc_atomtype]
    counter = 0
    invalid = True
    while (
        any(
            [idx > coefflim for tup in zip(reac.values(), prod.values()) for idx in tup]
        )
        or not reac
        or invalid
    ):
        if counter > lim - 1:
            print("Help compounds did not help. Extra reactant atoms")
            raise ValueError("Help compounds unable to balance")
        if hcid is not None:
            chempyp.remove(hc_atomtype[hcid]["formula"])
        hcid = keylist[counter]
        chempyp.add(hc_atomtype[hcid]["formula"])
        try:
            reac, prod = balance_stoichiometry(
                chempyr, chempyp, underdetermined=None, allow_duplicates=True
            )
            if any(idx < 0 for idx in reac.values()) or any(
                idx < 0 for idx in prod.values()
            ):  # Don't want negative stoich coefficients
                invalid = True
                raise Exception
            else:
                invalid = False
                counter += 1
                continue
        except Exception:
            counter += 1
            continue
    print("Reaction successfully balanced")
    return reac, prod, [hcid]



def balancerxn(
    Rdata,
    Pdata,
    Rgtdata={},
    Solvdata={},
    rxnsmiles0=None,
    first=True,
    usemapper=True,
    addedspecies=[],
    addedhc=[],
    hc_prod={},
    hc_react={},
    coefflim=6,
    msg="",
    mandrcts={},
    addrctonly=False,
    ignoreH=False,
):
    """
    Balances reactions given reactant species information (Rdata) and product species information (Pdata)

    Rgtdata is optional and refers to reagent species information
    Solvdata is optional and refers to solvent species information
    rxnsmiles0 refers to the reaction SMILES string as represented in Reaxys
    first is True/False depending on whether it is the first time running through the code
    usemapper is True if IBM RXN mapper is used to decide between possible reactants on LHS
    addedspecies refers to new species added to the LHS
    addedhc refers to small species (help reactants) added to the LHS
    hc_prod is optional and is a dictionary of small species (help compounds) to use for balancing the RHS
    hc_react is optional and is a dictionary of small species (help reactants) to use for balancing the LHS
    coefflim refers to the maximum allowed stoichiometric coefficient after balancing
    msg involves any warning messages or updates from the code
    mandrcts is a list of mandatory reactants (To avoid balancer from removing them)
    ignoreH refers to a boolean switch (True if all hydrogen species except H2 are ignored/not added)

    """

    if Solvdata:  # Add sovents
        Rgtdata = {**Rgtdata, **Solvdata}  # Newly added

    #%% Initialize added species/addedhc (help compound) and add small species
    if not mandrcts:
        mandrcts = copy.deepcopy(Rdata)
    if first:
        addedspecies = []
        addedhc = []
        # if not mandrcts:
        #     mandrcts = copy.deepcopy(Rdata)
        if Rgtdata:
            #             smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH'] if spec not in Solvdata]
            #             smallspecies=[spec for spec in Rgtdata if Rgtdata[spec]['formula'] in ['H2','O2','H2O','OH2','HOH']]
            smallspecies = [
                spec for spec in Rgtdata if Rgtdata[spec]["formula"] in ["H2", "O2"]
            ]
            if smallspecies:
                Rdata.update({spec: Rgtdata[spec] for spec in smallspecies})
                addedspecies += smallspecies

    #%% String output handling
    addedstr = ""
    if addedspecies:
        addedstr = ",".join(
            [str(species) for species in set(addedspecies) if species not in mandrcts]
        )
        if addedstr:
            addedstr = " with species: " + addedstr
    if addedhc:
        addedstr2 = ",".join(
            [
                hc_react[species]["formula"]
                for species in set(addedhc)
                if species not in mandrcts
            ]
        )
        if addedstr2:
            addedstr2 = " with help reactant(s): " + addedstr2
        if addedstr:
            addedstr = addedstr + ", " + addedstr2
        else:
            addedstr = addedstr2
    #     if 'Mandatory' in msg or 'Smiles discrepancy' in msg:
    if "Smiles discrepancy" in msg:
        msg = msg + addedstr
        return update_rxn(
            mandrcts,
            Pdata,
            hc_prod=hc_prod,
            hcrct=addedhc,
            rxnsmiles0=rxnsmiles0,
            msg=msg,
        )
    if "Hydrogen carriers" in msg:
        msg = msg + addedstr
        return update_rxn(
            Rdata, Pdata, hc_prod=hc_prod, hcrct=addedhc, rxnsmiles0=rxnsmiles0, msg=msg
        )

    Rcount = sum(
        [
            Counter(Rdata[ID]["atomdict"])
            for ID in Rdata
            for _ in range(Rdata[ID]["count"])
        ],
        start=Counter(),
    )  # Sum of atom counts/types on LHS
    Rcharge = sum(
        [Rdata[ID]["charge"] for ID in Rdata for _ in range(Rdata[ID]["count"])]
    )
    Pcount = sum(
        [
            Counter(Pdata[ID]["atomdict"])
            for ID in Pdata
            for _ in range(Pdata[ID]["count"])
        ],
        start=Counter(),
    )  # Sum of atom counts/types on RHS
    Pcharge = sum(
        [Pdata[ID]["charge"] for ID in Pdata for _ in range(Pdata[ID]["count"])]
    )

    #%% If reaction is balanced already
    if Rcount == Pcount and Rcharge == Pcharge:
        print("Reaction is fully balanced")
        if first:
            msg = "Already balanced"
        elif msg:
            msg = msg + ", " + "Balanced"
        else:
            msg = "Balanced"
        if addedstr:
            msg += addedstr
        return update_rxn(Rdata, Pdata, hcrct=addedhc, rxnsmiles0=rxnsmiles0, msg=msg)

    #%% Otherwise take difference between atom type/count RHS and LHS

    rem = Counter()  # Rem contains difference between product and reactant counters
    rem.update(
        Pcount
    )  # If atoms not balanced and same charge, attempt to balance. Can try if different charge but more tricky
    rem.subtract(
        Rcount
    )  # Subtracting reactant atom type index from product atom type index
    postype = {
        key: rem[key] for key in rem.keys() if rem[key] > 0
    }  # Finding only positive keys. Note that if counter is positive this means extra molecules need to be added on LHS (eg. reagent).
    negtype = {
        key: abs(rem[key]) for key in rem.keys() if rem[key] < 0
    }  # Finding only negative keys. If counter is negative this means extra molecules need to be added on RHS (eg. help compounds)

    #     breakpoint()

    if postype:  # Reactants, Reagents may be needed
        status = [
            mandrcts[ID0]["count"] > Rdata[ID0]["count"]
            for ID0 in mandrcts
            if ID0 in Rdata
        ]  # if ID0 in Rdata
        if any(status):
            if msg:
                msg = msg + ", " + "Mapping error" + addedstr
            else:
                msg = "Mapping error" + addedstr
            return update_rxn(
                Rdata,
                Pdata,
                hc_prod=hc_prod,
                hcrct=addedhc,
                rxnsmiles0=rxnsmiles0,
                msg=msg,
            )
        #%% Initializing variables
        candirxt = []
        candirgt = []
        candihc = []
        matches = []
        addedspecies_ = []
        addedhc_ = []
        #%% Get reactant match first
        candirxt = [
            rctid
            for rctid in Rdata
            if set(postype.keys()).issubset(set(Rdata[rctid]["atomdict"].keys()))
        ]
        #%% Get reagent match
        candirgt = [
            rgtid
            for rgtid in Rgtdata
            if set(postype.keys()).issubset(set(Rgtdata[rgtid]["atomdict"].keys()))
            and rgtid not in candirxt
        ]
        #%% Get help compound match
        if hc_react:
            candihc = [
                hcid
                for hcid in hc_react
                if set(postype.keys()).issubset(set(hc_react[hcid]["atomdict"].keys()))
            ]
        #         breakpoint()
        if candirxt and not candirgt:  # Only reactant matches
            Rdata, addedspecies_, msg_ = resolvecandidates(
                postype,
                Rdata,
                Rdata,
                candirxt,
                Pdata,
                validate=False,
                rctonly=addrctonly,
                ignoreH=ignoreH,
            )
        elif candirgt and not candirxt:
            Rdata, addedspecies_, msg_ = resolvecandidates(
                postype,
                Rdata,
                Rgtdata,
                candirgt,
                Pdata,
                validate=False,
                rctonly=False,
                ignoreH=ignoreH,
            )
        elif candirxt and candirgt:
            Rdata, addedspecies_, msg_ = resolvecandidates(
                postype,
                Rdata,
                copy.deepcopy({**Rdata, **Rgtdata}),
                candirxt + candirgt,
                Pdata,
                validate=False,
                rctonly=addrctonly,
                ignoreH=ignoreH,
            )
        elif not candirxt and not candirgt:
            combineddict = copy.deepcopy({**Rdata, **Rgtdata})
            candispec = [
                specid
                for specid in combineddict
                if set(combineddict[specid]["atomdict"].keys()).intersection(
                    set(postype.keys())
                )
            ]
            if not candispec:
                msg_ = "LHS species insufficient"
            else:
                Rdata, addedspecies_, msg_ = resolvecandidates(
                    postype,
                    Rdata,
                    combineddict,
                    candispec,
                    Pdata,
                    rctonly=addrctonly,
                    ignoreH=ignoreH,
                )
        elif candihc and not candirxt and not candirgt:
            Rdata, addedhc_, msg_ = resolvecandidates(
                postype,
                Rdata,
                hc_react,
                candihc,
                Pdata,
                validate=False,
                rctonly=addrctonly,
                ignoreH=ignoreH,
            )
        if (
            "Hydrogen carriers" in msg_
        ):  # Multiple candidates for hydrogen carriers (mapper won't help)
            if msg:
                msg = msg + ", " + msg_
            else:
                msg = msg_

        elif msg_ != "Valid":  # Atom surplus on RHS cannot be met by any LHS species
            if msg:
                msg = msg + ", " + msg_ + addedstr
            else:
                msg = msg_ + addedstr
            return update_rxn(
                Rdata,
                Pdata,
                hc_prod=hc_prod,
                hcrct=addedhc,
                rxnsmiles0=rxnsmiles0,
                msg=msg,
            )

        else:
            addedspecies += addedspecies_
            if addedhc_:
                adddedhc += addedhc_
            # New#
            if (
                len(postype) == 1
                and "H" in postype
                and "With hydrogen carriers" not in msg
            ):
                if msg:
                    msg = (
                        msg
                        + ", "
                        + "With hydrogen carriers: "
                        + ",".join([str(addedspec) for addedspec in addedspecies_])
                    )
                else:
                    msg = "With hydrogen carriers: " + ",".join(
                        [str(addedspec) for addedspec in addedspecies_]
                    )
            # New#
        #         breakpoint()
        return balancerxn(
            Rdata,
            Pdata,
            Rgtdata=Rgtdata,
            rxnsmiles0=rxnsmiles0,
            first=False,
            usemapper=usemapper,
            addedspecies=addedspecies,
            addedhc=addedhc,
            hc_prod=hc_prod,
            hc_react=hc_react,
            coefflim=coefflim,
            msg=msg,
            mandrcts=mandrcts,
            addrctonly=addrctonly,
            ignoreH=ignoreH,
        )

    elif negtype:
        #         breakpoint()
        if (
            usemapper and len(set(addedspecies)) > 1
        ):  # More than one choice or added small species, let the mapper decide
            rxnsmiles = buildrxn(Rdata, Pdata)
            mapped_rxn = maprxn([rxnsmiles])[0]
            if mapped_rxn == "Error":
                if not addrctonly:
                    addrctonly = True
                    candidates2 = [candi for candi in Rdata if candi in mandrcts]
                    if candidates2 and candidates2 != list(Rdata.keys()):
                        Rdata = {ID0: Rdata[ID0] for ID0 in candidates2}
                        addedspecies = [
                            addedspec
                            for addedspec in addedspecies
                            if addedspec in Rdata
                        ]
                        if addedhc:
                            addedhc = [addedh for addedh in addedhc if addedh in Rdata]
                        return balancerxn(
                            Rdata,
                            Pdata,
                            rxnsmiles0=rxnsmiles0,
                            first=False,
                            usemapper=usemapper,
                            addedspecies=addedspecies,
                            addedhc=addedhc,
                            hc_prod=hc_prod,
                            hc_react=hc_react,
                            coefflim=coefflim,
                            msg=msg,
                            mandrcts=mandrcts,
                            addrctonly=addrctonly,
                            ignoreH=ignoreH,
                        )

                if all([Rdata[ID0]["count"] == 1 for ID0 in Rdata]) or any(
                    [Rdata[ID0]["count"] >= 10 for ID0 in Rdata]
                ):
                    if msg:
                        msg = msg + ", " + "Mapping error" + addedstr
                    else:
                        msg = "Mapping error" + addedstr
                    return update_rxn(
                        Rdata,
                        Pdata,
                        hc_prod=hc_prod,
                        hcrct=addedhc,
                        rxnsmiles0=rxnsmiles0,
                        msg=msg,
                    )
                else:
                    status = [
                        mandrcts[ID0]["count"] >= Rdata[ID0]["count"]
                        for ID0 in mandrcts
                        if ID0 in Rdata
                    ]  # if ID0 in Rdata
                    if any(status):
                        for ID0 in Rdata:
                            Rdata[ID0]["count"] = Rdata[ID0]["count"] - 1
                    else:
                        mandrcts = copy.deepcopy(Rdata)
                        mincount = min([Rdata[ID0]["count"] for ID0 in Rdata])
                        for ID0 in Rdata:
                            Rdata[ID0]["count"] = mincount
                    return balancerxn(
                        Rdata,
                        Pdata,
                        rxnsmiles0=rxnsmiles0,
                        first=False,
                        usemapper=usemapper,
                        addedspecies=addedspecies,
                        addedhc=addedhc,
                        hc_prod=hc_prod,
                        hc_react=hc_react,
                        coefflim=coefflim,
                        msg=msg,
                        mandrcts=mandrcts,
                        addrctonly=addrctonly,
                        ignoreH=ignoreH,
                    )
            else:
                mappedrxn = mapped_rxn.get("mapped_rxn")
                if "With hydrogen carriers" in msg:
                    hcarriers = [
                        int(hcarrier)
                        for hcarrier in msg.split("With hydrogen carriers: ")[1]
                        .split(", ")[0]
                        .split(",")
                    ]
                else:
                    hcarriers = []
                #                 breakpoint()
                LHSdata, _, msg_ = checkrxn(
                    mappedrxn,
                    Rdata=Rdata,
                    updateall=False,
                    mandrcts=list(mandrcts.keys()),
                    hcarriers=hcarriers,
                )
                if (
                    "Mandatory" in msg_ and not addrctonly
                ):  # Reactants are unmapped (Try again with a smaller candidate selection)
                    return balancerxn(
                        mandrcts,
                        Pdata,
                        Rgtdata=Rgtdata,
                        rxnsmiles0=rxnsmiles0,
                        first=True,
                        usemapper=usemapper,
                        addedspecies=addedspecies,
                        addedhc=addedhc,
                        hc_prod=hc_prod,
                        hc_react=hc_react,
                        coefflim=coefflim,
                        msg=msg,
                        addrctonly=True,
                        ignoreH=ignoreH,
                    )
                if msg:
                    if msg_ != "Valid":
                        msg = msg_ + ", " + msg
                    msg = "Mapper used" + ", " + msg
                else:
                    if msg_ != "Valid":
                        msg = "Mapper used" + ", " + msg_
                    else:
                        msg = "Mapper used"
            if addedhc:
                addedhc = [addedh for addedh in addedhc if addedh in LHSdata]
            #             breakpoint()
            addedspecies = [
                addedspec for addedspec in addedspecies if addedspec in LHSdata
            ]
            return balancerxn(
                LHSdata,
                Pdata,
                Rgtdata=Rgtdata,
                rxnsmiles0=rxnsmiles0,
                first=False,
                usemapper=False,
                addedspecies=addedspecies,
                addedhc=addedhc,
                hc_prod=hc_prod,
                hc_react=hc_react,
                coefflim=coefflim,
                msg=msg,
                mandrcts=mandrcts,
                addrctonly=addrctonly,
                ignoreH=ignoreH,
            )
        else:
            #             breakpoint()
            # Find match with Rdata first
            candidates = [
                spec
                for spec in Rdata
                if Rdata[spec]["atomdict"] == negtype
                if spec in addedspecies
            ]
            if not candidates:
                candidates = [
                    spec
                    for spec in Rdata
                    if set(Rdata[spec]["atomdict"].keys()) == set(negtype.keys())
                    if spec in addedspecies
                ]
            if candidates:
                specremoved = False
                for candi in candidates:
                    matcharray = list(
                        set(
                            Counter(
                                {
                                    k: negtype[k] / Rdata[candi]["atomdict"][k]
                                    for k in negtype
                                }
                            ).values()
                        )
                    )
                    if all([mult.is_integer() for mult in matcharray]):
                        mult = min(matcharray)
                        newcount = int(Rdata[candi]["count"] - mult)
                        if newcount < 0:
                            newcount = 0
                        if newcount == 0:
                            del Rdata[candi]
                        else:
                            Rdata[candi]["count"] = newcount
                        specremoved = True
                        break
                if specremoved:
                    if addedhc:
                        addedhc = [addedh for addedh in addedhc if addedh in Rdata]
                    #                         refhc=Counter(addedhc)
                    #                         addedhc=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refhc[addedspec]])) if addedspec in addedhc]
                    addedspecies = [
                        addedspec for addedspec in addedspecies if addedspec in Rdata
                    ]
                    #                     refspec=Counter(addedspecies)
                    #                     addedspecies=[addedspec for addedspec in Rdata for _ in range(min([Rdata[addedspec]['count'],refspec[addedspec]])) if addedspec in addedspecies if addedspec not in addedhc]
                    if "Mixture" in msg:
                        msglist = msg.split(",")
                        msglist2 = copy.copy(msglist)
                        for i, msg_ in enumerate(msglist):
                            if "Mixture" in msg_:
                                mixmsg = msg_.split(":")
                                refmixtures = mixmsg[1].split(",")
                                rmixtures = {
                                    rmixture
                                    for rmixture in refmixtures
                                    if rmixture in Rdata
                                }
                                if rmixtures:
                                    msglist2[i] = mixmsg[0] + ",".join(rmixtures)
                                else:
                                    msglist2.remove(msglist2[i])
                                msg = ", ".join(msglist2)
                    return balancerxn(
                        Rdata,
                        Pdata,
                        Rgtdata=Rgtdata,
                        rxnsmiles0=rxnsmiles0,
                        first=False,
                        usemapper=False,
                        addedspecies=addedspecies,
                        addedhc=addedhc,
                        hc_prod=hc_prod,
                        hc_react=hc_react,
                        coefflim=coefflim,
                        msg=msg,
                        addrctonly=addrctonly,
                        ignoreH=ignoreH,
                    )

            # Then try help compounds
            if Rcharge == Pcharge:
                hc_prod = {
                    hcid: hc_prod[hcid]
                    for hcid in hc_prod
                    if hc_prod[hcid]["charge"] == 0
                }  # Atom type for help compounds

            hc_prod2 = {
                hcid: hc_prod[hcid]
                for hcid in hc_prod
                if hc_prod[hcid]["atomdict"] == negtype
            }  # Exact match for 1 compound
            if not hc_prod2:
                hc_prod2 = {
                    hcid: hc_prod[hcid]
                    for hcid in hc_prod
                    if hc_prod[hcid]["atomdict"].keys() == negtype.keys()
                }  # Key match for 1 compound
            if hc_prod2:
                balbefore = False
            else:
                balbefore = True
            reac = {}
            prod = {}
            hcid = []
            try:
                reac, prod, hcid, msg0 = balance(
                    Rdata,
                    Pdata,
                    hc_prod=hc_prod2,
                    balbefore=balbefore,
                    coefflim=coefflim,
                    addedspecies=[
                        addedspec
                        for addedspec in addedspecies
                        if addedspec not in mandrcts
                    ],
                    addedhc=[addedh for addedh in addedhc if addedh not in mandrcts],
                    hc_react=hc_react,
                )
                if msg:
                    msg = msg + ", " + msg0
                else:
                    msg = msg0
            except Exception:
                if not balbefore:
                    return balancerxn(
                        Rdata,
                        Pdata,
                        Rgtdata=Rgtdata,
                        rxnsmiles0=rxnsmiles0,
                        first=False,
                        usemapper=False,
                        addedspecies=addedspecies,
                        addedhc=addedhc,
                        hc_prod={},
                        hc_react=hc_react,
                        coefflim=coefflim,
                        msg=msg,
                        addrctonly=addrctonly,
                        ignoreH=ignoreH,
                    )
                if Rcharge != Pcharge:
                    if msg:
                        msg = msg + ", " + "Charge imbalance"
                    else:
                        msg = "Charge imbalance"
                if msg:
                    msg = msg + ", " + "RHS species insufficient" + addedstr
                else:
                    msg = "RHS species insufficient" + addedstr

                return update_rxn(
                    Rdata, Pdata, hcrct=addedhc, rxnsmiles0=rxnsmiles0, msg=msg
                )
            else:
                if Rcharge != Pcharge:
                    if msg:
                        msg = msg + ", " + "Charge imbalance" + addedstr
                    else:
                        msg = "Charge imbalance" + addedstr
                return update_rxn(
                    Rdata,
                    Pdata,
                    reac=reac,
                    prod=prod,
                    hcprod=hcid,
                    hc_prod=hc_prod2,
                    hcrct=addedhc,
                    rxnsmiles0=rxnsmiles0,
                    msg=msg,
                )
    else:
        if Rcharge != Pcharge:
            if msg:
                msg = msg + ", " + "Charge imbalance" + addedstr
            else:
                msg = "Charge imbalance" + addedstr
        return update_rxn(
            Rdata, Pdata, hc_prod=hc_prod, hcrct=addedhc, rxnsmiles0=rxnsmiles0, msg=msg
        )
        
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

# balance_rxn_dp('CCCCCCc1ccc(C(=O)Cl)cc1.CCCCc1ccc(C#Cc2ccc(CNc3ccc4c(c3)C(=O)OC(C)(C)O4)cc2)cc1.Cl>>CCCCCCc1ccc(C(=O)N(Cc2ccc(C#Cc3ccc(CCCC)cc3)cc2)c2ccc3c(c2)C(=O)OC(C)(C)O3)cc1')


