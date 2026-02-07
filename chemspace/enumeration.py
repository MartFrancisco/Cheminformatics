"""
Combinatorial library enumeration functions.

Tools for generating virtual libraries through chemical reactions,
particularly focused on amide coupling.
"""

import itertools
from rdkit import Chem
from rdkit.Chem import AllChem


def enumerate_library(rxn_mol, reagent_lol):
    """
    Enumerate all products from a combinatorial reaction.
    
    Takes a reaction SMARTS and lists of reagents, then generates all
    possible products from the combinatorial space.
    
    Parameters
    ----------
    rxn_mol : rdkit.Chem.rdChemReactions.ChemicalReaction
        RDKit reaction object (created from SMARTS)
    reagent_lol : list of lists
        List of lists, e.g., [amine_list, acid_list]
        Each inner list contains tuples of (mol, name/ID)
    
    Returns
    -------
    list
        List of [SMILES, name] for each product
    
    Notes
    -----
    Failed reactions (e.g., incompatible functional groups) are skipped
    and logged to console.
    
    Examples
    --------
    >>> # Define amide coupling reaction
    >>> rxn_smarts = "[#7H1,#7H2:7].[#6:1](=[#8:3])[#8H1]>>[#7:7]-[#6:1](=[#8:3])"
    >>> rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)
    >>> 
    >>> # Prepare reagents (mol objects with IDs)
    >>> amines = [(mol1, 'amine_A'), (mol2, 'amine_B')]
    >>> acids = [(mol3, 'acid_1'), (mol4, 'acid_2')]
    >>> 
    >>> # Enumerate library
    >>> products = enumerate_library(rxn_mol, [amines, acids])
    >>> print(f"Generated {len(products)} amides")
    Generated 4 amides
    >>> 
    >>> # Convert to DataFrame
    >>> import pandas as pd
    >>> df_products = pd.DataFrame(products, columns=['SMILES', 'Product_ID'])
    """
    prod_list = []
    
    # Use itertools.product explicitly to avoid conflicts
    for reagents in itertools.product(*reagent_lol):
        mol_list = [x[0] for x in reagents]  # Extract molecules
        name_list = [str(x[1]) for x in reagents]  # Extract names
        name = "_".join(name_list)
        
        try:
            # Run reaction
            prod = rxn_mol.RunReactants(mol_list)
            
            if prod and len(prod) > 0:
                product_mol = prod[0][0]
                Chem.SanitizeMol(product_mol)
                prod_list.append([Chem.MolToSmiles(product_mol), name])
        except Exception as e:
            print(f"Error with {name}: {e}")
            continue
    
    return prod_list