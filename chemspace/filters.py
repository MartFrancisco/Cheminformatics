"""
Drug-likeness and quality filtering functions.

Implements filters for medicinal chemistry including Lipinski's Rule of Five,
PAINS detection, and structural quality checks.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams


def lipinski_ro5(df, smiles_column='SMILES', keep_descriptors=False):
    """
    Filter molecules based on Lipinski's Rule of Five.
    
    Lipinski's Rule of Five criteria:
    - Molecular Weight < 500 Da
    - LogP ≤ 5
    - Hydrogen Bond Donors (HBD) ≤ 5
    - Hydrogen Bond Acceptors (HBA) ≤ 10
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with SMILES column
    smiles_column : str, default='SMILES'
        Name of column containing SMILES
    keep_descriptors : bool, default=False
        If True, keep MW, LogP, HBD, HBA columns in output
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame passing Lipinski's Ro5
    
    Notes
    -----
    Molecules passing all four criteria are considered "drug-like" and
    more likely to have favorable pharmacokinetic properties.
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'CC(C)CC1=CC=C(C=C1)C(C)C(=O)O']})
    >>> filtered = lipinski_ro5(df)
    Initial molecule count: 2
    
    LIPINSKI RULE OF 5 FILTERING RESULTS
    ========================================
    Initial molecules:        2
    Passed Ro5:               2 (100.0%)
    """
    if smiles_column not in df.columns:
        raise ValueError(f"DataFrame must have a '{smiles_column}' column")
    
    print(f"Initial molecule count: {len(df)}")
    
    df_copy = df.copy()
    
    # Convert SMILES to molecules
    df_copy['Molecule'] = df_copy[smiles_column].apply(Chem.MolFromSmiles)
    
    # Calculate Lipinski descriptors
    print("Calculating molecular descriptors...")
    df_copy['MolWt'] = df_copy['Molecule'].apply(lambda m: Descriptors.MolWt(m) if m else None)
    df_copy['LogP'] = df_copy['Molecule'].apply(lambda m: Descriptors.MolLogP(m) if m else None)
    df_copy['HBD'] = df_copy['Molecule'].apply(lambda m: Lipinski.NumHDonors(m) if m else None)
    df_copy['HBA'] = df_copy['Molecule'].apply(lambda m: Lipinski.NumHAcceptors(m) if m else None)
    
    # Remove molecules that couldn't be processed
    df_copy = df_copy[df_copy['Molecule'].notna()].copy()
    
    # Apply Lipinski's Rule of Five
    print("Applying Lipinski's Rule of Five filters...")
    conditions = (
        (df_copy['MolWt'] < 500) &
        (df_copy['LogP'] <= 5) &
        (df_copy['HBD'] <= 5) &
        (df_copy['HBA'] <= 10)
    )
    
    filtered_df = df_copy[conditions].copy()
    
    # Count violations for reporting
    violations = {
        'MolWt > 500': (df_copy['MolWt'] >= 500).sum(),
        'LogP > 5': (df_copy['LogP'] > 5).sum(),
        'HBD > 5': (df_copy['HBD'] > 5).sum(),
        'HBA > 10': (df_copy['HBA'] > 10).sum()
    }
    
    # Clean up temporary columns unless requested
    if not keep_descriptors:
        filtered_df = filtered_df.drop(columns=['Molecule', 'MolWt', 'LogP', 'HBD', 'HBA'])
    else:
        filtered_df = filtered_df.drop(columns=['Molecule'])
    
    # Print summary
    print(f"\n{'='*60}")
    print("LIPINSKI RULE OF 5 FILTERING RESULTS")
    print(f"{'='*60}")
    print(f"Initial molecules:     {len(df):>8,}")
    print(f"Passed Ro5:            {len(filtered_df):>8,} ({100*len(filtered_df)/len(df):.1f}%)")
    print(f"Filtered out:          {len(df) - len(filtered_df):>8,}")
    print(f"\nViolations breakdown:")
    for rule, count in violations.items():
        print(f"  {rule:20s} {count:>6,} molecules")
    print(f"{'='*60}\n")
    
    return filtered_df


def pains_filter(df, column_name='SMILES'):
    """
    Filter molecules based on PAINS (Pan-Assay Interference Compounds) filters.
    
    PAINS are compounds that frequently show activity in many assays due to
    non-specific reactivity or assay interference, rather than genuine biological
    activity. Removing them improves library quality.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with SMILES column
    column_name : str, default='SMILES'
        Name of column containing SMILES
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame without PAINS
    
    Notes
    -----
    Uses RDKit's built-in PAINS filter catalog based on:
    Baell, J. B.; Holloway, G. A. J. Med. Chem. 2010, 53, 2719–2740.
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'c1ccc(cc1)N=Nc2ccccc2']})
    >>> filtered = pains_filter(df)
    Initial molecule count: 2
    After filtering by PAINS: 1
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must have a '{column_name}' column")
    
    print(f"Initial molecule count: {len(df)}")

    df_copy = df.copy()
    
    # Initialize PAINS filter
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    catalog = FilterCatalog(params)
    
    # Function to check if a molecule hits PAINS filters
    def is_pains(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return catalog.HasMatch(mol)
    
    # Identify molecules that are PAINS
    df_copy['IsPAINS'] = df_copy[column_name].apply(is_pains)
    
    # Filter out PAINS
    filtered_df = df_copy[~df_copy['IsPAINS']].copy()
    
    # Drop the temporary column
    filtered_df.drop(columns=['IsPAINS'], inplace=True)
    
    print(f"After filtering by PAINS: {len(filtered_df)}")
    print(f"Check PAINS Done!\n")
    
    return filtered_df


def filter_by_atoms(df, column_name='SMILES'):
    """
    Filter molecules based on allowed atoms (removes unconventional atoms).
    
    Keeps only molecules containing conventional medicinal chemistry atoms.
    Allowed atoms: H, C, N, O, F, P, S, Cl, Br, I
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with SMILES column
    column_name : str, default='SMILES'
        Name of column containing SMILES
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame with only conventional atoms
    
    Notes
    -----
    Removes molecules containing metals, silicon, boron, and other
    non-standard atoms that may complicate synthesis or have unknown toxicity.
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'C[Si](C)(C)C', 'CCCl']})
    >>> filtered = filter_by_atoms(df)
    Initial molecule count: 3
    After filtering non-conventional atoms: 2
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must have a '{column_name}' column")
    
    print(f"Initial molecule count: {len(df)}")

    df_copy = df.copy()
    
    # Define allowed atomic numbers
    # H, C, N, O, F, P, S, Cl, Br, I
    allowed_atoms = {1, 6, 7, 8, 9, 15, 16, 17, 35, 53}
    
    # Check if molecule contains only allowed atoms
    def is_molecule_allowed(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return all(atom.GetAtomicNum() in allowed_atoms for atom in mol.GetAtoms())
    
    # Filter molecules
    df_copy['Allowed'] = df_copy[column_name].apply(is_molecule_allowed)
    filtered_df = df_copy[df_copy['Allowed']].copy()
    filtered_df.drop(columns=['Allowed'], inplace=True)
    
    print(f"After filtering non-conventional atoms: {len(filtered_df)}")
    print(f"Remove Unconventional Atoms Done!\n")
    
    return filtered_df


def filter_by_FGs(df, column_name='SMILES'):
    """
    Filter molecules based on unwanted (potentially harmful) functional groups.
    
    Removes molecules containing reactive or toxic functional groups that
    may cause issues in biological assays or have poor safety profiles.
    
    Filtered groups:
    - Anilines (metabolic toxicity concerns)
    - Nitro groups (mutagenicity)
    - Aldehydes (reactivity)
    - Acid halides (high reactivity)
    - Sulfonyl halides (high reactivity)
    - Peroxides (instability)
    - Epoxides (reactivity, potential toxicity)
    - Aziridines (reactivity, potential toxicity)
    - Michael acceptors (non-selective reactivity)
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with SMILES column
    column_name : str, default='SMILES'
        Name of column containing SMILES
    
    Returns
    -------
    pd.DataFrame
        Filtered DataFrame without unwanted functional groups
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'c1ccc(cc1)N', 'CCN']})
    >>> filtered = filter_by_FGs(df)
    Initial molecule count: 3
    After FG filtering: 2 (66.7%)
    Functional group filter done!
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must have a '{column_name}' column")
    
    print(f"Initial molecule count: {len(df)}")
    df_copy = df.copy()
    
    functional_groups = {
        'aniline': '[c][NH2]',
        'nitro': '[N+](=O)[O-]',
        'aldehyde': '[CX3H1](=O)[#6]',
        'acid_halide': 'C(=O)[Cl,Br,I]',
        'sulfonyl_halide': 'S(=O)(=O)[Cl,Br,I]',
        'peroxide': 'OO',
        'epoxide': '[OX2]1CC1',
        'aziridine': '[NX3]1CC1',
        'michael_acceptor': '[C,c]=[C,c][C,S](=O)'
    }
    
    patterns = {k: Chem.MolFromSmarts(v) for k, v in functional_groups.items()}
    
    def has_bad_fg(smiles):
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return False
        return any(mol.HasSubstructMatch(p) for p in patterns.values())
    
    df_copy['HasBadFG'] = df_copy[column_name].apply(has_bad_fg)
    filtered_df = df_copy[~df_copy['HasBadFG']].copy()
    filtered_df.drop(columns=['HasBadFG'], inplace=True)
    
    print(f"After FG filtering: {len(filtered_df)} "
          f"({100*len(filtered_df)/len(df):.1f}%)")
    print("Functional group filter done!\n")
    
    return filtered_df