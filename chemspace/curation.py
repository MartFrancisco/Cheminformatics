"""
Molecular curation and classification functions.

This module provides tools for cleaning, standardizing, and classifying
molecular datasets for cheminformatics workflows.
"""

import pandas as pd
from rdkit import Chem
from rdkit.Chem.SaltRemover import SaltRemover
from rdkit.Chem.MolStandardize import rdMolStandardize


def curate_dataframe(df, smiles_col='SMILES', verbose=True):
    """
    Curate a molecular dataset from a pandas DataFrame by removing salts,
    standardizing tautomers, and deduplicating based on InChIKey.
    
    Parameters
    ----------
    df : pd.DataFrame
        Input DataFrame containing a column with SMILES strings
    smiles_col : str, default='SMILES'
        Name of the column containing SMILES
    verbose : bool, default=True
        Print curation statistics
    
    Returns
    -------
    pd.DataFrame
        Curated DataFrame with columns: 'SMILES' (curated), 'InChIKey'
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CC(=O)O', 'CCN', 'CC(=O)O']})
    >>> curated = curate_dataframe(df)
    Dataset Curation Report:
      Input molecules:      3
      Duplicates removed:   1
      Final dataset size:   2
    """
    # Initialize tools
    salt_remover = SaltRemover()
    tautomer_enumerator = rdMolStandardize.TautomerEnumerator()
    tautomer_enumerator.SetMaxTautomers(1)
    
    # Track statistics
    stats = {'input': len(df), 'invalid': 0, 'desalted': 0, 'duplicates': 0, 'inchikey_fail': 0}
    
    seen_inchikeys = set()
    curated_rows = []
    
    for idx, row in df.iterrows():
        smi = row[smiles_col]
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            stats['invalid'] += 1
            continue
        
        # Remove salts
        mol_desalted = salt_remover.StripMol(mol)
        if Chem.MolToSmiles(mol_desalted) != Chem.MolToSmiles(mol):
            stats['desalted'] += 1
        
        # Standardize tautomer
        mol_std = tautomer_enumerator.Canonicalize(mol_desalted)
        
        # Generate InChIKey
        try:
            inchikey = Chem.MolToInchiKey(mol_std)
        except:
            stats['inchikey_fail'] += 1
            continue

        # Skip empty or invalid InChIKeys
        if not inchikey or not inchikey.startswith("InChIKey=") and len(inchikey) != 27:
            stats['inchikey_fail'] += 1
            continue
        
        # Deduplicate
        if inchikey in seen_inchikeys:
            stats['duplicates'] += 1
            continue
        seen_inchikeys.add(inchikey)
        
        # Add curated row
        curated_row = row.copy()
        curated_row[smiles_col] = Chem.MolToSmiles(mol_std)
        curated_row['InChIKey'] = inchikey
        curated_rows.append(curated_row)
    
    curated_df = pd.DataFrame(curated_rows)
    
    if verbose:
        print(f"Dataset Curation Report:")
        print(f"  Input molecules:      {stats['input']}")
        print(f"  Invalid SMILES:       {stats['invalid']}")
        print(f"  Salts removed:        {stats['desalted']}")
        print(f"  Duplicates removed:   {stats['duplicates']}")
        print(f"  InChiKey Failed:      {stats['inchikey_fail']}")
        print(f"  Final dataset size:   {len(curated_df)}")
        print(f"  Retention rate:       {len(curated_df)/stats['input']*100:.1f}%")
    
    return curated_df


def FG_classifier(df, column_name='SMILES'):
    """
    Classify molecules based on functional groups (FG).
    
    Separates acids, amines, and flags problematic functional groups
    specifically for amide production.
    
    Parameters
    ----------
    df : pd.DataFrame
        The DataFrame containing SMILES strings
    column_name : str, default='SMILES'
        The name of the column in df that contains the SMILES strings
    
    Returns
    -------
    pd.DataFrame
        A DataFrame with 'Classes' column. Optimized for amide coupling
    
    Notes
    -----
    Classification priorities for amide coupling:
    1. 'acid' - IDEAL for coupling (acyl donors)
    2. 'amine' - IDEAL for coupling (nucleophiles)
    3. 'acid_amine' - Can self-couple
    4. 'acid_aldehyde', 'amine_aldehyde' - PROBLEMATIC
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CC(=O)O', 'CCN']})
    >>> classified = FG_classifier(df)
    >>> print(classified['Classes'].tolist())
    ['acid', 'amine']
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must have a '{column_name}' column")
    
    print(f"Initial molecule count: {len(df)}")
    
    # Define SMARTS patterns
    functional_groups = {
        'acid': '[CX3](=O)[OX2H1]',           # Carboxylic acid
        'amine': '[NX3;H2,H1;!$(NC=O)]',      # Primary/secondary amine
        'aldehyde': '[CX3H1](=O)[#6]',        # Aldehyde
        'aryl_halide': '[c][F,Cl,Br,I]',      # Aryl halide
        'ester': '[CX3](=O)[OX2H0]',          # Ester
        'amide': '[NX3][CX3](=O)',            # Amide
    }
    
    patterns = {key: Chem.MolFromSmarts(value) for key, value in functional_groups.items()}
    
    def classify_for_coupling(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return 'invalid'
            
            # Check for all functional groups
            matches = [fg_name for fg_name, pattern in patterns.items() 
                      if mol.HasSubstructMatch(pattern)]
            
            if not matches:
                return 'other'
            
            # Priority for amide coupling:
            # 1. If it's ONLY acid → 'acid' (good for coupling)
            # 2. If it's ONLY amine → 'amine' (good for coupling)
            # 3. If it has BOTH → 'acid_amine' (can couple with itself!)
            # 4. If it has problematic FGs → flag them
            
            has_acid = 'acid' in matches
            has_amine = 'amine' in matches
            has_aldehyde = 'aldehyde' in matches
            has_ester = 'ester' in matches
            has_amide = 'amide' in matches
            
            # Classify by coupling utility
            if has_acid and has_amine:
                return 'acid_amine'  # Self-coupling possible
            elif has_acid and has_aldehyde:
                return 'acid_aldehyde'  # PROBLEMATIC
            elif has_amine and has_aldehyde:
                return 'amine_aldehyde'  # PROBLEMATIC (imine formation)
            elif has_acid and not has_aldehyde:
                if has_ester:
                    return 'acid_ester'  # May need protection
                return 'acid'  # IDEAL for coupling
            elif has_amine and not has_aldehyde:
                return 'amine'  # IDEAL for coupling
            elif has_amide:
                return 'amide'  # Already coupled
            else:
                return ', '.join(matches)
                
        except Exception as e:
            print(f"Error processing molecule {smiles}: {e}")
            return 'error'
    
    df_classified = df.copy()
    df_classified['Classes'] = df_classified[column_name].apply(classify_for_coupling)
    
    print(f"Classification complete!")
    print(f"\nClass distribution:")
    print(df_classified['Classes'].value_counts())
    print(f"\nCoupling-ready molecules:")
    acids = len(df_classified[df_classified['Classes'] == 'acid'])
    amines = len(df_classified[df_classified['Classes'] == 'amine'])
    print(f"  Acids (acyl donors): {acids}")
    print(f"  Amines (nucleophiles): {amines}")
    print(f"  Potential coupling pairs: {acids * amines:,}")
    
    return df_classified


def filter_single_functional_group(df, fg_type, column_name='SMILES'):
    """
    Filter molecules to keep only those with exactly ONE functional group.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with molecules classified according to their functional groups
    fg_type : str
        Either 'amine' or 'acid'
    column_name : str, default='SMILES'
        The name of the column containing SMILES strings
    
    Returns
    -------
    pd.DataFrame
        DataFrame with molecules containing exactly 1 type of functional group
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CC(=O)O', 'O=C(O)CC(=O)O']})
    >>> filtered = filter_single_functional_group(df, 'acid')
    Initial acid count: 2
    Molecules with exactly 1 acid: 1
    Removed: 1
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must have a '{column_name}' column")
    
    print(f"Initial {fg_type} count: {len(df)}")
    
    # SMARTS patterns
    patterns = {
        'amine': Chem.MolFromSmarts('[NX3;H2,H1;!$(NC=O)]'),
        'acid': Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    }
    
    if fg_type not in patterns:
        raise ValueError(f"fg_type must be 'amine' or 'acid', got '{fg_type}'")
    
    pattern = patterns[fg_type]
    
    def count_fg(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return -1
            return len(mol.GetSubstructMatches(pattern))
        except Exception as e:
            print(f"Error processing {smiles}: {e}")
            return -1
    
    # Count and filter
    df = df.copy()
    df['FG_Count'] = df[column_name].apply(count_fg)
    filtered_df = df[df['FG_Count'] == 1].copy()
    filtered_df.drop(columns=['FG_Count'], inplace=True)
    
    print(f"Molecules with exactly 1 {fg_type}: {len(filtered_df)}")
    print(f"Removed: {len(df) - len(filtered_df)}\n")
    
    return filtered_df