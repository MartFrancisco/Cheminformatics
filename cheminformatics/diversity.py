"""
Diversity selection functions.

Implements MaxMin algorithm for selecting diverse molecular subsets
using Tanimoto similarity.
"""

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.Chem import AllChem


def maxmin_from_fingerprints(df, n_select=100, fp_prefix='bit_', seed=None):
    """
    Optimized MaxMin selection using explicit fingerprint columns.
    
    Uses vectorized Tanimoto calculations for better performance.
    Best for n_select < 1000.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing fingerprints where each bit is a column
    n_select : int, default=100
        Number of molecules to select
    fp_prefix : str, default='bit_'
        Prefix of fingerprint columns
    seed : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    pd.DataFrame
        DataFrame with n_select diverse compounds, including 'Selection_Order' column
    
    Examples
    --------
    >>> df_fps = compute_morgan_fps_explicit(df, fpSize=1024)
    >>> diverse = maxmin_from_fingerprints(df_fps, n_select=100, seed=42)
    >>> print(diverse['Selection_Order'].tolist()[:5])
    [0, 1, 2, 3, 4]
    """
    if n_select > len(df):
        raise ValueError(f"Cannot select {n_select} molecules from {len(df)}")
    
    print(f"MaxMin selection: choosing {n_select} from {len(df)} molecules...")
    
    # Extract fingerprints as numpy array
    fp_cols = [col for col in df.columns if col.startswith(fp_prefix)]
    if not fp_cols:
        raise ValueError(f"No fingerprint columns found")
    
    fps = df[fp_cols].to_numpy(dtype=np.uint8)
    n_mols = len(fps)
    
    print(f"Using {len(fp_cols)} fingerprint bits")
    
    # Precompute fingerprint sums (for Tanimoto denominator)
    fp_sums = fps.sum(axis=1)
    
    # Vectorized Tanimoto distance function
    def tanimoto_distances_vectorized(fp_idx, selected_indices):
        """Compute Tanimoto distances from fp_idx to all selected molecules"""
        fp = fps[fp_idx]
        selected_fps = fps[selected_indices]
        
        # Intersection
        intersections = np.dot(selected_fps, fp)
        
        # Union
        unions = fp_sums[selected_indices] + fp_sums[fp_idx] - intersections
        
        # Similarity
        with np.errstate(divide='ignore', invalid='ignore'):
            similarities = intersections / unions
            similarities[unions == 0] = 0
        
        # Distance (1 - similarity)
        distances = 1 - similarities
        
        return distances.min()  # Return minimum distance
    
    # Initialize
    if seed is not None:
        np.random.seed(seed)
        first_idx = np.random.randint(0, n_mols)
    else:
        first_idx = 0
    
    selected_idx = [first_idx]
    unselected = set(range(n_mols)) - {first_idx}
    
    print("Selecting diverse molecules...")
    
    # MaxMin selection
    while len(selected_idx) < n_select:
        max_min_dist = -1
        next_idx = -1
        
        # Check each unselected molecule
        for i in unselected:
            min_dist = tanimoto_distances_vectorized(i, selected_idx)
            
            if min_dist > max_min_dist:
                max_min_dist = min_dist
                next_idx = i
        
        selected_idx.append(next_idx)
        unselected.remove(next_idx)
        
        if len(selected_idx) % 10 == 0 or len(selected_idx) == n_select:
            print(f"  Selected {len(selected_idx)}/{n_select} molecules... (diversity: {max_min_dist:.3f})")
    
    # Create result
    diverse_df = df.iloc[selected_idx].copy()
    diverse_df['Selection_Order'] = range(len(diverse_df))
    
    print(f"\nMaxMin selection complete!")
    
    return diverse_df


def maxmin_from_rdkit_fps(df, n_select=100, fp_col='FP', seed=None):
    """
    MaxMin selection using RDKit fingerprint objects.
    
    Suitable for n_select ~100–500. 
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing RDKit fingerprints as objects
    n_select : int, default=100
        Number of molecules to select
    fp_col : str, default='FP'
        Name of the column containing the fingerprint object
    seed : int, optional
        Random seed for reproducibility
    
    Returns
    -------
    pd.DataFrame
        DataFrame with n_select diverse compounds, including 'Selection_Order' column
    
    Examples
    --------
    >>> df_fps = compute_morgan_fps(df, fp_column='FP')
    >>> diverse = maxmin_from_rdkit_fps(df_fps, n_select=100, seed=42)
    MaxMin selection: choosing 100 from 5000 molecules
    Selecting diverse molecules...
      Selected 10/100 (diversity: 0.524)
      ...
    MaxMin selection complete!
    """
    from rdkit.Chem import AllChem
    from rdkit import DataStructs
    
    if n_select > len(df):
        raise ValueError(f"Cannot select {n_select} molecules from {len(df)}")        

    fps = df[fp_col].tolist()
    n_mols = len(fps)
    
    print(f"MaxMin selection: choosing {n_select} from {n_mols} molecules")
    
    # Initialize
    rng = np.random.default_rng(seed)
    first_idx = rng.integers(n_mols)
    
    selected = [first_idx]
    unselected = set(range(n_mols)) - {first_idx}
    
    print("Selecting diverse molecules...")
    
    while len(selected) < n_select:
        best_idx = None
        best_dist = -1.0
        
        for i in unselected:
            sims = DataStructs.BulkTanimotoSimilarity(
                fps[i],
                [fps[j] for j in selected]
            )
            min_dist = 1.0 - max(sims)
            
            if min_dist > best_dist:
                best_dist = min_dist
                best_idx = i
        
        selected.append(best_idx)
        unselected.remove(best_idx)
        
        if len(selected) % 10 == 0 or len(selected) == n_select:
            print(f"  Selected {len(selected)}/{n_select} (diversity: {best_dist:.3f})")
    
    result = df.iloc[selected].copy()
    result['Selection_Order'] = range(len(result))
    
    print("MaxMin selection complete!")
    
    return result