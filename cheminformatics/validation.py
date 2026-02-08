"""
Chemical space validation and diversity analysis functions.

Provides tools for comparing property distributions, analyzing scaffold diversity,
and computing pairwise similarity metrics.
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import stats
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


def compare_properties_distributions(original_df, selected_df, property_name='MW', smiles_col='SMILES'):
    """
    Compare distributions of molecular properties between original and selected datasets.
    
    Uses Kolmogorov-Smirnov test to assess statistical similarity and plots
    overlapping histograms for visual comparison.
    
    Parameters
    ----------
    original_df : pd.DataFrame
        Original dataset before selection
    selected_df : pd.DataFrame
        Dataset containing the selected compounds
    property_name : str, default='MW'
        Name of the molecular property to compare.
        Options: 'MW', 'LogP', 'TPSA', 'HBD', 'HBA', 'NumRings', 'NumRotBonds'
    smiles_col : str, default='SMILES'
        Name of the column containing SMILES
    
    Returns
    -------
    None
        Displays plot and prints statistics
    
    Examples
    --------
    >>> compare_properties_distributions(df_all, df_diverse, property_name='MW')
    Calculating MW...
    
    MW Coverage:
      Original range: [120.15, 498.67]
      Selected range: [125.23, 495.34]
      Coverage: 94.5%
    """
    prop_functions = {
        'MW': Descriptors.MolWt,
        'LogP': Descriptors.MolLogP,
        'TPSA': Descriptors.TPSA,
        'HBD': Descriptors.NumHDonors,
        'HBA': Descriptors.NumHAcceptors,
        'NumRings': Descriptors.RingCount,
        'NumRotBonds': Descriptors.NumRotatableBonds,
    }
    
    if property_name not in prop_functions:
        raise ValueError(f"Unknown property: {property_name}. Choose from {list(prop_functions.keys())}")
    
    prop_func = prop_functions[property_name]
    
    # Convert SMILES to mol if needed
    if 'mol' not in original_df.columns:
        original_df['mol'] = original_df[smiles_col].apply(Chem.MolFromSmiles)
    
    if 'mol' not in selected_df.columns:
        selected_df['mol'] = selected_df[smiles_col].apply(Chem.MolFromSmiles)
    
    # Calculate property for both sets
    print(f"Calculating {property_name}...")
    original_vals = original_df['mol'].apply(prop_func)
    selected_vals = selected_df['mol'].apply(prop_func)
    
    # Plot overlapping histograms
    plt.figure(figsize=(10, 6))
    plt.hist(original_vals, bins=50, alpha=0.5, label=f'Original (n={len(original_df)})', 
             density=True, edgecolor='black')
    plt.hist(selected_vals, bins=30, alpha=0.7, label=f'Selected (n={len(selected_df)})', 
             density=True, edgecolor='black', color='red')
    plt.xlabel(property_name)
    plt.ylabel('Density')
    plt.title(f'{property_name} Distribution: Original vs Selected')
    plt.legend()
    plt.grid(alpha=0.3)
    
    # Statistical test
    ks_stat, p_value = stats.ks_2samp(original_vals, selected_vals)
    plt.text(0.02, 0.98, f'KS test: p={p_value:.4f}', 
             transform=plt.gca().transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    plt.show()
    
    print(f"\n{property_name} Coverage:")
    print(f"  Original range: [{original_vals.min():.2f}, {original_vals.max():.2f}]")
    print(f"  Selected range: [{selected_vals.min():.2f}, {selected_vals.max():.2f}]")
    print(f"  Coverage: {(selected_vals.max()-selected_vals.min())/(original_vals.max()-original_vals.min())*100:.1f}%")


def scaffold_diversity_plot(original_df, selected_df, smiles_col='SMILES'):
    """
    Compare scaffold diversity between original and selected datasets.
    
    Uses Murcko scaffolds to assess structural diversity. Displays side-by-side
    histograms showing the distribution of molecules per scaffold.
    
    Parameters
    ----------
    original_df : pd.DataFrame
        Original dataset before selection
    selected_df : pd.DataFrame
        Dataset containing the selected compounds
    smiles_col : str, default='SMILES'
        Name of the column containing SMILES
    
    Returns
    -------
    None
        Displays plot and prints scaffold diversity metrics
    
    Notes
    -----
    Higher scaffold diversity (ratio of unique scaffolds to total molecules)
    indicates better structural coverage of chemical space.
    
    Examples
    --------
    >>> scaffold_diversity_plot(df_all, df_diverse)
    Calculating scaffolds...
    
    Scaffold Diversity Metrics:
      Original: 2500 scaffolds / 10000 molecules = 25.00%
      Selected: 85 scaffolds / 100 molecules = 85.00%
    """
    # Convert SMILES to mol if needed
    if 'mol' not in original_df.columns:
        original_df['mol'] = original_df[smiles_col].apply(Chem.MolFromSmiles)
    
    if 'mol' not in selected_df.columns:
        selected_df['mol'] = selected_df[smiles_col].apply(Chem.MolFromSmiles)
    
    def get_scaffold_counts(df):
        print("Calculating scaffolds...")
        scaffolds = df['mol'].apply(
            lambda m: MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None
        )
        return scaffolds.value_counts()
    
    orig_scaffolds = get_scaffold_counts(original_df)
    sel_scaffolds = get_scaffold_counts(selected_df)
    
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    # Original dataset scaffold frequency
    axes[0].hist(orig_scaffolds.values, bins=50, edgecolor='black', alpha=0.7)
    axes[0].set_xlabel('Molecules per Scaffold')
    axes[0].set_ylabel('Number of Scaffolds')
    axes[0].set_title(f'Original Dataset\n{len(orig_scaffolds)} unique scaffolds from {len(original_df)} molecules')
    axes[0].set_yscale('log')
    
    # Selected dataset scaffold frequency
    axes[1].hist(sel_scaffolds.values, bins=20, edgecolor='black', alpha=0.7, color='red')
    axes[1].set_xlabel('Molecules per Scaffold')
    axes[1].set_ylabel('Number of Scaffolds')
    axes[1].set_title(f'Selected Dataset\n{len(sel_scaffolds)} unique scaffolds from {len(selected_df)} molecules')
    
    plt.tight_layout()
    plt.show()
    
    print(f"\nScaffold Diversity Metrics:")
    print(f"  Original: {len(orig_scaffolds)} scaffolds / {len(original_df)} molecules = {len(orig_scaffolds)/len(original_df):.2%}")
    print(f"  Selected: {len(sel_scaffolds)} scaffolds / {len(selected_df)} molecules = {len(sel_scaffolds)/len(selected_df):.2%}")


def pairwise_similarity_histogram(selected_df, smiles_col='SMILES', sample_size=100, fp_col='FP'):
    """
    Plot distribution of pairwise Tanimoto similarities.
    
    Directly shows how diverse the selection is by computing all pairwise
    similarities. Lower average similarity indicates better diversity.
    
    Parameters
    ----------
    selected_df : pd.DataFrame
        Dataset containing the selected compounds
    smiles_col : str, default='SMILES'
        Name of the column containing SMILES
    sample_size : int, default=100
        Maximum number of molecules to use (for computational efficiency)
    fp_col : str, default='FP'
        Name of the fingerprint column (will generate if not present)
    
    Returns
    -------
    None
        Displays histogram and prints similarity statistics
    
    Notes
    -----
    Interpretation guide:
    - Mean < 0.4: Excellent diversity (low similarity)
    - Mean 0.4-0.6: Good diversity (moderate similarity)
    - Mean > 0.6: Poor diversity (high similarity)
    
    Examples
    --------
    >>> pairwise_similarity_histogram(df_diverse, sample_size=100)
    Calculating 4950 pairwise similarities...
    
    Pairwise Similarity Statistics:
      Mean: 0.325
      Median: 0.310
    
    Interpretation:
      ✓ Excellent diversity! Low average similarity.
    """
    # Generate fingerprints if needed
    if fp_col not in selected_df.columns:
        selected_df['mol'] = selected_df[smiles_col].apply(Chem.MolFromSmiles)
        morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        selected_df[fp_col] = selected_df['mol'].apply(lambda m: morgan_gen.GetFingerprint(m) if m is not None else None)
        
        selected_df = selected_df[selected_df[fp_col].notna()].copy()
    
    # Limit sample size for computational efficiency
    if len(selected_df) > sample_size:
        print(f"Sampling {sample_size} molecules from {len(selected_df)} for pairwise comparison...")
        df_sample = selected_df.sample(sample_size, random_state=42)
    else:
        df_sample = selected_df
    
    fps = df_sample[fp_col].tolist()
    n = len(fps)
    
    # Calculate all pairwise similarities
    print(f"Calculating {n*(n-1)//2} pairwise similarities...")
    similarities = []
    for i in range(n):
        for j in range(i+1, n):
            sim = DataStructs.TanimotoSimilarity(fps[i], fps[j])
            similarities.append(sim)
    
    plt.figure(figsize=(10, 6))
    plt.hist(similarities, bins=50, edgecolor='black', alpha=0.7)
    plt.xlabel('Tanimoto Similarity')
    plt.ylabel('Count')
    plt.title(f'Pairwise Similarity Distribution ({len(similarities)} pairs from {n} molecules)')
    plt.axvline(np.mean(similarities), color='red', linestyle='--', 
                linewidth=2, label=f'Mean: {np.mean(similarities):.3f}')
    plt.axvline(np.median(similarities), color='green', linestyle='--', 
                linewidth=2, label=f'Median: {np.median(similarities):.3f}')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.show()
    
    print(f"\nPairwise Similarity Statistics:")
    print(f"  Mean: {np.mean(similarities):.3f}")
    print(f"  Median: {np.median(similarities):.3f}")
    print(f"  Min: {np.min(similarities):.3f}")
    print(f"  Max: {np.max(similarities):.3f}")
    print(f"  Std: {np.std(similarities):.3f}")
    
    # Interpretation guide
    print(f"\nInterpretation:")
    if np.mean(similarities) < 0.4:
        print("  ✓ Excellent diversity! Low average similarity.")
    elif np.mean(similarities) < 0.6:
        print("  ○ Good diversity. Moderate similarity.")
    else:
        print("  ✗ Poor diversity. High average similarity suggests molecules are very similar.")