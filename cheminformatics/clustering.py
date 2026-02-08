"""
Molecular clustering functions.

Implements Butina clustering algorithm for both explicit fingerprints
and RDKit fingerprint objects.
"""

import numpy as np
import pandas as pd
from rdkit import DataStructs
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.ML.Cluster import Butina


def butina_cluster_expandedFP(df, cutoff=0.35, fp_prefix="bit_", print_stats=True):
    """
    Cluster molecules using Butina algorithm with explicit fingerprint bits.
    
    Requires DataFrame with fingerprint columns where each bit is a separate column.
    Uses RDKit's lazy Tanimoto distance calculations.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with fingerprint columns (binary bits)
    cutoff : float, default=0.35
        Tanimoto distance cutoff (0–1). Lower values = tighter clusters
    fp_prefix : str, default="bit_"
        Prefix of fingerprint columns
    print_stats : bool, default=True
        Print clustering statistics
    
    Returns
    -------
    tuple
        (df_result, clusters) where:
        - df_result: DataFrame with ClusterID column added
        - clusters: list of clusters (each cluster is a list of indices)
    
    Examples
    --------
    >>> # Assuming df has columns: bit_0, bit_1, ..., bit_1023
    >>> df_clustered, clusters = butina_cluster_expandedFP(df, cutoff=0.35)
    >>> print(f"Found {len(clusters)} clusters")
    """
    # Identify fingerprint columns
    fingerprint_cols = [c for c in df.columns if c.startswith(fp_prefix)]
    if not fingerprint_cols:
        raise ValueError(f"No fingerprint columns found with prefix '{fp_prefix}'")

    n_mols = len(df)
    n_bits = len(fingerprint_cols)

    print(f"Clustering {n_mols} molecules with Butina algorithm...")
    print(f"Cutoff: {cutoff} (Tanimoto distance)")
    print(f"Using {n_bits} fingerprint bits")

    print("Converting fingerprints to RDKit BitVects...")

    fps = []
    for _, row in df[fingerprint_cols].iterrows():
        bv = DataStructs.ExplicitBitVect(n_bits)
        on_bits = np.flatnonzero(row.values)
        for bit in on_bits:
            bv.SetBit(int(bit))
        fps.append(bv)

    print("Clustering (lazy Tanimoto distance)...")

    def tanimoto_distance(fp1, fp2):
        return 1.0 - TanimotoSimilarity(fp1, fp2)

    clusters = Butina.ClusterData(
        fps,
        n_mols,
        cutoff,
        isDistData=False,
        distFunc=tanimoto_distance,
    )

    clusters = sorted(clusters, key=len, reverse=True)

    cluster_ids = np.full(n_mols, -1, dtype=int)
    for cid, members in enumerate(clusters):
        for idx in members:
            cluster_ids[idx] = cid

    df_result = df.copy()
    df_result["ClusterID"] = cluster_ids

    if print_stats:
        sizes = np.array([len(c) for c in clusters])
        singletons = np.sum(sizes == 1)

        print("\n" + "=" * 60)
        print("CLUSTERING RESULTS")
        print("=" * 60)
        print(f"Total clusters: {len(clusters)}")
        print(f"Largest cluster: {sizes.max()}")
        print(f"Smallest cluster: {sizes.min()}")
        print(f"Average cluster size: {sizes.mean():.1f}")
        print(f"Median cluster size: {np.median(sizes):.0f}")
        print(f"Singletons: {singletons} ({100*singletons/len(clusters):.1f}%)")
        print("=" * 60)

        print("\nTop 10 largest clusters:")
        for i, c in enumerate(clusters[:10], 1):
            print(f"  Cluster {i:3d}: {len(c):4d} molecules")

    return df_result, clusters


def butina_cluster_fp(df, fp_column='FP', cutoff=0.35, print_stats=True):
    """
    Cluster molecules using Butina algorithm with RDKit fingerprint objects.
    
    More memory-efficient than explicit fingerprints. Suitable for large datasets.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing RDKit fingerprints
    fp_column : str, default='FP'
        Column containing fingerprints
    cutoff : float, default=0.35
        Tanimoto distance cutoff (1 - similarity)
    print_stats : bool, default=True
        Print clustering statistics
    
    Returns
    -------
    tuple
        (df_result, clusters) where:
        - df_result: DataFrame with ClusterID column
        - clusters: list of clusters (lists of indices)
    
    Notes
    -----
    Cutoff interpretation:
    - cutoff=0.2 → similarity threshold of 0.8 (very tight clusters)
    - cutoff=0.35 → similarity threshold of 0.65 (moderate)
    - cutoff=0.5 → similarity threshold of 0.5 (loose clusters)
    
    Examples
    --------
    >>> df_fps = compute_morgan_fps(df)
    >>> df_clustered, clusters = butina_cluster_fp(df_fps, cutoff=0.35)
    """
    if fp_column not in df.columns:
        raise ValueError(f"DataFrame must contain '{fp_column}' column")

    fps = df[fp_column].tolist()
    n_mols = len(fps)

    print(f"Clustering {n_mols} molecules with Butina algorithm...")
    print(f"Cutoff: {cutoff} (Tanimoto distance)")
    print(f"Equivalent similarity threshold: {1 - cutoff:.2f}")

    def tanimoto_distance(fp1, fp2):
        return 1.0 - TanimotoSimilarity(fp1, fp2)

    clusters = Butina.ClusterData(
        fps,
        n_mols,
        cutoff,
        isDistData=False,
        distFunc=tanimoto_distance,
    )

    clusters = sorted(clusters, key=len, reverse=True)

    cluster_ids = np.full(n_mols, -1, dtype=int)
    for cid, members in enumerate(clusters):
        for idx in members:
            cluster_ids[idx] = cid

    df_result = df.copy()
    df_result['ClusterID'] = cluster_ids

    if print_stats:
        sizes = np.array([len(c) for c in clusters])
        singletons = np.sum(sizes == 1)

        print("\n" + "=" * 60)
        print("CLUSTERING RESULTS")
        print("=" * 60)
        print(f"Total clusters: {len(clusters)}")
        print(f"Largest cluster: {sizes.max()}")
        print(f"Smallest cluster: {sizes.min()}")
        print(f"Average cluster size: {sizes.mean():.1f}")
        print(f"Median cluster size: {np.median(sizes):.0f}")
        print(f"Singletons: {singletons} ({100*singletons/len(clusters):.1f}%)")
        print("=" * 60)

        print("\nTop 10 largest clusters:")
        for i, c in enumerate(clusters[:10], 1):
            print(f"  Cluster {i:3d}: {len(c):4d} molecules")

    return df_result, clusters