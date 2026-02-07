"""
Molecular fingerprint generation functions.

Provides utilities for computing Morgan (ECFP-like) fingerprints in both
explicit bit vector and RDKit object formats.
"""

import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator


def compute_morgan_fps_explicit(df, column_name='SMILES', radius=2, fpSize=1024):
    """
    Compute Morgan fingerprints using the modern MorganGenerator API.
    
    Generates columns containing each bit as a separate column (bit_0, bit_1, ...).
    Useful for machine learning applications that require explicit feature vectors.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing molecules' SMILES
    column_name : str, default='SMILES'
        Name of the column containing SMILES strings
    radius : int, default=2
        The radius parameter for Morgan fingerprints (radius=2 ≈ ECFP4)
    fpSize : int, default=1024
        The size of the fingerprint in bits
    
    Returns
    -------
    pd.DataFrame
        DataFrame with Morgan fingerprint bits as new columns (bit_0 to bit_{fpSize-1})
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'CCN']})
    >>> df_fps = compute_morgan_fps_explicit(df, fpSize=512)
    >>> print(df_fps.shape)
    (2, 513)  # Original column + 512 bit columns
    """
    # Ensure the DataFrame contains the specified column
    if column_name not in df.columns:
        raise ValueError(f"The DataFrame must contain the column '{column_name}'.")
    
    print(f"Computing Morgan fingerprints for {len(df)} molecules...")
    print(f"Parameters: radius={radius}, fpSize={fpSize}")
    
    # Create MorganGenerator
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize)
    
    def smiles_to_fp(smiles):
        """Convert SMILES to Morgan fingerprint bit vector."""
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol:
                fp = morgan_gen.GetFingerprint(mol)
                return np.array(fp, dtype=np.uint8)
            else:
                return np.zeros(fpSize, dtype=np.uint8)
        except Exception as e:
            print(f"Error processing {smiles}: {e}")
            return np.zeros(fpSize, dtype=np.uint8)
    
    # Compute fingerprints for all molecules
    fps = df[column_name].apply(smiles_to_fp)
    
    # Convert to DataFrame with proper column names
    fp_df = pd.DataFrame(
        fps.tolist(),
        columns=[f'bit_{i}' for i in range(fpSize)],
        index=df.index
    )
    
    # Concatenate with original DataFrame
    result_df = pd.concat([df, fp_df], axis=1)
    
    print(f"Fingerprint computation complete!")
    print(f"DataFrame shape: {result_df.shape}")
    
    return result_df


def compute_morgan_fps(df, column_name='SMILES', radius=2, fpSize=1024, fp_column='FP'):
    """
    Compute Morgan fingerprints and store them as RDKit fingerprint objects.
    
    Stores fingerprints as RDKit objects in a single DataFrame column, which is
    more memory-efficient and compatible with RDKit's similarity functions.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing molecules' SMILES
    column_name : str, default='SMILES'
        Column containing SMILES
    radius : int, default=2
        Morgan radius (radius=2 ≈ ECFP4)
    fpSize : int, default=1024
        Fingerprint size in bits
    fp_column : str, default='FP'
        Name of the output fingerprint column
    
    Returns
    -------
    pd.DataFrame
        DataFrame with an added column containing RDKit fingerprint objects
    
    Examples
    --------
    >>> df = pd.DataFrame({'SMILES': ['CCO', 'CCN']})
    >>> df_fps = compute_morgan_fps(df)
    >>> print(df_fps.columns.tolist())
    ['SMILES', 'FP']
    """
    if column_name not in df.columns:
        raise ValueError(f"DataFrame must contain '{column_name}' column")

    print(f"Computing Morgan fingerprints for {len(df)} molecules...")
    print(f"Parameters: radius={radius}, fpSize={fpSize}")

    # Create Morgan fingerprint generator
    morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fpSize)

    def smiles_to_fp(smiles):
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            return morgan_gen.GetFingerprint(mol)
        except Exception as e:
            print(f"Error processing {smiles}: {e}")
            return None

    df_fp = df.copy()
    df_fp[fp_column] = df_fp[column_name].apply(smiles_to_fp)

    n_invalid = df_fp[fp_column].isna().sum()
    print(f"Fingerprint computation complete!")
    print(f"Invalid molecules (no FP): {n_invalid}")

    return df_fp