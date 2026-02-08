"""
ChemSpace Explorer: Diverse Amide Library Design Pipeline

A comprehensive cheminformatics toolkit for curating molecular datasets,
selecting diverse building blocks, and designing drug-like combinatorial libraries.

Main modules:
    - curation: Dataset cleaning, standardization, and functional group classification
    - fingerprints: Morgan fingerprint generation (explicit and object formats)
    - clustering: Butina clustering for diversity analysis
    - diversity: MaxMin selection for diverse subset picking
    - validation: Chemical space coverage and diversity metrics
    - enumeration: Combinatorial library generation
    - filters: Drug-likeness filters (Ro5, PAINS, etc.)

Example workflow:
    >>> from chemspace import *
    >>> import pandas as pd
    >>> 
    >>> # Load and curate dataset
    >>> df = pd.read_csv('molecules.csv')
    >>> df_curated = curate_dataframe(df)
    >>> 
    >>> # Classify and filter
    >>> df_classified = FG_classifier(df_curated)
    >>> acids = df_classified[df_classified['Classes'] == 'acid']
    >>> acids_single = filter_single_functional_group(acids, 'acid')
    >>> 
    >>> # Select diverse subset
    >>> acids_fps = compute_morgan_fps(acids_single)
    >>> diverse_acids = maxmin_from_rdkit_fps(acids_fps, n_select=100, seed=42)
    >>> 
    >>> # Enumerate and filter library
    >>> products = enumerate_library(rxn_mol, [amines, acids])
    >>> df_products = pd.DataFrame(products, columns=['SMILES', 'Product_ID'])
    >>> df_filtered = lipinski_ro5(df_products)
"""

__version__ = "0.1.0"
__author__ = "Francisco Martins"
__email__ = "francisco.qui.martins@gmail.com"

# Curation functions
from .curation import (
    curate_dataframe,
    FG_classifier,
    filter_single_functional_group
)

# Fingerprint functions
from .fingerprints import (
    compute_morgan_fps_explicit,
    compute_morgan_fps
)

# Clustering functions
from .clustering import (
    butina_cluster_expandedFP,
    butina_cluster_fp
)

# Diversity selection functions
from .diversity import (
    maxmin_from_fingerprints,
    maxmin_from_rdkit_fps
)

# Validation functions
from .validation import (
    compare_properties_distributions,
    scaffold_diversity_plot,
    pairwise_similarity_histogram
)

# Enumeration functions
from .enumeration import (
    enumerate_library
)

# Filter functions
from .filters import (
    lipinski_ro5,
    pains_filter,
    filter_by_atoms,
    filter_by_FGs
)

__all__ = [
    # Curation
    'curate_dataframe',
    'FG_classifier',
    'filter_single_functional_group',
    
    # Fingerprints
    'compute_morgan_fps_explicit',
    'compute_morgan_fps',
    
    # Clustering
    'butina_cluster_expandedFP',
    'butina_cluster_fp',
    
    # Diversity
    'maxmin_from_fingerprints',
    'maxmin_from_rdkit_fps',
    
    # Validation
    'compare_properties_distributions',
    'scaffold_diversity_plot',
    'pairwise_similarity_histogram',
    
    # Enumeration
    'enumerate_library',
    
    # Filters
    'lipinski_ro5',
    'pains_filter',
    'filter_by_atoms',
    'filter_by_FGs',
]