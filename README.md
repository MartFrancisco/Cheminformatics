# ChemSpace Explorer: Diverse Amide Library Design Pipeline

A comprehensive cheminformatics pipeline for curating large molecular datasets, selecting diverse building blocks, and designing drug-like amide libraries through virtual combinatorial chemistry.

## 🎯 Project Overview

This pipeline processes large molecular databases (tested on 500k+ molecules) to:
1. **Curate** and standardize molecular structures
2. **Classify** molecules by functional groups (acids, amines, etc.)
3. **Select** diverse representative subsets using cheminformatics algorithms
4. **Design** combinatorial amide libraries via virtual coupling reactions
5. **Filter** products for drug-likeness (Lipinski's Ro5, PAINS, etc.)

## 🔬 Scientific Workflow
```
Raw Dataset (500k+ molecules)
         ↓
    Curation (desalt, tautomer standardization, deduplication)
         ↓
    Functional Group Classification
         ↓
    Filter for mono-functional acids & amines
         ↓
    Diversity Selection (MaxMin algorithm)
         ↓
    Virtual Amide Coupling (100 acids × 100 amines = 10k products)
         ↓
    Drug-likeness Filtering (Ro5, PAINS, atom filters)
         ↓
    Final Focused Library
```

## 🛠️ Key Features

### Molecular Curation
- **Salt removal** using RDKit's `SaltRemover`
- **Tautomer standardization** for consistent representation
- **InChIKey-based deduplication** to remove redundancy

### Intelligent Classification
- SMARTS-based functional group detection
- Priority-based classification for amide coupling reactions
- Identifies problematic functional groups (aldehydes, esters, etc.)

### Diversity Selection
- **Butina clustering** for scaffold diversity analysis
- **MaxMin algorithm** for optimal chemical space coverage
- **Morgan fingerprints** (ECFP4-like) for similarity calculations

### Quality Control
- Lipinski's Rule of Five compliance
- PAINS (Pan-Assay Interference Compounds) filtering
- Unwanted functional group removal
- Conventional atom filtering (H, C, N, O, F, P, S, Cl, Br, I)

## 📦 Installation

### Prerequisites

- Python 3.8 or higher
- pip package manager
- Git (for cloning the repository)

### Step 1: Clone the Repository
```bash
# Clone from GitHub
git clone https://github.com/MartFrancisco/cheminformatics.git

# Navigate into the project directory
cd cheminformatics
```

### Step 2: Install the Package
```bash
# Install in editable mode (recommended for development)
pip install -e .

# This will automatically install all required dependencies:
# - rdkit>=2023.3.1
# - pandas>=1.5.0
# - numpy>=1.23.0
# - scipy>=1.9.0
# - matplotlib>=3.6.0
# - mols2grid>=2.0.0
```

### Step 3: Verify Installation
```bash
# Test that the package was installed correctly
python -c "from cheminformatics import curate_dataframe; print('✅ Installation successful!')"
```

### Alternative: Install Dependencies Separately

If you prefer to install dependencies first:
```bash
# Install dependencies from requirements.txt
pip install -r requirements.txt

# Then install the package
pip install -e .
```

### Optional Dependencies
```bash
# For development tools (testing, linting)
pip install -e ".[dev]"

# For running Jupyter notebooks
pip install -e ".[notebooks]"
```

## 🚀 Quick Start

### Basic Example
```python
from chemspace import *
import pandas as pd

# 1. Load your dataset
df = pd.read_csv('your_molecules.csv')  # Must have 'SMILES' column

# 2. Curate the dataset
df_curated = curate_dataframe(df, smiles_col='SMILES')

# 3. Classify functional groups
df_classified = FG_classifier(df_curated)

# 4. Filter for single functional groups
acids = df_classified[df_classified['Classes'] == 'acid']
amines = df_classified[df_classified['Classes'] == 'amine']

acids_single = filter_single_functional_group(acids, 'acid')
amines_single = filter_single_functional_group(amines, 'amine')

# 5. Compute fingerprints
acids_fps = compute_morgan_fps(acids_single, radius=2, fpSize=1024)
amines_fps = compute_morgan_fps(amines_single, radius=2, fpSize=1024)

# 6. Select diverse subset
diverse_acids = maxmin_from_rdkit_fps(acids_fps, n_select=100, seed=42)
diverse_amines = maxmin_from_rdkit_fps(amines_fps, n_select=100, seed=42)

# 7. Enumerate amide library
from rdkit.Chem import AllChem

rxn_smarts = "[#7H1,#7H2:7].[#6:1](=[#8:3])[#8H1]>>[#7:7]-[#6:1](=[#8:3])"
rxn_mol = AllChem.ReactionFromSmarts(rxn_smarts)

acid_list = diverse_acids[['mol','ID']].values
amine_list = diverse_amines[['mol','ID']].values

products = enumerate_library(rxn_mol, [amine_list, acid_list])
df_products = pd.DataFrame(products, columns=['SMILES', 'Product_ID'])

# 8. Filter for drug-likeness
df_filtered = lipinski_ro5(df_products)
df_filtered = pains_filter(df_filtered)
df_filtered = filter_by_atoms(df_filtered)
df_filtered = filter_by_FGs(df_filtered)

print(f"Final library size: {len(df_filtered)}")
```

### Using Sample Data

The repository includes sample data for testing:
```python
from chemspace import *
import pandas as pd

# Load sample molecules
df = pd.read_csv('sample_data/example_molecules.csv')
print(f"Loaded {len(df)} sample molecules")

# Run the pipeline
df_curated = curate_dataframe(df)
df_classified = FG_classifier(df_curated)
# ... continue with your analysis
```

## 📊 Validation & Quality Metrics

The pipeline includes built-in validation functions:
```python
# Compare property distributions (original vs. selected)
compare_properties_distributions(original_df, selected_df, 
                                property_name='MW')

# Analyze scaffold diversity
scaffold_diversity_plot(original_df, selected_df)

# Check pairwise similarity distribution
pairwise_similarity_histogram(selected_df, sample_size=100)
```

## 📁 Project Structure
```
cheminformatics/
├── README.md
├── requirements.txt
├── setup.py
├── pyproject.toml
├── LICENSE
├── .gitignore
│
├── cheminformatics/                # Core package (renamed from chemspace/)
│   ├── __init__.py
│   ├── curation.py                 # Dataset curation
│   ├── fingerprints.py             # Fingerprint generation
│   ├── clustering.py               # Butina clustering
│   ├── diversity.py                # MaxMin selection
│   ├── validation.py               # Quality metrics
│   ├── enumeration.py              # Library enumeration
│   └── filters.py                  # Drug-likeness filters
│
├── notebooks/
│   └── example_workflow.ipynb      # Full workflow demonstration
│
└── sample_data/
    └── example_molecules.csv       # Small test dataset
```

## 🔍 Function Reference

### Core Functions

| Function | Purpose |
|----------|---------|
| `curate_dataframe()` | Standardize, desalt, and deduplicate molecules |
| `FG_classifier()` | Classify molecules by functional groups |
| `filter_single_functional_group()` | Keep only mono-functional molecules |
| `compute_morgan_fps()` | Generate Morgan fingerprints |
| `butina_cluster_fp()` | Perform Butina clustering |
| `maxmin_from_rdkit_fps()` | Select diverse subset via MaxMin |
| `enumerate_library()` | Generate combinatorial product library |
| `lipinski_ro5()` | Filter by Lipinski's Rule of Five |
| `pains_filter()` | Remove PAINS compounds |
| `filter_by_atoms()` | Keep only conventional atoms |
| `filter_by_FGs()` | Remove unwanted functional groups |

## 📖 Example Results

From a 500k molecule database:
- **Curated**: 450k unique molecules (after deduplication)
- **Acids**: 50k → 100 diverse representatives
- **Amines**: 130k → 100 diverse representatives
- **Virtual library**: 10,000 amides generated
- **Drug-like library**: ~4,000 compounds (after filtering)

## 🐛 Troubleshooting

### RDKit Installation Issues

If you encounter issues installing RDKit via pip, try using conda:
```bash
# Create a conda environment
conda create -n cheminformatics python=3.10
conda activate cheminformatics

# Install RDKit via conda
conda install -c conda-forge rdkit

# Then install the package
cd cheminformatics
pip install -e .
```

### Import Errors

If you get `ModuleNotFoundError: No module named 'cheminformatics'`:
```bash
# Make sure you're in the project directory
cd /path/to/cheminformatics

# Reinstall in editable mode
pip install -e .
```

### Permission Errors

If you get permission errors during installation:
```bash
# Install with user flag
pip install -e . --user
```

## ⚠️ Note on Data

This repository contains **code only**. The example dataset is a small synthetic set for demonstration purposes. Users should supply their own molecular databases (SDF, CSV with SMILES, etc.).

Commercial datasets like **eMolecules**, **ZINC**, or **ChEMBL** can be used with this pipeline, subject to their respective licenses.

## 🤝 Contributing

Contributions welcome! Areas for improvement:
- [ ] Add support for other reaction types
- [ ] Implement GPU-accelerated fingerprint calculations
- [ ] Add more diversity selection algorithms (sphere exclusion, OptiSim)
- [ ] Create visualization dashboard

## 📄 License

MIT License - see [LICENSE](LICENSE) file for details

## 📧 Contact

Questions? Open an issue or contact [francisco.qui.martins@gmail.com]

## 🙏 Acknowledgments

This project is built using the following open-source libraries:

- [RDKit](https://www.rdkit.org/) - Cheminformatics and machine learning toolkit (BSD-3-Clause)
- [Pandas](https://pandas.pydata.org/) - Data manipulation and analysis (BSD-3-Clause)
- [NumPy](https://numpy.org/) - Numerical computing (BSD-3-Clause)
- [SciPy](https://scipy.org/) - Scientific computing (BSD-3-Clause)
- [Matplotlib](https://matplotlib.org/) - Plotting and visualization (PSF-based)
- [mols2grid](https://github.com/cbouy/mols2grid) - Interactive molecular visualization (Apache-2.0)

Special thanks to the open-source cheminformatics community!

---

**Note**: This is an educational/research tool. Always validate computational results experimentally.