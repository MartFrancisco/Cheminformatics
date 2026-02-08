from setuptools import setup, find_packages

# Read the README for the long description
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="chemspace-explorer",
    version="0.1.0",
    author="Francisco Martins",
    author_email="francisco.qui.martins@gmail.com",
    description="Diverse amide library design pipeline for cheminformatics",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/chemspace-explorer",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Chemistry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
    ],
    python_requires=">=3.8",
    install_requires=[
        "rdkit>=2023.3.1",
        "pandas>=1.5.0",
        "numpy>=1.23.0",
        "scipy>=1.9.0",
        "matplotlib>=3.6.0",
        "mols2grid>=2.0.0",
    ],
    extras_require={
        "dev": ["pytest>=7.0", "black", "flake8"],
        "notebooks": ["jupyter>=1.0.0", "ipykernel>=6.0.0"],
    },
)