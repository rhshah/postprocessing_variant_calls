# Post-processing of variant calls

This package provides a variety of commands for manipulating different types of common outputs (e.g. mafs, vcf and txt files) from different bioinformatic variant callers such as mutect and vardict. 

Supported File Types:
- [maf](docs/MAF.md) 
- [vardict](docs/VARDICT.md)

# Installation  

For general use you can run: `pip install postprocessing_variant_calls`
or a tagged version with `pip install git+https://github.com/msk-access/postprocessing_variant_calls.git@<version>`

For setting up a development environment please see the [Setting up a Dev Environment](#Setting-up-a-Dev-Environment) section.

# Usage

See [CLI](docs/CLI.md) for commmand line usage of the package.

# Setting up a Dev Environment 

## Install External Dependencies
Have an environment with python >= 3.8 installed. 

Install poetry: 

```bash
pip install poetry
```

## Install Package Dependencies

Then install project dependencies with Poetry.

```bash
cd /path/to/postprocessing_variant_calls
poetry install .
```

## Accessing Environment

To access the environment after initial setup up run: 

```bash
poetry shell
```
