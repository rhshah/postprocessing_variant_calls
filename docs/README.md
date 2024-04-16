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

# Contributing to Documentation

The Gitbook for this repository is configured so changes are written in Gitbook and synced with the `docs` branch. 

To contribute to the documentation, you can write your changes in [Gitbook](https://app.gitbook.com/o/-LhMNgvjydB3TFWAUMVb/s/VBp8SqbRAs28AQCVNoIS/), request a review, and merge the changes. Keep in mind, you will need access to the organization to contribute. 

Each file-type supported should have a section in the Gitbook detailing the implementation of the file-type and a justification of it's operations. For example, the `maf` file-type has it's own section, which includes a description of how a maf is defined internally in the package and a justification of it's operations and how to use them.

Beyond file-type sections, you will also notice a section called `cli`, which lists all commands in the `postprocessing_variant_calls` package. Do not manually edit this section. This section is created using the `typer-cli` package, which uses the typer `help` parameters specified in typer commands to generate documentation. It is automatically updated by the git-action, `.github/workflows/document_package.yml` upon a push to the `main` branch. To make sure the `cli.md` document updates to include newly added commands, specify all relevant typer `help` parameters.

If you'd like to see a mock of your typer commands as they'll appear in the `cli.md` document, you can run: `poetry run typer postprocessing_variant_calls.main utils docs > docs/cli.md` in a properly configured dev environment. Note that this file is in the git-ignore and should not be included in your PRs. The `cli.md` should be only updated through the git-action, `.github/workflows/document_package.yml`.
