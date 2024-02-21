# Main Import
import typer
from .vardict import vardict_process
from .mutect import mutect_process

# from .maf import main
from .maf import main
from .maf.annotate import annotate_process
import logging
import time

# app
app = typer.Typer()

# Vardict filter App
app.add_typer(
    vardict_process.app,
    name="vardict",
    help="post-processing commands for VarDict version 1.4.6 VCFs.",
)

# muTect filter App
app.add_typer(
    mutect_process.app,
    name="mutect1",
    help="post-processing commands for MuTect version 1.1.5 VCFs.",
)

# Add Annote Maf
app.add_typer(
    main.app,
    name="maf",
    help="operations for manipulating maf files based on a given input.",
)

# versioning
__version__ = "0.2.4"


def version_callback(value: bool):
    if value:
        typer.echo(f"{__version__}")
        raise typer.Exit()


@app.callback()
def main(
    version: bool = typer.Option(
        None, "--version", callback=version_callback, is_eager=True
    ),
):
    return


# Main App
if __name__ == "__main__":
    app()
