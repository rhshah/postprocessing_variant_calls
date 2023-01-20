# Main Import
import typer
from .vardict import vardict_process
from .annotate import annotate_process
from .concat import concat_process
import logging
import time
# setup logger
logger = logging.getLogger("filter")
app = typer.Typer()

# Vardict filter App 
app.add_typer(vardict_process.app, name="vardict", help="post-processing commands for VarDict version 1.4.6 VCFs.")

# Add Annote App 
app.add_typer(annotate_process.app, name="annotate", help="annotate maf files based on a given input. Currently supports bed and maf files as references.")

# Add Concat App 
app.add_typer(concat_process.app, name="concat", help="concat common variant file formats")

# Haplo filter

# Main App
if __name__ == "__main__":
    app()
